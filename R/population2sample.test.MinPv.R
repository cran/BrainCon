#' Identify differences of partial correlations between two populations using Genovese and Wasserman's method
#'
#' Identify differences of partial correlations between two populations
#' in two groups of time series data,
#' based on controlling the rate of the false discovery proportion (FDP) exceeding \eqn{c0}
#' at \eqn{\alpha}. The method is based on the minimum of the p-values.
#' Input two \code{popEst} class objects returned by \code{\link{population.est}}
#' (the number of individuals in two groups can be different).
#' \cr
#' \cr
#'
#'@param popEst1 A \code{popEst} class object.
#'@param popEst2 A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of FDP,
#'default value is \code{0.1}.
#'@param targetSet a two-column matrix. Each row contains two index corresponding to a pair of variables of interest.
#'If \code{NULL}, any pair of two variables is considered to be of interest.
#'@param simplify a logical indicating whether results should be simplified if possible.
#'
#'@return If \code{simplify} is \code{FALSE}, a \eqn{p*p} matrix with values 0 or 1 is returned, and 1 means unequal.
#'
#'And if \code{simplify} is \code{TRUE}, a two-column matrix is returned,
#'indicating the row index and the column index of recovered unequal partial correlations.
#'Those with lower p values are sorted first.
#'
#'@examples
#' ## Quick example for the two-sample case inference
#' data(popsimA)
#' data(popsimB)
#' # estimating partial correlation coefficients by lasso (scaled lasso does the same)
#' pc1 = population.est(popsimA, type = 'l')
#' pc2 = population.est(popsimB, type = 'l')
#' # conducting hypothesis test
#' Res = population2sample.test.MinPv(pc1, pc2)
#'
#' @references
#' Genovese C., and Wasserman L. (2006).
#' Exceedance Control of the False Discovery Proportion,
#' \emph{Journal of the American Statistical Association}, 101, 1408-1417
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

population2sample.test.MinPv <- function(popEst1, popEst2, alpha = 0.05, c0 = 0.1, targetSet = NULL, simplify = !is.null(targetSet)){
  force(simplify)
  if (!inherits(popEst1, 'popEst') | !inherits(popEst1, 'popEst'))
    stop("The arguments popEst1 and popEst2 require 'popEst' class inputs!\n")
  EstAll1 = popEst1$coef
  EstAll2 = popEst2$coef
  p = nrow(EstAll1)
  MC1 = length(popEst1[['ind.est']])
  MC2 = length(popEst2[['ind.est']])
  if (is.null(targetSet)){
    targetSet = which(upper.tri(EstAll1), arr.ind = T)
    Mp = p * (p - 1) / 2
  } else {
    simplify = TRUE
    targetSet = normalize.set(targetSet, p)
    Mp = nrow(targetSet)
  }
  CAll1 = array(dim = c(p, p, MC1))
  CAll2 = array(dim = c(p, p, MC2))
  for (sub in 1 : MC1)
    CAll1[, , sub] = popEst1[['ind.est']][[sub]][['coef']]
  for (sub in 1 : MC2)
    CAll2[, , sub] = popEst2[['ind.est']][[sub]][['coef']]
  SdAll1 = apply(CAll1, c(1, 2), sd)
  SdAll2 = apply(CAll2, c(1, 2), sd)
  EstT = (EstAll1 - EstAll2) / sqrt(SdAll1^2 / MC1 + SdAll2^2 / MC2 + 1e-6)
  pv0 = 2 * (1 - pnorm(abs(EstT)))
  pv1 = sort(pv0[targetSet])
  Beta = qbeta(alpha, 1, Mp:1)
  a0 = which(pv1 > Beta)[1]
  a1 = ceiling(a0 / (1 - c0))
  pvThreshold = ifelse(is.na(a1), pv1[1]/2, pv1[a1])
  if (simplify){
    index = which(pv0[targetSet] < pvThreshold)
    ord.index = index[order(pv0[targetSet][index])]
    return(matrix(targetSet[ord.index, ], ncol=2,
                  dimnames = list(NULL, c("row", "col"))))
  }
  MinPv = 1 * (pv0 < pvThreshold)
  return(MinPv)
}

