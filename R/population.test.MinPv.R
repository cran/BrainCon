#' The one-sample population inference using Genovese and Wasserman's method
#'
#' Identify the nonzero partial correlations in one-sample population,
#' based on controlling the rate of the false discovery proportion (FDP) exceeding \eqn{c0}
#' at \eqn{\alpha}. The method is based on the minimum of the p-values.
#' Input a \code{popEst} class object returned by \code{\link{population.est}}.
#' \cr
#' \cr
#'
#'@param popEst A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of FDP,
#'default value is \code{0.1}.
#'@param targetSet a two-column matrix. Each row contains two index corresponding to a pair of variables of interest.
#'If \code{NULL}, any pair of two variables is considered to be of interest.
#'@param simplify a logical indicating whether results should be simplified if possible.
#'
#'@return If \code{simplify} is \code{FALSE}, a \eqn{p*p} matrix with values 0 or 1 is returned, and 1 means nonzero.
#'
#'And if \code{simplify} is \code{TRUE}, a two-column matrix is returned,
#'indicating the row index and the column index of recovered nonzero partial correlations.
#'Those with lower p values are sorted first.
#'
#'@seealso \code{\link{population.test}}.
#'
#'@examples
#' ## Quick example for the one-sample population inference
#' data(popsimA)
#' # estimating partial correlation coefficients
#' pc = population.est(popsimA)
#' # conducting hypothesis test
#' Res  = population.test.MinPv(pc)
#'
#' @references
#' Genovese C. and Wasserman L. (2006).
#' Exceedance Control of the False Discovery Proportion,
#' \emph{Journal of the American Statistical Association}, 101, 1408-1417.
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

population.test.MinPv <- function(popEst, alpha = 0.05, c0 = 0.1, targetSet = NULL, simplify = !is.null(targetSet)){
  force(simplify)
  if (class(popEst) != 'popEst')
    stop("The argument popEst requires a 'popEst' class input!\n")
  EstAll = popEst$coef
  p = nrow(EstAll)
  MC = length(popEst[['ind.est']])
  if (is.null(targetSet)){
    targetSet = which(upper.tri(EstAll), arr.ind = T)
    Mp = p * (p - 1) / 2
  } else {
    simplify = TRUE
    targetSet = normalize.set(targetSet, p)
    Mp = nrow(targetSet)
  }
  CAll = array(dim = c(p, p, MC))
  for (sub in 1 : MC){
    CAll[, , sub] = popEst[['ind.est']][[sub]][['coef']]
  }
  SdAll = apply(CAll, c(1, 2), sd)
  EstT = sqrt(MC) * EstAll / (SdAll + 1e-6)
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
