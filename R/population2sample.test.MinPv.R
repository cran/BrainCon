#' Identify differences of partial correlations between two populations using Genovese and Wasserman's method
#'
#' Identify differences of partial correlations between two populations
#' in two groups of time series data
#' by controlling the exceedance rate of the false discovery proportion (FDP)
#' at \eqn{\alpha=0.05}. The method is based on the minimum of the p-values.
#' Input two groups of data \eqn{Z_1} and \eqn{Z_2}, each contains values of p interested
#' variables of individuals
#' (the number of individuals in two groups can be different) over n periods.
#' \cr
#' \cr
#'
#'@param popEst1 A \code{popEst} class object.
#'@param popEst2 A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of the false discovery proportion (FDP),
#'default value is \code{0.1}.
#'The choice of \code{c0} depends on the empirical problem. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'
#'@return A \eqn{p*p} matrix with values 0 or 1.
#'If the j-th row and k-th column of the matrix is 1,
#'then the partial correlation coefficients between
#'the j-th variable and the k-th variable in two populations
#'are identified to be unequal.
#'
#'@examples
#' ## Quick example for the two-sample case inference
#' data(popsimA)
#' data(popsimB)
#' pc1 = population.est(popsimA)
#' pc2 = population.est(popsimB)
#' Res = population2sample.test.MinPv(pc1, pc2)         # conducting hypothesis test
#'
#' @references
#' Genovese C., and Wasserman L. (2006).
#' Exceedance Control of the False Discovery Proportion,
#' \emph{Journal of the American Statistical Association}, 101, 1408-1417
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

population2sample.test.MinPv <- function(popEst1, popEst2, alpha = 0.05, c0 = 0.1){
  if (class(popEst1) != 'popEst' | class(popEst2) != 'popEst')
    stop("The arguments popEst1 and popEst2 require 'popEst' class inputs!\n")
  EstAll1 = popEst1$coef
  EstAll2 = popEst2$coef
  n = nrow(popEst1[['ind.est']][[1]][['asym.ex']])
  p = nrow(EstAll1)
  MC1 = length(popEst1[['ind.est']])
  MC2 = length(popEst2[['ind.est']])
  Mp = p * (p - 1) / 2
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
  pv1 = sort(pv0[upper.tri(pv0, diag = FALSE)])
  Beta = qbeta(alpha, 1, Mp:1)
  a0 = which(pv1 > Beta)[1]
  a1 = ceiling(a0 / (1 - c0))
  pvThreshold = ifelse(is.na(a1), pv1[1]/2, pv1[a1])
  MinPv = 1 * (pv0 < pvThreshold)
  return(MinPv)
}
