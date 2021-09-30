#' The one-sample population inference using Genovese and Wasserman's method
#'
#' Identify the nonzero partial correlations in one-sample population,
#' based on controlling the exceedance rate of the false discovery proportion (FDP)
#' at \eqn{\alpha=0.05}. The method is based on the minimum of the p-values.
#' Input  data \eqn{Z} , contains values of p interested variables.
#' \cr
#' \cr
#'
#'@param popEst A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of the false discovery proportion (FDP),
#'default value is \code{0.1}.
#'@return A  \eqn{p*p} matrix with values 0 or 1.
#'
#'@examples
#' ## Quick example for the one-sample population inference
#' data(popsimA)
#' pc = population.est(popsimA)            # estimating partial correlation coefficients
#' Res  = population.test.MinPv(pc)        # conducting hypothesis test
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

population.test.MinPv <- function(popEst, alpha = 0.05, c0 = 0.1){
  if (class(popEst) != 'popEst')
    stop("The argument popEst requires a 'popEst' class input!\n")
  EstAll = popEst$coef
  n = nrow(popEst[['ind.est']][[1]][['asym.ex']])
  MC = length(popEst[['ind.est']])
  p = nrow(EstAll)
  Mp = p * (p - 1) / 2
  CAll = array(dim = c(p, p, MC))
  for (sub in 1 : MC){
    CAll[, , sub] = popEst[['ind.est']][[sub]][['coef']]
  }
  SdAll = apply(CAll, c(1, 2), sd)
  EstT = sqrt(MC) * EstAll / (SdAll + 1e-6)
  pv0 = 2 * (1 - pnorm(abs(EstT)))
  pv1 = sort(pv0[upper.tri(pv0, diag = FALSE)])
  Beta = qbeta(alpha, 1, Mp:1)
  a0 = which(pv1 > Beta)[1]
  a1 = ceiling(a0 / (1 - c0))
  pvThreshold = ifelse(is.na(a1), pv1[1]/2, pv1[a1])
  MinPv = 1 * (pv0 < pvThreshold)
  return(MinPv)
}

