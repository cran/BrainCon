#' Estimate population-level partial correlation coefficients
#'
#' Estimate population-level partial correlation coefficients in time series data.
#' And also return each individual-level coefficients.
#' \cr
#' \cr
#'
#'@param Z is a \eqn{n*p*m} dimensional array, where \eqn{m} is number of subjects.
#'@param alpha significance level, default value is \code{0.05}.
#'@param lambda a penalty parameter used in lasso of order \code{sqrt(log(p)/n)}, if \code{NULL}, \code{2*sqrt(log(p)/n)} will be used.
#'@param ind.ci  a logical indicating whether to compute \eqn{1-\alpha} confidence interval of each subject, default value is \code{FALSE}.
#'
#'@return A \code{popEst} class object containing two components.
#'
#' \code{coef} a \eqn{p*p} partial correlation coefficients matrix.
#'
#' \code{ind.est} a \eqn{m}-length list, containing estimates for each individuals.
#'
#'
#'@examples
#' ## Quick example for the individual-level estimates
#' data(popsimA)
#' pc = population.est(popsimA)        # estimating partial correlation coefficients
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

population.est <- function(Z, alpha = 0.05, lambda = NULL, ind.ci = FALSE){
  n = dim(Z)[1]
  p = dim(Z)[2]
  MC = dim(Z)[3]
  if (is.null(lambda)) lambda = 2 * sqrt(log(p) / n)
  ind.est = list()
  CAll = array(dim = c(p, p, MC))
  for (sub in 1 : MC){
    ind.est[[sub]] = individual.est(Z[, , sub], alpha = alpha, lambda = lambda, ci = ind.ci)
    CAll[, , sub] = ind.est[[sub]][['coef']]
  }
  Est = apply(CAll, c(1, 2), mean)
  return(structure(list(coef = Est, ind.est = ind.est), class = 'popEst'))
}
