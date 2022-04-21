#' Estimate population-level partial correlation coefficients
#'
#' Estimate population-level partial correlation coefficients in time series data.
#' And also return coefficients for each individual.
#' Input time series data for population as a 3-dimensional array or a list.
#' \cr
#' \cr
#'
#'@param Z If each individual shares the same number of periods of time, \code{Z} can be a \eqn{n*p*m} dimensional array, where \eqn{m} is number of individuals.
#'In general, \code{Z} should be a m-length list, and each element in the list is a \eqn{n_i*p} matrix, where \eqn{n_i} stands for the number of periods of time of the i-th individual.
#'@param lambda a scalar or a m-length vector, representing the penalty parameters of order \eqn{\sqrt{\log(p)/n_i}} for each individual.
#'If a scalar, the penalty parameters used in each individual are the same.
#'If a m-length vector, the penalty parameters for each individual are specified in order.
#'And if \code{NULL}, penalty parameters are specified by \code{type}.
#'More details about the penalty parameters are in \code{\link{individual.est}}.
#'@param type  a character string representing the method of estimation. \code{"slasso"} means scaled lasso, and \code{"lasso"} means lasso. Default value is \code{"slasso"}.
#'@param alpha a numeric scalar, default value is \code{0.05}. It is used when \code{ind.ci} is \code{TRUE}.
#'@param ind.ci  a logical indicating whether to compute \eqn{1-\alpha} confidence intervals of each subject, default value is \code{FALSE}.
#'
#'@return A \code{popEst} class object containing two components.
#'
#' \code{coef} a \eqn{p*p} partial correlation coefficients matrix.
#'
#' \code{ind.est} a \eqn{m}-length list, containing estimates for each individuals.
#'
#' \code{type} regression type in estimation.
#'
#'@examples
#' ## Quick example for the population-level estimates
#' data(popsimA)
#' # estimating partial correlation coefficients by scaled lasso
#' pc = population.est(popsimA)
#'
#' ## Inference on the first subject in population
#' Res_1 = individual.test(pc$ind.est[[1]])
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

population.est <- function(Z, lambda = NULL, type = c("slasso", "lasso"), alpha = 0.05, ind.ci = FALSE){
  if (!is.array(Z) & !is.list(Z))
    stop("The argument Z requires an 3-D array or a list!")
  if (is.array(Z))
    Z = lapply(apply(Z, 3, list), '[[', 1)
  n = sapply(Z, nrow)
  p = unique(sapply(Z, ncol))
  if (length(p)>1)
    stop("Each individual has to have the same number of variables!")
  MC = length(Z)
  type = match.arg(type)
  if (length(lambda) == 0){
    if (type == "slasso"){
      lambda = sqrt(2 * 2.01 * log(p * (log(p))^(1.5) / sqrt(n)) / n)
    } else if (type == "lasso"){
      lambda = sqrt(2 * log(p) / n)
    }
  } else if (length(lambda) == 1){
    lambda = rep(lambda, MC)
  } else if (length(lambda) != MC){
    stop("The argument lambda requires a scalar or a m-length vector!")
  }
  ind.est = list()
  CAll = array(dim = c(p, p, MC))
  for (sub in 1 : MC){
    ind.est[[sub]] = individual.est(Z[[sub]], lambda = lambda[sub], type = type, alpha = alpha, ci = ind.ci)
    CAll[, , sub] = ind.est[[sub]][['coef']]
  }
  Est = apply(CAll, c(1, 2), mean)
  return(structure(list(coef = Est, ind.est = ind.est, type = type), class = 'popEst'))
}
