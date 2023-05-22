#' Estimate individual-level partial correlation coefficients
#'
#' Estimate individual-level partial correlation coefficients in time series data
#' with \eqn{1-\alpha} confidence intervals.
#' Note that these are confidence intervals for single parameters, not simultaneous confidence intervals.
#' \cr
#' \cr
#'
#'@param X time series data of an individual which is a \eqn{n*p} numeric matrix, where \eqn{n} is the number of periods of time and \eqn{p} is the number of variables.
#'@param lambda a penalty parameter of order \eqn{\sqrt{\log(p)/n}}.
#'If \code{NULL}, \eqn{\sqrt{2*2.01/n*\log(p*(\log(p))^{1.5}/n^{0.5})}} is used in scaled lasso, and \eqn{\sqrt{2*\log(p)/n}} is used in lasso.
#'Increasing the penalty parameter may lead to larger residuals in the node-wise regression,
#'causing larger absolute values of estimates of partial correlation coefficients, which may cause more false positives in subsequent tests.
#'@param type  a character string representing the method of estimation. \code{"slasso"} means scaled lasso, and \code{"lasso"} means lasso. Default value is \code{"slasso"}.
#'@param alpha significance level, default value is \code{0.05}.
#'@param ci  a logical indicating whether to compute \eqn{1-\alpha} confidence interval, default value is \code{TRUE}.
#'
#'@return An \code{indEst} class object containing two or four components.
#'
#' \code{coef} a \eqn{p*p} partial correlation coefficients matrix.
#'
#' \code{ci.lower} a \eqn{p*p} numeric matrix containing the lower bound of \eqn{1-\alpha} confidence interval,
#' returned if \code{ci} is \code{TRUE}.
#'
#' \code{ci.upper} a \eqn{p*p} numeric matrix containing the upper bound of \eqn{1-\alpha} confidence interval,
#' returned if \code{ci} is \code{TRUE}.
#'
#' \code{asym.ex} a matrix measuring the asymptotic expansion of estimates, which will be used for multiple tests.
#'
#' \code{type} regression type in estimation.
#'
#'@seealso \code{\link{population.est}}.
#'
#'@examples
#' ## Quick example for the individual-level estimates
#' data(indsim)
#' # estimating partial correlation coefficients by scaled lasso
#' pc = individual.est(indsim)
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.
#' @references
#' Sun T. and Zhang C. (2012).
#' Scaled Sparse Linear Regression,
#' \emph{Biometrika}, 99, 879–898.
#' @references
#' Liu W. (2013).
#' Gaussian Graphical Model Estimation With False Discovery Rate Control,
#' \emph{The Annals of Statistics}, 41, 2948–2978.
#' @references
#' Ren Z., Sun T., Zhang C. and Zhou H. (2015).
#' Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models,
#' \emph{The Annals of Statistics}, 43, 991–1026.

individual.est <- function(X, lambda = NULL, type = c("slasso", "lasso"), alpha = 0.05, ci = TRUE){
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]
  Mp = p * (p - 1) / 2
  X = scale(X, scale = FALSE)
  XS = scale(X, center = FALSE)
  sdX = apply(X, 2, sd)
  if (min(sdX) == 0)
    stop("The argument X should not have any constant column!\n")
  Eresidual = matrix(0, n, p)
  CoefMatrix = matrix(0, p, p - 1)
  type = match.arg(type)
  if (is.null(lambda)){
    if (type == "slasso"){
      lambda = sqrt(2 * 2.01 * log(p * (log(p))^(1.5) / sqrt(n)) / n)
    } else if (type == "lasso"){
      lambda = sqrt(2 * log(p) / n)
    }
  }
  if (type == "slasso"){
    for (i in 1 : p){
      slasso = scaledlasso(X = XS[, -i], y = X[, i], lam0 = lambda)
      Eresidual[, i] = slasso$residuals
      CoefMatrix[i, ] = slasso$coefficients / sdX[-i]
    }
  } else if (type == "lasso"){
    for (i in 1 : p){
      lasso = glmnet(x = XS[,-i], y = X[,i], intercept = FALSE, standardize = FALSE)
      Coef = coef.glmnet(lasso, s = lambda * sdX[i])
      CoefMatrix[i, ] = as.vector(Coef)[-1] / sdX[-i]
      Predict = predict.glmnet(lasso, s = lambda* sdX[i], newx = XS[,-i])
      Eresidual[, i] = X[,i] - Predict[,1]
    }
  }
  CovRes = t(Eresidual) %*% Eresidual / n
  m = 1
  Est = matrix(1, p, p)
  BTAll = matrix(0, n, Mp)
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      temp = CovRes[i, j] + diag(CovRes)[i] * CoefMatrix[j, i] + diag(CovRes)[j] * CoefMatrix[i, j - 1]
      Est[j, i] = Est[i, j] = pmin(pmax(-1, temp / sqrt(diag(CovRes)[i] * diag(CovRes)[j])), 1)
      omegaHat = - temp / (diag(CovRes)[i] * diag(CovRes)[j])
      BTAll[, m] = ( Eresidual[, i] * Eresidual[, j] + temp ) / sqrt(diag(CovRes)[i] * diag(CovRes)[j]) - omegaHat * sqrt(diag(CovRes)[j]) * ( Eresidual[, i]^2 - CovRes[i, i] ) / (2 * sqrt(diag(CovRes)[i]))  - omegaHat * sqrt(diag(CovRes)[i]) * ( Eresidual[, j]^2 - CovRes[j, j] ) / (2 * sqrt(diag(CovRes)[j]))
      m = m + 1
    }
  }
  BTAllcenter = scale(BTAll, scale = FALSE)
  if (!ci) return(structure(list(coef = Est, asym.ex = BTAllcenter, type = type), class='indEst'))
  NumAll = c()
  DenAll = c()
  for(i in 1 : Mp){
    AR1 = ar(BTAllcenter[, i], aic = FALSE, order.max = 1)
    rhoEst = AR1$ar
    sigma2Est = AR1$var.pred
    NumAll[i] = 4 * (rhoEst * sigma2Est)^2 / (1 - rhoEst)^8
    DenAll[i] = sigma2Est^2 / (1 - rhoEst)^4
  }
  a2All = sum(NumAll) / sum(DenAll)
  bandwidthAll = 1.3221 * (a2All * n)^(0.2)
  diagW1 = colSums(BTAll^2) / n
  for (h in 1 : (n - 1)){
    gammah = colSums(matrix(BTAll[(1 + h):n,] * BTAll[1:(n - h),], ncol=Mp))
    diagW1 = diagW1 + 2 * QS(h / bandwidthAll) * gammah / n
  }
  m = 1
  ci.upper = ci.lower = diag(rep(1, p))
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      ci.upper[j, i] = ci.upper[i, j] = min(1, Est[i, j] + qnorm(1 - alpha / 2) * sqrt(diagW1[m] / n))
      ci.lower[j, i] = ci.lower[i, j] = max(-1, Est[i, j] - qnorm(1 - alpha / 2) * sqrt(diagW1[m] / n))
      m = m + 1
    }
  }
  return(structure(list(coef = Est, ci.lower = ci.lower, ci.upper = ci.upper, asym.ex = BTAllcenter, type = type), class = 'indEst'))
}
