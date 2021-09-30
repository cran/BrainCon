#' Estimate individual-level partial correlation coefficients
#'
#' Estimate individual-level partial correlation coefficients in time series data
#' with \eqn{1-\alpha} confidence interval.
#' It's not a joint confidence interval for multiple tests.
#' \cr
#' \cr
#'
#'@param X time series data of an individual which is a \eqn{n*p} numeric matrix.
#'@param alpha significance level, default value is \code{0.05}.
#'@param lambda a penalty parameter used in lasso of order \code{sqrt(log(p)/n)}, if \code{NULL}, \code{2*sqrt(log(p)/n)} will be used.
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
#' \code{asym.ex} a matrix measuring the asymptotical expansion of estimates, which will be used for multiple tests.
#'
#'@examples
#' ## Quick example for the individual-level estimates
#' data(indsim)
#' pc = individual.est(indsim)       # estimating partial correlation coefficients
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

individual.est <- function(X, alpha = 0.05, lambda = NULL, ci = TRUE){
  n = dim(X)[1]
  p = dim(X)[2]
  Mp = p * (p - 1) / 2
  X = scale(X, scale = FALSE)
  XS = scale(X, center = FALSE)
  sdX = apply(X, 2, sd)
  Eresidual = matrix(0, n, p)
  CoefMatrix = matrix(0, p, p - 1)
  if (is.null(lambda)) lambda = 2 * sqrt(log(p) / n)
  for (i in 1 : p){
    lasso = glmnet(x = XS[,-i], y = X[,i], intercept = FALSE, standardize = FALSE)
    Coef = coef.glmnet(lasso, s = lambda)
    CoefMatrix[i, Coef@i] = Coef@x / sdX[-i][Coef@i]
    Predict = predict.glmnet(lasso, s = lambda, newx = XS[,-i])
    Eresidual[, i] = X[,i] - Predict[,1]
  }
  CovRes = t(Eresidual) %*% Eresidual / n
  m = 1
  Est = matrix(1, p, p)
  BTAll = matrix(0, n, Mp)
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      temp = CovRes[i, j] + diag(CovRes)[i] * CoefMatrix[j, i] + diag(CovRes)[j] * CoefMatrix[i, j - 1]
      Est[j, i] = Est[i, j] = temp / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
      omegaHat = - temp / (diag(CovRes)[i] * diag(CovRes)[j])
      BTAll[, m] = ( Eresidual[, i] * Eresidual[, j] + temp ) / sqrt(diag(CovRes)[i] * diag(CovRes)[j]) - omegaHat * sqrt(diag(CovRes)[j]) * ( Eresidual[, i]^2 - CovRes[i, i] ) / (2 * sqrt(diag(CovRes)[i]))  - omegaHat * sqrt(diag(CovRes)[i]) * ( Eresidual[, j]^2 - CovRes[j, j] ) / (2 * sqrt(diag(CovRes)[j]))
      m = m + 1
    }
  }
  BTAllcenter = scale(BTAll, scale = FALSE)
  if (!ci) return(structure(list(coef = Est, asym.ex=BTAllcenter),class='indEst'))
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
  return(structure(list(coef = Est, ci.lower = ci.lower, ci.upper = ci.upper, asym.ex=BTAllcenter), class='indEst'))
}

