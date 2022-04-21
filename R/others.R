#' tool functions
#'
#' @param u numeric value.
#' @param set two-column numeric matrix.
#' @param p the number of variables.
#' @param X the input matrix of scaled lasso.
#' @param y response variable of scaled lasso.
#' @param lam0 numeric value, the penalty parameter of scaled lasso.
#'
#' @return Intermediate results.
#'
#' @name tool
#' @keywords internal
NULL

#' @rdname tool
QS = function(u){
  if (u == 0) ker = 1
  else ker = 25 * ( sin(6 * pi * u / 5) / (6 * pi * u / 5) - cos(6 * pi * u / 5) ) / (12 * pi^2 * u^2)
  return(ker)
}

#' @rdname tool
normalize.set <- function(set, p){
  set = as.matrix(set)
  if (ncol(set) != 2)
    stop('The argument targetSet requires a two-column matrix!\n')
  colnames(set) = c('row', 'col')
  set = rbind(set, set[,2:1])
  set = set[set[,1] < set[,2],]
  set = set[set[,2] <= p,]
  return(set[!duplicated(set),])
}

#' @rdname tool
scaledlasso <- function (X, y, lam0 = NULL){
  objlasso = glmnet(x = X, y = y, intercept = FALSE, standardize = FALSE)
  sigmaint = 0.1
  sigmanew = 5
  flag = 0
  while (abs(sigmaint - sigmanew) > 1e-04 & flag <= 100) {
    flag = flag + 1
    sigmaint = sigmanew
    lam = lam0 * sigmaint
    hy = predict.glmnet(objlasso, s = lam, newx = X)[, 1]
    sigmanew = sqrt(mean((y - hy)^2))
  }
  hbeta = as.vector(coef.glmnet(objlasso, s = lam))[-1]
  return(list(hsigma = sigmanew, coefficients = hbeta,
              fitted.values = hy, residuals = y - hy))
}
