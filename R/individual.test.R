#' Identify nonzero individual-level partial correlations
#'
#' Identify nonzero individual-level partial correlations in time series data
#' by controlling the rate of the false discovery proportion (FDP) exceeding \eqn{c0}
#' at \eqn{\alpha}, considering time dependence.
#' Input an \code{indEst} class object returned by \code{\link{individual.est}} or \code{\link{population.est}}.
#' \cr
#' \cr
#'
#'@param indEst An \code{indEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of FDP,
#'default value is \code{0.1}.
#'The choice of \code{c0} depends on the empirical problem. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'@param targetSet a two-column matrix. Each row contains two index corresponding to a pair of variables of interest.
#'If \code{NULL}, any pair of two variables is considered to be of interest.
#'@param MBT times of multiplier bootstrap, default value is \code{3000}.
#'@param simplify a logical indicating whether results should be simplified if possible.
#'
#'@return If \code{simplify} is \code{FALSE}, a \eqn{p*p} matrix with values 0 or 1 is returned.
#'If the j-th row and k-th column of the matrix is 1,
#'then the partial correlation coefficient between
#'the j-th variable and the k-th variable is identified to be nonzero.
#'
#'And if \code{simplify} is \code{TRUE}, a two-column matrix is returned,
#'indicating the row index and the column index of recovered nonzero partial correlations.
#'We only retain the results which the row index is less than the column index.
#'Those with larger test statistics are sorted first.
#'
#'@seealso \code{\link{population.est}} for making inferences on one individual in the population.
#'
#'@examples
#' ## Quick example for the individual-level inference
#' data(indsim)
#' # estimating partial correlation coefficients by scaled lasso
#' pc = individual.est(indsim)
#' # conducting hypothesis test
#' Res = individual.test(pc)
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

individual.test <- function(indEst, alpha = 0.05, c0 = 0.1, targetSet = NULL, MBT = 3000, simplify = !is.null(targetSet)){
  force(simplify)
  if (!inherits(indEst, 'indEst'))
    stop("The argument indEst requires an 'indEst' class input!\n")
  Est=indEst$coef
  BTAllcenter = indEst$asym.ex
  n = nrow(BTAllcenter)
  p = nrow(Est)
  Mp = p * (p - 1) / 2
  if (is.null(targetSet)){
    targetSet = lower.tri(Est)
    index = 1 : Mp
  } else {
    simplify = TRUE
    targetSet = normalize.set(targetSet, p)
    index = (2 * p - targetSet[, 1]) * (targetSet[, 1] - 1) / 2 + targetSet[, 2] - targetSet[, 1]
  }
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
  BTcovAll = matrix(0, n, n)
  for (i in 1 : n){
    for (j in 1 : n){
      BTcovAll[i, j] = QS(abs(i - j) / bandwidthAll)
    }
  }
  BTAllsim = matrix(0, Mp, MBT)
  for (i in 1 : MBT){
    temp = mvrnorm(1, rep(0, n), BTcovAll)
    BTAllsim[, i] = (n)^(-0.5) * colSums(temp * BTAllcenter)
  }
  WdiagAllEmp = colSums(BTAllcenter ^ 2) / n
  TestAllstandard = WdiagAllEmp^(-1/2) * Est[lower.tri(Est)]
  BTAllsim0 = WdiagAllEmp^(-1/2) * BTAllsim
  SignalID=c()
  TestPro = Est[targetSet]
  TestProstandard = TestAllstandard[index]
  BTPro = abs(BTAllsim0)[index, ]
  repeat{
    PCmaxIndex = which.max(abs(TestProstandard))
    SignalIDtemp = which(Est == TestPro[PCmaxIndex], arr.ind = T)
    SignalID = rbind(SignalID, SignalIDtemp)
    TestPro = TestPro[-PCmaxIndex]
    BTPro = BTPro[-PCmaxIndex, ]
    TestProstandard = TestProstandard[-PCmaxIndex]
    TestStatPro = sqrt(n) * max(abs(TestProstandard))
    BTAllsimPro = apply(BTPro, 2, max)
    QPro = sort(BTAllsimPro)[(1 - alpha) * MBT]
    if (TestStatPro < QPro)	break
  }
  aug = floor(c0 * dim(SignalID)[1] / (2 * (1 - c0)))
  if (aug > 0){
    PCmaxIndex = order(-abs(TestProstandard))[1 : aug]
    for (q in 1 : length(PCmaxIndex)){
      SignalIDtemp = which(Est == TestPro[PCmaxIndex[q]], arr.ind = TRUE)
      SignalID = rbind(SignalID, SignalIDtemp)
    }
  }
  if (simplify) return(subset(SignalID, SignalID[,1] < SignalID[,2]))
  recovery = diag(rep(1, p))
  recovery[SignalID[,1]+(SignalID[,2]-1)*p]=1
  return(recovery)
}
