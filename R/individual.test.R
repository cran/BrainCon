#' Identify nonzero individual-level partial correlations
#'
#' Identify nonzero individual-level partial correlations in time series data
#' by controlling the exceedance rate of the false discovery proportion (FDP)
#' at \eqn{\alpha=0.05}, considering time dependence.
#' Input data \eqn{X} contains values of p interested
#' variables over n periods.
#' \cr
#' \cr
#'
#'@param indEst An \code{indEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of the false discovery proportion (FDP),
#'default value is \code{0.1}.
#'The choice of \code{c0} depends on the empirical problem. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'@param MBT times of multiplier bootstrap, default value is \code{3000}.
#'
#'@return A \eqn{p*p} matrix with values 0 or 1.
#'If the j-th row and k-th column of the matrix is 1,
#'then the partial correlation coefficient between
#'the j-th variable and the k-th variable is identified to be nonzero.
#'
#'@seealso \code{\link{population.test}} for making inferences on one individual in the population.
#'
#'@examples
#' ## Quick example for the individual-level inference
#' data(indsim)
#' pc = individual.est(indsim)       # estimating partial correlation coefficients
#' Res = individual.test(pc)         # conducting hypothesis test
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

individual.test <- function(indEst, alpha = 0.05, c0 = 0.1, MBT = 3000){
  if (class(indEst) != 'indEst')
    stop("The argument indEst requires an 'indEst' class input!\n")
  Est=indEst$coef
  BTAllcenter = indEst$asym.ex
  TestAllmean = Est[lower.tri(Est)]
  n = nrow(BTAllcenter)
  p = nrow(Est)
  Mp = p * (p - 1) / 2
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
  TestAllstandard = WdiagAllEmp^(-1/2) * TestAllmean
  BTAllsim0 = WdiagAllEmp^(-1/2) * BTAllsim
  SignalID=c()
  TestPro = TestAllmean; TestProstandard = TestAllstandard; BTPro = BTAllsim0
  repeat{
    PCtemp = TestPro
    TPStemp = round(TestProstandard, 5)
    PCmaxIndex0 = order(-abs(TPStemp))[1]
    PCmaxIndex = which(TPStemp %in% c(TPStemp[PCmaxIndex0], -TPStemp[PCmaxIndex0]))
    SignalIDtemp = c()
    for (q in 1 : length(PCmaxIndex)){
      SignalIDtemp1 = which(abs(Est - PCtemp[PCmaxIndex[q]]) == min(abs(Est - PCtemp[PCmaxIndex[q]])), arr.ind = TRUE)
      SignalIDtemp = rbind(SignalIDtemp, SignalIDtemp1)
    }
    SignalID = rbind(SignalID, SignalIDtemp)
    TestPro = TestPro[-PCmaxIndex]
    BTPro = BTPro[-PCmaxIndex, ]
    TestProstandard = TestProstandard[-PCmaxIndex]
    TestStatPro = sqrt(n) * max(abs(TestProstandard))
    BTAllsimPro = c()
    for (i in 1 : MBT){
      BTAllsimPro[i] = max(abs(BTPro[, i]))
    }
    QPro = sort(BTAllsimPro)[(1 - alpha) * MBT]
    if (TestStatPro < QPro)	break
  }
  aug = ceiling(c0 * dim(SignalID)[1] / (2 * (1 - c0)))
  PCtemp = TestPro
  TPStemp = round(TestProstandard, 5)
  PCmaxIndex0 = order(-abs(TPStemp))[1 : aug]
  PCmaxIndex = which(TPStemp %in% c(TPStemp[PCmaxIndex0], -TPStemp[PCmaxIndex0]))
  SignalIDtemp = c()
  for (q in 1 : length(PCmaxIndex)){
    SignalIDtemp1 = which(abs(Est - PCtemp[PCmaxIndex[q]]) == min(abs(Est - PCtemp[PCmaxIndex[q]])), arr.ind = TRUE)
    SignalIDtemp = rbind(SignalIDtemp, SignalIDtemp1)
  }
  SignalID = rbind(SignalID, SignalIDtemp)
  recovery = diag(rep(1, p))
  recovery[SignalID[,1]+(SignalID[,2]-1)*p]=1
  return(recovery)
}
