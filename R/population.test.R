#' The one-sample population inference
#'
#' Identify the nonzero partial correlations in one-sample population,
#' based on false discovery proportion controlling.
#' at \eqn{\alpha=0.05}, considering time dependence.
#' Input  data \eqn{Z} , contains values of p interested variables.
#' \cr
#' \cr
#'
#'@param popEst A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of the false discovery proportion (FDP),
#'default value is \code{0.1}.
#'The choice of \code{c0} depends on the empirical problem. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'@param MBT times of multiplier bootstrap, default value is \code{3000}.
#'
#'@return A \eqn{p*p} matrix with values 0 or 1.
#'
#'@seealso \code{\link{individual.test}}.
#'
#'@examples
#' ## Quick example for the one-sample population inference
#' data(popsimA)
#' pc = population.est(popsimA)     # estimating partial correlation coefficients
#' Res = population.test(pc)        # conducting hypothesis test
#'
#' ## Inference on the first subject in population
#' Res1 = individual.test(pc$ind.est[[1]])
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

population.test <- function(popEst, alpha = 0.05, c0 = 0.1, MBT = 3000){
  if (class(popEst) != 'popEst')
    stop("The argument popEst requires a 'popEst' class input!\n")
  EstAll = popEst$coef
  n = nrow(popEst[['ind.est']][[1]][['asym.ex']])
  p = nrow(EstAll)
  MC = length(popEst[['ind.est']])
  Mp = p * (p - 1) / 2
  EstVec = matrix(0, MC, Mp)
  for (i in 1 : MC){
    Est = popEst[['ind.est']][[i]][['coef']]
    EstVec[i,] = Est[lower.tri(Est)]
  }
  EstVecCenter = scale(EstVec, scale = FALSE)
  BTAllsim = matrix(0, Mp, MBT)
  for (i in 1 : MBT){
    temp = rnorm(MC)
    BTAllsim[, i] = (MC)^(-0.5) * colSums(temp * EstVecCenter)
  }
  SignalID = c()
  TestPro = TestProstandard = EstAll[lower.tri(EstAll)]
  BTPro = BTAllsim
  repeat{
    PCtemp = TestPro
    TPStemp = round(TestProstandard, 5)
    PCmaxIndex0 = order(-abs(TPStemp))[1]
    PCmaxIndex = which(TPStemp %in% c(TPStemp[PCmaxIndex0], -TPStemp[PCmaxIndex0]))
    SignalIDtemp = c()
    for (q in 1 : length(PCmaxIndex)){
      SignalIDtemp1 = which(abs(EstAll - PCtemp[PCmaxIndex[q]]) == min(abs(EstAll - PCtemp[PCmaxIndex[q]])), arr.ind = TRUE)
      SignalIDtemp = rbind(SignalIDtemp, SignalIDtemp1)
    }
    SignalID = rbind(SignalID, SignalIDtemp)
    TestPro = TestPro[-PCmaxIndex]
    BTPro = BTPro[-PCmaxIndex, ]
    TestProstandard = TestProstandard[-PCmaxIndex]
    TestStatPro = sqrt(MC) * max(abs(TestProstandard))
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
    SignalIDtemp1 = which(abs(EstAll - PCtemp[PCmaxIndex[q]]) == min(abs(EstAll - PCtemp[PCmaxIndex[q]])), arr.ind = TRUE)
    SignalIDtemp = rbind(SignalIDtemp, SignalIDtemp1)
  }
  SignalID = rbind(SignalID, SignalIDtemp)
  recovery = diag(rep(1, p))
  recovery[SignalID[,1]+(SignalID[,2]-1)*p]=1
  return(recovery)
}
