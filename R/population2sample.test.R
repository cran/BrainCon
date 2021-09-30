#' Identify differences of partial correlations between two populations
#'
#' Identify differences of partial correlations between two populations
#' in two groups of time series data
#' by controlling the exceedance rate of the false discovery proportion (FDP)
#' at \eqn{\alpha=0.05}, considering time dependence.
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
#'@param MBT times of multiplier bootstrap, default value is \code{3000}.
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
#' Res = population2sample.test(pc1, pc2)          # conducting hypothesis test
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15

population2sample.test <- function(popEst1, popEst2, alpha = 0.05, c0 = 0.1, MBT = 3000){
  if (class(popEst1) != 'popEst' | class(popEst2) != 'popEst')
    stop("The arguments popEst1 and popEst2 require 'popEst' class inputs!\n")
  EstAll1 = popEst1$coef
  EstAll2 = popEst2$coef
  n = nrow(popEst1[['ind.est']][[1]][['asym.ex']])
  p = nrow(EstAll1)
  MC1 = length(popEst1[['ind.est']])
  MC2 = length(popEst2[['ind.est']])
  Mp = p * (p - 1) / 2
  EstVec1 = matrix(0, MC1, Mp)
  EstVec2 = matrix(0, MC2, Mp)
  for (i in 1 : MC1){
    Est = popEst1[['ind.est']][[i]][['coef']]
    EstVec1[i,] = Est[lower.tri(Est)]
  }
  for (i in 1 : MC2){
    Est = popEst2[['ind.est']][[i]][['coef']]
    EstVec2[i,] = Est[lower.tri(Est)]
  }
  EstVecCenter1 = scale(EstVec1, scale = FALSE)
  EstVecCenter2 = scale(EstVec2, scale = FALSE)
  TestAllstandard1 = EstAll1[lower.tri(EstAll1)]
  TestAllstandard2 = EstAll2[lower.tri(EstAll2)]
  EstAll = EstAll1 - EstAll2
  BTAllsim = matrix(0, Mp, MBT)
  for (i in 1 : MBT){
    temp1 = rnorm(MC1)
    temp2 = rnorm(MC2)
    BTAllsim[, i] = (MC1)^(-0.5) * colSums(temp1 * EstVecCenter1) - (MC1)^(0.5) * colMeans(temp2 * EstVecCenter2)
  }
  SignalID=c()
  TestPro = TestProstandard = TestAllstandard1 - TestAllstandard2
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
    TestStatPro = sqrt(MC1) * max(abs(TestProstandard))
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
  recovery = matrix(0, p, p)
  recovery[SignalID[,1]+(SignalID[,2]-1)*p]=1
  return(recovery)
}
