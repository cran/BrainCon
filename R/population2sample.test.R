#' Identify differences of partial correlations between two populations
#'
#' Identify differences of partial correlations between two populations
#' in two groups of time series data by
#' controlling the rate of the false discovery proportion (FDP) exceeding \eqn{c0}
#' at \eqn{\alpha}, considering time dependence.
#' Input two \code{popEst} class objects returned by \code{\link{population.est}}
#' (the number of individuals in two groups can be different).
#' \cr
#' \cr
#'
#'@param popEst1 A \code{popEst} class object.
#'@param popEst2 A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of FDP,
#'default value is \code{0.1}. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'@param targetSet a two-column matrix. Each row contains two index corresponding to a pair of variables of interest.
#'If \code{NULL}, any pair of two variables is considered to be of interest.
#'@param MBT times of multiplier bootstrap, default value is \code{5000}.
#'@param simplify a logical indicating whether results should be simplified if possible.
#'
#'@return If \code{simplify} is \code{FALSE}, a \eqn{p*p} matrix with values 0 or 1 is returned.
#'If the j-th row and k-th column of the matrix is 1,
#'then the partial correlation coefficients between
#'the j-th variable and the k-th variable in two populations
#'are identified to be unequal.
#'
#'And if \code{simplify} is \code{TRUE}, a two-column matrix is returned,
#'indicating the row index and the column index of recovered unequal partial correlations.
#'We only retain the results which the row index is less than the column index.
#'Those with larger test statistics are sorted first.
#'
#'@examples
#' ## Quick example for the two-sample case inference
#' data(popsimA)
#' data(popsimB)
#' # estimating partial correlation coefficients by lasso (scaled lasso does the same)
#' pc1 = population.est(popsimA, type = 'l')
#' pc2 = population.est(popsimB, type = 'l')
#' # conducting hypothesis test
#' Res = population2sample.test(pc1, pc2)
#' # conducting hypothesis test and returning simplified results
#' Res_s = population2sample.test(pc1, pc2, simplify = TRUE)
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

population2sample.test <- function(popEst1, popEst2, alpha = 0.05, c0 = 0.1, targetSet = NULL, MBT = 5000, simplify = !is.null(targetSet)){
  force(simplify)
  if (class(popEst1) != 'popEst' | class(popEst2) != 'popEst')
    stop("The arguments popEst1 and popEst2 require 'popEst' class inputs!\n")
  EstAll1 = popEst1$coef
  EstAll2 = popEst2$coef
  p = nrow(EstAll1)
  MC1 = length(popEst1[['ind.est']])
  MC2 = length(popEst2[['ind.est']])
  if (is.null(targetSet)){
    targetSet = upper.tri(EstAll1)
    Mp = p * (p - 1) / 2
  } else {
    simplify = TRUE
    targetSet = normalize.set(targetSet, p)
    Mp = nrow(targetSet)
  }
  EstVec1 = matrix(0, MC1, Mp)
  EstVec2 = matrix(0, MC2, Mp)
  for (i in 1 : MC1){
    Est = popEst1[['ind.est']][[i]][['coef']]
    EstVec1[i,] = Est[targetSet]
  }
  for (i in 1 : MC2){
    Est = popEst2[['ind.est']][[i]][['coef']]
    EstVec2[i,] = Est[targetSet]
  }
  EstVecCenter1 = scale(EstVec1, scale = FALSE)
  EstVecCenter2 = scale(EstVec2, scale = FALSE)
  TestAllstandard1 = EstAll1[targetSet]
  TestAllstandard2 = EstAll2[targetSet]
  EstAll = EstAll1 - EstAll2
  BTAllsim = matrix(0, Mp, MBT)
  for (i in 1 : MBT){
    temp1 = rnorm(MC1)
    temp2 = rnorm(MC2)
    BTAllsim[, i] = (MC1)^(-0.5) * colSums(temp1 * EstVecCenter1) - (MC1)^(0.5) * colMeans(temp2 * EstVecCenter2)
  }
  SignalID=c()
  TestPro = TestAllstandard1 - TestAllstandard2
  BTPro = abs(BTAllsim)
  repeat{
    PCmaxIndex = which.max(abs(TestPro))
    SignalIDtemp = which(EstAll == TestPro[PCmaxIndex], arr.ind = T)
    SignalID = rbind(SignalID, SignalIDtemp)
    TestPro = TestPro[-PCmaxIndex]
    BTPro = BTPro[-PCmaxIndex, ]
    TestStatPro = sqrt(MC1) * max(abs(TestPro))
    BTAllsimPro = apply(BTPro, 2, max)
    QPro = sort(BTAllsimPro)[(1 - alpha) * MBT]
    if (TestStatPro < QPro)	break
  }
  aug = round(c0 * dim(SignalID)[1] / (2 * (1 - c0)) + 1e-3)
  if (aug > 0){
    PCmaxIndex = order(-abs(TestPro))[1 : aug]
    for (q in 1 : length(PCmaxIndex)){
      SignalIDtemp = which(EstAll == TestPro[PCmaxIndex[q]], arr.ind = TRUE)
      SignalID = rbind(SignalID, SignalIDtemp)
    }
  }
  if (simplify) return(subset(SignalID, SignalID[,1] < SignalID[,2]))
  recovery = matrix(0, p, p)
  recovery[SignalID[,1]+(SignalID[,2]-1)*p]=1
  return(recovery)
}
