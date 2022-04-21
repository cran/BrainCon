#' The one-sample population inference
#'
#' Identify the nonzero partial correlations in one-sample population,
#' based on controlling the rate of the false discovery proportion (FDP) exceeding \eqn{c0}
#' at \eqn{\alpha}, considering time dependence.
#' Input a \code{popEst} class object returned by \code{\link{population.est}}.
#' \cr
#' \cr
#'
#'@param popEst A \code{popEst} class object.
#'@param alpha significance level, default value is \code{0.05}.
#'@param c0  threshold of the exceedance rate of FDP,
#'default value is \code{0.1}. A smaller value of \code{c0} will
#'reduce false positives, but it may also cost more false negatives.
#'@param targetSet a two-column matrix. Each row contains two index corresponding to a pair of variables of interest.
#'If \code{NULL}, any pair of two variables is considered to be of interest.
#'@param MBT times of multiplier bootstrap, default value is \code{5000}.
#'@param simplify a logical indicating whether results should be simplified if possible.
#'
#'@return If \code{simplify} is \code{FALSE}, a \eqn{p*p} matrix with values 0 or 1 is returned, and 1 means nonzero.
#'
#'And if \code{simplify} is \code{TRUE}, a two-column matrix is returned,
#'indicating the row index and the column index of recovered nonzero partial correlations.
#'We only retain the results which the row index is less than the column index.
#'Those with larger test statistics are sorted first.
#'
#'@seealso \code{\link{individual.test}}.
#'
#'@examples
#' ## Quick example for the one-sample population inference
#' data(popsimA)
#' # estimating partial correlation coefficients by scaled lasso
#' pc = population.est(popsimA)
#' # conducting hypothesis test
#' Res = population.test(pc)
#' # conducting hypothesis test in variables of interest
#' set = cbind(rep(7:9, each = 10), 1:10)
#' Res_like = population.test(pc, targetSet = set)
#'
#' @references
#' Qiu Y. and Zhou X. (2021).
#' Inference on multi-level partial correlations
#' based on multi-subject time series data,
#' \emph{Journal of the American Statistical Association}, 00, 1-15.

population.test <- function(popEst, alpha = 0.05, c0 = 0.1, targetSet = NULL, MBT = 5000, simplify = !is.null(targetSet)){
  force(simplify)
  if (class(popEst) != 'popEst')
    stop("The argument popEst requires a 'popEst' class input!\n")
  EstAll = popEst$coef
  p = nrow(EstAll)
  MC = length(popEst[['ind.est']])
  if (is.null(targetSet)){
    targetSet = upper.tri(EstAll)
    Mp = p * (p - 1) / 2
  } else {
    simplify = TRUE
    targetSet = normalize.set(targetSet, p)
    Mp = nrow(targetSet)
  }
  EstVec = matrix(0, MC, Mp)
  for (i in 1 : MC){
    Est = popEst[['ind.est']][[i]][['coef']]
    EstVec[i,] = Est[targetSet]
  }
  EstVecCenter = scale(EstVec, scale = FALSE)
  BTAllsim = matrix(0, Mp, MBT)
  for (i in 1 : MBT){
    temp = rnorm(MC)
    BTAllsim[, i] = (MC)^(-0.5) * colSums(temp * EstVecCenter)
  }
  SignalID = c()
  TestPro = EstAll[targetSet]
  BTPro = abs(BTAllsim)
  repeat{
    PCmaxIndex = which.max(abs(TestPro))
    SignalIDtemp = which(EstAll == TestPro[PCmaxIndex], arr.ind = T)
    SignalID = rbind(SignalID, SignalIDtemp)
    TestPro = TestPro[-PCmaxIndex]
    BTPro = BTPro[-PCmaxIndex, ]
    TestStatPro = sqrt(MC) * max(abs(TestPro))
    BTAllsimPro = apply(BTPro, 2, max)
    QPro = sort(BTAllsimPro)[(1 - alpha) * MBT]
    if (TestStatPro < QPro)	break
  }
  aug = floor(c0 * dim(SignalID)[1] / (2 * (1 - c0)))
  if (aug > 0){
    PCmaxIndex = order(-abs(TestPro))[1 : aug]
    for (q in 1 : length(PCmaxIndex)){
      SignalIDtemp = which(EstAll == TestPro[PCmaxIndex[q]], arr.ind = TRUE)
      SignalID = rbind(SignalID, SignalIDtemp)
    }
  }
  if (simplify) return(subset(SignalID, SignalID[,1] < SignalID[,2]))
  recovery = diag(rep(1, p))
  recovery[SignalID[,1] + (SignalID[,2] - 1) * p] = 1
  return(recovery)
}
