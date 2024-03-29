% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/individual.est.R
\name{individual.est}
\alias{individual.est}
\title{Estimate individual-level partial correlation coefficients}
\usage{
individual.est(
  X,
  lambda = NULL,
  type = c("slasso", "lasso"),
  alpha = 0.05,
  ci = TRUE
)
}
\arguments{
\item{X}{time series data of an individual which is a \eqn{n*p} numeric matrix, where \eqn{n} is the number of periods of time and \eqn{p} is the number of variables.}

\item{lambda}{a penalty parameter of order \eqn{\sqrt{\log(p)/n}}.
If \code{NULL}, \eqn{\sqrt{2*2.01/n*\log(p*(\log(p))^{1.5}/n^{0.5})}} is used in scaled lasso, and \eqn{\sqrt{2*\log(p)/n}} is used in lasso.
Increasing the penalty parameter may lead to larger residuals in the node-wise regression,
causing larger absolute values of estimates of partial correlation coefficients, which may cause more false positives in subsequent tests.}

\item{type}{a character string representing the method of estimation. \code{"slasso"} means scaled lasso, and \code{"lasso"} means lasso. Default value is \code{"slasso"}.}

\item{alpha}{significance level, default value is \code{0.05}.}

\item{ci}{a logical indicating whether to compute \eqn{1-\alpha} confidence interval, default value is \code{TRUE}.}
}
\value{
An \code{indEst} class object containing two or four components.

\code{coef} a \eqn{p*p} partial correlation coefficients matrix.

\code{ci.lower} a \eqn{p*p} numeric matrix containing the lower bound of \eqn{1-\alpha} confidence interval,
returned if \code{ci} is \code{TRUE}.

\code{ci.upper} a \eqn{p*p} numeric matrix containing the upper bound of \eqn{1-\alpha} confidence interval,
returned if \code{ci} is \code{TRUE}.

\code{asym.ex} a matrix measuring the asymptotic expansion of estimates, which will be used for multiple tests.

\code{type} regression type in estimation.
}
\description{
Estimate individual-level partial correlation coefficients in time series data
with \eqn{1-\alpha} confidence intervals.
Note that these are confidence intervals for single parameters, not simultaneous confidence intervals.
\cr
\cr
}
\examples{
## Quick example for the individual-level estimates
data(indsim)
# estimating partial correlation coefficients by scaled lasso
pc = individual.est(indsim)

}
\references{
Qiu Y. and Zhou X. (2021).
Inference on multi-level partial correlations
based on multi-subject time series data,
\emph{Journal of the American Statistical Association}, 00, 1-15.

Sun T. and Zhang C. (2012).
Scaled Sparse Linear Regression,
\emph{Biometrika}, 99, 879–898.

Liu W. (2013).
Gaussian Graphical Model Estimation With False Discovery Rate Control,
\emph{The Annals of Statistics}, 41, 2948–2978.

Ren Z., Sun T., Zhang C. and Zhou H. (2015).
Asymptotic Normality and Optimalities in Estimation of Large Gaussian Graphical Models,
\emph{The Annals of Statistics}, 43, 991–1026.
}
\seealso{
\code{\link{population.est}}.
}
