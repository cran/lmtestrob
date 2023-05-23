\name{robfmtest}
\alias{robfmtest}
\title{
Robust Functional Specification Test
}
\description{
Outlier robust test for functional miss-specification. It can be used to test linearity. The test is based on the robust Wald-type test by Heritier and Ronchetti (1994). 
}
\usage{
robfmtest(formula, power = 2:3, type = c("regressor"), data,
          x.weights = c("HAT", "MCD"), testtype = "Wald", ...)
}
\arguments{
  \item{formula}{
a symbolic description of the model to be tested.
}
  \item{power}{
integer(s). A vector of positive integers specifying the powers of the variables that should be tested. The default option tests second and third powers.
}
  \item{type}{
currently, only powers of regressors can be used.
}
  \item{data}{
an optional data frame containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which {robfmtest} is called.
}
  \item{x.weights}{
a string, indicating how the robustness weights on the covariates should be computed. The default option uses hat-matrix-based weights, second option allows to use robust Mahalanobis distance-based weights, where the Minimum Covariance Determinant is used to estimate location and scatter.
}
  \item{testtype}{
currently, the robust version of Wald test is implemented.
}
  \item{\dots}{
currently not used.
}
}
\details{
Since the classical tests including {resettest}, {raintest} and {harvtest} implemented in {lmtest} are not resistant to outliers and can become misleading even in the presence of one outlier, we provide a test which is resistant to outliers. The price to pay for robustness is a small loss of power, when the model holds exactly.
}
\value{
A list with class robfmtest containing the following components:
\item{statistic}{the value of the test statistic.}
\item{dof}{the number of degrees of freedom.}
\item{method}{a character string indicating what type of test was performed.}
\item{p.value}{the p-value of the test.}
\item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Heritier, S., and Ronchetti, E. (1994) Robust Bounded-Influence Tests in General Parametric Models. \emph{Journal of the American Statistical Association}, 89, p. 897-904.
}
\author{
Mikhail Zhelonkin
}
\examples{
set.seed(123)
n <- 50
x = runif(n, -3, 3)
y = rnorm(n)
example.dat <- data.frame(x, y)
robfmtest(y ~ x, data = example.dat)
library(lmtest)
resettest(y ~ x, data = example.dat, type = "fitted")
x[50] <- -3
y[50] <- -10
example.dat <- data.frame(x, y)
robfmtest(y ~ x, data = example.dat)
resettest(y ~ x, data = example.dat, type = "fitted")
}
\keyword{ Robust test }