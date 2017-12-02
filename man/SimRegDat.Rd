\name{SimRegDat}
\alias{SimRegDat}
\title{ Simulate Incomplete Data for High-Dimensional Linear Regression.}
\description{Simulate incomplete data for high-dimensional linear regression with dependent or independent covariates\code{RegICC(x,y...)}.}

\usage{
SimRegDat(n = 100, p = 200, type = "indep", rate = 0.1)
}

\arguments{
  \item{ n }{ Number of observations, default of 100.}
  \item{ p }{ Number of covariates, default of 200.}
  \item{ type }{ When \code{type=="indep"}, it simulates the data with independent covariates, or \code{type=="dep"}, it simulates the data with dependent covariates, the default type is "indep".}
  \item{ rate }{ Missing rate, the default value is 0.1.}
}

\value{
   \item{ x }{ \emph{\eqn{n}}x\emph{p} covariates matrix.}
  \item{ y }{  \emph{\eqn{n}}x\emph{1} responses.}
  \item{ coef }{ \emph{\eqn{p}}x\emph{1} coefficients for generating responses from the covariates matrix.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
library(ICmiss)
SimRegDat(n = 100, p = 200, type = "dep", rate = 0.1)
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{SimRegDat}
