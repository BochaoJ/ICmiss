\name{EyeICC}
\alias{EyeICC}
\title{ Variable selection for Bardet-Biedl syndrome data with missing observations.}
\description{The imputation-conditional consistency (ICC) algorithm is used to select variables for the Bardet-Biedl syndrome data with missing observations: We first randomly delete a specified percentage of observations and then apply the ICC algorithm for variable selection.}

\usage{
EyeICC(x, y, rate = 0.05, alpha1 = 0.1, alpha2 = 0.1, iteration = 30, warm = 20)
}

\arguments{
  \item{ x }{ a \emph{\eqn{n}}x\emph{p} covariates matrix.}
  \item{ y }{ a \emph{\eqn{n}}x\emph{1} responses.}
  \item{ rate }{ Missing rate, the default value is 0.05 .}
  \item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see \pkg{equSA}. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in  
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.1.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.1.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}

}

\value{
  \item{topVar}{Variables ranked by the frequency of appearance in the last few iterations.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
## long time ##
#library(ICmiss)
#data(eye_norm)
#EyeICC(eye_norm$x, eye_norm$y, rate = 0.05, alpha1 = 0.1, alpha2 = 0.1)
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{EyeICC}
