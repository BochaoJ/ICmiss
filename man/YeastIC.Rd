\name{YeastIC}
\alias{YeastIC}
\title{ Learning gene regulatory networks for Yeast Cell Expression Data.}
\description{An Imputation Consistency (IC) algorithm for learning gene regulatory networks with missing data. The dataset is collected from the yeast Saccharomyces cerevisiae responding to diverse environmental changes and is available at http://genome-www.stanford.edu/yeast-stress/.}

\usage{
YeastIC(data, alpha1 = 0.05, alpha2 = 0.01, alpha3 = 0.01, iteration = 30, warm = 20)
}

\arguments{
  \item{ data }{ \emph{\eqn{n}}x\emph{p} Yeast Cell expression data.}
 \item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see \pkg{equSA}. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.05.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.01.}
  \item{ alpha3  }{ The significance level of integrative \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix of IC_Ave method, the default value is 0.01.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}
}

\value{
  \item{ A }{  \emph{\eqn{p}}x\emph{p} Estimated adjacency matrix for network construction.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
## long time ##
#library(ICmiss)
#library(huge)
#data(yeast)
#A <- YeastIC(yeast, alpha1 = 0.05, alpha2 = 0.01, alpha3 = 0.01, iteration = 30, warm = 20)
#huge.plot(A)  ## plot gene regulatory network by our estimated adjacency matrix. 
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{YeastIC}
