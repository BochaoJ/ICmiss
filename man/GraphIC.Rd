\name{GraphIC}
\alias{GraphIC}
\title{ Learning high-dimensional Gaussian Graphical Models with Missing Observations.}
\description{The imputation-consistency (IC) algorithm for learning high-dimensional Gaussian Graphical Models with simulated incomplete data.}

\usage{
GraphIC(data, A, alpha1 = 0.05, alpha2 = 0.05, alpha3 = 0.05, iteration = 30, warm = 20)
}

\arguments{
  \item{ data }{\emph{\eqn{n}}x\emph{p} Dataset with missing values.}
  \item{ A }{ True adjacency matrix for evaluating the performance of the IC algorithm.}
\item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see \pkg{equSA}. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.05.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.05.}
\item{ alpha3  }{ The significance level of integrative \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix of IC_Ave method, the default value is 0.05.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}
}

\value{
   \item{ RecPre }{ The output of Recall and Precision values for the IC algorithm.}
  \item{ Adj }{  \emph{\eqn{p}}x\emph{p} Estimated adjacency matrix by our IC algorithm.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
#library(ICmiss)
#library(huge)
#result <- SimGraDat(n = 200, p = 100, type = "band", rate = 0.1)
#Est <- GraphIC(result$data, result$A, alpha1 = 0.05, alpha2 = 0.05, alpha3 = 0.05, iteration = 10, warm = 5)
#huge.plot(Est$Adj)  ## plot network by our estimated adjacency matrix.
#plot(Est$RecPre[,1], Est$RecPre[,2], type="l", xlab="Recall", ylab="Precision")  ## plot the Recall-Precision curve.
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{GraphIC}
