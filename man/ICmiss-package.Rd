\name{ICmiss-package}
\alias{ICmiss-package}
\docType{package}
\title{

An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond

}
\description{

Missing data are frequently encountered in high-dimensional data analysis, but they are usually difficult to deal with using standard algorithms, such as the EM algorithm and its variants. This package provides a general algorithm, the so-called imputation-consistency (IC) algorithm, for treating high-dimensional missing data problems. A variant of the IC algorithm, the so-called imputation-conditional consistency (ICC) algorithm, has also been provided in the package.

}
\details{
\tabular{rl}{
Package: \tab ICmiss\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
Date: \tab 2017-10-03\cr
License: \tab GPL-2\cr
}



This package illustrates the use of the IC/ICC algorithms in three modules:
\bold{}

The first module is to apply the IC algorithm to learning high-dimensional Gaussian Graphical Models (GGMs) in presence
of missing data with a simulated dataset \code{SimGraDat(n,p,...)} and
Yeast cell example \code{YeastIC(data,...)}.
\bold{}

The second module is to apply the ICC algorithm to variable selection for high-dimensional linear regression in presence of missing data. The simulation study covers both cases, the covariates are mutually independent and generally dependent, with the code \code{SimRegDat(n,p,...)}. The real data example is for Bardet-Biedl syndrome (Scheetz et al., 2006) with the dataset available in the R package \emph{flare}.
\bold{}

The third module is to apply the ICC algorithm to random coefficient models, where the random coefficients are treated as missing data. A simulated dataset \code{data(RCDat)} is included in the package, which can be used in \code{RCLM(RCDat)}.

}
\author{
Bochao Jia, Faming Liang
Maintainer: Bochao Jia<jbc409@ufl.edu>
}
\references{


Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.<doi:10.1080/01621459.2015.1012391>

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.<doi:10.1093/biomet/asn036>

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation-Consistency Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.


Jia, B., Xu, S., Xiao, G., Lamba, V., Liang, F. (2017) Inference of Genetic Networks from Next Generation Sequencing Data. Biometrics.


}
\keyword{ package }
\examples{
#library(ICmiss)
#result <- SimRegDat(n = 100, p = 200, type = "dep", rate = 0.1)
#RegICC(result$x, result$y, result$coef, type = "dep", iteration = 30, warm = 20)
}
