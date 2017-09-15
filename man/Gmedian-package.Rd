\name{Gmedian-package}
\alias{Gmedian-package}
%\alias{Gmedian}
\docType{package}
\title{
\packageTitle{Gmedian}
}
\description{
The geometric median (also called spatial median or L1 median) is a robust multivariate indicator of central position. This library provides fast estimation procedures  that can handle rapidly large samples of high dimensional data. Function \code{\link{Gmedian}} computes the geometric median of a numerical data set with averaged stochastic gradient algorithms, whereas \code{\link{GmedianCov}} computes the median covariation matrix, a useful indicator for robust PCA. Robust clustering, based on the geometric k-medians, can also be performed with the same type of recursive algorithm thanks to \code{\link{kGmedian}}.
Less fast estimation procedures  based on Weiszfeld's algorithm are also available : function  \code{\link{Weiszfeld}} computes the geometric median whereas \code{\link{WeiszfeldCov}} computes the median covariation matrix. These procedures may be preferred for small and moderate sample sizes. Note that weighting statistical units (for example with survey sampling weights) is allowed. 
%\packageDescription{Gmedian}
}
\details{

\packageDESCRIPTION{Gmedian}
\packageIndices{Gmedian}

}
\author{
\packageAuthor{Gmedian}

Maintainer: \packageMaintainer{Gmedian}
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Component7s Analysis. \emph{TEST}, 26, 461-480.

Cardot, H., Cenac, P. and Monnez, J-M. (2012). A fast and recursive algorithm for clustering large datasets with k-medians.  \emph{Computational Statistics and Data Analysis}, 56, 1434-1449.

Lardin, P., Cardot, H. and Goga, C. (2014). Analyzing large datasets of functional data : a survey sampling point of view.  \emph{Journal de la SFdS}, 155, 70-94.

Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.

}
\keyword{Gmedian}
