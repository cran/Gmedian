\name{WeiszfeldCov}
\alias{WeiszfeldCov}
\title{WeiszfeldCov}

\description{Estimation of the Geometric median covariation matrix  with Weiszfeld's algorithm. Weights (such as sampling weights) for statistical units are allowed.
}
\usage{
WeiszfeldCov(X, weights=NULL, scores=2, epsilon=1e-08, nitermax = 100) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{Data matrix, with n (rows) observations in dimension d (columns).}
\item{weights}{When \code{NULL}, all observations have the same weight, say 1/n. Else, the user can provide a size n vector of weights (such as sampling weights). These weights are used in the estimating equation (see details).} 
\item{scores}{An integer \code{q}, by default \code{q=2}. The function computes the eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues and the corresponding principal component scores. No output if \code{scores=NULL}.}
\item{epsilon}{Numerical tolerance. By defaut 1e-08.}
\item{nitermax}{Maxium number of iterations of the algorithm. By default set to 100.}
}
\details{
This fast and accurate iterative algorithm can deal with moderate size datasets (see . For large datasets use preferably \code{\link{GmedianCov}}, if fast estimations are required.
Weights can be given for statistical units, allowing to deal with data drawn from  unequal probability sampling designs (see Lardin-Puech, Cardot and Goga, 2014). 
}
\value{ 
\item{median}{Vector of the geometric median}
\item{covmedian}{Median covariation matrix}
\item{vectors}{The \code{scores=q} eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues}
\item{scores}{Principal component scores corresponding to the \code{scores=q} eigenvectors.}
\item{iterm}{Number of iterations needed to estimate the median}
\item{itercov}{Number of iterations needed to estimate the median covariation matrix}
}

\references{
Cardot, H. and Godichon-Baggioni, A. (2015). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \url{http://arxiv.org/abs/1504.02852}

Lardin-Puech, P., Cardot, H.  and Goga, C. (2014). Analysing large datasets of functional data: a survey sampling point of view,   Journal de la Soc. Fr. de Statis., 155(4), 70-94. 
}
\seealso{
See also \code{\link{Weiszfeld}} and \code{\link{GmedianCov}}.
}
\examples{
## Simulated data - Brownian paths
n <- 1e3
d <- 50
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))

## Estimation
median.est <- WeiszfeldCov(x)

par(mfrow=c(1,2))
image(median.est$covmedian) ## median covariation function
plot(c(1:d)/d,median.est$vectors[,1]*sqrt(d),type="l",xlab="Time",
ylab="Eigenvectors",ylim=c(-1.4,1.4))
lines(c(1:d)/d,median.est$vectors[,2]*sqrt(d),lty=2)
}
\keyword{WeiszfeldCov}
