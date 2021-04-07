\name{GmedianCov}
\alias{GmedianCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GmedianCov}
\description{Computes recursively the Geometric median and the (geometric) median covariation matrix  with fast averaged stochastic gradient algorithms. The estimation of the Geometric median is performed first and then the median covariation matrix is estimated, as well as its leading eigenvectors. The original recursive estimator of the median covariation matrix may not be a non negative matrix. A fast projected estimator onto the convex closed cone of the non negative matrices allows to get a non negative solution.}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\usage{
GmedianCov(X, init=NULL, nn=TRUE, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{Data matrix, with n observations (rows) in dimension d (columns).}
\item{init}{When \code{NULL} the starting point of the algorithm estimating the median is the first observation.} 
\item{nn}{When \code{TRUE} the algorithm provides a non negative estimates of the median covariation matrix. When \code{nn=FALSE}, the original algorithm is performed, with no guaranty that all the eigenvalues of the estimates are non negative}
\item{scores}{An integer \code{q}, by default \code{q=2}. The function computes the eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues and the corresponding principal component scores. No output if \code{scores=0}.}
\item{gamma}{Value (positive) of the constant controling the descent steps (see details) for the algorithm computing median.}
\item{gc}{Value (positive) of the constant controling the descent steps (see details) for algorithm computing the median covariation matrix}
\item{alpha}{Rate of decrease of the descent steps, \eqn{ 1/2 < alpha <=1}.}
\item{nstart}{Number of time the algorithms are ran.}
}
\details{
The (fast) computation of the eigenvectors is performed by  \code{\link[RSpectra]{eigs_sym}} of package \code{RSpectra}.
See Cardot, H. and Godichon-Baggioni (2017) for more details on the recursive algorithm. See also \code{\link{Gmedian}}.
When \code{nn=TRUE}, the descent step is bounded above so that the solution remains non negative at each iteration.
The principal components standard deviation is estimed robustly thanks to function \code{\link[robustbase]{scaleTau2}} from package \code{robustbase}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
\item{median}{Vector of the geometric median}
\item{covmedian}{Median covariation matrix}
\item{vectors}{The \code{scores=q} eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues}
\item{scores}{Principal component scores corresponding to the \code{scores=q} eigenvectors}
\item{sdev}{The \code{scores=q} estimates of the standard deviation of the \code{scores=q} principal components.}
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

Cardot, H. and Godichon-Baggioni, A. (2017). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis. TEST, 26, 461-480.
%% ~put references to the literature/web site here ~
}

\seealso{
See also  \code{\link{Gmedian}} and \code{\link{WeiszfeldCov}}.
}
\examples{
## Simulated data - Brownian paths
n <- 1e3
d <- 20
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))

## Estimation
median.est <- GmedianCov(x)

par(mfrow=c(1,2))
image(median.est$covmedian) ## median covariation function
plot(c(1:d)/d,median.est$vectors[,1]*sqrt(d),type="l",xlab="Time",
ylab="Eigenvectors",ylim=c(-1.4,1.4))
lines(c(1:d)/d,median.est$vectors[,2]*sqrt(d),lty=2)
}
\keyword{Gmedian}

