\name{GmedianCov}
\alias{GmedianCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GmedianCov}
\description{Computes recursively the Geometric median and the (geometric) median covariation matrix  with fast averaged stochastic gradient algorithms. The estimation of the Geometric median is performed first and then the median covariation matrix is estimated.}
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\usage{
GmedianCov(X, init=NULL, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{data matrix, with n observations (rows) in dimension d (columns).}
\item{init}{When \code{NULL} the starting point of the algorithm estimating the median is the first observation.}  
\item{scores}{an integer \code{q}, by default \code{q=2}. The function computes the eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues and the corresponding principal component scores. No output if \code{scores=NULL}.}
\item{gamma}{Value (positive) of the constant controling the descent steps (see details) for the algorithm computing median.}
\item{gc}{Value (positive) of the constant controling the descent steps (see details) for algorithm computing the median covariation matrix}
\item{alpha}{Rate of decrease of the descent steps, \eqn{ 1/2 < alpha <=1}.}
\item{nstart}{Number of time the algorithms are ran.}
}
\details{
The (fast) computation of the eigenvectors is performed by  \code{\link[rARPACK]{eigs_sym}} of package \code{rARPACK}.
See the technical report Cardot, H. and Godichon-Baggioni (2015) for more details on the recursive algorithm. See also \code{\link{Gmedian}}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ 
\item{median}{Vector of the geometric median}
\item{covmedian}{Median covariation matrix}
\item{vectors}{The \code{scores=q} eigenvectors of the median covariation matrix associated to the \code{q} largest eigenvalues}
\item{scores}{Principal component scores corresponding to the \code{scores=q} eigenvectors.}
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

Cardot, H. and Godichon-Baggioni, A. (2015). Fast Estimation of the Median Covariation Matrix with Application to Online Robust Principal Components Analysis.  \url{http://arxiv.org/abs/1504.02852}
%% ~put references to the literature/web site here ~
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}


\seealso{
See also  \code{\link{Gmedian}}.
}
\examples{
## Simulated data - Brownian paths
n <- 1e3
p <- 50
x <- matrix(rnorm(n*p,sd=1/sqrt(p)), n, p)
x <- t(apply(x,1,cumsum))

## Estimation
median.est <- GmedianCov(x)

par(mfrow=c(1,2))
image(median.est$covmedian) ## median covariation function
plot(c(1:p)/p,median.est$vectors[,1]*sqrt(p),type="l",xlab="Time",
ylab="Eigenvectors",ylim=c(-1.4,1.4))
lines(c(1:p)/p,median.est$vectors[,2]*sqrt(p),lty=2)
}
\keyword{Gmedian}

