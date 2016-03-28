\name{Gmedian}
\alias{Gmedian}
\title{Gmedian}

\description{Computes recursively the Geometric median (also named spatial median or L1-median) with a fast averaged stochastic gradient algorithms that can deal rapidly with large samples of high dimensional data.
}
\usage{
Gmedian(X, init = NULL, gamma = 2, alpha = 0.75, nstart=2, epsilon=1e-08) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{data matrix, with n (rows) observations in dimension d (columns).}
\item{init}{When \code{NULL} the starting point of the algorithm is the first observation. Else the starting point of the algorithm is provided by \code{init}.}  
\item{gamma}{Value (positive) of the constant controling the descent steps (see details).}
\item{alpha}{Rate of decrease of the descent steps (see details). Should satisfy \eqn{1/2< alpha <= 1}.}
\item{nstart}{Number of times the algorithm is ran over all the data set.}
\item{epsilon}{numerical tolerance.}
}
\details{
The recursive averaged algorithm is described in Cardot, Cenac, Zitt (2013), with descent steps defined as \eqn{\alpha_n = gamma/n^{alpha}}.
}
\value{
Vector of the geometric median.
}
\references{
Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.
}
\seealso{
See also \code{\link{GmedianCov}} and \code{\link{kGmedian}}.
}
\examples{
## Simulated data - Brownian paths
n <- 1e4
d <- 500
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))

## Computation speed
system.time(replicate(10, {
  median.est = Gmedian(x)})) 
system.time(replicate(10, {
  mean.est = apply(x,2,mean)}))
##

## Accuracy with contaminated data
n <- 1e03
d <- 100
n.contaminated <- 50 ## 5% of contaminated observations
n.experiment <- 100
err.L2 <- matrix(NA,ncol=3,nrow=n.experiment)
colnames(err.L2) = c("mean (no contam.)", "mean (contam.)","Gmedian")

for (n.sim in 1:n.experiment){
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))
err.L2[n.sim,1] <- sum((apply(x,2,mean))^2/d)
ind.contaminated <- sample(1:n,n.contaminated) ## contam. units
x[ind.contaminated,] <- 5 
err.L2[n.sim,2] <- sum((apply(x,2,mean))^2/d)
err.L2[n.sim,3] <- sum(Gmedian(x)^2/d)
}
boxplot(err.L2,main="L2 error")
}
\keyword{Gmedian}
