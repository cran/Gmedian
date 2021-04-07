\name{Weiszfeld}
\alias{Weiszfeld}
\title{Weiszfeld}

\description{Computes the Geometric median (also named spatial median or L1-median) with Weiszfeld's algorithm.
}
\usage{
Weiszfeld(X, weights = NULL, epsilon=1e-08, nitermax = 100) 
}
\arguments{
\item{X}{Data matrix, with n (rows) observations in dimension d (columns).}
\item{weights}{When \code{NULL}, all observations have the same weight, say 1/n. Else, the user can provide a size n vector of weights (such as sampling weights). These weights are used in the estimating equation (see details).} 
\item{epsilon}{Numerical tolerance. By defaut 1e-08.}
\item{nitermax}{Maximum number of iterations of the algorithm. By default set to 100.}
}
\details{
Weizfeld's algorithm (see Vardi and Zhang, 2000) is fast and accurate and can deal with large samples of high dimension data. However it is not as fast as the recursive approach proposed in \code{\link{Gmedian}}, which may be preferred for very large samples in high dimension. 
Weights can be given for statistical units, allowing to deal with data drawn from  unequal probability sampling designs (see Lardin-Puech, Cardot and Goga, 2014). 
}
\value{
\item{median}{Vector of the geometric median.}
\item{iter}{Number of iterations}
}
\references{
Lardin-Puech, P., Cardot, H.  and Goga, C. (2014). Analysing large datasets of functional data: a survey sampling point of view,   \emph{J. de la SFdS}, 155(4), 70-94. 

Vardi, Y. and Zhang, C.-H. (2000). The multivariate L1-median and associated data depth. \emph{Proc. Natl. Acad. Sci. USA}, 97(4):1423-1426.
}
\seealso{
See also \code{\link{Gmedian}} and \code{\link{WeiszfeldCov}}.
}
\examples{
## Robustness of the geometric median of n=3 points in dimension d=2.
a1 <- c(-1,0); a2 <- c(1,0); a3 <-c(0,1)
data.mat <- rbind(a1,a2,a3)
plot(data.mat,xlab="x",ylab="y")
med.est <- Weiszfeld(data.mat)
points(med.est$median,pch=19)

 ### weighted units
poids = c(3/2,1,1)
plot(data.mat,xlab="x",ylab="y")
med.est <- Weiszfeld(data.mat,weights=poids)
plot(data.mat,xlab="x",ylab="y")
points(med.est$median,pch=19)

## outlier
data.mat[3,] <- c(0,10) 
plot(data.mat,xlab="x",ylab="y")
med.est <- Weiszfeld(data.mat)
points(med.est$median,pch=19)

## Computation speed
## Simulated data - Brownian paths
n <- 1e2 ## choose n <- 1e5 for better evaluation
d <- 20
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))

system.time(replicate(10, {
  median.est = Weiszfeld(x)}))
system.time(replicate(10, {
  median.est = Gmedian(x)})) 
system.time(replicate(10, {
  mean.est = apply(x,2,mean)}))
}
\keyword{Weiszfeld}
