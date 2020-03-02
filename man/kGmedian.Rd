\name{kGmedian}
\alias{kGmedian}
\title{kGmedian}
\description{Fast k-medians clustering based on recursive averaged stochastic gradient algorithms.
The procedure is similar to the \code{\link[stats]{kmeans}} clustering technique performed recursively with the \code{MacQueen} algorithm. The advantage of the kGmedian algorithm compared to MacQueen strategy is that it deals with sum of norms instead of sum of squared norms, ensuring a more robust behaviour against outlying values. 
}
\usage{
kGmedian(X, ncenters=2, gamma=1, alpha=0.75, nstart = 10, nstartkmeans = 10, 
iter.max = 20)
}
\arguments{
\item{X}{matrix, with n observations (rows) in dimension d (columns).}
\item{ncenters}{Either the number of clusters, say k, or a set of initial (distinct) cluster centres. If a number, the initial centres are chosen as the output of the \code{\link[stats]{kmeans}} function computed with the \code{MacQueen} algorithm.}
\item{gamma}{Value of the constant controling the descent steps (see details).}
\item{alpha}{Rate of decrease of the descent steps.}
\item{nstart}{Number of times the algorithm is ran, with  random sets of initialization centers chosen among the observations.}
\item{nstartkmeans}{Number of initialization points in the \code{\link[stats]{kmeans}} function for choosing the starting point of \code{kGmedian}.}
\item{iter.max}{Maximum number of iterations considered in the \code{\link[stats]{kmeans}} function for choosing the starting point of \code{kGmedian}.}
}
\details{
See Cardot, Cenac and Monnez (2012).
}
\value{
\item{cluster}{A vector of integers (from 1:k) indicating the cluster to which each point is allocated.}
\item{centers}{A matrix of cluster centres.}
\item{withinsrs}{Vector of within-cluster sum of norms, one component per cluster.}
\item{size}{The number of points in each cluster.}
}
\references{
Cardot, H., Cenac, P. and Monnez, J-M. (2012). A fast and recursive algorithm for clustering large datasets with k-medians.  \emph{Computational Statistics and Data Analysis}, 56, 1434-1449.

Cardot, H., Cenac, P. and Zitt, P-A. (2013). Efficient and fast estimation of the geometric median in Hilbert spaces with an averaged stochastic gradient algorithm. \emph{Bernoulli}, 19, 18-43.

MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. In Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability, eds L. M. Le Cam and J. Neyman, 1, pp. 281-297. Berkeley, CA: University of California Press.
}

\seealso{
See also \code{\link{Gmedian}} and \code{\link[stats]{kmeans}}.
}

\examples{
# a 2-dimensional example 
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
           matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")

cl.kmeans <- kmeans(x, 2)
cl.kmedian <- kGmedian(x)

par(mfrow=c(1,2))
plot(x, col = cl.kmeans$cluster, main="kmeans")
points(cl.kmeans$centers, col = 1:2, pch = 8, cex = 2)

plot(x, col = cl.kmedian$cluster, main="kmedian")
points(cl.kmedian$centers, col = 1:2, pch = 8, cex = 2)
}
\keyword{Gmedian}

