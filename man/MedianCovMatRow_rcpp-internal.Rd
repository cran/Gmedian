\name{MedianCovMatRow_rcpp}
\alias{MedianCovMatRow_rcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MedianCovMatRow_rcpp
%%  ~~function to do ... ~~
}
\description{Internal function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MedianCovMatRow_rcpp(X, Gmedian, gamma = 2, alpha = 0.75, nstart = 1L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Gmedian}{
%%     ~~Describe \code{Gmedian} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{nstart}{
%%     ~~Describe \code{nstart} here~~
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Gmedian, gamma = 2, alpha = 0.75, nstart = 1L) 
{
    .Call("Gmedian_MedianCovMatRow_rcpp", PACKAGE = "Gmedian", 
        X, Gmedian, gamma, alpha, nstart)
  }
}
