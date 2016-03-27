\name{stoKmed_rcpp}
\alias{stoKmed_rcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{stoKmed_rcpp
%%  ~~function to do ... ~~
}
\description{Internal function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stoKmed_rcpp(X, Xtot, centers, gamma = 2, alpha = 0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Xtot}{
%%     ~~Describe \code{Xtot} here~~
}
  \item{centers}{
%%     ~~Describe \code{centers} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, Xtot, centers, gamma = 2, alpha = 0.75, nstart = 1L) 
{
    .Call("Gmedian_stoKmed_rcpp", PACKAGE = "Gmedian", X, Xtot, 
        centers, gamma, alpha, nstart)
  }
}

