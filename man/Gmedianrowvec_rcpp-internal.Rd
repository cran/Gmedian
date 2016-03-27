\name{Gmedianrowvec_rcpp}
\alias{Gmedianrowvec_rcpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gmedianrowvec_rcpp
%%  ~~function to do ... ~~
}
\description{Internal function
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Gmedianrowvec_rcpp(X, gamma = 2, alpha = 0.75, nstart = 1L, epsilon = 1e-08)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
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
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, gamma = 2, alpha = 0.75, nstart = 1L, epsilon = 1e-08) 
{
    .Call("Gmedian_Gmedianrowvec_rcpp", PACKAGE = "Gmedian", 
        X, gamma, alpha, nstart, epsilon)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
