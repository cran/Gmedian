GmedianCov <- function(X, init=NULL, scores=2, gamma=2, gc=2, alpha=0.75, nstart=1){
  ### Computation of the Geometric covariation matrix 
  ### with averaged stochastic gradient algorithms
  ### input : X   (n x p matrix, n observations, in dimension p)
  ### output : (geometric) median (1 x p numeric vector) and (geometric) median covariation matrix (p x p)  
  ### require library(rARPACK)
  Gmed.est = Gmedian(X,init=init,gamma=gamma,alpha=alpha,nstart=nstart)
  GMCM.est = MedianCovMatRow_rcpp(X,Gmedian=Gmed.est,gamma=gc,alpha=alpha,nstart=nstart)
  if (scores==FALSE){ 
  return(list(median = Gmed.est,covmedian=GMCM.est))
  }
  else {
    ### Computation of the eigenvectors and scores 
    vectors <- eigs_sym(GMCM.est, scores)$vectors
    scores = sweep(X,2,Gmed.est)%*%vectors
    return(list(median=Gmed.est,covmedian=GMCM.est,scores=scores,vectors=vectors))
  }
}
