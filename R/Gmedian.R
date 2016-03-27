Gmedian <- function (X, init = NULL, gamma = 2, alpha = 0.75, nstart=2, epsilon=1e-08) 
{
  X <- as.matrix(X)
  if (!is.null(init)) X = rbind(init,X) ### initialisation   
  med.X = Gmedianrowvec_rcpp(X,gamma=gamma,alpha=alpha,nstart=nstart,epsilon=epsilon)
  return(med.X)
}
