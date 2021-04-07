#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*using namespace Rcpp;*/


 // [[Rcpp::export]] 
Rcpp::List Weiszfeld_rcpp(const arma::mat X, const arma::rowvec weights, const double epsilon = 1e-08, int nitermax = 100)
{
// X : n * p  matrix  
// Inputs
const int n = X.n_rows ;
const int p = X.n_cols ;
// Containers and intialisation of the algorithm
arma::rowvec meanvec(p);
arma::rowvec medvec(p);
arma::rowvec poids(n);
double diffxn, diffmax = 0, normxm = 1;
int iter = 0;
/* Initialisation avec la moyenne */
meanvec.fill(0.0);
for (int it=0 ; it < n ; it++)
{
meanvec += (X.row(it)-meanvec)/(it+1);  
}
medvec = meanvec; /* to deal with constant value */
/* Boucle Weiszfeld */
while (iter < nitermax and normxm > epsilon )
{
  for (int it=0 ; it < n ; it++)
  {
    diffxn = arma::norm(X.row(it)-meanvec);
    if (diffxn > 0) { 
      poids(it) = weights(it)/diffxn;
      diffmax = 2; 
    }
    else poids(it)=0;
  }  
  if (diffmax > 0) {
    poids = poids/arma::sum(poids); /* normalisation */
    medvec = poids*X; /* mise a jour de la mediane */
  }
  else medvec = meanvec;
  normxm = arma::norm(medvec-meanvec)/sqrt(double(p));
  meanvec = medvec;
  iter++;
}
// Returns ;
    Rcpp::List ret ;
    ret["median"] = medvec ;
    ret["iter"] = iter ;
    return Rcpp::wrap(ret);
}  

// [[Rcpp::export]]
Rcpp::List MedianCovMatW_rcpp(const arma::mat  X, const arma::rowvec median_est, const arma::rowvec weights, const double epsilon = 1e-08, int nitermax = 100)
{
// Estimation of the Median Covariation Matrix with Weiszfeld's algorithm  
    // X : n * p  matrix
    // Inputs
    const int n = X.n_rows ;
    const int p = X.n_cols ;
    // Containers
    arma::mat Xcent(n,p);
 for (int it=0 ; it < n ; it++)
  {
    Xcent.row(it) = X.row(it)-median_est;
  }
arma::mat medinit = arma::trans(Xcent) * Xcent/n; 
arma::mat  medest(p,p);
arma::rowvec poids(n);
double diffxn, normxm = 1;
int iter = 0;

    
 while (iter < nitermax and normxm > epsilon )
{
medest.fill(0.0);
  for (int it=0 ; it < n ; it++)
  {
    diffxn = arma::norm(arma::trans(Xcent.row(it))*Xcent.row(it)-medinit,"fro");
    if (diffxn > 0) poids(it) = weights(it)/diffxn;
    else poids(it) = 0;
    medest +=poids(it)*arma::trans(Xcent.row(it))*Xcent.row(it);
    }  
  medest = medest/arma::sum(poids); /* normalisation */
  normxm = arma::norm(medest-medinit,"fro")/p;
  medinit = medest;
  iter++;
}
// Returns ;
    Rcpp::List ret ;
    ret["median"] = medest ;
    ret["iter"] = iter ;
    return Rcpp::wrap(ret);
}  

