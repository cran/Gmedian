#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*using namespace Rcpp;*/

 // [[Rcpp::export]] 
Rcpp::NumericVector Gmedianrowvec_rcpp(const arma::mat X, const double gamma = 2, const double alpha = 0.75, const int nstart = 1, const double epsilon = 1e-8)
{
// X : n * p  matrix  
// Inputs
const int n = X.n_rows ;
const int p = X.n_cols ;
// Containers and intialisation of the algorithm
arma::rowvec medvec = X.row(0) ;
arma::rowvec medrm = X.row(0);
double poids, normxm ;
// Number of replications of the algorithm
for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
// Stochastic gradient algorithms
  for (int it = 1 ; it < n ; it++)
    {
      normxm = arma::norm(X.row(it)-medrm);
      if (normxm > epsilon) {
        poids = sqrt(p)*gamma*pow(it+1,-alpha)/normxm;
        medrm += poids * (X.row(it)-medrm) ;
      }
    medvec += (medrm-medvec)/(it+1);
    }
}  
// Returns ;
return Rcpp::wrap(medvec);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix MedianCovMatRow_rcpp(const arma::mat  X, const arma::rowvec  Gmedian, const double gamma = 2, const double alpha = 0.75, const int nstart = 1){
    // X : n * p  matrix
    // Inputs
    const int n = X.n_rows ;
    const int p = X.n_cols ;
    // Containers
    arma::rowvec diffmed = X.row(0)-Gmedian ;
    arma::mat medav = arma::trans(diffmed) * diffmed ; 
    arma::mat  diffmat(p,p), medrm(p,p) ;
    double nrmrm , poids ; 

    
// Initialization of the algorithm
//    diffmed = X.row(0)-Gmedian ;
//    medav =   arma::trans(diffmed) * diffmed ;
    medrm = medav ;
 
    
    // Number of replications of the algorithm
    for (int nbcomp = 0 ; nbcomp < nstart ; nbcomp++){
        // Stochastic gradient algorithms
        for (int it = 1 ; it < n ; it++)
        {
            diffmed = X.row(it)-Gmedian ;
            diffmat = arma::trans(diffmed)*diffmed;
            diffmat -= medrm ;
            nrmrm = arma::norm(diffmat,"fro") ; // Frobenius norm divided by the dimension
 //           nrmrm = arma::sqrt(arma::sum(arma::square(diffmat)))/p ;
            poids = p  * gamma * pow(it+1,-alpha) * pow(nrmrm,-1);
            medrm +=   poids*diffmat ;
            medav += (medrm-medav)/(it+1);
        }
    }
    return Rcpp::wrap(medav);
}



// [[Rcpp::export()]]  
    Rcpp::List stoKmed_rcpp(const arma::mat X, const arma::mat Xtot, const arma::mat centers, const double gamma=2, const double alpha = 0.75)
{
    // Inputs
    // X : n * p  matrix
    // centers : k * p matrix
    int n = X.n_rows,  p = X.n_cols, ntot = Xtot.n_rows;
    int k = centers.n_rows; /* Nombre de classes */
    int i, j,  it, inew = 0;
    arma::mat  centersrm = centers;
    arma::mat centersav = centers;
    arma::vec nc(k);
    arma::vec wss(k);
//    arma::vec nc(k);
    nc.fill(1); /* initialisation des tailles de classe */
    arma::vec clind(ntot);
    clind.fill(0.0);
    double alph = alpha;
  
	double poids, best, dd;
    /*Rboolean updated;*/
    
	/*for(j = 0; j < k; j++) nc(j) = 1;*/

    for(i = 0; i < n; i++) { /* boucle sur les individus */
        best = R_PosInf;
        for(j = 0; j < k; j++) {  /*calcul du centre le plus proche*/
            dd = arma::norm(X.row(i)-centersrm.row(j))/sqrt(p);
          if(dd < best) {
		best = dd;
		inew = j+1;
	    }
	 } 
        clind(i) = inew; /* affectation de l'ind. dans la nouvelle classe */
        it = clind(i) - 1; /* numero de la classe */
    /* mise a jour de RM*/
    /*for(c = 0; c < p; c++)*/
        if (best>0){
        poids = gamma/(pow(nc[it],alph)*sqrt(best));
        centersrm.row(it) +=   (X.row(i) - centersrm.row(it))*poids;
        }
        nc(it)++; /* incrementation de la taille de la classe */
	  /* mise a jour de AVeraged*/
	/*  for(c = 0; c < p; c++)*/
        centersav.row(it) +=  (centersrm.row(it)-centersav.row(it))/nc[it];
    
	} /* fin boucle sur ind */
  
	/* Calcul des distances intra-classes */
	/* initialisation des classes, des tailles de classe et des erreurs */
	wss.fill(0.0);
	nc.fill(0);
	clind.fill(0);

    /* boucle sur tous les n+k individus de la table xtot*/
    for(i = 0; i < ntot; i++) {
	best = R_PosInf;
	for(j = 0; j < k; j++) { /**/ /*calcul du centre j le plus proche*/
        dd = arma::norm(Xtot.row(i)-centersav.row(j))/sqrt(p);
        
	    /*for(c = 0; c < p; c++) {
		tmp = xtot[i+ntot*c] - centav[j+k*c];
		dd += tmp * tmp;
	    }
         */
        if(dd < best) {
		best = dd;
		inew = j+1;
	    }
	 } 
	 clind(i) = inew; /* affectation de l'ind. dans la nouvelle classe */
	 it = clind(i) - 1; /* numero de la classe */
	 nc(it)++;
	wss(it) += sqrt(best);
	}
// Returns ;
    Rcpp::List ret ;
    ret["wss"] = wss ;
    ret["centers"] = centersav ;
    ret["cl"] = clind ;
    ret["nc"] = nc ;
    return Rcpp::wrap(ret);
	
    }

