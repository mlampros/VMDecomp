# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]

# include "VMDheader.h"

using namespace vmdR;


/**
 * Variational Mode Decomposition (1-dimensional)
 */

// [[Rcpp::export]]
Rcpp::List vmd_1d(arma::vec signal,
                  double alpha,
                  double tau,
                  arma::uword K,
                  bool DC,
                  arma::uword init,
                  double tol,
                  bool verbose = false) {

  VarModeDecomp vmdecomp;
  return vmdecomp.VMD_1D(signal,
                         alpha,
                         tau,
                         K,
                         DC,
                         init,
                         tol,
                         verbose);
}


/**
 * Variational Mode Decomposition (2-dimensional)
 */

// [[Rcpp::export]]
Rcpp::List vmd_2d(arma::mat signal,
                  double alpha,
                  double tau,
                  arma::uword K,
                  bool DC,
                  arma::uword init,
                  double tol,
                  bool verbose = false) {

  VarModeDecomp vmdecomp;
  return vmdecomp.VMD_2D(signal,
                         alpha,
                         tau,
                         K,
                         DC,
                         init,
                         tol,
                         verbose);
}
