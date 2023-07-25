#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXi;                  // variable size matrix, integers
using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
void computeDerativeCpp(Map<VectorXd> p, Map<MatrixXi> M, Map<VectorXd> r,
                                 double lambda0, 
                                 double lambda1, 
                                 double gamma01, 
                                 double gamma10, 
                                 double phi) {
  
  // get the size of the array
  size_t num_states = p.size();
  
  // initialize the return value
  // VectorXd r(num_states);
  // VectorXd r = p;
  
  // loop over states
  for(size_t i = 0; i < num_states; ++i) {
    
    // make zero
    double dp = 0.0;
    size_t instate;
    
    // get the flow in from lambda0
    // if ( M(i,3) > 0 ) {
      instate = M(i,2);
      dp += lambda0 * p(instate) * (double)M(i,3);  
    // }
    
    // get the flow in from lambda1
    // if ( M(i,5) > 0 ) {
      instate = M(i,4);
      dp += lambda1 * p(instate) * (double)M(i,5);
    // }
    
    // get the flow in from gamma01
    // if ( M(i,7) > 0 ) {
      instate = M(i,6);
      dp += gamma01 * p(instate) * (double)M(i,7);
    // }
    
    // get the flow in from gamma10
    // if ( M(i,9) > 0) {
      instate = M(i,8);
      dp += gamma10 * p(instate) * (double)M(i,9);
    // }
    
    // get the flow out
    double this_i = (double)M(i,0);
    double this_j = (double)M(i,1);
    dp -= (lambda0 * this_i + lambda1 * this_j + this_i * gamma01 + this_j * gamma10 + this_i * phi + this_j * phi) * p(i);
    
    // update vector
    r(i) = dp;
    
  }
  
  // return r;
  
}