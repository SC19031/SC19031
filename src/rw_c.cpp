#include <Rcpp.h>
using namespace Rcpp;
//' @title A random walk Metropolis sampler for generating the standard Laplace distribution
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution
//' @param sigma The sample variance
//' @param x0 Initial value
//' @param N Cycles
//' @return A set of data about x
//' @examples
//' \dontrun{
//' rw_c(1,20,200)
//' @export
// [[Rcpp::export]]

NumericVector rw_C (double sigma,double x0,int N){
  NumericVector x(N); 
  x[1] = x0; 
  NumericVector u = runif(N); 
  double k = 0;
  for (int i = 2; i < N+1 ;i++){ 
    double y = as<double>(rnorm(1, x[i-1], sigma));
    double t = exp(-abs(y))/exp(-abs(x[i-1]));
    if (u[i] <= t){
      x[i] = y; 
      k = k + 1;
    }
    else{ 
      x[i] = x[i-1];} 
  };
  return x;}