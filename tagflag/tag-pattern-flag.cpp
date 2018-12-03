#include <Rcpp.h>
using namespace Rcpp;    

// [[Rcpp::export]]
int countatleastonemaf(NumericMatrix x, NumericVector p) {
  int f = 0; // flagged
  int n=x.nrow();
  int m=x.ncol();
  double tmp12, tmp21, tmp;
  for(int i=0; i<n; i++)
    p[i]=sqrt(p[i]*(1-p[i]));
  for(int i=0; i < n-1 ; i++) {
    for(int j=i+1; j<n; j++) {
      for(int k=0; k<m; k++) {
	if(k==i || k==j)
	  continue;
	tmp=std::abs(x(i,k)/p(i) + x(j,k)/p(j));
	tmp12=std::abs(1.0/p(i) + x(i,j)/p(j));
	tmp21=std::abs(1.0/p(j) + x(i,j)/p(i));
	if (tmp > tmp12 && tmp > tmp21) {
	  // Rcpp::Rcout << i+1 << ',' << j+1 << ' ' << k+1 << ' ' << tmp << ' ' << tmp12 << std::endl;
	    f++;
	    break;
	}
      }
    }
  }
  return f;
}

// [[Rcpp::export]]
int countatleastone(NumericMatrix x) {
  int f = 0; // flagged
  int n=x.nrow();
  int m=x.ncol();
  double tmp12, tmp;
  for(int i=0; i < n-1; i++) {
    for(int j=i+1; j<n; j++) {
      for(int k=0; k<m; k++) {
	if(k==i || k==j)
	  continue;
	tmp=std::abs(x(i,k) + x(j,k));
	tmp12=std::abs(1.0+x(i,j));
	if (tmp > tmp12) {
	  // Rcpp::Rcout << i << ',' << j << ' ' << k << ' ' << tmp << ' ' << tmp12 << std::endl;
	  f++;
	  break;
	}
      }
    }
  }
  return f;
}

// [[Rcpp::export]]
int countpatt(NumericMatrix x) {
  int f = 0; // flagged
  int n=x.nrow();
  int m=x.ncol();
  double tmp, tmp12;
  for(int i=0; i < n-1; i++) {
    for(int j=i+1; j<n; j++) {
      for(int k=0; k<m; k++) {
	if(k==i || k==j)
	  continue;
	tmp=std::abs(x(i,k) + x(j,k));
	tmp12=std::abs(1.0+x(i,j));
	/* Rcpp::Rcout << i+1 << ',' << j+1 << ' ' << k+1 << ' ' << tmp << ' ' << tmp12 << std::endl; */
	if (tmp > tmp12) 
	  f++;
      }
    }
  }
  return f;
}

// [[Rcpp::export]]
int countpattmaf(NumericMatrix x, NumericVector p) {
  int f = 0; // flagged
  int n=x.nrow();
  int m=x.ncol();
  double tmp12, tmp21, tmp;
  for(int i=0; i<n; i++)
    p[i]=sqrt(p[i]*(1-p[i]));
  for(int i=0; i < n-1 ; i++) {
    for(int j=i+1; j<n; j++) {
      for(int k=0; k<m; k++) {
	if(k==i || k==j)
	  continue;
	tmp=std::abs(x(i,k)/p(i) + x(j,k)/p(j));
	tmp12=std::abs(1.0/p(i) + x(i,j)/p(j));
	tmp21=std::abs(1.0/p(j) + x(i,j)/p(i));
	/* Rcpp::Rcout << i+1 << ',' << j+1 << ' ' << k+1 << ' ' << tmp << ' ' << tmp12 << std::endl; */
	if (tmp > tmp12 && tmp > tmp21) {
	    f++;
	}
      }
    }
  }
  return f;
}

