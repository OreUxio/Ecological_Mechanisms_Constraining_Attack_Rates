// -*- mode: c++ -*-
// $Id: NewMatrix.h 3 2005-12-01 07:13:32Z cvsrep $
#ifndef _NEWMATRIX_H_
#define _NEWMATRIX_H_

#define ARMA_DONT_USE_WRAPPER
//#include "/home/axel/Downloads/armadillo/armadillo-6.700.6-ordqz/include/armadillo"
#include "armadillo"

#define NewMatrix arma::mat
#define NewIdentityMatrix arma::eye
#define NewZeroMatrix(X,Y) NewMatrix(X,Y,arma::fill::zeros)
#define NewVector arma::vec


// trick to absorbe () behind .size1() and .size2()
inline int zero_int_func(){
  return 0;
}
#define SIZE1 n_rows+zero_int_func
#define SIZE2 n_cols+zero_int_func
#define prod(X,Y) X*Y
#define inner_prod(X,Y) as_scalar((X).t() * (Y))
#define outer_prod(X,Y) ((X) * (Y).t())
#define COLUMN(M,i) M.col(i)
#define ROW(M,i) M.row(i)

inline void NewEigen(const NewMatrix &m, NewVector &er, NewVector &ei){
  arma::cx_vec ec;
#if 0 && (ARMA_VERSION_MAJOR > 7 || (ARMA_VERSION_MAJOR == 7 && ARMA_VERSION_MINOR >= 600))
  eig_gen(ec,m,"dc"); // use divide and concer method
#else
  eig_gen(ec,m); 
#endif
  er=real(ec);
  ei=imag(ec);
  return;
}
inline void NewEigen(const NewMatrix &m, NewVector &er){
  arma::cx_vec ec;
  eig_gen(ec,m);
  er=real(ec);
  return;
}

inline double NewVectorAbs(NewVector & v){
  return arma::norm(v,2);
}

void eigen_report(const arma::mat & M,const char * filename);
void eigen_report_cmplx(const arma::cx_mat & M,const char * filename);
NewMatrix Schur_complement(const NewMatrix & m, arma::uvec sel);
NewMatrix Schur_complement(const NewMatrix & m, int mini, int maxi);


#endif // _NEWMATRIX_H_
