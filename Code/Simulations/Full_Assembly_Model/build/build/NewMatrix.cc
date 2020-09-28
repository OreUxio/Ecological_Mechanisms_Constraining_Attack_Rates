// -*- mode: c++ -*-
// $Id: NewMatrix.cc 3 2005-12-01 07:13:32Z cvsrep $

#include "NewMatrix.h"
#include "error.h"

void eigen_report(const arma::mat & M,const char * filename){
  ASSERT(M.SIZE1() == M.SIZE2());
  const int S = M.SIZE1();
  
  NewVector ddr(S),ddi(S);
  NewEigen(M,ddr,ddi);    
  
  std::ofstream os(filename);
  for(int i=0;i<S;i++){
    os << ddr[i] << " " << ddi[i] << std::endl;
  }
  return;
}

void eigen_report_cmplx(const arma::cx_mat & M,const char * filename){
  ASSERT(M.SIZE1() == M.SIZE2());
  const int S = M.SIZE1();
  
  arma::cx_vec ec;
  eig_gen(ec,M);
  
  std::ofstream os(filename);
  for(int i=0;i<S;i++){
    os << real(ec[i]) << " " << imag(ec[i]) << std::endl;
  }
  return;
}

NewMatrix Schur_complement(const NewMatrix & m, arma::uvec sel){
  if(sel.size()==0) return NewMatrix(0,0);
  arma::ivec which(m.SIZE1());
  which.ones();
  which(sel).zeros();
  arma::uvec nonsel = arma::find(which);
  if(nonsel.size()==0) return m;

  return m(sel,sel) - m(sel,nonsel)*solve(m(nonsel,nonsel),m(nonsel,sel));
}

NewMatrix Schur_complement(const NewMatrix & m, int mini, int maxi){
  ALWAYS_ASSERT(maxi >= mini);
  arma::uvec sel(maxi-mini+1);
  for(int i=sel.size();i-->0;){
    sel[i]=mini+i;
  }
  return Schur_complement(m,sel);
}

