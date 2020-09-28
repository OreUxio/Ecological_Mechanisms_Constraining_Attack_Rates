// $Id: binomial_dist.cc 838 2007-01-09 06:46:53Z cvsrep $

#include "binomial_dist.h"
#include "random.h"
#include "error.h"

double binomial_dist(int n, int k, double p){
  if(k>n/2){
    k=n-k;
    p=1-p;
  }
  double q=1-p;
  double pq=p*q;
  double pBinomial=1;
  for(int i=1;i<=k;i++){
    pBinomial*=double(n-k+i)/i*pq;
  }
  for(int i=0;i<n-2*k;i++){
    pBinomial*=q;
  }
  return pBinomial;
}

void test_BinomialDist(){
  double one=0;
  double p=unirand();
  int n=int(unirand()*20);
  for(int k=0;k<=n;k++)
    one+=binomial_dist(n,k,p);
  REPORT(one);
}
