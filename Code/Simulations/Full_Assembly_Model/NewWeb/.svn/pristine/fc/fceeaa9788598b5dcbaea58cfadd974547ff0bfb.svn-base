// -*- mode: c++ -*-
// $Id: CompMatrix.cc 3 2005-12-01 07:13:32Z cvsrep $

#include "error.h"
#include "Statistics.h"
#include "CompMatrix.h"

using namespace std;

int CompMatrix(NewMatrix & hatC, NewMatrix & ihatC,
	       NewMatrix & epsAT, NewMatrix & A, NewMatrix C){
  // Iteration to get effective competition matrix hatC and overlaps alpha

  const int S=C.SIZE1();
  
  NewMatrix alpha=NewZeroMatrix(0,0);
  NewMatrix norm(S,S);
  norm = arma::diagmat(norm);

  for(int rep=1; rep <= 1000; rep++){
    cerr << rep << "      \r";
    cerr << endl;
    cerr.flush();
    
    int negative_count=0;
    for(int j=S;j-->0;){
      if(hatC(j,j) < 0) negative_count++;
      double absolute=abs(hatC(j,j));
      if(absolute==0){
	REPORT(j);
	WARNING("zero self-competition, isolating species!!");
	hatC.row(j)*=0;
	hatC.col(j)*=0;
	A.row(j)*=0;
	A.col(j)*=0;
	epsAT.row(j)*=0;
	epsAT.col(j)*=0;
	C.row(j)*=0;
	C.col(j)*=0;
	C(j,j)=1;
	hatC(j,j)=1;
	norm(j,j)=1;
      }else{
	norm(j,j)=1/sqrt(absolute);
      }
    }
    norm = arma::diagmat(norm);

    REPORT(negative_count);
    WARNING("Getting alpha");

    NewMatrix alpha_old = alpha;
    alpha = norm * hatC * norm;

    WARNING("Checking convergence");

    // Check convergence
    if(S==0) goto competition_convergence;
    if(alpha_old.SIZE1() == alpha.SIZE1() ){
      double max_diff=0;
      for(int i=S;i-->0;){
	for(int j=S;j-->0;){
	  if(isnan(alpha(i,j)) or isinf(alpha(i,j))){
	    REPORT(i);
	    REPORT(j);
	    REPORT(alpha(i,j));
	    FATAL_ERROR("isnan(alpha(i,j))");
	  }
	  double diff=alpha(i,j)-alpha_old(i,j);
	  max_diff=max(max_diff,abs(diff));
	}
      }
      REPORT(max_diff);
      if(max_diff<1e-5){
	goto competition_convergence;
      }
    }

    WARNING("Inverting alpha");

    NewMatrix ialpha(S,S);
    ialpha = arma::inv(alpha);

    WARNING("Computing inverse competition matrix");
     
    ihatC = norm * ialpha * norm;

    WARNING("Computing new competition matrix");

    hatC = epsAT * ihatC * A  + C;

    WARNING("End of iteration");
  }
  WARNING("Computation of competition matrix did not converge");

 competition_convergence:

  return 0;
}

using namespace arma;

#define JOIN_HORIZ(A,B,C) join_horiz(A,join_horiz(B, C))
#define JOIN_VERT(A,B,C) join_vert(A,join_vert(B, C))

template <typename T> 
class smaller_alpha_beta_ev {
  const T & A;
  const T & B;
public:
  smaller_alpha_beta_ev(const T & AA, const T & BB):A(AA),B(BB){};
  bool operator()(int n,int m) const{
    return abs(A(n,n)*B(m,m)) < abs(A(m,m)*B(n,n));
  }
};

void write_Schur_spectrum(const arma::cx_mat& A, const arma::cx_mat& B,
			  const char * filename){
  std::ofstream os(filename);
      
  for(int i=0;i<A.SIZE1();i++){
    if(abs(B(i,i))>1e-20*abs(A(i,i))){
      std::complex< double > ev = A(i,i)/B(i,i);
      os << ev.real() << " " << ev.imag() << std::endl;
    }else{
      os << "1e20" << " " << 0 << std::endl;
    }
  }
  return;
}


int CompMatrixSchur(NewMatrix & hatC, NewMatrix & ihatC,
		    NewMatrix & epsAT, NewMatrix & A, NewMatrix C){

  const int S=C.SIZE1();

  double nC = max(norm(A),norm(C));
  
  mat LL0 =
    JOIN_VERT(JOIN_HORIZ(zeros(S,S), zeros(S,S), eye(S,S)   ),
	      JOIN_HORIZ(-C/nC     , eye(S,S)  , epsAT/nC   ),
	      JOIN_HORIZ(A/nC      , zeros(S,S), zeros(S,S) ) );
	      
  mat KK0 =
    JOIN_VERT(JOIN_HORIZ(eye(S,S)  , zeros(S,S), zeros(S,S) ),
	      JOIN_HORIZ(zeros(S,S), zeros(S,S), zeros(S,S) ),
	      JOIN_HORIZ(zeros(S,S),-eye(S,S)  , zeros(S,S) ) );

  cx_mat LL=conv_to<cx_mat>::from(LL0);

  cx_mat KK=conv_to<cx_mat>::from(KK0);

  cx_mat KK1, LL1, QQT1, ZZ1;
  WARNING("Starting qz factorisation");
  REPORT(qz(LL1,KK1,QQT1,ZZ1,KK,LL));

  // Prepares sorting
  std::vector<double> index(3*S);
  for(int i=index.size();i-->0;){
    index[i]=i;
  }
  sort(index.begin(),index.end(),smaller_alpha_beta_ev<cx_mat>(KK1, LL1));

  uvec selection(3*S,fill::zeros);
  for(int i=S;i-->0;) selection(index[i])=1;

  int replace_n=0;
  for(int i=0; i < replace_n; i++){
    selection[index[S+i]]=1;
    selection[index[S-1-i]]=0;
  }
  
  write_Schur_spectrum(KK1, LL1, "Schur_spectrum.dat");

  double pl, pr, difu, difl;
  
  WARNING("Extracting deflating subspace");
  FATAL_ERROR("ordqz not available, code needs updating");
  //REPORT(ordqz(LL1,KK1,ZZ1,LL1,KK1,ZZ1,selection,pl,pr,difu,difl));

  REPORT(pl);
  REPORT(pr);
  REPORT(difu);
  REPORT(difl);
  
  hatC = nC      * real( ZZ1.submat(S,0,2*S-1,S-1) * inv(ZZ1.submat(0,0,S-1,S-1)) );
  ihatC = (1/nC) * real( ZZ1.submat(0,0,S-1,S-1) * inv(ZZ1.submat(S,0,2*S-1,S-1)) );

  
  REPORT(norm(imag( ZZ1.submat(S,0,2*S-1,S-1) * inv(ZZ1.submat(0,0,S-1,S-1)))));
  REPORT(norm(epsAT*ihatC*A + C - hatC)/norm(hatC));
  
  return 0;
}
