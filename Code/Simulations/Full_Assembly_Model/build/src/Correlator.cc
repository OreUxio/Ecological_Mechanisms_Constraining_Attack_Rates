// $Id: Correlator.cc 1762 2010-03-22 22:50:07Z axel $
#include <string>
const std::string version("$Revision: 1762 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/Correlator.cc,v $ $Date: 2010-03-22 22:50:07 +0000 (Mon, 22 Mar 2010) $");

#include <math.h>
#include <algorithm>
#include <CLHEP/Matrix/DiagMatrix.h>

#include "simple_vector.h"
#include "sequence.h"
#include "random.h"
#include "linpack_eigen.h"
#include "Statistics.h"
#include "niche_width_finder.h"
#include "gsl/gsl_sf.h"


#include <set>
#include <CLHEP/Matrix/Matrix.h>

class Correlator 
{
  class degree_t:public  sequence< sequence < double > >{
  public:
    double & operator()(int i,int j){
      if(i>=j)
	return (*this)[i][j];
      else
	return (*this)[j][i];
    }
  };
  
  typedef std::set < int > descendants_set_t;
  typedef descendants_set_t::iterator diter;

  degree_t get_degrees_of_separation(const int S0);
  CLHEP::HepMatrix get_spread_covariance(const int clade_size,degree_t degrees);
  inline double variance(const int S0,double time_since_speciaton);
  int _tree_uses;
  double _r;
  double _sigma;
  int _S0;
protected:
  CLHEP::
  HepMatrix correlated_normals();
public:
  Correlator(int S0, double r, double sigma, int tree_uses);
  void do_stats(int n_repetitions);
};


using namespace std;

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    //    CFGDOUBLE(),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);

Correlator::degree_t Correlator::get_degrees_of_separation(const int S0){
  
  degree_t degree;
  simple_vector< descendants_set_t > descendants(S0);

  for(int i=S0;i-->0;){
    descendants[i].insert(i);
  }

  int open_degrees=(S0*(S0-1))/2;
  double time_before_present=0;

  while(open_degrees){
    // link two species (speciation backward in time)
    int p,c; // parent and child indices
    
    // assume a speciation rate s=1:
    time_before_present+=random_exponential(1.0/S0);

    // get parent
    p=random_integer(S0);

    if(!descendants[p].empty()){

      // get child
      do{
	c=random_integer(S0);
      }while(p==c);
      
      if(!descendants[c].empty()){
	for(diter c_desc=descendants[c].begin();
	    c_desc!=descendants[c].end();
	    c_desc++ ){
	  for(diter p_desc=descendants[p].begin();
	      p_desc!=descendants[p].end();
	      p_desc++ ){
	    degree(*p_desc,*c_desc)=time_before_present;
	    open_degrees--;
	  }
	}
	for(diter c_desc=descendants[c].begin();
	    c_desc!=descendants[c].end();
	    c_desc++ ){
	  descendants[p].insert(*c_desc);
	}
	descendants[c].clear();
      }
    }
  }
//   {
//     HepMatrix deg(S0,S0);
//     for(int i=S0;i-->0;)
//       for(int j=S0;j-->0;)
// 	deg[i][j]=degree(i,j);
//     REPORT(deg);
//   }
  return degree;
}

inline double Correlator::variance(const int S0,double time_since_speciaton){
  return 2*_sigma*_sigma*(1-exp(-2*_r*time_since_speciaton));
}

CLHEP::HepMatrix Correlator::get_spread_covariance(const int S0,degree_t degrees){
  const int size_m=S0-1;
  CLHEP::HepMatrix cov(size_m,size_m,0);
  for(int i=size_m;i-->0;){
    cov[i][i]=variance(S0,degrees(i,size_m));
    for(int j=i;j-->0;){
      cov[i][j]=cov[j][i]=
	(variance(S0,degrees(i,size_m))+
	 variance(S0,degrees(j,size_m))-
	 variance(S0,degrees(i,j)) 
	 )*0.5;
    }
  }
  return cov;
}

Correlator::Correlator(int S0, double r, double sigma, int tree_uses):
  _tree_uses(tree_uses),
  _r(r),
  _sigma(sigma),
  _S0(S0)
{
}

static double smax(sequence< double > x){
  if(x.size()==0)
    return 0;
  double xmax=x[0];
  for(int i=x.size();i-->0;){
    if(x[i]>xmax) xmax=x[i];
  }
  return xmax;
}


CLHEP::HepMatrix Correlator::correlated_normals(){
    // compute phylogenetic tree:
  const degree_t degrees=get_degrees_of_separation(_S0);

  // compute time since first speciation:
  const double max_degree=smax(Map(::smax,degrees));

  // compute phylogenetic correlations:
  CLHEP::HepMatrix cholesky=
    get_spread_covariance(_S0,degrees);

  // compute Cholesky factorization of covariance:
  CholeskyL(cholesky);
  // alternatively, this might be faster: rCholeskyL(cholesky);
  
  // compute positions:
  CLHEP::HepMatrix normals(_S0-1,_tree_uses),ri(_S0,_tree_uses);
  for(int i=_S0-1;i-->0;){
    for(int d=_tree_uses;d-->0;){
      normals[i][d]=gaussian(0,1);
    }
  }
  ri.sub(1,1,cholesky^normals);

  return ri;
}

void Correlator::do_stats(int n_samples){
  int tree_uses_left=0;
  CLHEP::HepMatrix y;
  vector< double > amount(_S0);
  sequence< double > frac;
  ofstream maxfg("tmp.dat");
  double threshold=1e-5;
  int k=0;
  int last_k=0;
  int rounds=0;

  for(int nth_sample=0;nth_sample<n_samples;nth_sample++){
    if(nth_sample%100==0){
      cerr << 1-double(nth_sample)/n_samples << "\r";
      cerr.flush();
    }

    if(tree_uses_left==0){
      if(_tree_uses){
	y=correlated_normals();
	tree_uses_left=_tree_uses;
      }else{
	tree_uses_left=1;
      }
    }
    tree_uses_left--;

    //// first sweep: generation of data and basic stats
    double sum=0;// sum of all A
    double sum2=0; // sum of all A^2
    double maxA=0; // largest A
    double sub_maxA=0; // second largest A

    for(int i=_S0;i-->0;){
      double A=exp(_tree_uses ? y[i][tree_uses_left] : gaussian(0,_sigma));
      if(A>sub_maxA){
	if(A>maxA){
	  sub_maxA=maxA;
	  maxA=A;
	}else{
	  sub_maxA=A;
	}
      }
      sum+=
	amount[i]=A;
      sum2+=A*A;
    }
    
    double f_normalizer=1/sum;
    double fg_normalizer=sum/sum2;
    maxfg << (maxA * fg_normalizer-1) << " " 
	  << sub_maxA/maxA<< endl;

    //// Second sweep: collect data for DPF
    for(int i=_S0;i-->0;){
      double f=amount[i]*f_normalizer;
      if(f>threshold)
	frac[k++]=f;
    }
    last_k=k;
    rounds++;

  }
  
  //// Compute and save DPF:
  cerr << rounds << " rounds" << endl;
  sort(frac.begin(),frac.end());
  for(int i=0;i<frac.size();i++){
    double f=frac[i];
    if(1-f>0){
      double Zc=double(frac.size()-i)/rounds;
      cout << f/(1-f)
	   << " " 
	   << Zc 
	   << std::endl;
      if(f>=0.01 && (i==0 || frac[i-1]<0.01) ){
	cerr << "Zc01 " << Zc << std::endl;
      }
    }
  }
}

int main(void){
  set_random_seed(getpid());
  //double nu=0.35;//gives actual nu=0.5 with r=0.1, S=100
  //double nu=0.5;//gives actual nu=~0.7 with r=0.1, S=100
  //double nu=0.2;//gives actual nu=~0.5 with r=0.01, S=1000
  double nu=0.2;
  int S=1000;
  double sigma=sqrt(2*log(double(S)))/nu;
  cerr << "sigma = " << sigma << endl;
  Correlator c(/*S0=*/S,/*r=*/0.01,/*sigma=*/sigma,/*tree_uses=*/100);
  c.do_stats(1000);
}

