// $Id: OUWeb.cc 1404 2009-03-15 15:05:30Z axel $
#include <string>
const std::string version("$Revision: 1404 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/OUWeb.cc,v $ $Date: 2009-03-15 15:05:30 +0000 (Sun, 15 Mar 2009) $");

#include "OUWeb.h"
#include "simple_vector.h"
#include "sequence.h"
#include "random.h"
#include "linpack_eigen.h"
#include "Statistics.h"
#include "gsl/gsl_sf.h"
#include <gsl/gsl_cdf.h>

using namespace std;

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);



inline double OUWeb::variance(const int S0,double time_since_speciaton){
  if(_saturation){
    double V=1.0/_saturation;
    return 2*V*(1-exp(-2/V*time_since_speciaton));
  }else
    return 4*time_since_speciaton;
}



OUWeb::OUWeb(int ndim,double D_bodymass,double C0,
	       double lambda,double saturation,double q,int co_co):
  _ndim(ndim),
  _D_bodymass(D_bodymass),_C0(C0),
  _lambda(lambda),
  _saturation(saturation),
  _q(q),_correct_correlations(co_co)
{
  if(_q<1){
    // select a subset of species
    FATAL_ERROR("q<1 not implemented");
  }
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


CLHEP::HepMatrix OUWeb::neutral_resource_traits(const int S0,
					  double & total_variance){
  // compute clade total variance
  total_variance=2/_saturation;


  // compute positions:
  CLHEP::HepMatrix ri(S0,_ndim);
  double total_std=sqrt(1.0/_saturation);
  for(int i=S0;i-->0;){
    for(int d=_ndim;d-->0;){
      ri[i][d]=gaussian(0,total_std);
    }
  }
  return ri;
}


Interaction_Matrix OUWeb::draw_sample(const double rS0){

  const int S0=int(rS0+unirand());
  
  double total_variance;
  CLHEP::HepMatrix ri=neutral_resource_traits(S0,total_variance);

  // center first coordinate (body mass):
  double mean=0;
  for(int i=S0;i-->0;) mean+=ri[i][0];
  mean/=S0;
  for(int i=S0;i-->0;) ri[i][0]-=mean;

  // compute phenotypic distances between species:
#if 1 // this version is faster if your BLAS are fast:
  CLHEP::HepMatrix d2=-2*ri^(ri.T());
  for(int i=S0;i-->0;){
    double r2i=-0.5*(d2[i][i]);
    for(int j=S0;j-->0;){
      d2[i][j]+=r2i;
      d2[j][i]+=r2i;
    }
  }
#else // this version is potentially more accurate
  CLHEP::HepMatrix d2(S0,S0);
  for(int i=S0;i-->0;){
    d2[i][i]=0;
    for(int j=i;j-->0;){
      double sum=0;
      for(int n=_ndim;n-->0;){
	double delta=ri[i][n]-ri[j][n];
	sum+=delta*delta;
      }
      d2[i][j]=d2[j][i]=sum;
    }
  }
#endif

  // sample target species:
  simple_vector< int > target(S0);
  for(int i=S0;i-->0;){
    target[i]=random_integer(S0);
  }

  // compute niche radius:
  double niche_radius_squared;
  {
    static sequence< double > radius;
    niche_radius_squared=radius[S0];
    if(!niche_radius_squared){
      radius[S0]=
	niche_radius_squared=
	gsl_cdf_chisq_Pinv((_C0*S0-1)/(S0-1),_ndim)*total_variance;
    }
  }


  // compute links:
  Interaction_Matrix m(S0);

  for(int fi=S0;fi-->0;){
    for(int vi=S0;vi-->0;){
      if(d2[target[fi]][vi]<niche_radius_squared 
	 &&
	 gsl_sf_erf_Q(ri[fi][0]/total_variance)>
	 gsl_sf_erf_Q(ri[vi][0]/total_variance)-_lambda
	 ){
      //if(unirand()<_C0){
	m[fi][vi]=
	  NetworkAnalysis::eats;
      }else{
	m[fi][vi]=
	  NetworkAnalysis::none;
      }
    }
  }
  return m;
}
