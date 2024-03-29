// $Id: MoranIpp.cc 1404 2009-03-15 15:05:30Z axel $
#include <string>
const std::string version("$Revision: 1404 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/MoranIpp.cc,v $ $Date: 2009-03-15 15:05:30 +0000 (Sun, 15 Mar 2009) $");

#include "MoranIpp.h"
#include "simple_vector.h"
//#include "sequence.h"
#include "random.h"
//#include "linpack_eigen.h"
#include "Statistics.h"
#include "niche_width_finder.h"
#include "gsl/gsl_sf.h"
#include <algorithm>

using namespace std;

MoranIpp::MoranIpp(int ndim,double D,double C0,double lambda,
		   double saturation,double q,int co_co,
		   double T_foraging, 
		   double D_foraging, double forager_over_variability, 
		   double resource_sensitivity,double granularity):
  MoranI(ndim,D,C0,lambda,saturation,q,co_co,0),
  _T_foraging(T_foraging),
  _D_foraging(D_foraging),
  _forager_over_variability(forager_over_variability),
  _resource_sensitivity(resource_sensitivity),
  _granularity(granularity)
{
  if(_q<1){
    // select a subset of species
    FATAL_ERROR("q<1 not implemented");
  }
}

inline double MoranIpp::niche_potential(double r2){
  return 
    -(_D_foraging/_T_foraging)*exp(-r2*1/(2*_niche_radius_squared));
}

inline double MoranIpp::radial_niche_force(double r2){
  return 
    -sqrt(r2)*(1/_niche_radius_squared)*
    niche_potential(r2);
}

inline double MoranIpp::abs2(const sequence< double > &v){
  double sum=0;
  for(int i=v.size();i-->0;){
    sum+=(v(i)*v(i));
  }
  return sum;
}

static const double niche_over_width=2;

Interaction_Matrix MoranIpp::draw_sample(const double rS0){
  
  const int S0=int(rS0+unirand());
  
  // find some initial values for vulnerability traits:
  double total_variance;
  CLHEP::HepMatrix V_mat=neutral_resource_traits(S0,total_variance);
  
  // find some initial values for foraging traits:
  sequence< sequence< double > > V(S0,sequence< double >(_ndim));
  sequence< sequence< double > > F(S0,sequence< double >(_ndim));
  for(int i=S0;i-->0;){
    int target=random_integer(S0);
    for(int j=_ndim;j-->0;){
      V[i][j]=V_mat[i][j];
      F[i][j]=V_mat[target][j];
    }
  }

//   REPORT(F);
//   REPORT(V);

  // go niche_over_width standard deviations beyond mean of chi^2
  // distribution.
  const double radius_multiplier_squared=
    (1+niche_over_width*sqrt(2.0/_ndim));
  REPORT(radius_multiplier_squared);

  // compute niche radius:
  {
    static sequence< double > radius;
    _niche_radius_squared=radius[S0];
    if(!_niche_radius_squared){
      radius[S0]=
	_niche_radius_squared=
	niche_width_finder(_saturation,S0,_C0,_ndim)/
	radius_multiplier_squared;
    }
  }
  REPORT(_niche_radius_squared);

  const double forager_distribution_width_squared=
    (_T_foraging<1?
     _T_foraging*_niche_radius_squared :
     _niche_radius_squared );
  REPORT(forager_distribution_width_squared);
  
  const double forager_diffusion_stepsize=
    forager_distribution_width_squared*
    _granularity/(2*_D_foraging);
  REPORT(forager_diffusion_stepsize);

  const double forager_drift_stepsize=
    _granularity*sqrt(_niche_radius_squared)/
    radial_niche_force(_niche_radius_squared);
  REPORT(forager_drift_stepsize);

  const double stepsize=
    min(forager_diffusion_stepsize,forager_drift_stepsize);
  REPORT(stepsize);
  
  const double relaxation_time=4;
    
  const double sigma_F=sqrt(2*_D_foraging*stepsize);
  const double sigma_V=sqrt(2.0*stepsize);
  REPORT(sigma_F);
  REPORT(sigma_V);

  long int steps_left=(long int)(relaxation_time/stepsize);
  REPORT(steps_left);
  
  sequence< sequence< double > > delta_F(S0,sequence< double >(_ndim));
  sequence< sequence< double > > delta_V(S0,sequence< double >(_ndim));
  sequence< double > energy_list(S0);
  sequence< double > diff(_ndim);
  
  /***************************************************************/
  /*********************** main loop start ***********************/
  /***************************************************************/

//     REPORT(F);
//     REPORT(V);

  while(steps_left){
    
    double total_energy=0;

    // compute speciations;
    for(int n=poisson(stepsize*S0);n-->0;){
      // get parent
      int p=random_integer(S0);
      // get child
      int c;
      do{
	c=random_integer(S0);
      }while(p==c);
      F(c)=F(p);
      V(c)=V(p);
    }
    

    // compute physiological limitations:
    delta_V=V;
    delta_V*=-2*_saturation*stepsize;
    delta_F=F;
    delta_F*=-2*_saturation*(_D_foraging/_forager_over_variability)*stepsize;


    //     REPORT(delta_F);
    
    for(int i=S0;i-->0;){
      double energy=0;
      
      energy+=abs2(F(i))*_saturation*(_D_foraging/_forager_over_variability);


      // compute drift:
      for(int j=S0;j-->0;){
	diff=F(i)-V(j);
// 	if(i==0){
// 	  REPORT(diff);
// 	  REPORT(niche_potential(abs2(diff)));
// 	  REPORT(diff*((1.0/_niche_radius_squared)*niche_potential(abs2(diff))));
// 	}
	double U=niche_potential(abs2(diff));
	diff*=stepsize*(1.0/_niche_radius_squared)*U;
	delta_F(i)+=diff;
	delta_V(j)+=_resource_sensitivity*diff;
	  
	energy+=U;
      }
      
//     REPORT(delta_F);
      // compute diffusion;
      for(int k=_ndim;k-->0;){
	delta_F(i)(k)+=gaussian(0,sigma_F);
	delta_V(i)(k)+=gaussian(0,sigma_V);
      }
//     REPORT(delta_F);

      energy_list(i)=energy;
      total_energy+=energy;
    }
    
    F+=delta_F;
    V+=delta_V;
    
//     REPORT(F(0));
//     REPORT(V(3));
    

    steps_left--;


    if(steps_left%100==0){
      REPORT(steps_left);
      sort(&energy_list(0),&energy_list(S0));
      REPORT(energy_list);
      REPORT(total_energy);
    }
    
    if(steps_left%1000==0){
      ofstream ff("F.dat");
      ofstream vv("V.dat");
      for(int i=0;i<S0;i++){
	ff << F(i) << endl;
	vv << V(i) << endl;
      }
      cout << "wrote positions" << endl;
    }
  }
  
  /***************************************************************/
  /*********************** main loop end ***********************/
  /***************************************************************/

  ofstream ff("F.dat");
  ofstream vv("V.dat");
  for(int i=0;i<S0;i++){
    ff << F(i) << endl;
    vv << V(i) << endl;
  }



//   REPORT(F);
//   REPORT(V);

//   for(int i=S0;i-->0;){
//     F[i]=V[random_integer(S0)];
//   }
  


  sequence< double > d2_list(S0*S0);
  int pos=0;
  for(int fi=S0;fi-->0;){
    for(int vi=S0;vi-->0;){
      d2_list(pos++)=abs2(F(fi)-V(vi));
    }
  }
  sort(&d2_list(0),&d2_list(S0*S0));
  double cutoff_d2=d2_list(int(_C0*pos));
  REPORT(cutoff_d2);
  REPORT(_niche_radius_squared*radius_multiplier_squared);

  // compute links:
  Interaction_Matrix m(S0);

  for(int fi=S0;fi-->0;){
    for(int vi=S0;vi-->0;){
      diff=F(fi)-V(vi);
      if(abs2(diff) <= cutoff_d2//
	 //_niche_radius_squared*radius_multiplier_squared
  	 &&
  	 gsl_sf_erf_Q(V[fi][0]/total_variance)>
  	 gsl_sf_erf_Q(V[vi][0]/total_variance)-_lambda
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
