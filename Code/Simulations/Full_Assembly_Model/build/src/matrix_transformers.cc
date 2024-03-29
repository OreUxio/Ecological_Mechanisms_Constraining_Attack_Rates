// -*- mode: c++ -*-
// $Id: matrix_transformers.cc 2036 2010-12-21 23:02:18Z axel $

#include <fstream>
#include <algorithm>
#include <float.h>
#include <cmath>
#include "matrix_transformers.h"
#include "sequence.h"
#include "polyfit.h"
#include "evaluate.h"

static my_evaluator_t eval_here;

const double gram=eval_here("1*gram"); //get a dimensionless gram
const double unit_mass=eval_here("1*kilogram"); //for output files
const double meter2=eval_here("1*meter^2");
const double unit_area=eval_here("1*meter^2"); //for output files
const double year=eval_here("1*year");

const double log_DBL_MAX=5*log(DBL_MAX)/7.0; // we added some safety margin
const double log_DBL_MIN=5*log(DBL_MIN)/7.0; // we added some safety margin
static double strength_distribution_lowerthreshold=-1; // small enough to include all species

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(strength_distribution_lowerthreshold),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

Interaction_Matrix 
threshold(const link_strength_matrix & l,double th){
  Interaction_Matrix im(l.size());
  for(int i=l.size();i-->0;){
    for(int j=l.size();j-->0;){
      if(l[i][j]>th)
	im[i][j]=NetworkAnalysis::eats;
      else
	im[i][j]=NetworkAnalysis::none;
    }
  }
  return im; 
}

void
strength_distribution(const link_strength_matrix & l,double th, 
		      const char * filename){
  double Zc01,nu;
  strength_distribution(l,th,filename,Zc01,nu);
  REPORT(Zc01);
  REPORT(nu);
}

void
strength_distribution(const link_strength_matrix & l,double th, 
		      const char * filename,double &Zc01,double &nu){
  std::ofstream os(filename);
  sequence<double> strength;
  average_meter mean_c_star;
  int animals=0;
  for(int i=l.size();i-->0;){
    bool is_big_animal=false;
    double mean_c_star_sum=0;
    // A species is only considered if it is large enough 
    for(int j=l.size();j-->0;){
      double f=l[i][j];
      mean_c_star_sum+=f*f;
      if(f>th){
	strength[strength.size()]=f;
	is_big_animal=true;
      }
    }
    if(is_big_animal){
      animals++;
      mean_c_star.sample(mean_c_star_sum);
    }
  }
  std::sort(strength.begin(),strength.end());
  
  sequence<double> logr;
  sequence<average_meter> logZc;
  for(int i=strength.size();i-->0;){
    double Zc=double(strength.size()-i)/animals;
    if(strength[i]<1){
      double r= strength[i]/(1-strength[i]);
      os << r << " " << Zc
	 << std::endl;
      if(0.01<r && r<10){
	int k=logr.size();
	logr[k]=log(r);
	logZc[k].sample(log(Zc)+0.5);
	logZc[k].sample(log(Zc)-0.5);
      }
    }
    if(strength[i]>=0.01 && (i==0 || strength[i-1]<0.01) ){
      Zc01=Zc;
    }
  }
  REPORT(mean_c_star);
  if(logr.size()>=3){
    fitted_function f(logr,logZc,2);
    nu = -f[1];
  }else{
    nu = -1;
  }
}

void
strength_distribution_new(const link_strength_matrix & l,double th, 
	const char * filename,const sequence<double> &M){
  double Zc01,nu;
  strength_distribution_new(l,th,filename,M,Zc01,nu);
  REPORT(Zc01);
  REPORT(nu);
}

void
strength_distribution_new(const link_strength_matrix & l,double th, 
			  const char * filename,const sequence<double> &M,
			  double &Zc01,double &nu){
  std::ofstream os(filename);
  sequence<double> strength;
  average_meter mean_c_star;
  int animals=0;
  for(int i=l.size();i-->0;){
    bool is_big_animal=false;
    double mean_c_star_sum=0;
    // A species is only considered if it is large enough 
    if(M[i] > strength_distribution_lowerthreshold){
      for(int j=l.size();j-->0;){
	double f=l[i][j];
	mean_c_star_sum+=f*f;
	if(f>th){
	  strength[strength.size()]=f;
	  is_big_animal=true;
	}
      }
    }
    if(is_big_animal){
      animals++;
      mean_c_star.sample(mean_c_star_sum);
    }
  }
  std::sort(strength.begin(),strength.end());
  
  sequence<double> logr;
  sequence<average_meter> logZc;
  for(int i=strength.size();i-->0;){
    double Zc=double(strength.size()-i)/animals;
    if(strength[i]<1){
      double r= strength[i]/(1-strength[i]);
      os << r << " " << Zc
	 << std::endl;
      if(0.01<r && r<10){
	int k=logr.size();
	logr[k]=log(r);
	logZc[k].sample(log(Zc)+0.5);
	logZc[k].sample(log(Zc)-0.5);
      }
    }
    if(strength[i]>=0.01 && (i==0 || strength[i-1]<0.01) ){
      Zc01=Zc;
    }
  }
  REPORT(mean_c_star);
  if(logr.size()>=3){
    fitted_function f(logr,logZc,2);
    nu = -f[1];
  }else{
    nu = -1;
  }
}

link_strength_matrix in_fraction(link_strength_matrix ff){
  double s;
  for(int i=ff.size();i-->0;){
    s=0;
    for(int j=ff.size();j-->0;){
      s+=ff[i][j];
    }
    if(s>1e5*DBL_MIN){
      for(int j=ff.size();j-->0;){
	ff[i][j]*=1.0/s;
      }
    }else{
      if(s)
	WARNING("COMPUTATION OF IN FRACTIONS FAILED!!");
    }
  }
  return ff;
}
