// -*- mode: c++ -*-
// $Id: mass_integrator.cc,v 1.1 2005/12/01 07:13:32 cvsrep Exp $

#include <fstream>
#include <sstream>
#include "error.h"
#include "mass_integrator.h"
#include "evaluate.h"

using namespace std;
static int mass_integrator_version=1;
static double conversion_factor_to_Mmat=1/0.571;
static double LFI_L_threshold_CS=50;
static double LFI_L_threshold_NS=40;
static double logLmax_logMmat_reg_a_CS=1.81817;
static double logLmax_logMmat_reg_b_CS=0.30819;
static double proplargefishB_Lmax_reg_a_CS=1.39089;
static double proplargefishB_Lmax_reg_b_CS=17.89367;
static double logLmax_logMmat_reg_a_NS=1.83214;
static double logLmax_logMmat_reg_b_NS=0.23904;
static double proplargefishB_Lmax_reg_a_NS_noqcorrection=1.39300;
static double proplargefishB_Lmax_reg_b_NS_noqcorrection=23.98780;
static double proplargefishB_Lmax_reg_a_NS_qcorrection=1.35103;
static double proplargefishB_Lmax_reg_b_NS_qcorrection=24.40813;
static double logLmax_logMmat_reg_a_CS_RMA=1.876;
static double logLmax_logMmat_reg_b_CS_RMA=0.3852;
static double proplargefishB_logLmax_reg_a_CS_RMA=0.008698;
static double proplargefishB_logLmax_reg_b_CS_RMA=2.552;
static double logLmax_logMmat_reg_a_NS_RMA=1.926;
static double logLmax_logMmat_reg_b_NS_RMA=0.3421;
static double proplargefishB_logLmax_reg_a_NS_noqcorrection_RMA=-0.09929;
static double proplargefishB_logLmax_reg_b_NS_noqcorrection_RMA=2.452;
static double proplargefishB_logLmax_reg_a_NS_qcorrection_RMA=-0.06042;
static double proplargefishB_logLmax_reg_b_NS_qcorrection_RMA=2.386;


// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGINT(mass_integrator_version),
  CFGDOUBLE(conversion_factor_to_Mmat),
  CFGDOUBLE(LFI_L_threshold_CS),
  CFGDOUBLE(LFI_L_threshold_NS),
  CFGDOUBLE(logLmax_logMmat_reg_a_CS),
  CFGDOUBLE(logLmax_logMmat_reg_b_CS),
  CFGDOUBLE(logLmax_logMmat_reg_a_NS),
  CFGDOUBLE(logLmax_logMmat_reg_b_NS),
  CFGDOUBLE(proplargefishB_Lmax_reg_a_CS),
  CFGDOUBLE(proplargefishB_Lmax_reg_b_CS),
  CFGDOUBLE(proplargefishB_Lmax_reg_a_NS_noqcorrection),
  CFGDOUBLE(proplargefishB_Lmax_reg_b_NS_noqcorrection),
  CFGDOUBLE(proplargefishB_Lmax_reg_a_NS_qcorrection),
  CFGDOUBLE(proplargefishB_Lmax_reg_b_NS_qcorrection),
  CFGDOUBLE(logLmax_logMmat_reg_a_CS_RMA),
  CFGDOUBLE(logLmax_logMmat_reg_b_CS_RMA),
  CFGDOUBLE(logLmax_logMmat_reg_a_NS_RMA),
  CFGDOUBLE(logLmax_logMmat_reg_b_NS_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_a_CS_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_b_CS_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_a_NS_noqcorrection_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_b_NS_noqcorrection_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_a_NS_qcorrection_RMA),
  CFGDOUBLE(proplargefishB_logLmax_reg_b_NS_qcorrection_RMA),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

Mass_Integrator::Mass_Integrator()
{
  switch(mass_integrator_version){
  case 1:
    {
      init("/usr/local/etc/NewWeb/massSum.dat");
    }
    break;
  case 2:
    {
      init("/home/tfung/NewWeb/massSum.dat");
    }
    break;
  case 3:
    {
      init("/home/tfung/NewWeb/massSum.dat");
    }
  break;
  case 4:
    {
    }
  break;
  case 5:
    {
    }
  break;
  case 6:
    {
    }
  break;
  case 7:
    {
    }
  break;
  case 8:
    {
    }
  break;
  case 9:
    {
    }
  break;
  default:
    FATAL_ERROR("mass_integrator_version " << mass_integrator_version << " unknown");
  }
}

void Mass_Integrator::init(const char * filename){
  ifstream sum_file(filename);
  if(sum_file.is_open()){
    while(!sum_file.eof()){
      double x,s;
      sum_file >> x >> s;
      _standard_sum[_standard_sum.size()]=std::pair<double,double>(x,s);
    }
  }else{
    WARNING("Cound not open """ << filename << """, using built-in table");
    for(int i=0;i< sum_data_length;i++){
      double x=sum_data[i][0];
      double s=sum_data[i][1];
      _standard_sum[_standard_sum.size()]=std::pair<double,double>(x,s);
    }
  }
  REPORT(_standard_sum.size());
}

/* Example:

_standard_sum:
0 1
1 0.5
2 0.25
3 0.12
4 0.06

n_points = 5
xmax=4
indexer=4/4=1
*/

/// Computes the biomass of individuals larger than startM
/** Now using linear interpolation of _standard_sum **/
double 
Mass_Integrator::integrate(const sequence< double > & B,
			   const sequence< double > & M,
			   double startM){
  const int n_points=_standard_sum.size();
  const double xmax=_standard_sum[n_points-1].first;
  const double indexer=(n_points-1)/xmax;
  
  ALWAYS_ASSERT(B.size()==M.size());
  double sum=0;

  for(int i=B.size();i-->0;){// for each species
    double index=startM*indexer/M(i);//index must be double to avoid
    //conversion errors!

    if(mass_integrator_version < 2)
      index/=conversion_factor_to_Mmat;
    
    // add up biomass above startM from all species
    if(index<n_points){
      if(index >= 1){
	sum+=B(i)*
	  (_standard_sum[index].second+
	   (index-int(index))*
	   ((index+1<n_points ? _standard_sum[index+1].second : 0)-
	    _standard_sum[index].second ))/_standard_sum[1].second;
      }else{
	sum+=B(i);
      }
    }
  }
  return sum;
}

// integrateThreshold is same as integrate except that only species with M above and below 2 thresholds are considered.
// used to calculate LFI where we only want to consider fish with M above and below 2 thresholds. 
double 
Mass_Integrator::integrateThreshold(const sequence< double > & B,
			   const sequence< double > & M,
				    double startM, double Mlowerthreshold, double Mupperthreshold){
  ALWAYS_ASSERT(B.size()==M.size());
  double sum=0;

  switch(mass_integrator_version){
  case 1: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav above and below 2 thresholds
	  double index=(startM*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid conversion errors!
	  // add up biomass above startM from all species
	  if(index<n_points){
	    sum+=B(i)*
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ));
	  }
	}
     }
   }
  case 2: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav above and below 2 thresholds
	  double index=(startM*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid conversion errors!
	  // add up biomass above startM from all species
	  if(index<n_points){
	    sum+=B(i)*
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ));
	  }
	}
     }
   }
   break;
  case 3: // No need to convert M(i) to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
      if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	double index=startM*indexer/M(i);//index must be double to avoid conversion errors!
	// add up biomass above startM from all species
	if(index<n_points){
	  sum+=B(i)*
	    (_standard_sum[index].second+
	     (index-int(index))*
	     ((index+1<n_points ? _standard_sum[index+1].second : 0)-
	      _standard_sum[index].second ));
	}
      }
     }
   }
   break;
  case 4: // There IS conversion of M(i) to maturation body mass. Based on Celtic Sea data.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
   {
     double Mkg=0;
     double log10Lmax=0;
     double Lmax=0;
     double proplargefishB=0;
     for(int i=B.size();i-->0;){// for each species
       if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	 // First convert M(i) into maturation body mass and then Lmax
	 Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	 log10Lmax=logLmax_logMmat_reg_a_CS+logLmax_logMmat_reg_b_CS*(log(Mkg)/log(10));
	 Lmax=pow(10,log10Lmax);
	 // Then convert Lmax into proportion
	 if(Lmax<LFI_L_threshold_CS){
	   proplargefishB=0;
	 }else{
	   proplargefishB=pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)/
	     (pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)+pow(proplargefishB_Lmax_reg_b_CS,proplargefishB_Lmax_reg_a_CS));
	 }
	 // Add biomass of large fish to sum
	 sum+=B(i)*proplargefishB;
       }
     }
   }
   break;
  case 5: // There IS conversion of M(i) to maturation body mass. Based on North Sea data, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)/
	      (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)
	       +pow(proplargefishB_Lmax_reg_b_NS_noqcorrection,proplargefishB_Lmax_reg_a_NS_noqcorrection));
	  }
	  // Add biomass of large fish to sum
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 6: // There IS conversion of M(i) to maturation body mass. Based on North Sea data, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)/
	      (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)
	       +pow(proplargefishB_Lmax_reg_b_NS_qcorrection,proplargefishB_Lmax_reg_a_NS_qcorrection));
	  }
	  // Add biomass of large fish to sum
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 7: // There IS conversion of M(i) to maturation body mass. Based on Celtic Sea data with RMA regression.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_CS_RMA+logLmax_logMmat_reg_b_CS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_CS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_CS_RMA+proplargefishB_logLmax_reg_b_CS_RMA*(log10Lmax-(log(LFI_L_threshold_CS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add biomass of large fish to sum
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 8: // There IS conversion of M(i) to maturation body mass. Based on North Sea data, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_noqcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_noqcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add biomass of large fish to sum
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 9: // There IS conversion of M(i) to maturation body mass. Based on North Sea data, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_qcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_qcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add biomass of large fish to sum
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  default:
    FATAL_ERROR("mass_integrator_version " << mass_integrator_version << " unknown");
  }

  return sum;
}
  
void
Mass_Integrator::spectrum(const sequence< double > & B,
			  const sequence< double > & M,
			  const char * filename,
			  const bool take_log10){
  const double step_factor=1.1;
  const double bin_center_factor=sqrt(step_factor);
  const int n_points=_standard_sum.size();
  const double xmax=_standard_sum[n_points-1].first;
  stringstream Bname;
  Bname << "B" << filename ;
  ofstream spectrum(filename);
  ofstream Bspectrum(Bname.str().c_str());
  
  ALWAYS_ASSERT(M.size()>0);

  double minM=M(0),maxM=M(0);
  for(int i=M.size();i-->0;){
    minM=min<double>(M(i),minM);
    maxM=max<double>(M(i),maxM);
  }

  double unit_mass=eval("1*kilogram");

  if(mass_integrator_version>2){
    for(double m=minM;m<maxM*xmax;m*=step_factor){
      double b=(integrate(B,M,m)-integrate(B,M,m*step_factor))/
	(m*(step_factor-1));
      spectrum << m*bin_center_factor/unit_mass << " "
	       << b
	       << endl;
      Bspectrum << m*bin_center_factor/unit_mass << " "
		<< m*sqrt(step_factor)*b/unit_mass
		<< endl;
    }
  }else{ // This is for case where M for each species is Mav rather than the maturation body mass. 
    // Mav is multiplied by a factor conversion_factor_to_Mmat to convert to the maturation body mass. 
    for(double m=minM*conversion_factor_to_Mmat;m<maxM*xmax*conversion_factor_to_Mmat;m*=step_factor){
      double b=(integrate(B,M,m)-integrate(B,M,m*step_factor))/
	(m*(step_factor-1));
      if(take_log10){
	spectrum << log10(m*bin_center_factor/unit_mass) << " "
		 << log10(b+1e-100)
		 << endl;
	Bspectrum << log10(m*bin_center_factor/unit_mass) << " "
		  << log10(m*sqrt(step_factor)*b/unit_mass+1e-100)
		  << endl;
      }else{
	spectrum << m*bin_center_factor/unit_mass << " "
		 << b
		 << endl;
	Bspectrum << m*bin_center_factor/unit_mass << " "
		  << m*sqrt(step_factor)*b/unit_mass
		  << endl;
      }	
    }
  }
}

// Computes the biomass of all species with M above and below 2 thresholds. 
// No integration is involved, but function is here
// because it is used to calculate LSI, and other functions
// here are used to calcualte the LFI. 
double 
Mass_Integrator::integrateLSI(const sequence< double > & B,
			   const sequence< double > & M,
			      double Mlowerthreshold, double Mupperthreshold){
  ALWAYS_ASSERT(B.size()==M.size());
  double sum=0;
  for(int i=B.size();i-->0;){
    if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){
      sum+=B(i);
    }
  }
  return sum;
}

// integrateProportion is same as integrate except that it works out the proportion of the biomass of species 
// within a certain range of body masses that is above the LFI body mass threshold.
double 
Mass_Integrator::integrateProportion(const sequence< double > & B,
			   const sequence< double > & M,
				     double Mlowerthreshold, double Mupperthreshold, double LFIMthreshold){
  ALWAYS_ASSERT(B.size()==M.size());
  double sum=0;
  double sumtotalbiomass=0;
  double proportion=0;

  switch(mass_integrator_version){
  case 1: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav > some threshold
	  double index=(LFIMthreshold*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid conversion errors!
	  sumtotalbiomass+=B(i);
	  // add up biomass above LFIMthreshold from all species
	  if(index<n_points){
	    sum+=B(i)*
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ));
	  }
	}
     }
   }
  case 2: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav > some threshold
	  double index=(LFIMthreshold*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid conversion errors!
	  sumtotalbiomass+=B(i);
	  // add up biomass above LFIMthreshold from all species
	  if(index<n_points){
	    sum+=B(i)*
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ));
	  }
	}
     }
   }
   break;
  case 3: // No need to convert M(i) to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=B.size();i-->0;){// for each species
       if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass > some threshold
	 double index=LFIMthreshold*indexer/M(i);//index must be double to avoid conversion errors!
	 sumtotalbiomass+=B(i);
	 // add up biomass above LFIMthreshold from all species
	 if(index<n_points){
	   sum+=B(i)*
	     (_standard_sum[index].second+
	      (index-int(index))*
	      ((index+1<n_points ? _standard_sum[index+1].second : 0)-
	       _standard_sum[index].second ));
	 }
       }
     }
   }
   break;
  case 4: // M(i) IS converted to maturation body mass. Based on Celtic Sea data.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_CS+logLmax_logMmat_reg_b_CS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_CS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)/
	      (pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)+pow(proplargefishB_Lmax_reg_b_CS,proplargefishB_Lmax_reg_a_CS));
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 5: // M(i) IS converted to maturation body mass. Based on North Sea data, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)/
	      (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)
	       +pow(proplargefishB_Lmax_reg_b_NS_noqcorrection,proplargefishB_Lmax_reg_a_NS_noqcorrection));
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 6: // M(i) IS converted to maturation body mass. Based on North Sea data, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)/
	      (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)
	       +pow(proplargefishB_Lmax_reg_b_NS_qcorrection,proplargefishB_Lmax_reg_a_NS_qcorrection));
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 7: // M(i) IS converted to maturation body mass. Based on Celtic Sea data.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_CS_RMA+logLmax_logMmat_reg_b_CS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_CS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_CS_RMA+proplargefishB_logLmax_reg_b_CS_RMA*(log10Lmax-(log(LFI_L_threshold_CS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 8: // M(i) IS converted to maturation body mass. Based on North Sea data, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_noqcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_noqcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  case 9: // M(i) IS converted to maturation body mass. Based on North Sea data, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=B.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_qcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_qcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add total biomss to sumtotalbiomss and biomass of large fish to sum
	  sumtotalbiomass+=B(i);
	  sum+=B(i)*proplargefishB;
	}
      }
    }
    break;
  default:
    FATAL_ERROR("mass_integrator_version " << mass_integrator_version << " unknown");
  }

  if(sumtotalbiomass==0){
    proportion=-1;
  } else{
    proportion=sum/sumtotalbiomass;
  }
  return proportion;
}

// integrateThresholdSpecies outputs a list showing the proportion of biomass above a bodymass threshold
// for each species with M above and below 2 thresholds. 
sequence< double > 
Mass_Integrator::integrateThresholdSpecies(const sequence< double > & M, 
					   double startM, double Mlowerthreshold, double Mupperthreshold){
  sequence< double > LFIalphas;

  switch(mass_integrator_version){
  case 1: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav above and below 2 thresholds
	  double index=(startM*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid
	  //conversion errors!
  
	  // calculate proportion of biomass above startM for a species and append to LFIalphas
	  // if(index>=n_points), then this is 0
	  if(index<n_points){
	    LFIalphas.push_back(
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ))
				);
	    //cout << "M(" << i << ") in IntegrateThresholdSpecies = " << M(i) << endl;
	  }else{
	    LFIalphas.push_back(0);
	  }
	}
      }
   }
  case 2: // M(i) is Mav rather than the maturation body mass.
    // Hence, M(i) is multiplied by factor conversion_factor_to_Mmat to convert to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with Mav above and below 2 thresholds
	  double index=(startM*indexer)/(M(i)*conversion_factor_to_Mmat);//index must be double to avoid
	  //conversion errors!
  
	  // calculate proportion of biomass above startM for a species and append to LFIalphas
	  // if(index>=n_points), then this is 0
	  if(index<n_points){
	    LFIalphas.push_back(
	      (_standard_sum[index].second+
	       (index-int(index))*
	       ((index+1<n_points ? _standard_sum[index+1].second : 0)-
		_standard_sum[index].second ))
				);
	    //cout << "M(" << i << ") in IntegrateThresholdSpecies = " << M(i) << endl;
	  }else{
	    LFIalphas.push_back(0);
	  }
	}
      }     
   }
   break;
  case 3: // No need to convert M(i) to maturation body mass.
   {
     const int n_points=_standard_sum.size();
     const double xmax=_standard_sum[n_points-1].first;
     const double indexer=(n_points-1)/xmax;
     for(int i=M.size();i-->0;){// for each species
      if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	double index=startM*indexer/M(i);//index must be double to avoid
	//conversion errors!
  
	// calculate proportion of biomass above startM for a species and append to LFIalphas
	// if(index>=n_points), then this is 0
	if(index<n_points){
	  LFIalphas.push_back(
	    (_standard_sum[index].second+
	     (index-int(index))*
	     ((index+1<n_points ? _standard_sum[index+1].second : 0)-
	      _standard_sum[index].second ))
			      );
	}else{
	  LFIalphas.push_back(0);
	}
      }
     }
   }
   break;
  case 4: // M(i) IS converted to maturation body mass. Based on data for Celtic Sea.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_CS+logLmax_logMmat_reg_b_CS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_CS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)/
	     (pow(Lmax-LFI_L_threshold_CS,proplargefishB_Lmax_reg_a_CS)+pow(proplargefishB_Lmax_reg_b_CS,proplargefishB_Lmax_reg_a_CS));
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  case 5: // M(i) IS converted to maturation body mass. Based on data for North Sea, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)/
	     (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_noqcorrection)
	      +pow(proplargefishB_Lmax_reg_b_NS_noqcorrection,proplargefishB_Lmax_reg_a_NS_noqcorrection));
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  case 6: // M(i) IS converted to maturation body mass. Based on data for North Sea, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted curve relating this proportion to Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS+logLmax_logMmat_reg_b_NS*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)/
	     (pow(Lmax-LFI_L_threshold_NS,proplargefishB_Lmax_reg_a_NS_qcorrection)
	      +pow(proplargefishB_Lmax_reg_b_NS_qcorrection,proplargefishB_Lmax_reg_a_NS_qcorrection));
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  case 7: // M(i) IS converted to maturation body mass. Based on data for Celtic Sea.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_CS_RMA+logLmax_logMmat_reg_b_CS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_CS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_CS_RMA+proplargefishB_logLmax_reg_b_CS_RMA*(log10Lmax-(log(LFI_L_threshold_CS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  case 8: // M(i) IS converted to maturation body mass. Based on data for North Sea, not corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_noqcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_noqcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  case 9: // M(i) IS converted to maturation body mass. Based on data for North Sea, corrected for q.
    // For each species, the proportion of biomass due to large fish is taken from a fitted straight line relating this proportion to log10Lmax
    // M(i) is first converted to Lmax using a regression line.
    {
      double Mkg=0;
      double log10Lmax=0;
      double Lmax=0;
      double proplargefishB=0;
      for(int i=M.size();i-->0;){// for each species
	if((M(i)>Mlowerthreshold)&&(M(i)<Mupperthreshold)){// only consider those species with maturation body mass above and below 2 thresholds
	  // First convert M(i) into maturation body mass and then Lmax
	  Mkg=(M(i)*conversion_factor_to_Mmat)/1000;
	  log10Lmax=logLmax_logMmat_reg_a_NS_RMA+logLmax_logMmat_reg_b_NS_RMA*(log(Mkg)/log(10));
	  Lmax=pow(10,log10Lmax);
	  // Then convert log10Lmax into proportion
	  if(Lmax<LFI_L_threshold_NS){
	    proplargefishB=0;
	  }else{
	    proplargefishB=proplargefishB_logLmax_reg_a_NS_qcorrection_RMA
	      +proplargefishB_logLmax_reg_b_NS_qcorrection_RMA*(log10Lmax-(log(LFI_L_threshold_NS)/log(10)));
	    // If >1, set to 1
	    if(proplargefishB>1){proplargefishB=1;}
	    // If <0, set to 0
	    if(proplargefishB<0){proplargefishB=0;}
	  }
	  // Add proportion of large fish to LFIalphas
	  LFIalphas.push_back(proplargefishB);
	}
      }
    }
    break;
  default:
    FATAL_ERROR("mass_integrator_version " << mass_integrator_version << " unknown");
  }

  return LFIalphas;
}


const int Mass_Integrator::sum_data_length=401;
const double Mass_Integrator::sum_data[][2]={
  {0,   1},
  {0.01,   0.7900922626281807},
  {0.02,   0.7455683129398812},
  {0.02999999999999999,   0.7157141924088227},
  {0.04,   0.6926201378338531},
  {0.05,   0.6735275302842935},
  {0.05999999999999999,   0.6571174126288373},
  {0.07,   0.6426472299717093},
  {0.08,   0.6296538230022664},
  {0.08999999999999999,   0.6178271945175133},
  {0.1,   0.6069487788303948},
  {0.11,   0.5968580996175924},
  {0.12,   0.5874337491167599},
  {0.13,   0.5785812721840545},
  {0.14,   0.5702257021603396},
  {0.1499999999999999,   0.562306474039871},
  {0.16,   0.5547738611720035},
  {0.17,   0.5475864836762642},
  {0.1799999999999999,   0.5407095634149987},
  {0.19,   0.5341134848301302},
  {0.2,   0.5277728515877809},
  {0.2099999999999999,   0.521665725436624},
  {0.22,   0.5157729592273012},
  {0.23,   0.5100778206211684},
  {0.2399999999999999,   0.5045654474094557},
  {0.25,   0.4992227432467749},
  {0.26,   0.4940380131724699},
  {0.27,   0.4890007405228247},
  {0.28,   0.4841014968770232},
  {0.2899999999999999,   0.4793318004508849},
  {0.2999999999999999,   0.4746838809408141},
  {0.31,   0.470150773323879},
  {0.32,   0.4657260107090854},
  {0.33,   0.461403791630859},
  {0.34,   0.4571786828126178},
  {0.35,   0.4530457838628146},
  {0.3599999999999999,   0.4490004668838966},
  {0.37,   0.4450385941849708},
  {0.38,   0.4411561781230792},
  {0.39,   0.4373496323144558},
  {0.4,   0.4336155625956116},
  {0.41,   0.4299507966417111},
  {0.4199999999999999,   0.4263524303335386},
  {0.4299999999999999,   0.4228176592009721},
  {0.44,   0.4193439435655309},
  {0.45,   0.4159288101594134},
  {0.46,   0.4125699902374538},
  {0.47,   0.4092653619878295},
  {0.4799999999999999,   0.4060128305281181},
  {0.4899999999999999,   0.4028104950315428},
  {0.5,   0.3996565481533047},
  {0.51,   0.3965492094501639},
  {0.52,   0.3934868457893873},
  {0.53,   0.3904678725221881},
  {0.54,   0.3874907687662739},
  {0.55,   0.3845540969598932},
  {0.56,   0.3816564437600102},
  {0.57,   0.3787964737315122},
  {0.5799999999999999,   0.3759728740489636},
  {0.5899999999999999,   0.3731843760321692},
  {0.5999999999999999,   0.3704297466728583},
  {0.6099999999999999,   0.3677077710632741},
  {0.62,   0.3650172608406364},
  {0.63,   0.3623570309741427},
  {0.64,   0.359725916664612},
  {0.65,   0.3571227576563361},
  {0.66,   0.3545463871137856},
  {0.67,   0.3519956343440533},
  {0.68,   0.3494693115067749},
  {0.69,   0.3469662167977828},
  {0.7,   0.3444851174151882},
  {0.7099999999999999,   0.3420247563234138},
  {0.7199999999999999,   0.339583833362786},
  {0.7299999999999999,   0.3371610096618546},
  {0.7399999999999999,   0.3347548955908347},
  {0.75,   0.3323640473505051},
  {0.76,   0.3299869583161712},
  {0.77,   0.3276220560633538},
  {0.78,   0.3252676956702573},
  {0.79,   0.3229221550770961},
  {0.8,   0.3205836559999081},
  {0.81,   0.3182502873502238},
  {0.82,   0.3159200221409738},
  {0.83,   0.3135908320987847},
  {0.8399999999999999,   0.3112605884514878},
  {0.8499999999999999,   0.3089269983093888},
  {0.8599999999999999,   0.3065877457215142},
  {0.87,   0.3042403936678463},
  {0.88,   0.3018824404754522},
  {0.89,   0.2995113252187974},
  {0.9,   0.2971244058177229},
  {0.91,   0.2947190201833527},
  {0.92,   0.2922924633887011},
  {0.93,   0.2898420288658483},
  {0.94,   0.2873650195754731},
  {0.95,   0.284858796009157},
  {0.9599999999999999,   0.2823207453066692},
  {0.9699999999999999,   0.2797484126003904},
  {0.9799999999999999,   0.2771394096891017},
  {0.9899999999999999,   0.2744915442973346},
  {1.,   0.2718028394901176},
  {1.01,   0.2690714632753657},
  {1.02,   0.2662959602951487},
  {1.03,   0.2634750356687274},
  {1.04,   0.2606078493038269},
  {1.05,   0.2576937613620397},
  {1.06,   0.2547326244790618},
  {1.07,   0.2517245518458485},
  {1.08,   0.2486701257692813},
  {1.09,   0.2455702622398257},
  {1.1,   0.2424262734516164},
  {1.11,   0.2392398678488733},
  {1.12,   0.2360130535970648},
  {1.13,   0.2327482555822421},
  {1.14,   0.2294481066187802},
  {1.15,   0.2261156283763675},
  {1.159999999999999,   0.2227539718895366},
  {1.169999999999999,   0.2193666128181927},
  {1.179999999999999,   0.2159571139576319},
  {1.189999999999999,   0.2125291649939169},
  {1.2,   0.2090866012849112},
  {1.21,   0.2056332665604093},
  {1.22,   0.2021730335557937},
  {1.23,   0.1987097450055109},
  {1.24,   0.1952472302028992},
  {1.25,   0.1917891896057431},
  {1.26,   0.1883392303889769},
  {1.27,   0.1849009003808382},
  {1.28,   0.181477484994502},
  {1.29,   0.1780721950997978},
  {1.3,   0.1746880901479451},
  {1.31,   0.171327930549113},
  {1.32,   0.1679944299877997},
  {1.33,   0.1646900490536796},
  {1.34,   0.1614170037984831},
  {1.35,   0.158177461169107},
  {1.36,   0.1549732716471956},
  {1.37,   0.1518061327762833},
  {1.38,   0.1486776733652268},
  {1.39,   0.1455891775684632},
  {1.4,   0.1425418256649517},
  {1.409999999999999,   0.1395367745016039},
  {1.419999999999999,   0.1365749003963651},
  {1.429999999999999,   0.1336568975460718},
  {1.439999999999999,   0.130783456038564},
  {1.45,   0.1279551352699365},
  {1.46,   0.1251722579152436},
  {1.47,   0.1224351224271231},
  {1.48,   0.1197439923655623},
  {1.49,   0.1170989322023125},
  {1.5,   0.1144999377411368},
  {1.51,   0.1119470002953051},
  {1.52,   0.1094400036407489},
  {1.53,   0.1069787250716325},
  {1.54,   0.1045629376624174},
  {1.55,   0.1021923845008498},
  {1.56,   0.09986669355313396},
  {1.57,   0.09758546532276728},
  {1.58,   0.09534830031324743},
  {1.59,   0.09315477711964599},
  {1.6,   0.091004393515325},
  {1.61,   0.08889663084048654},
  {1.62,   0.08683097046753567},
  {1.63,   0.08480688465679811},
  {1.64,   0.08282380066835655},
  {1.65,   0.08088113394838344},
  {1.66,   0.07897830004951347},
  {1.669999999999999,   0.07711471135138534},
  {1.679999999999999,   0.0752897602467943},
  {1.689999999999999,   0.07350283332381983},
  {1.7,   0.07175331733819151},
  {1.71,   0.07004059818441048},
  {1.72,   0.06836405750735039},
  {1.73,   0.06672307861286634},
  {1.74,   0.0651170451732312},
  {1.75,   0.06354534045290397},
  {1.76,   0.06200735060344292},
  {1.77,   0.06050247297660111},
  {1.78,   0.05903010619292305},
  {1.79,   0.05758964887295327},
  {1.8,   0.05618049987728011},
  {1.81,   0.05480207017524264},
  {1.82,   0.05345378758573021},
  {1.83,   0.05213508113135065},
  {1.84,   0.05084537983471187},
  {1.85,   0.0495841130919906},
  {1.86,   0.04835072587204182},
  {1.87,   0.04714468389885166},
  {1.88,   0.04596545432996464},
  {1.89,   0.0448125043229253},
  {1.9,   0.04368530070482718},
  {1.91,   0.0425833222783822},
  {1.919999999999999,   0.0415060786284063},
  {1.929999999999999,   0.04045308300215866},
  {1.939999999999999,   0.03942384864689854},
  {1.95,   0.03841788880988518},
  {1.96,   0.03743471673888003},
  {1.97,   0.03647385571344126},
  {1.98,   0.03553486346850241},
  {1.99,   0.03461730507689867},
  {2.,   0.03372074561146527},
  {2.01,   0.03284475014503736},
  {2.02,   0.03198888375045021},
  {2.03,   0.03115271142597061},
  {2.04,   0.03033582302381576},
  {2.049999999999999,   0.0295378372990688},
  {2.06,   0.02875837410969898},
  {2.069999999999999,   0.02799705331367566},
  {2.08,   0.02725349476896807},
  {2.089999999999999,   0.02652731833354557},
  {2.1,   0.0258181438653774},
  {2.109999999999999,   0.02512559455446672},
  {2.12,   0.0244493267536247},
  {2.129999999999999,   0.02378901520509857},
  {2.14,   0.02314433479797067},
  {2.149999999999999,   0.02251496042132341},
  {2.16,   0.02190056696423911},
  {2.169999999999999,   0.02130082931580023},
  {2.18,   0.02071542236508907},
  {2.189999999999999,   0.0201440237899837},
  {2.2,   0.01958634001118242},
  {2.21,   0.01904209381694553},
  {2.22,   0.01851100813996137},
  {2.23,   0.01799280591291836},
  {2.24,   0.01748721006850485},
  {2.25,   0.01699394353940928},
  {2.26,   0.01651272925831998},
  {2.27,   0.01604329250207586},
  {2.28,   0.01558538316565289},
  {2.29,   0.01513876527748285},
  {2.3,   0.01470320298600347},
  {2.31,   0.01427846043965253},
  {2.319999999999999,   0.01386430178686778},
  {2.33,   0.01346049117608694},
  {2.339999999999999,   0.01306679275574781},
  {2.35,   0.01268297027030631},
  {2.359999999999999,   0.01230880355414076},
  {2.37,   0.01194409540726601},
  {2.379999999999999,   0.01158864978157816},
  {2.39,   0.01124227062897323},
  {2.399999999999999,   0.01090476190134731},
  {2.41,   0.01057592755059643},
  {2.419999999999999,   0.01025557152861669},
  {2.43,   0.009943497787304123},
  {2.439999999999999,   0.009639510278554806},
  {2.45,   0.009343414557771899},
  {2.46,   0.009055039859720854},
  {2.47,   0.008774225072280152},
  {2.48,   0.008500807705872961},
  {2.49,   0.008234625270922399},
  {2.5,   0.007975515277851636},
  {2.51,   0.007723315237083791},
  {2.52,   0.007477862659042031},
  {2.53,   0.00723899511876773},
  {2.54,   0.007006558446688978},
  {2.55,   0.0067804137212534},
  {2.56,   0.006560423553798731},
  {2.569999999999999,   0.006346450555662675},
  {2.58,   0.006138357338182938},
  {2.589999999999999,   0.005936006512697251},
  {2.6,   0.005739260690543312},
  {2.609999999999999,   0.005547982483058852},
  {2.62,   0.005362038487449716},
  {2.629999999999999,   0.005181309667224106},
  {2.64,   0.005005680035070529},
  {2.649999999999999,   0.004835033602600032},
  {2.66,   0.00466925438142363},
  {2.669999999999999,   0.004508226383152373},
  {2.68,   0.004351833619397275},
  {2.689999999999999,   0.004199960101769384},
  {2.7,   0.004052491470775737},
  {2.71,   0.003909324993562715},
  {2.72,   0.003770362771400424},
  {2.73,   0.003635506907188176},
  {2.74,   0.003504659503825254},
  {2.75,   0.003377722664210968},
  {2.76,   0.003254598491244603},
  {2.77,   0.003135189087825467},
  {2.78,   0.003019397079710264},
  {2.79,   0.002907133447883749},
  {2.8,   0.002798315623787628},
  {2.81,   0.002692861180337578},
  {2.819999999999999,   0.00259068769044926},
  {2.83,   0.002491712727038338},
  {2.839999999999999,   0.002395853863020486},
  {2.85,   0.002303028671311364},
  {2.859999999999999,   0.002213154828110786},
  {2.87,   0.002126155289814357},
  {2.879999999999999,   0.002041960301379671},
  {2.89,   0.00196050059411347},
  {2.899999999999999,   0.001881706899322511},
  {2.91,   0.001805509948313539},
  {2.919999999999999,   0.001731840472393309},
  {2.93,   0.001660629202868567},
  {2.939999999999999,   0.001591806873664225},
  {2.95,   0.001525307028042598},
  {2.96,   0.001461070563807542},
  {2.97,   0.001399039494877415},
  {2.98,   0.001339155835170592},
  {2.99,   0.001281361598605433},
  {3.,   0.00122559879910031},
  {3.01,   0.001171809450573585},
  {3.02,   0.001119935566943628},
  {3.03,   0.001069920762100524},
  {3.04,   0.001021715079655693},
  {3.05,   0.0009752696988383236},
  {3.06,   0.0009305357806096034},
  {3.07,   0.0008874644859307138},
  {3.08,   0.0008460069757628434},
  {3.089999999999999,   0.0008061144110671743},
  {3.1,   0.0007677382198571249},
  {3.109999999999999,   0.0007308338464508743},
  {3.12,   0.0006953598088182133},
  {3.129999999999999,   0.0006612747006141005},
  {3.14,   0.0006285371154934881},
  {3.149999999999999,   0.0005971056471113339},
  {3.16,   0.0005669388891225906},
  {3.169999999999999,   0.0005379954702951382},
  {3.18,   0.0005102363788637357},
  {3.189999999999999,   0.0004836263636656396},
  {3.2,   0.0004581305046013722},
  {3.21,   0.0004337138815714604},
  {3.22,   0.0004103415744764269},
  {3.23,   0.0003879786632167984},
  {3.24,   0.00036659022769637},
  {3.25,   0.0003461424299213639},
  {3.26,   0.0003266052068696228},
  {3.27,   0.0003079493048225441},
  {3.28,   0.0002901454700615217},
  {3.29,   0.0002731644488679525},
  {3.3,   0.0002569769875232306},
  {3.31,   0.000241553832308753},
  {3.32,   0.0002268661589508481},
  {3.33,   0.0002128882887683201},
  {3.339999999999999,   0.0001995958744531829},
  {3.35,   0.0001869645689390867},
  {3.359999999999999,   0.0001749700251596833},
  {3.37,   0.0001635878960486223},
  {3.379999999999999,   0.0001527938345395556},
  {3.39,   0.0001425635867335079},
  {3.399999999999999,   0.0001328750524440808},
  {3.41,   0.000123708242914454},
  {3.419999999999999,   0.000115043253240417},
  {3.43,   0.0001068601785177576},
  {3.439999999999999,   0.00009913911384226513},
  {3.45,   0.00009186015430972743},
  {3.46,   0.00008500341025529345},
  {3.47,   0.00007855042615973278},
  {3.48,   0.00007248462445438934},
  {3.49,   0.00006678943963467043},
  {3.5,   0.0000614483061959845},
  {3.51,   0.00005644465863373897},
  {3.52,   0.00005176211552402465},
  {3.53,   0.00004738587884370828},
  {3.54,   0.00004330187677882647},
  {3.55,   0.00003949603952379541},
  {3.56,   0.00003595429727303197},
  {3.57,   0.00003266262135512502},
  {3.58,   0.00002960809431030573},
  {3.589999999999999,   0.00002677887560085258},
  {3.6,   0.00002416316310743469},
  {3.609999999999999,   0.00002174915471072154},
  {3.62,   0.00001952507009393599},
  {3.629999999999999,   0.00001747994897262958},
  {3.64,   0.00001560385164822105},
  {3.649999999999999,   0.00001388689929402965},
  {3.66,   0.0000123192130833743},
  {3.669999999999999,   0.00001089091505461895},
  {3.68,   9.592664761949843e-6},
  {3.689999999999999,   8.416400503805155e-6},
  {3.7,   7.354228339796448e-6},
  {3.71,   6.39825432953548e-6},
  {3.72,   5.540584532633831e-6},
  {3.73,   4.773425383659846e-6},
  {3.74,   4.090046438552644e-6},
  {3.75,   3.484313718986277e-6},
  {3.76,   2.949696142662469e-6},
  {3.77,   2.479982657921501e-6},
  {3.78,   2.069906772025203e-6},
  {3.79,   1.713762615413793e-6},
  {3.8,   1.406339416130424e-6},
  {3.81,   1.143206496869819e-6},
  {3.82,   9.201362588984711e-7},
  {3.83,   7.326254079389377e-7},
  {3.839999999999999,   5.758880083115572e-7},
  {3.85,   4.448342062867608e-7},
  {3.859999999999999,   3.380760351198248e-7},
  {3.87,   2.521707875854227e-7},
  {3.879999999999999,   1.839809829125983e-7},
  {3.89,   1.312709681369503e-7},
  {3.899999999999999,   9.082770567756031e-8},
  {3.91,   6.110193090287373e-8},
  {3.919999999999999,   3.877415793922973e-8},
  {3.93,   2.372779300192175e-8},
  {3.939999999999999,   1.386871858013727e-8},
  {3.95,   7.266974170588816e-9},
  {3.96,   3.745029133634972e-9},
  {3.97,   1.258512259687735e-9},
  {3.98,   2.46215233032655e-10},
  {3.99,   1.97948978454335e-11},
  {4.,   2.130221377363245e-11}
};