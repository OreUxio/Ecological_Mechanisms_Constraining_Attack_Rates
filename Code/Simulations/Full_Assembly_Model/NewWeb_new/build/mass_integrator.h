// -*- mode: c++ -*-
// $Id: mass_integrator.h,v 1.1 2005/12/01 07:13:32 cvsrep Exp $
#ifndef _MASS_INTEGRATOR_H_
#define _MASS_INTEGRATOR_H_

#include <utility>
#include "sequence.h"

/// Computes and outputs individual-level size spectra.
class Mass_Integrator 
{
  sequence< std::pair< double, double > > _standard_sum;
  void init(const char * filename);
  static const int sum_data_length; 
  static const double sum_data[][2]; 
 public:
  Mass_Integrator();
  double integrate(const sequence< double > & B,const sequence< double > & M,
		   double startM);
  double integrateThreshold(const sequence< double > & B,const sequence< double > & M,
			    double startM, double Mlowerthreshold, double Mupperthreshold);
  void spectrum(const sequence< double > & B,const sequence< double > & M,
		const char * filename, const bool take_log10=false);
  double integrateLSI(const sequence< double > & B,const sequence< double > & M,
		      double Mlowerthreshold, double Mupperthreshold);
  double integrateProportion(const sequence< double > & B,const sequence< double > & M,
			     double Mlowerthreshold, double Mupperthreshold, double LFIMthreshold);
  sequence< double > integrateThresholdSpecies(const sequence< double > & M,
					       double startM, double Mlowerthreshold, double Mupperthreshold);
};

#endif // _MASS_INTEGRATOR_H_
