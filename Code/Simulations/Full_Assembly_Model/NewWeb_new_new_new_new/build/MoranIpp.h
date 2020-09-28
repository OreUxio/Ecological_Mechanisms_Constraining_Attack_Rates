// -*- mode: c++ -*-
// $Id: MoranIpp.h 914 2007-07-02 13:36:03Z cvsrep $
#ifndef _MORANIPP_H_
#define _MORANIPP_H_

#include "MoranI.h"

class MoranIpp : public MoranI
{
protected:
  double _T_foraging;
  double _foraging_speedup;
  double _forager_over_variability;
  double _resource_sensitivity;
  double _granularity;

  double _niche_radius_squared;
  double _D_foraging;

  inline double niche_potential(double r2);  
  inline double radial_niche_force(double r2);
  inline double abs2(const sequence< double > &v);

public:
  MoranIpp(int ndim,double D,double C0,double lambda,
	   double saturation,double q,int co_co,
	   double T_foraging, 
	   double D_foraging, double forager_over_variability, 
	   double resource_sensitivity,double granularity);
  virtual Interaction_Matrix draw_sample(const double rS0);
};

#endif // _MORANIPP_H_
