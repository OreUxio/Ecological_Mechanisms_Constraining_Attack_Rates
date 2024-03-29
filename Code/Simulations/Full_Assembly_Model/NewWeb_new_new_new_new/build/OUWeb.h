// -*- mode: c++ -*-
// $Id: OUWeb.h 1404 2009-03-15 15:05:30Z axel $
#ifndef _OUWEB_H_
#define _OUWEB_H_

#include "topology_generator.h"
#include "memory_power.h"

#include <set>
#include <CLHEP/Matrix/Matrix.h>

class OUWeb : public Topology_Generator
{
  inline double variance(const int S0,double time_since_speciaton);
protected:
  CLHEP::
  HepMatrix neutral_resource_traits(const int S0,double & total_variance);
  const int _ndim;
  const double _D_bodymass;
  const double _C0;
  const double _lambda;
  const double _saturation;
  const double _q;
  const int _correct_correlations;
public:
  OUWeb(int ndim,double D_bodymass,double C0,double lambda,
	 double saturation,double q,int co_co);
  virtual Interaction_Matrix draw_sample(const double rS0);
};

#endif // _OUWEB_H_
