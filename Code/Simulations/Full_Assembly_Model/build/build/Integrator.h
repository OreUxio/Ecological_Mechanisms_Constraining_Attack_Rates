// -*- mode: c++ -*-
// $Id: Integrator.h 2067 2011-01-26 19:22:44Z tak $
#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "sequence.h"

#define USE_SIMPLE_INTEGRATOR


/// Compute time integral of variables given upon construction.
class Integrator 
{
private:
  const sequence< double * > & the_locations;
  sequence< double > the_integrals;

  sequence< double > y1;
#ifdef USE_SIMPLE_INTEGRATOR
  bool write_back,started;
  double start_t,last_t;
#else // use complex integrator
  sequence< double > y2;
  bool write_back;
  double start_t,last_t,x1,x2;
  int num;
#endif

public:
  /// \a l list is of pointers to variables to be integrated.
  Integrator(const sequence< double * > & l);
  /// Writes integral back to variables if not disabled.
  ~Integrator();
  /// Disables write-back upon destruction of Integrator.
  void disable_write_back();
  /// Sample variables for integration.
  void sample(double t);
  /// Multiply all integrals by \a x.
  void operator*=(double x);
  // Return the_integrals.
  sequence< double > get_integrals();
  // Resets integrals to 0 and started to false.
  void reset();
};





#endif // _INTEGRATOR_H_
