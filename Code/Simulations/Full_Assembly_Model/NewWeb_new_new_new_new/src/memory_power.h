// -*- mode: c++ -*-
// $Id: memory_power.h 893 2007-06-25 14:23:00Z cvsrep $
#ifndef _MEMORY_POWER_H_
#define _MEMORY_POWER_H_

#include "sequence.h"

class memory_power{
  const double base;
  sequence < double > pow;
  memory_power():base(0){};//no default constructor;
public:
  memory_power(double b);
  double operator()(int i);
};

#endif // _TEMPLATE_H_
