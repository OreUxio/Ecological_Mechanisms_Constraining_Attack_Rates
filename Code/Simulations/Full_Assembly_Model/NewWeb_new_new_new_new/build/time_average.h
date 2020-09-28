// -*- mode: c++ -*-
// $Id: time_average.h 1966 2010-11-10 19:35:42Z axel $
#ifndef _TIME_AVERAGE_H_
#define _TIME_AVERAGE_H_

#include "sequence.h"
#include "NewWeb.h"

/// Computes time averages of abundances and diagnostics.
class time_average 
{
public:
  time_average(NewWeb & web,
	       double time_average_max_t=10); // this may take a long
                                              // time!
  ~time_average();
};

#endif // _TIME_AVERAGE_H_
