// -*- mode: c++ -*-
// $Id: topology_generator.h 891 2007-06-25 14:22:15Z cvsrep $
#ifndef _TOPOLOGY_GENERATOR_H_
#define _TOPOLOGY_GENERATOR_H_

#include "NetworkAnalysis.h"

class Topology_Generator 
{
 public:
  Topology_Generator();
  virtual ~Topology_Generator();
  virtual Interaction_Matrix draw_sample(const double scale_parameter)=0;
};

#endif // _TOPOLOGY_GENERATOR_H_
