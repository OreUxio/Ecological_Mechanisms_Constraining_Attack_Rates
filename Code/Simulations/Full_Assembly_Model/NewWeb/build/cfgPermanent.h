// -*- mode: c++ -*-
// $Id: cfgPermanent.h 1431 2009-05-04 12:22:40Z axel $
#ifndef _CFGPERMANENT_H_
#define _CFGPERMANENT_H_

#include "remember.h"

/// Encapsulates all configuration parameters for use by XMLStore.
class cfgPermanent: public permanent
{
 public:
  virtual void data_mapping();
};

#endif // _CFGPERMANENT_H_
