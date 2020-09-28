// -*- mode: c++ -*-
// $Id: cfgPermanent.h 2502 2017-02-27 17:35:59Z axel $
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
