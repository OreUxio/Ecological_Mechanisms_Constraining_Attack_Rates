// -*- mode: c++ -*-
// $Id: memory_power.cc 893 2007-06-25 14:23:00Z cvsrep $

#include "memory_power.h"

memory_power::memory_power(double b):
  base(b),pow(1)
{
  pow[0]=1;
};

double memory_power::operator()(int i){
  ASSERT(i>=0);
  if(i>=pow.size()){
    for(int j=pow.size();j<=i;j++){
      pow[j]=pow(j-1)*base;
    }
  }
  return pow[i];
}


