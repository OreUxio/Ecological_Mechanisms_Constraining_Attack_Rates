// -*- mode: c++ -*-
// $Id: cfgPermanent.cc 2502 2017-02-27 17:35:59Z axel $

#include "cfgPermanent.h"
#include "cfgList.h"
#include "error.h"

void cfgPermanent::data_mapping(){
  cfgStruct * cfg=full_cfg_list();
  for(cfgStruct* i=cfg; i->type!=CFG_END; i++){
    switch(i->type){
    case CFG_BOOL:
      the_remember->sync(i->parameterName,*(bool *)(i->value));
      break;
    case CFG_STRING:
      FATAL_ERROR("string parameters won't work yet");
      break;
    case CFG_INT:
      the_remember->sync(i->parameterName,*(int *)(i->value));
      break;
    case CFG_UINT:
      the_remember->sync(i->parameterName,*(unsigned int *)(i->value));
      break;
    case CFG_LONG:
      the_remember->sync(i->parameterName,*(long *)(i->value));
      break;
    case CFG_ULONG:
      the_remember->sync(i->parameterName,*(unsigned long *)(i->value));
      break;
    case CFG_STRING_LIST:
      FATAL_ERROR("string_list parameters won't work yet");
      break;
    case CFG_FLOAT:
      the_remember->sync(i->parameterName,*(float *)(i->value));
      break;
    case CFG_DOUBLE:
      the_remember->sync(i->parameterName,*(double *)(i->value));
      break;
    }
  }
  delete[] cfg;
}


