// $Id: standardize.cc 1431 2009-05-04 12:22:40Z axel $

#include "NetworkAnalysis.h"

// adjustable parameters:
static int lump_lowest_level=1;
static int never_standardize=0;
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGINT(lump_lowest_level),
    CFGINT(never_standardize),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);


Interaction_Matrix standardize(const Interaction_Matrix & rim){
  if(never_standardize) return rim;
  if(rim.Number_of_Species_S()>0){
    if(lump_lowest_level){
      return rim.largest_connected_subweb().lump_lowest_level().trophic();
    }else{
      return rim.largest_connected_subweb().trophic();
    }
  }else{
    return rim;
  }
}


