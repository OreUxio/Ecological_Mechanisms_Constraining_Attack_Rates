// -*- mode: c++ -*-
// $Id: cell.cc 674 2006-10-17 13:08:49Z cvsrep $

#include <math.h>
#include "cell.h"
#include "NewSpecies.h"

static double cell_grid_resolution=1000;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(cell_grid_resolution),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

cellid_t::cellid_t(const NewSpecies & pt):
  std::vector<int>(niche_space_dimensions_D){
  for(int i=niche_space_dimensions_D;i-->0;){
    (*this)[i]=(int)floor(pt.vulnerability_V()[i]*cell_grid_resolution);
  }
}

bool cellid_t::operator<(const cellid_t & other) const{
  for(int i=niche_space_dimensions_D;i-->0;){
    if((*this)[i] < other[i]){
      return true;
    }else if((*this)[i] > other[i]){
      return false;
    }
  }
  return false;
};


int Cells::check_internal_consistency(){
  int fullsize=0;
  for(std::map<cellid_t,Cell>::iterator i=std::map<cellid_t,Cell>::begin();
      i!=std::map<cellid_t,Cell>::end();
      i++){
    ASSERT(i->second.size()>0);
    fullsize+=i->second.size();
  }
  std::vector<int> indices(fullsize);
  int k=0;
  for(std::map<cellid_t,Cell>::iterator i=std::map<cellid_t,Cell>::begin();
      i!=std::map<cellid_t,Cell>::end();
      i++){
    for(Cell::iterator j=i->second.begin();j!=i->second.end();j++){
      indices[k++]=*j;
    }
  }
  std::sort(indices.begin(),indices.end());
  for(int k=fullsize;k-->0;){
    if(indices[k]!=k){
      for(int l=0;l<fullsize;l++){
	std::cout << l << "|" << indices[l] <<" ";
      }
      std::cout << std::endl;
      FATAL_ERROR("cells garbled");
    }
  }
};
