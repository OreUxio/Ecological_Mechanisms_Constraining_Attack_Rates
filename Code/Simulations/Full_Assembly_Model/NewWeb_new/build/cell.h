// -*- mode: c++ -*-
// $Id: cell.h 2502 2017-02-27 17:35:59Z axel $
#ifndef _CELL_H_
#define _CELL_H_

#include <vector>
#include <map>

#include "error.h"

class NewSpecies;
typedef const NewSpecies * species_pt;

class cellid_t : private std::vector<int>{
public:
  cellid_t(const NewSpecies & pt);
  bool operator<(const cellid_t & other) const; 
};

typedef std::vector<int> Cell_base;

class Cell : public Cell_base{ //public because to lazy to define gateways
public:
  void insert(int i){
    push_back(i);
  }
  void move(int old_i,int new_i){
    for(iterator i=begin();i!=end();i++){
      if(*i==old_i){
	*i=new_i;
	break;
      }
    }
  }
  void remove(int j){
    for(iterator i=begin();i!=end();i++){
      if(*i==j){
	erase(i);
	break;
      }
    }
  }
  bool empty() const{
    return Cell_base::empty();
  }
  bool contains(int j) const{
    for(const_iterator i=begin();i!=end();i++){
      if(*i==j){
	return true;
      }
    }
    return false;
  }
};

class Cells : public std::map<cellid_t,Cell>{ //public because to lazy
					      //to define gateways
 public:
  void insert(const NewSpecies & pt,int i)  __attribute__((always_inline))
  {
    (*this)[pt].insert(i);
  }
  void move(const NewSpecies & pt, int old_i, int new_i){
    std::map<cellid_t,Cell>::iterator i=find(pt);
    ALWAYS_ASSERT(i!=(this->end()));
    i->second.move(old_i,new_i);
  }
  void remove(const NewSpecies & pt,int j){
    std::map<cellid_t,Cell>::iterator i=find(pt);
    ALWAYS_ASSERT(i!=(this->end()));
    i->second.remove(j);
    if(i->second.empty()) erase(i);
  }
  int check_internal_consistency(const NewSpecies & pt, int i){
    std::map<cellid_t,Cell>::iterator loc=find(pt);
    ALWAYS_ASSERT(loc!=(this->end()));
    ALWAYS_ASSERT(loc->second.contains(i));
  };
  int check_internal_consistency();
};



#endif // _CELL_H_
