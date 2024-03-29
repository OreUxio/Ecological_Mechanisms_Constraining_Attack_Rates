// -*- mode: c++ -*-
// $Id: index_set.h 429 2006-04-28 04:16:22Z cvsrep $
#ifndef _INDEX_SET_H_
#define _INDEX_SET_H_

// An index set is a set of indices.  Add, remove, and move take such
// indices as arguments and are all O(1) operations.  The set of
// indices is represented by a list of indices that is efficiently
// traversed in vectorized loops.  It is the responsability of the
// user that operations on index sets are logically consistent
// (e.g. not double add of the same index).

#include "sequence.h"
typedef sequence<int> Index_Set_Base;

class Index_Set : private Index_Set_Base
{
  sequence<int> the_position_of;
 public:
  Index_Set(){};
  ~Index_Set(){};
  void add(int i){
    int s=size();
    resize(s+1);
    Index_Set_Base::operator()(s)=i;
    the_position_of[i]=s;  // automatic allocation if required
    return;
  }
  void remove(int i){
    int s=size();
    int last=Index_Set_Base::operator()(s-1);
    Index_Set_Base::operator()(the_position_of(i))=last;
    the_position_of(last)=the_position_of(i);
    resize(s-1);
    return;
  }
  void move(int old_i,int new_i){
    int pos=the_position_of(old_i);
    Index_Set_Base::operator()(pos)=new_i;
    the_position_of[new_i]=pos;
    return;
  };
  bool internally_consistent(){
    bool ok=true;
    for(int i=size();i-->0;){
      if(i!=the_position_of((*this)(i)))
	ok=false;
    }
    return ok;
  }
  inline int operator[](int i){
    return Index_Set_Base::operator()(i);
  }
  inline int size(){
    return Index_Set_Base::size();
  }
  inline void clear(){
    resize(0);
  }
};
#endif // _INDEX_SET_H_
