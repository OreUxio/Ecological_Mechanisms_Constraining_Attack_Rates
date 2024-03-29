// $Id: indices.cc 12 2005-12-03 04:22:02Z cvsrep $

#include "indices.h"
#include "error.h"


int index_space::search_backward_for_beginning_of_free_range(int i){
  ASSERT(0<=i && i<=(signed int) the_size.size());
  if(i==0) return 0;
  int j=i;
  while(--j >= 0 && the_size[j]==0);
  if(filled[j])
    return i;
  else
    return j;
}

int index_space::alloc(int length){
  ASSERT(length > 0);
  for(unsigned int i=0;i<the_size.size();i+=the_size[i]){
    if(!filled[i] && the_size[i]>=length){
      // we found a usable pice of index space, now reserve it for us:
      set_filled(i,i+length,true);
      if(the_size[i]>length){
	filled[i+length]=false;
	the_size[i+length]=the_size[i]-length;
	the_size[i]=length;
      }
      return i;
    }
  }
  // we did not find sufficient space, so we have to extend the index
  // space:
  int startpoint=
    search_backward_for_beginning_of_free_range(the_size.size());
  the_size.resize(startpoint+length);
  filled.resize(startpoint+length);
  the_size[startpoint]=length;
  set_filled(startpoint,startpoint+length,true);
  return startpoint;
}

void index_space::free(int i){
  if(!filled[i] || the_size[i]<=0)
    FATAL_ERROR("Cannot free this index.");
  set_filled(i,i+the_size[i],false);
  int startpoint=
    search_backward_for_beginning_of_free_range(i);
  int full_length=
    i-startpoint+the_size[i];
  if(i+the_size[i] < (signed int)the_size.size() && !filled[i+the_size[i]]){
    // the block block behind was also free:
    full_length+=the_size[i+the_size[i]];
    the_size[i+the_size[i]]=0;
  }
  the_size[i]=0;
  the_size[startpoint]=full_length;
  return;
}

int index_space::max(){
  return
    search_backward_for_beginning_of_free_range(the_size.size());
}
