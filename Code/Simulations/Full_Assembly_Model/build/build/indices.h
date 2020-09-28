//$Id: indices.h 152 2006-01-16 10:02:54Z cvsrep $

#ifndef __INDICES_H__
#define __INDICES_H__

#include <vector>

typedef std::vector<int> v_int;
typedef std::vector<int> v_bool;

// allocate and deallocate indices, keeping the total index range as
// compact as possible

class index_space {
  v_int the_size;
  v_bool filled;
  int search_backward_for_beginning_of_free_range(int i);
  inline void set_filled(int start, int stop, bool value);
 public:
  // allocate a consecutive range of indices:
  int alloc(int length=1);
  // free a consecutive range
  void free(int position);
  int max();
  inline bool uses(int i){return filled[i];};
};

inline void index_space::set_filled(int start, int stop, bool value){
  for(int i=start;i<stop;i++)
    filled[i]=value;
}

      



#endif // __INDICES_H__