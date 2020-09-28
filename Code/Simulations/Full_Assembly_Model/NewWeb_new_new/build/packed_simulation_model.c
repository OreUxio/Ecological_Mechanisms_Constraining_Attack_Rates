#include <assert.h>
#include <stdlib.h>

#define value_block_alignment (2*sizeof(double))

typedef double value_t;

int const bs=8; //int const packed_simulation_block_size;
typedef value_t aligned_value_t 
__attribute__ ((__aligned__(value_block_alignment)));


typedef long int size_T;

inline char * aligned_pointer(char * ptr){
  if(sizeof(size_T)!=sizeof(char *)) exit(-1);
  size_T mask=value_block_alignment-1;
  return (char *) (((size_T)(ptr+mask)) & ~mask);
}

inline
aligned_value_t * aligned_value_pointer(char * ptr){
  const size_T mask=value_block_alignment-1;
  return (aligned_value_t *) (((size_T)(ptr+mask)) & ~mask);
}

struct location_struct {
  short int row, column;
};
typedef struct location_struct location_t;

double asm_fast_sandwich_product(const double * v1,const double * v2,double factor,  char * __restrict__ tsk){
  double sum=0,last;
  double vec1[bs];
  char * matrix_end = *(char **) tsk;
  tsk += sizeof(char *);
  aligned_value_t * atsk=aligned_value_pointer(tsk);
  tsk = matrix_end;
  do{
    int i;
    for(i=0; i<bs; ++i){
      // Alternatively, try multiplying while reading here already;
      // multiplication then cannot be parallel, but it is streched in
      // time.
      vec1[i]=v1[((location_t *)(atsk))[i].row]*
    	v2[((location_t *)(atsk))[i].column];
    }
    atsk+=bs*sizeof(location_t)/sizeof(double);
    double * tt=atsk;
    for(i=0; i<bs; ++i){
      sum+=vec1[i] * *tt++;
    }
    atsk+=bs;
    last = *(((value_t *)atsk)-1);
  }while(last*factor>sum);
  return sum;
}
