// -*- mode: c++ -*-
// $Id: packed_simulation.h 3 2005-12-01 07:13:32Z cvsrep $
#ifndef _PACKED_SIMULATION_H_
#define _PACKED_SIMULATION_H_

//////////////////////////////////////////
// experimental packed simulation state //
//////////////////////////////////////////

#include "NewWeb.h"

#ifndef __APPLE__  // Apples gcc does not understand INTEL's assembly code !?
#define TRY_ASSEMBLER_CODE
#endif

#define value_block_alignment (2*sizeof(double))
const int packed_simulation_block_size=8;
const int bs=packed_simulation_block_size; //abbreviation

class packed_simulation {

public:
  typedef double value_t;
private:
  typedef std::vector<int> partition_t;
  typedef value_t aligned_value_t 
  __attribute__ ((__aligned__(value_block_alignment)));

  typedef SortedMatrix::location_t  location_t; 
  typedef SortedVector::location_t  vector_index_t; 

  static partition_t chunkify(int size,int n_threads);


  inline char * aligned_pointer(char * ptr);
  inline aligned_value_t * aligned_value_pointer(char * ptr);
  inline size_t alignement_expand(size_t s);

  // Helper classes to assign locations in allocated memory chunk:
  template < typename T >
  class array_following {
    T * const the_start;
    const size_t the_size;
    array_following(char * mem, const size_t s):
      the_start((T *)mem),
      the_size(s){};
    template <typename U>
    array_following(U & u, const size_t s):
      the_start((T *)u.end()),
      the_size(s){}; 
    template <typename U>
    array_following(U & u, const int s):
      the_start((T *)u.end()),
      the_size(s){}; 
    template <typename U, typename V>
    array_following(U & u, V v):
      the_start((T *)u.end()),
      the_size(v.size()){
      copy(v.begin(),v.end(),the_start);
    };
    T & operator[](short int i) const {return *(the_start+i);}
    friend class packed_simulation;
  public:
    const T * begin() const {return the_start;}
    const T * end() const {return the_start + the_size;}
    size_t size() const {return the_size;}
  };

  template < typename T >
  class const_array_following {
    T * const the_start;
    const size_t the_size;
    template <typename U, typename V>
    const_array_following(U & u, V v):
      the_start((T *)u.end()),
      the_size(v.size()){
      copy(v.begin(),v.end(),the_start);
    };
    const T & operator[](short int i) const {return *(the_start+i);}
    friend class packed_simulation;
  public:
    const T * begin() const {return the_start;}
    const T * end() const {return the_start + the_size;}
    size_t size() const {return the_size;}
  };

  // Read data from chunc of memory:
  template< typename T > T postinc(char * __restrict__ & ptr);
  // Write data to memory chunc:
  template< typename T > void push(T val,char * & ptr);
  template< typename T > void push_value_t(T val,char * & ptr);
  
  inline value_t 
  fast_sandwich_product(const value_t * v1,const value_t * v2,value_t factor,
			char * __restrict__ & tsk);
  inline value_t 
  fast_dot_product(const value_t * v1,value_t factor,
		   char * __restrict__  &tsk);
  
  char * push_sorted_matrix(const SortedMatrix & m,char *& tsk);
  char * push_sorted_vector(const SortedVector & v,char *& tsk);

  class do_dynamics_manager_t{
    typedef packed_simulation base_class_t;
  public:
    do_dynamics_manager_t();
    do_dynamics_manager_t(const do_dynamics_manager_t & other);
    do_dynamics_manager_t& operator=(const do_dynamics_manager_t& other);
    ~do_dynamics_manager_t();
    void initialize_threads_maybe(base_class_t * base_class);
    void start_threads(int nthreads,base_class_t * base_class);
    void stop_threads();
    int get_num_threads(){return _num_threads;}
    boost::thread_group  threads;
    struct task_t {
      base_class_t * dispatcher; 
      int index;
      void operator()();
    };
    std::vector< task_t > task;
    boost::mutex mutex;
    boost::mutex B_max_mutex;
    boost::condition condition;
    boost::barrier * barrier;
    std::vector< value_t > common_factor_max;
    static const value_t unset; // some value < 0;
    bool max_not_set(int i){return common_factor_max[i]<0;}
    const ODE_vector * state;
    ODE_vector * time_derivative;
    bool stop_now;
    int threads_ready;
    std::vector< pthread_t > pthread_id;
  private:
    int _num_threads;
  };
  do_dynamics_manager_t do_dynamics_manager;

  inline bool I_am_last_of(int num_threads);

  //data:
  int arrival_counter;
  const int n_threads;
  const size_t memory_size;
  char * const memory;
  const NewWeb & web;
  double * turnover_rate;
  const int S;
  const int A;
  const value_t matrix_accuracy;
  const value_t vector_accuracy;
  const int plant_phys_vers;
public:
  value_t B_max;
  value_t common_factor_max;
  array_following<value_t> biomass_B;
  const_array_following<int> biomass_B_chunks;
  array_following<value_t> common_factor;
private:
  array_following<char *> eating_start;
  array_following<int> eating_chunks;
  array_following<char *> plant_start;
  array_following<int> plant_chunks;
  array_following<char *> being_eaten_start;  
  array_following<int> being_eaten_chunks;
  int data_end;

  inline size_t 
  required_size_of_memory(species_list_t & s, NewWeb::precomputed_t & pre);
  packed_simulation(const packed_simulation & other); //forbidden

public:
  packed_simulation(NewWeb & w);
  ~packed_simulation();
  int dynamics(ODE_vector const & state, 
	       ODE_vector & time_derivative);
  void do_dynamics(int t,   // thread id
		   ODE_vector const & state, 
		   ODE_vector & time_derivative);
  void do_dynamics1(int t,   // thread id
		    ODE_vector const & state, 
		    ODE_vector & time_derivative);
  void do_dynamics2(int t,   // thread id
		    ODE_vector const & state, 
		    ODE_vector & time_derivative);
  void do_dynamics3(int t,   // thread id
		    ODE_vector const & state, 
		    ODE_vector & time_derivative);
  void precondition(ODE_vector const & state,
		    ODE_vector const & in,
		    ODE_vector & out,
		    realtype gamma );
};

packed_simulation::aligned_value_t * packed_simulation::aligned_value_pointer(char * ptr){
  const size_t mask=value_block_alignment-1;
  return (aligned_value_t *) ((size_t(ptr)+mask) & ~mask);
}
char * packed_simulation::aligned_pointer(char * ptr){
  const size_t mask=value_block_alignment-1;
  return (char *) ((size_t(ptr)+mask) & ~mask);
}

size_t packed_simulation::alignement_expand(size_t s){
  const size_t mask=value_block_alignment-1;
  return (s+mask) & ~mask;
}

template< typename T >
T packed_simulation::postinc(char * __restrict__ & ptr){
  T val= *(T *) ptr;
  ptr+=sizeof(T);
  return val;
};

// Write data to memory chunk:
template< typename T >
void packed_simulation::push(T val,char * & ptr){
  *(T *) ptr=val;
  ptr+=sizeof(T);
};

// Write data to memory chunk:
template< typename T >
void packed_simulation::push_value_t(T val,char * & ptr){
  *(value_t *) ptr=val;
  ptr+=sizeof(value_t);
};
  

#endif // _PACKED_SIMULATION_H_
