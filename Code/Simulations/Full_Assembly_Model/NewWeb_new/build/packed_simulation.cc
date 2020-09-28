// -*- mode: c++ -*-
// $Id: packed_simulation.cc 3 2005-12-01 07:13:32Z cvsrep $

#include "packed_simulation.h"
#include "packed_simulation_asm.h"

// stolen and modified from kernel header arch/x86/include/asm/system.h
static inline void clflush(volatile char *p)
{
  asm volatile("clflush %0" : "+m" (*p));
}
// stolen from kernel header arch/x86/include/asm/system.h
#define mb()    asm volatile("mfence":::"memory")
// stolen and modified from linux kernel code arch/x86/mm/pageattr.c
/**
 * clflush_cache_range - flush a cache range with clflush
 *
 * clflush is an unordered instruction which needs fencing with mfence
 * to avoid ordering issues.
 */
const int cache_line_size=64; //should be adjustable
void clflush_cache_range(char *ptr, char * end)
{
  end--;
 
  mb();
	 
  for (; ptr < end; ptr += cache_line_size)
    clflush(ptr);
  /*
   * Flush any possible final partial cacheline:
   */
  clflush(end);
  
  mb();
}

//////////////////////////////////////////
// Experimental interlocking scope.     //
// Use to sync threads                  //
//////////////////////////////////////////

#include <boost/thread.hpp>
#include "error.h" // for debugging
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

#ifdef _GNU_SOURCE 
#define USE_SPINLOCK
#endif
#ifdef USE_SPINLOCK
extern pthread_spinlock_t spinlock;
#else
extern pthread_mutex_t spinlock;
#endif

static volatile int interlocking_scope_counter[3]={0,0,0};
static volatile int interlocking_scope_current=0;

class interlocking_scope {
  inline volatile int & _current_counter(){
    return interlocking_scope_counter[interlocking_scope_current];
  }
  inline void _switch_counter(){
    interlocking_scope_current+=1;
    if(interlocking_scope_current==3) interlocking_scope_current=0;
  }

  bool _I_am_last;
  int _my_counter;

public:
  interlocking_scope (const std::vector< pthread_t > & pthread_id):
    _my_counter(interlocking_scope_current)
  {
    int value_now=__sync_add_and_fetch(interlocking_scope_counter+_my_counter,1);
    
    _I_am_last=(value_now == pthread_id.size() );

    if(_I_am_last){
      _switch_counter();
      interlocking_scope_counter[_my_counter]=0;  
      // From here other threads potentially start running.
    }
  }
  ~interlocking_scope(){
    if (interlocking_scope_counter[_my_counter] > 0) {
      int b = 100;
      do {
    	for (int i = b; i; i--) 
    	  asm volatile ("");
    	b =  int(b * 1.01);
    	// spin with reads
      } while (interlocking_scope_counter[_my_counter]>0);
    }
  }

  bool I_am_last() const {return _I_am_last;}
};
    
//////////////////////////////////////////
// experimental packed simulation state //
//////////////////////////////////////////

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

const packed_simulation::value_t 
log_DBL_MAX=5*log(DBL_MAX)/7.0; ///< we added some safety margin
const packed_simulation::
value_t log_DBL_MIN=5*log(DBL_MIN)/7.0; ///< we added some safety margin

#ifdef USE_SPINLOCK
pthread_spinlock_t spinlock;
#else
pthread_mutex_t spinlock;
#endif

	
const packed_simulation::
value_t packed_simulation::do_dynamics_manager_t::unset=-1;

packed_simulation::do_dynamics_manager_t::do_dynamics_manager_t():
  _num_threads(1),
  threads_ready(0),
  common_factor_max(1,0),
  barrier(new boost::barrier(1)),
  pthread_id(1)
{
  pthread_id[0]=pthread_self();
}

/// Don't copy do_dynamics_manager_t when copying packed_simulation:
packed_simulation::do_dynamics_manager_t::
do_dynamics_manager_t(const do_dynamics_manager_t & other):
  _num_threads(1),
  threads_ready(0),
  common_factor_max(1,0)
{
}

/// Don't assign to do_dynamics_manager_t when assigning to packed_simulation:
packed_simulation::do_dynamics_manager_t& packed_simulation::do_dynamics_manager_t::
operator=(const packed_simulation::do_dynamics_manager_t& other){
}

void packed_simulation::do_dynamics_manager_t::
start_threads(int n,packed_simulation * base_web){
  if(n!=_num_threads){
    stop_threads();
    _num_threads=n;
    task.resize(_num_threads-1);
    pthread_id.resize(_num_threads);
    delete barrier;
    barrier=new boost::barrier(_num_threads);
    common_factor_max.resize(_num_threads);
    stop_now=false;
    for(int i=_num_threads-1;i-->0;){
      task[i].dispatcher=base_web;
      task[i].index=i+1;
      threads.create_thread(task[i]);
    }
    // while(threads_ready<_num_threads-1){
    //   boost::this_thread::yield();
    // }
  }
}

void packed_simulation::do_dynamics_manager_t::stop_threads(){
  if(_num_threads>1){
    stop_now=true;

    {interlocking_scope main_stop(pthread_id);
    }
    // {boost::mutex::scoped_lock lock(mutex);
    //   condition.notify_all();
    // }
    threads.join_all();

    _num_threads=1;

    delete barrier;
    barrier=new boost::barrier(_num_threads);
    pthread_id.resize(_num_threads);
    common_factor_max.resize(_num_threads);
    task.resize(_num_threads-1);
  }
}

packed_simulation::do_dynamics_manager_t::~do_dynamics_manager_t(){
  stop_threads();
  delete barrier;
}

inline bool packed_simulation::
I_am_last_of(int num_threads){
  int arrival_number=__sync_add_and_fetch(&arrival_counter,1);
  bool am_last=(arrival_number == num_threads );
  if(am_last){
    arrival_counter=0;
  }
  return am_last;
  //return false;
}


inline packed_simulation::value_t packed_simulation::
fast_sandwich_product(const value_t * __restrict__ v1,const value_t * __restrict__ v2,value_t factor,
		      char * __restrict__ & tsk){
#ifdef TRY_ASSEMBLER_CODE
  value_t result=asm_fast_sandwich_product(v1,v2,factor,tsk);
  tsk = *(char **)tsk;
  return result;
#else
  value_t sum=0,last;
  value_t vec1[bs];
  char * matrix_end = postinc<char *>(tsk);
  aligned_value_t * atsk=aligned_value_pointer(tsk);
  tsk = matrix_end;
  do{
    for(int i=0; i<bs; ++i){
      // Alternatively, try multiplying while reading here already;
      // multiplication then cannot be parallel, but it is streched in
      // time.
      vec1[i]=v1[((location_t *)(atsk))[i].row]*
    	v2[((location_t *)(atsk))[i].column];
    }
    atsk+=bs*sizeof(location_t)/sizeof(value_t);
    value_t * tt=atsk;
    for(int i=0; i<bs; ++i){
      sum+=vec1[i] * *tt++;
    }
    atsk+=bs;
    last = *(((value_t *)atsk)-1);
  }while(last*factor>sum);
  return sum;
#endif
}

   

inline packed_simulation::value_t packed_simulation::
fast_dot_product(const value_t * v1,value_t factor,
		 char * __restrict__  &tsk){
  value_t sum=0,last;
  value_t vec1[bs];
  char * vector_end = postinc<char *>(tsk);
  aligned_value_t * atsk=aligned_value_pointer(tsk);
  tsk = vector_end;
  do{
    for(int i=0; i<bs; ++i){
      vec1[i]=v1[((vector_index_t *)(atsk))[i].row];
    }
    atsk+=bs*sizeof(vector_index_t)/sizeof(value_t);
    value_t * tt=atsk;
    for(int i=0; i<bs; ++i){
      sum+=vec1[i]* *tt++;
    }
    atsk+=bs;
    last = *(((value_t *)atsk)-1);
  }while(last*factor>sum);
  return sum;
}
    
char * packed_simulation::push_sorted_matrix(const SortedMatrix & m,char *& tsk){
  SortedMatrix::_container & e=m._entries;

  // matrix start:
  char * matrix_start=tsk;
  // matrix end:
  char * matrix_end_loc=tsk;
  push(aligned_pointer(tsk+sizeof(char *))+
       bs*(1+e.size()/bs)*
       (sizeof(SortedMatrix::location_t)+
	sizeof(value_t) ),
       tsk );
    
  tsk=aligned_pointer(tsk);
  for(int i=0; i<1+e.size(); i+=bs){
    int j;
    for(j=0;j<bs;++j){
      if(i+j<e.size()){
	push(SortedMatrix::location_t(e[i+j]),tsk);
      }else{
	push(SortedMatrix::location_t(0,0),tsk);
      }
    }
    for(j=0;j<bs;++j){
      if(i+j<e.size()){
	push_value_t(e[i+j].value,tsk);
      }else{
	push_value_t(0,tsk);
      }
    }
  }
  ASSERT(tsk==*(char **)matrix_end_loc);
  return matrix_start;
}
  
char * packed_simulation::
push_sorted_vector(const SortedVector & v,char *& tsk){
  SortedVector::_container & e=v._entries;

  // vector end:
  char * vector_start=tsk;

  char * vector_end_loc=tsk;
  push(aligned_pointer(tsk+sizeof(char *))+
       bs*(1+e.size()/bs)*
       (sizeof(SortedVector::location_t)+
	sizeof(value_t) ),
       tsk );
    
  tsk=aligned_pointer(tsk);
  for(int i=0; i<1+e.size(); i+=bs){
    int j;
    for(j=0;j<bs;++j){
      if(i+j<e.size()){
	push(SortedVector::location_t(e[i+j]),tsk);
      }else{
	push(SortedVector::location_t(0),tsk);
      }
    }
    for(j=0;j<bs;++j){
      if(i+j<e.size()){
	push_value_t(e[i+j].value,tsk);
      }else{
	push_value_t(0,tsk);
      }
    }
  }
  ASSERT(tsk==*(char **)vector_end_loc);
  return vector_start;
}




size_t packed_simulation::
required_size_of_memory(species_list_t & s, NewWeb::precomputed_t & pre){
  const size_t S=s.size();
  const size_t A=s.n_animals;
  size_t sum=0;
  sum+=sizeof(packed_simulation);
  sum+=S*sizeof(value_t); // biomass_B
  sum+=(n_threads+1)*sizeof(int); // biomass_B_chunks
  sum+=A*sizeof(value_t); // precomputed_t
  sum+=n_threads*sizeof(int *); // task for each thread

  for(int k=0;k<3*n_threads;++k){
    sum+=sizeof(int);  // for the starting index
  }
  sum+=A*sizeof(char *); // eating starts
  sum+=(S-A)*sizeof(char *); // plant starts
  sum+=S*sizeof(char *); // being eaten starts

  sum=alignement_expand(sum);

  for(int k=0;k<A;++k){
    sum=alignement_expand(sum+sizeof(char *));
    pre[k].csc_eating.cleanup();
    sum+=bs * ( 1 + (pre[k].csc_eating._entries.size()) / bs ) * 
      (sizeof(SortedMatrix::location_t)+sizeof(value_t));
    sum=alignement_expand(sum+sizeof(char *));
    pre[k].c.cleanup();
    sum+=bs * ( 1 + (pre[k].c._entries.size()) / bs ) * 
      (sizeof(SortedVector::location_t)+sizeof(value_t));
    sum+=4*sizeof(value_t);
  }

  for(int k=A;k<S;++k){
    sum=alignement_expand(sum+sizeof(char *));
    pre[k].c.cleanup();
    sum+=bs * ( 1 + (pre[k].c._entries.size()) / bs ) * 
      (sizeof(SortedVector::location_t)+sizeof(value_t));
    sum+=sizeof(value_t);  // for plant_growth_rate_sigma
    sum+=sizeof(value_t);  // for second parameter
  }

  for(int k=0;k<S;++k){
    sum=alignement_expand(sum+sizeof(char *));
    pre[k].csc_being_eaten.cleanup();
    sum+=bs * ( 1 + (pre[k].csc_being_eaten._entries.size()) / bs ) * 
      (sizeof(SortedMatrix::location_t)+sizeof(value_t));
  }
  
  sum+=S*sizeof(double); // for turnover rates

  sum+=25; // The computation above is not perfectly precise, because
	   // we do not do all the alignments at the right places. (??)

#ifdef DEBUGGING
  size_t size_of_packed_memory_in_k=sum/value_t(1<<10);
  REPORT(size_of_packed_memory_in_k);
#endif

  return sum;
}

packed_simulation::partition_t 
packed_simulation::chunkify(int size,int n_threads){
  partition_t ch(n_threads+1,0);
  int k=0;
  int small_chunk_size=size/n_threads;
  int number_of_large=size % n_threads;
  for(int i=0;i<=n_threads;++i){
    ch[i]=k;
    k+=small_chunk_size;
    if( i >= n_threads-number_of_large){
      k+=1;
    }
  }
  //for(int i=0;i<n_threads;++i) ch[i]=i;
  return ch;
};

void handle_alrm(int sig) {
}

packed_simulation::
packed_simulation(NewWeb & w):
  arrival_counter(0),
  n_threads(max_num_threads),
  memory_size(required_size_of_memory(w.s,w.precomputed)),
  memory((char *)malloc(memory_size)),
  web(w),
  S(web.s.size()),
  A(web.s.n_animals),
  matrix_accuracy(web.precomputed[0].csc_eating.accuracy()),
  vector_accuracy(web.precomputed[0].c.accuracy()),
  plant_phys_vers(plant_physiology_version),
  biomass_B(aligned_pointer(memory+sizeof(packed_simulation)),S),
  biomass_B_chunks(biomass_B,chunkify(S,n_threads)),
  common_factor(biomass_B_chunks,A),
  eating_start(common_factor,A),
  eating_chunks(eating_start,chunkify(A,n_threads)),
  plant_start(eating_chunks,S-A),
  plant_chunks(plant_start,chunkify(S-A,n_threads)),
  being_eaten_start(plant_chunks,S),
  being_eaten_chunks(being_eaten_start,chunkify(S,n_threads))
{
  if(!memory){
    FATAL_ERROR("could not allocate memory for packed simulation state");
  }

  signal(SIGALRM, handle_alrm);

  ALWAYS_ASSERT(A>=n_threads);
  ALWAYS_ASSERT(S-A>=n_threads);
  ALWAYS_ASSERT(do_switching);
    
  char * tsk=(char *)being_eaten_chunks.end();
        
  for(int k=0;k<A;k++){

    ALWAYS_ASSERT(aligned_pointer((char *)(bs*sizeof(SortedMatrix::location_t)))==
		  (char *)(bs*sizeof(SortedMatrix::location_t)));
    eating_start[k]=
      push_sorted_matrix(web.precomputed[k].csc_eating,tsk);
    push_sorted_vector(web.precomputed[k].c,tsk);
    
    push_value_t(web.s(k).handling_time_T()*web.s(k).attack_rate_a(),tsk);
    push_value_t(web.s(k).attack_rate_a(),tsk);
    push_value_t(web.s(k).conversion_efficiency_eps(),tsk);
    push_value_t(web.s(k).turnover_rate_r(),tsk);
  }

  for(int k=0;k<S-A;k++){
    plant_start[k]=
      push_sorted_vector(web.precomputed[k+A].c,tsk);
    push_value_t(web.s(k+A).plant_growth_rate_sigma(),tsk);
    if(plant_phys_vers==3){
      push_value_t(web.s(k+A).env_effect(),tsk);
    }else{
      push_value_t(web.s(k+A).get_loss_rate_over_max_production_rate_r(),tsk);
    }
  }

  for(int k=0;k < S;k++){
    being_eaten_start[k]=
      push_sorted_matrix(web.precomputed[k].csc_being_eaten,tsk);
  }

  tsk=aligned_pointer(tsk);

  turnover_rate=(double *)tsk;
  for(int i=0;i<S;++i){
    push(web.s(i).turnover_rate_r(),tsk);
  }
  
  if((tsk-memory)>=memory_size){
    ALWAYS_ASSERT((tsk-memory)>=memory_size);
  }

  memcpy(memory,this,sizeof(this));

#ifdef USE_SPINLOCK
  pthread_spin_init(&spinlock, 0);
#else
  pthread_mutex_init(&spinlock, NULL);
#endif
#ifndef USE_OPENMP
  do_dynamics_manager.start_threads(max_num_threads,this);
#endif
#ifdef SET_CPU_AFFINITY
  sched_setaffinity(0,CPU_SETSIZE, &CPU_set_list[0] );
#endif
};
  
packed_simulation::~packed_simulation(){
#ifndef USE_OPENMP
  do_dynamics_manager.stop_threads();
#endif
  free((void *)memory);
#ifdef USE_SPINLOCK
  pthread_spin_destroy(&spinlock);
#else
  pthread_mutex_destroy(&spinlock);
#endif
};
  

/// Main function for spawned threads:
void packed_simulation::do_dynamics_manager_t::task_t::operator()(){

  packed_simulation::do_dynamics_manager_t & mgr=
    dispatcher->do_dynamics_manager;
  
#ifdef SET_CPU_AFFINITY
  sched_setaffinity(0,CPU_SETSIZE, &CPU_set_list[index] );
#endif
  mgr.pthread_id[index]=pthread_self();

  while(true){
    
    // {boost::mutex::scoped_lock lock(mgr.mutex);
    //   mgr.common_factor_max[index]=do_dynamics_manager_t::unset;  

    //   mgr.threads_ready++;
    //   mgr.condition.wait(mgr.mutex);
    //   mgr.threads_ready--;
    // }
    {interlocking_scope main_stop(mgr.pthread_id);
    }
    
    if(mgr.stop_now){
      return;
    }else{
      dispatcher->do_dynamics(index,
			      *mgr.state,
			      *mgr.time_derivative);
    }
  }
}  

int next_e;
int next_p;
int next_be;

int packed_simulation::dynamics(ODE_vector const & state, 
				ODE_vector & time_derivative){

#ifndef USE_OPENMP

  do_dynamics_manager_t & mgr = do_dynamics_manager;

  // Pass data to other threads:
  mgr.state=&state;
  mgr.time_derivative=&time_derivative;

  B_max=0;
  common_factor_max=0;
  mgr.common_factor_max[0] = do_dynamics_manager_t::unset;

  // Run parallel computation:
  {interlocking_scope main_stop(mgr.pthread_id);
  }
  // {boost::mutex::scoped_lock lock(mgr.mutex);
  //   mgr.condition.notify_all();
  // }
  do_dynamics(0,state,time_derivative);


  if(newton_failure){
    const int i=newton_failure;
    WARNING("state[" << i << "] = " << state[i]);
    REPORT(web.assigned_column[i]);
    REPORT(web.s(i).turnover_rate_r());
    REPORT(web.s(i).get_trade_off_multiplier());
    REPORT(web.s(i).plant_growth_rate_sigma());
    newton_failure=true;
    return 1; // recoverable failure calculating dynamics
  }

  // Wait until all threads are ready to start the next iteration.
  // while(mgr.threads_ready<
  // 	mgr.get_num_threads()-1){
  //   boost::this_thread::yield();
  // }

  return 0; // success calculating dynamics
}


void packed_simulation::do_dynamics(int t, // thread number
				    ODE_vector const & state, 
				    ODE_vector & time_derivative){

  //   do_dynamics1(t,state,time_derivative);
  //   do_dynamics2(t,state,time_derivative);
  //   do_dynamics3(t,state,time_derivative);
  // }
  // void packed_simulation::do_dynamics1(int t, // thread number
  // 				    ODE_vector const & state, 
  // 				    ODE_vector & time_derivative){
  const int num_threads=do_dynamics_manager.get_num_threads();
  
#else // use OPENMP:
#pragma omp parallel num_threads(n_threads)
  {
    int t=omp_get_thread_num();
#ifdef SET_CPU_AFFINITY
    sched_setaffinity(0,CPU_SETSIZE, &CPU_set_list[t] );
#endif
    const int num_threads=n_threads;
#if 0
  } //balance superfluous brace
#endif
#endif // use OPENMP

  // Maximum of biomass_B computed in this thread.
  value_t B_thread_max = 0;
  
#ifdef USE_OPENMP
#pragma omp for nowait
  for(short int i=0; i<S ; ++i){
#if 0
  } //balance superfluous brace
#endif
#else
  for(short int i=biomass_B_chunks[t]; i<biomass_B_chunks[t+1]; ++i){
#endif
    if(state[i] > log_DBL_MAX/2 /*the '2' here is the assumed
				  switching exponent!*/  ||
       state[i] < log_DBL_MIN/2){
      newton_failure=i;
    }else{
      // In principle, these exps could also be computed in parallel,
      // but the cost of meeting at a barrier afterwards is too
      // high (and it is linear in problem size!).
      value_t B=exp(state[i]);
      B_thread_max = std::max(B, B_thread_max);
      biomass_B[i]=B;
    }
  }
  
#ifdef USE_OPENMP
#pragma omp critical 
#endif
  {
#ifndef USE_OPENMP
    boost::mutex::scoped_lock lock(do_dynamics_manager.B_max_mutex);
#endif
    B_max = std::max(B_max, B_thread_max);
  }
#ifdef USE_OPENMP
#pragma omp barrier
#else
  { interlocking_scope wait_B_max(do_dynamics_manager.pthread_id);
  }
#endif
  //do_dynamics_manager.barrier->wait();

  // }
  // void packed_simulation::do_dynamics2(int t, // thread number
  // 				    ODE_vector const & state, 
  // 				    ODE_vector & time_derivative){
  //   const int num_threads=do_dynamics_manager.get_num_threads();

  value_t common_factor_thread_max=0;
  {
    int k=eating_chunks[t];
    int & k_end=eating_chunks[(t+1)%num_threads];
    
    char * tsk=(char *)eating_start[k];
    
    const value_t sand_factor=B_max*B_max/matrix_accuracy;
    const value_t halfsat_factor=B_max/vector_accuracy;
    const value_t shadowing_factor=B_max/vector_accuracy;
    
    do{
      //
      value_t sand=fast_sandwich_product(biomass_B.begin(),biomass_B.begin(),sand_factor,tsk);
      value_t halfsat=fast_dot_product(biomass_B.begin(),halfsat_factor,tsk);
      //
      const value_t handling_time_X_attack_rate=postinc<value_t>(tsk);
      const value_t attack_rate=postinc<value_t>(tsk);
      const value_t conversion_efficiency=postinc<value_t>(tsk);
      const value_t turnover_rate=postinc<value_t>(tsk);
      const value_t fact0=
	( halfsat>0 ? 1/(halfsat+handling_time_X_attack_rate*sand) : 0 );
      const value_t cf=biomass_B[k]*(attack_rate*fact0);
      // REPORT(sand);
      // REPORT(cf);
      common_factor[k]=cf;
      common_factor_thread_max=std::max(common_factor_thread_max,cf);
      time_derivative[k]=
	conversion_efficiency*sand*(attack_rate*fact0)
	-turnover_rate;	
      //
      k++;
      if(k==A){
	k=0;
	tsk=eating_start[0];
      }
    }while(k!=k_end);// end while eating
#ifdef USE_OPENMP
#pragma omp critical 
    {
      common_factor_max = 
	std::max(common_factor_max, common_factor_thread_max);
    }
#else
    if(t & 1){
      //boost::mutex::scoped_lock lock(do_dynamics_manager.B_max_mutex);
#ifdef USE_SPINLOCK
      pthread_spin_lock(&spinlock);
#else
      pthread_mutex_lock(&spinlock);
#endif
      common_factor_max = std::max(common_factor_max, common_factor_thread_max);
#ifdef USE_SPINLOCK
      pthread_spin_unlock(&spinlock);
#else
      pthread_mutex_unlock(&spinlock);
#endif
    }
#endif // not USE_OPENMP
  } // eating block
 
  {
    // Plant population growth:
    const value_t shadowing_factor=B_max/vector_accuracy;
    int k=plant_chunks[t];
    int k_end=plant_chunks[(t+1)%num_threads];
    char * tsk=(char *)plant_start[k];

    do{
      //
      value_t shadowing=
	fast_dot_product(biomass_B.begin(),shadowing_factor,tsk);
      // REPORT(shadowing);
      //
      // saturation effects:
      const value_t plant_growth_rate = postinc<value_t>(tsk);
      // REPORT(plant_growth_rate);
      if(plant_phys_vers==3){
	const value_t env_eff = postinc<value_t>(tsk);
	time_derivative[A+k]=plant_growth_rate*(1-shadowing+env_eff);
      }else{
	const value_t loss_rate_over_max_production = postinc<value_t>(tsk);
	time_derivative[A+k]=plant_growth_rate*
	  (exp(-shadowing)-loss_rate_over_max_production);
      }
      k++;
      if(k==(S-A)){
	k=0;
	tsk=plant_start[0];
      }
    }while(k!=k_end);// end while population growth

  } // plant block

#ifndef USE_OPENMP
  // Now even threads declare their maximum:
  if((t & 1)==0){
    //boost::mutex::scoped_lock lock(do_dynamics_manager.B_max_mutex);
#ifdef USE_SPINLOCK
    pthread_spin_lock(&spinlock);
#else
    pthread_mutex_lock(&spinlock);
#endif
    common_factor_max = std::max(common_factor_max, common_factor_thread_max);
#ifdef USE_SPINLOCK
    pthread_spin_unlock(&spinlock);
#else
    pthread_mutex_unlock(&spinlock);
#endif
    
  }
#endif
  
#ifdef USE_OPENMP
#pragma omp barrier
#else
  { interlocking_scope set_common_factor_max(do_dynamics_manager.pthread_id);
    if(set_common_factor_max.I_am_last()){
      int next_t=(t+1)%num_threads;
      int up=(plant_chunks[t]+1)%(S-A);
      if(up!=plant_chunks[next_t]){
	plant_chunks[t]=up;
	// std::cout << "P ";
	// for(int i=0;i<num_threads;i++){
	//   std::cout << (-plant_chunks[i]+plant_chunks[(i+1)%num_threads]+S-A)%(S-A)  << " ";
	// }
	// std::cout << std::endl;
      }else{
	int up2=(eating_chunks[t]+1)%A;
	if(up2!=eating_chunks[next_t]){
	  eating_chunks[t]=up2;
	}
	// std::cout << "E ";
	// for(int i=0;i<num_threads;i++){
	//   std::cout << (-eating_chunks[i]+eating_chunks[(i+1)%num_threads]+A)%A  << " ";
	// }
	// std::cout << std::endl;
      }
    }
  }// End of interlocking_scope for setting common_factor_max
#endif // not USE_OPENMP


  // }
  // void packed_simulation::do_dynamics3(int t, // thread number
  // 				    ODE_vector const & state, 
  // 				    ODE_vector & time_derivative){
  //   const int num_threads=do_dynamics_manager.get_num_threads();


  // Predation/consumption mortality for all species:
  {
    // Consumption mortality:
    int k=being_eaten_chunks[t];
    int k_end=being_eaten_chunks[(t+1)%num_threads];
    char * tsk=(char *)being_eaten_start[k];
    const value_t being_eaten_factor=
      B_max*common_factor_max/matrix_accuracy;

    do{
      value_t being_eaten=
	fast_sandwich_product(biomass_B.begin(),
			      common_factor.begin(),
			      being_eaten_factor,tsk);
      // REPORT(being_eaten);
      time_derivative[k]-=being_eaten;
      k++;
      if(k==S){
	k=0;
	tsk=being_eaten_start[0];
      }
    }while(k!=k_end); // end while being eaten
  }
  // being eaten block
    
#ifndef USE_OPENMP
  {interlocking_scope exit_barrier(do_dynamics_manager.pthread_id);
    if(exit_barrier.I_am_last()){
      int next_t=(t+1)%num_threads;
      int up=(being_eaten_chunks[t]+1)%S;
      if(up!=being_eaten_chunks[next_t]){
	being_eaten_chunks[t]=up;
      }
      // for(int i=0;i<num_threads;i++){
      // 	std::cout << (-being_eaten_chunks[i]+being_eaten_chunks[(i+1)%num_threads]+S)%S  << " ";
      // }
      // std::cout << std::endl;
    }
  } 
#endif

  //clflush_cache_range(memory,memory+memory_size);
  
#ifdef USE_OPENMP
#if 0 
  { // balance superfluous brace
#endif
  } // omp parallel
  if(newton_failure){
    const int i=newton_failure;
    WARNING("state[" << i << "] = " << state[i]);
    REPORT(web.assigned_column[i]);
    REPORT(web.s(i).turnover_rate_r());
    REPORT(web.s(i).get_trade_off_multiplier());
    REPORT(web.s(i).plant_growth_rate_sigma());
    newton_failure=true;
    return 1; // recoverable failure calculating dynamics
  }
  return 0;
#endif
}

#ifdef USE_OPENMP
void packed_simulation::do_dynamics(int t, // thread number
				    ODE_vector const & state, 
				    ODE_vector & time_derivative){
  // Dummy.  This function is not used with Open MP.
}
#endif

void packed_simulation::
precondition(ODE_vector const & state,
	     ODE_vector const & in,
	     ODE_vector & out,
	     realtype gamma ){
#ifdef USE_OPENMP
#pragma omp parallel num_threads(n_threads)
  {
    int t = omp_get_thread_num();
    sched_setaffinity(0,CPU_SETSIZE, &CPU_set_list[t] );
    short int i_end=biomass_B_chunks[t+1];
    for(short int i=biomass_B_chunks[t]; i<i_end; ++i){
      out[i]=in[i]/(1+gamma*turnover_rate[i]);
    }
  }
#else
  for(int i=0;i<S;++i){
    out[i]=in[i]/(1+gamma*turnover_rate[i]);
  }
#endif
}

