// $Id: props_t.cc 888 2007-06-22 15:23:00Z cvsrep $

#include <time.h>
#include "random.h"
#include "sequence.h"
#include "props_t.h"
#include "error.h"

// adjustable parameters:
static int force_dull_mode_selection=1;
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGINT(force_dull_mode_selection),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);


using namespace std;


const int props_t::nwords=props_t_::nwords;
const int props_t::accuracy_in_bits[]=
  {props_t_::wbits,props_t_::wbits,8,10,12};
/* at accuracy 12 parallel can still compete for some p */
const int props_t::precision_goal_in_bits=7;


void props_t::set_random(){
  for(int i=0;i<nwords;i++){
    word[i]=random_int();
  }
}

void props_t::swap_a_fraction(const double p_swap, 
			      const swapping_mode mode){
  switch(mode){
  case logarithmic:
    {
      static double p_hold=-1;
      static double logCp;
      if(p_swap!=p_hold){
	p_hold=p_swap;
	logCp=log(1-p_swap);
      }
      double n_logCp=n_property_bits*logCp;
      double log_rand=log(unirand());
      if(__builtin_expect(log_rand<n_logCp,0))
	return;
      int bi=int(log_rand/logCp);
      while(__builtin_expect(bi<n_property_bits,1)){
	word[bi/props_t_::wbits]^=(1<<(bi%props_t_::wbits));
	log_rand=log(unirand());
	if(__builtin_expect(log_rand<n_logCp,0))
	  return;
	bi++;
	bi+=int(log_rand/logCp);
      }
    }
    return;
  case multiplying:
    {
      //p_swap is small:
      const double p=1-p_swap;
      register double random=unirand();
      register double threshold=p;
      for(int i=0;i<props_t_::nwords;i++){
	register unsigned int w=word[i];
	for(register unsigned int b=1;
	    __builtin_expect (b!=0, 1);
	    b<<=1){
	  if(__builtin_expect (random<threshold,1)){
	    threshold*=p;
	  }else{
	    w^=b;
	    threshold=p;
	    random=unirand();
	  }
	}
	word[i]=w;
      }
    }
    return;
#ifdef PARALLEL
#undef PARALLEL
#endif
#define PARALLEL(N)						\
  case parallel##N:						\
    {								\
      if(p_swap==0.5){						\
	set_random();						\
	return;							\
      }								\
      const int bit_accuracy=accuracy_in_bits[parallel##N];	\
      unsigned int prob_bits=unsigned(4*p_swap*(1<<(props_t_::wbits-2)));	\
      prob_bits&=((-1)^((1<<(props_t_::wbits-bit_accuracy))-1));		\
      if(prob_bits){						\
	/*normal case*/						\
	for(int i=0;__builtin_expect(i<nwords,1);i++){		\
	  unsigned int swap_these=0;				\
	  unsigned int b=(unsigned(1)<<(props_t_::wbits-bit_accuracy));	\
	  while(!(b&prob_bits)) b<<=1;				\
	  while(__builtin_expect(b,1)){				\
	    if(b&prob_bits){					\
	      swap_these|=random_int();				\
	    }else{						\
	      swap_these&=random_int();				\
	    }							\
	    b<<=1;						\
	  }							\
	  word[i]^=swap_these;					\
	}							\
      }else{							\
	/*swap all (p_swap~1)*/					\
	for(int i=0;i<nwords;i++){				\
	  word[i]^=unsigned(-1);				\
	}							\
      }								\
    }								\
    return;
    PARALLEL(1);
    PARALLEL(2);
    PARALLEL(3);
#undef PARALLEL
  }
}

props_t::swapping_mode props_t::fastest_mode(double p){
  {
#ifdef CLOCK_PROCESS_CPUTIME_ID
    timespec tp={0,0};
    if(force_dull_mode_selection ||
       0!=clock_settime(CLOCK_PROCESS_CPUTIME_ID, &tp))
#else
    if(true)
#endif
      {
	// looks like clock does not work
	WARNING("dull selection of swapping mode");
	if(p==0.5)
	  return parallel3;
	else if(p > 0.05){
	  // select by accuray:
	  for(swapping_mode m=parallel1;
	      m<n_swapping_modes;
	      m=swapping_mode(int(m)+1) ){
	    if(int(p*(2<<(accuracy_in_bits[m]-precision_goal_in_bits))))
	      return m;
	  }
	  FATAL_ERROR("could not find appropriate swapping mode");
	} else
	  return logarithmic;
      }
  }
  DEBUG("getting fastest swapping_mode for p=" << p);
#ifdef CLOCK_PROCESS_CPUTIME_ID
  props_t props;
  const int n_reps=int(1e5);
  sequence<double> extime;

  for(swapping_mode m=swapping_mode(0);m<n_swapping_modes;
      m=swapping_mode(int(m)+1) ){
    timespec tp={0,0};
    SYSCALL(clock_settime(CLOCK_PROCESS_CPUTIME_ID, &tp));
    double p_hold=p;
    for(int n=0;n<n_reps;n++){
      props.swap_a_fraction(p_hold,m);
    }
    SYSCALL(clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp));
    extime[m]=(tp.tv_sec+tp.tv_nsec*1e-9)/n_reps*1e6;
  }
  swapping_mode fastest_mode=swapping_mode(0);
  for(swapping_mode m=swapping_mode(0);m<n_swapping_modes;
      m=swapping_mode(m+1) ){
    if(extime[m]<extime[fastest_mode]){
      if(int(p*(2<<(accuracy_in_bits[m]-precision_goal_in_bits))))
	fastest_mode=m;
    }
  }
  cout << "found mode " << fastest_mode << " fastest at p=" << p 
       << " (" << extime[fastest_mode] << "ms)" << endl;
  return fastest_mode;
#else
  return parallel3; //this point should never be reached
#endif
}

void props_t::measure_swapping_methods(){
#ifdef CLOCK_PROCESS_CPUTIME_ID
  props_t props;
  double p;
  swapping_mode m;
  const int n_reps=int(1e5);
  sequence<double> pp;
  sequence< sequence<double> > extime;
  int i;
  
  for(m=swapping_mode(0);m<n_swapping_modes;
      m=swapping_mode(int(m)+1) ){
    cout << "mode " << m << endl;
    for(p=0.5,i=0;p>0.001;p*=(1/1.3),i++ ){
      cout << "p=" << p << ": ";
      pp[i]=p;
      timespec tp={0,0};
      SYSCALL(clock_settime(CLOCK_PROCESS_CPUTIME_ID, &tp));
      for(int n=0;n<n_reps;n++){
	props.swap_a_fraction(p,m);
      }
      SYSCALL(clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tp));
      extime[m][i]=(tp.tv_sec+tp.tv_nsec*1e-9)/n_reps*1e6;
      cout << extime[m][i] << " ms per call" << endl;
    }
  }
  for(i=0;i<int(pp.size());i++){
    swapping_mode fastest_mode=swapping_mode(0);
    for(m=swapping_mode(0);m<n_swapping_modes;
	m=swapping_mode(int(m)+1) ){
      if(extime[m][i]<extime[fastest_mode][i]){
	if(int(pp[i]*(2<<(accuracy_in_bits[m]-precision_goal_in_bits))))
	  fastest_mode=m;
      }
    }
    cout << "mode " << fastest_mode 
	 << " is fastest at p=" << pp[i] 
	 << " (" << extime[fastest_mode][i] << "ms)" << endl;
    cout << "mode " << props_t::fastest_mode(pp[i]) 
	 << " recommended" << endl;
  }
#endif 
  exit(0);
}

#include "Statistics.h"

void props_t::test_modes(){
  props_t props,props0;
  average_meter p_hit;
  const int n_reps=100;
  double p=unirand();
  swapping_mode m;
  for(m=swapping_mode(0);m<n_swapping_modes;
      m=swapping_mode(int(m)+1) ){
    for(int i=n_reps;i>0;i--){
      props0.set_random();
      props=props0;
      props.swap_a_fraction(p,m);
      p_hit.sample(1-props.matches_with(props0)/double(n_property_bits));
    }
    cout << "mode " << m 
	      << " swappes " << p_hit.readout() 
	      << " +/- " << p_hit.error()
	      << ", expected " << p 
	 << "(t=" << (p-p_hit.readout())/p_hit.error() << ")"
	      << endl;
  }
}

    
