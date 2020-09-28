// $Id: hotspots.cc 2504 2017-02-28 10:54:54Z axel $

/// \file: The computation intensive parts of NewWeb simulation.

//**** fast exponential *****//
/* Source: A Fast, Compact Approximation of the Exponential Function
to appear in Neural Computation 11(4) Technical Report IDSIA-07-98
Nicol N. Schraudolph A, Corso Elveznic@idsia.ch
*/
#include <math.h>
#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN
#endif
static union eco
{
  double d;
  struct 
  {
#ifdef LITTLE_ENDIAN
    int j, i;
#else
    int i, j;
#endif
  } n;
} _eco;
#define EXP_A (1048576/M_LN2)
#define EXP_C 60801
/* use 1512775 for integer version */
/* see text for choice of c values */
inline double quickEXP(double y){
  eco _eco;
  _eco.n.i = EXP_A*(y) + (1072693248 - EXP_C);
  return _eco.d;
}
//**** end of fast exponential *****//
  

#define myEXP exp  // Make this quickEXP to give it a try.

#include "NewWeb.h"
#include "SortedMatrix.h"
#include "Integrator.h"
#include "packed_simulation.h"

static double relative_yearly_modulation=0;///< not implemented
static int multithreading_only_fit=0;
static double averaging_time=10; // in years FIXME: use 'eval("10*year")'
static int max_invasion_attempts=9999;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(relative_yearly_modulation),
  CFGINT(multithreading_only_fit),
  CFGINT(max_invasion_attempts),
  CFGDOUBLE(averaging_time),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

const double log_DBL_MAX=5*log(DBL_MAX)/7.0; ///< we added some safety margin
const double log_DBL_MIN=5*log(DBL_MIN)/7.0; ///< we added some safety margin
int max_num_threads=1; ///< maximal number of allowed threads (changed
		       ///by main())

const double NewWeb::do_dynamics_manager_t::unset=-1;


/// This is just a diagnostic helper function, not needed for
/// population dynamics.
int NewWeb::compute_flows(const ODE_vector & state,bool initialize){
  ODE_vector time_derivative(number_of_variables());
  bool initialized;

  // sets common_factor and biomass_B as a side effect:
  if(dynamics_for_side_effects(state,time_derivative)!=0) return 1;

  
  if(do_switching){

    if(initialize){
      compute_flows_csc=
	std::vector< compute_flows_csc_matrix_t >
	(s.n_animals,
	 compute_flows_csc_matrix_t(s.size()) );
      
      for(int k=s.n_animals;k-->0;){
	for(int i=s.size();i-->0;){
	  double sum=0;
	  for(int j=s.size();j-->0;){
	    compute_flows_csc[k][i][j]=
	      precomputed[k].c[i]*
	      s(i).switching_similarity_to(s(j),
					   s(k).prey_similarity_width())*
	      precomputed[k].c[j];
	  }
	}
      }
      initialized = true;
    }else{//User claims this is initialized, let's double check!
      if(compute_flows_csc.size()!=s.n_animals or
	 (s.n_animals and compute_flows_csc[0].size()!=s.size()) ){
	WARNING("compute_flows is not initialized.");
	WARNING("This computation will take less space but up to O(S^3) time!");
	initialized = false;
      }else{
	initialized = true;
      }
    }
    
    if(initialized){
      for(int k=s.n_animals;k-->0;){
	for(int i=s.size();i-->0;){
	  fx(i,k)=biomass_B(i)*common_factor(k)*
	    compute_flows_csc[k][i].dot(biomass_B);
	}
      }
    }else{
#if 0 //// I need to remove prior computation of the switching matrix
      // here and below to accommodate variable switching
      // similarities.  If speed is a problem we have to do combined
      // the new and the old version of the code here.  We did not
      // initialize.  Try to be as fast as possible anyway.
      link_strength_matrix similarity;
      similarity.resize(s.size());
      for(int i=s.size();i-->0;){
	similarity[i][i]=s(i).switching_similarity_to(s(i));
	for(int j=i;j-->0;){
	  similarity[i][j]=similarity[j][i]=s(i).switching_similarity_to(s(j));
	}
      }
#endif
      for(int k=s.n_animals;k-->0;){
	for(int i=s.size();i-->0;){
	  fx(i,k)=0;
	}
      }
      std::vector< double > availability(s.size());
      std::vector< int > available_prey(s.size());
      for(int k=s.n_animals;k-->0;){
	int n_prey=0;
	for(int i=s.size();i-->0;){
	  if(precomputed[k].c[i]>0){
	    availability[i]=biomass_B(i)*precomputed[k].c[i];
	    available_prey[n_prey]=i;
	    n_prey++;
	  }
	}
	for(int ii=n_prey;ii-->0;){
	  int i=available_prey[ii];
	  double sa_sum=0;
	  for(int jj=n_prey;jj-->0;){
	    int j=available_prey[jj];
	    sa_sum+=
	      s(i).switching_similarity_to(s(j),s(k).prey_similarity_width())*
	      availability[j];
	  }
	  //	  sa_sum*=  // What the hell was this line good for?
	  fx(i,k)=availability[i]*sa_sum*common_factor(k);
	}
      }
    }
  }else{// no switching

    for(int k=s.n_animals;k-->0;){
      for(int i=s.size();i-->0;){
	fx(i,k)=precomputed[k].c[i]*biomass_B(i)*common_factor(k);
      }
    }

  }
  fix_fx();
  return 0; // OK
}

/// This is just a diagnostic helper function, not needed for
/// population dynamics.
void NewWeb::get_diet(int k /*predator*/,sequence< double > & diet){
  // This code closely follows correspondig code in compute_flows
  std::vector< double > availability(s.size());
  std::vector< int > available_prey(s.size());
  int n_prey=0;
  for(int i=s.size();i-->0;){
    diet[i]=0;
    if(precomputed[k].c[i]>0){
      availability[i]=biomass_B(i)*precomputed[k].c[i];
      available_prey[n_prey]=i;
      n_prey++;
    }
  }
  for(int ii=n_prey;ii-->0;){
    int i=available_prey[ii];
    double sa_sum=0;
    for(int jj=n_prey;jj-->0;){
      int j=available_prey[jj];
      sa_sum+=
	s(i).switching_similarity_to(s(j),s(k).prey_similarity_width())*
	availability[j];
    }
    //	  sa_sum*=  // What the hell was this line good for?
    diet[i]=availability[i]*sa_sum;
  }
  diet/=sum(diet);
  return;
}  


void NewWeb::fix_fx(){
  for(int k=s.size();k-->s.n_animals;){
    for(int i=s.size();i-->0;){
      fx(i,k)=0;
    }
  }
}

//#define NAN_DIAGNOSTIC //comment this out to hunt nan and inf
#ifdef NAN_DIAGNOSTIC
static void nantest(ODE_vector & dt){
  for(int i=dt.size();i-->0;){
    if(dt[i]==dt[i]+1){
      WARNING(dt);
      abort();
    }
  }
}
#endif


/// helper thread to compute dynamics()
void NewWeb::do_dynamics(int offset,
			 ODE_vector const & state, 
			 ODE_vector & time_derivative){

  const int num_threads=do_dynamics_manager.get_num_threads();

  // Maximum common_factor computed in this thread.
  double common_factor_thread_max=0;

#ifdef NAN_DIAGNOSTIC
  if(time_derivative.size()==0){
    abort();
  }
  for(int i=time_derivative.size();i-->0;){
    time_derivative[i]=0;
  }
  nantest(time_derivative);
#endif
  
  // animal population growth:
  if(do_switching){
    for(int k=offset;k<s.n_animals;k+=num_threads){
      // This is an efficint implementation of the functional response
      // (and the resulting consumer dynamics) described in the draft
      // vanLeeuwen2008.pdf
      const double sand=
	precomputed[k].csc_eating.sandwich_product(biomass_B,biomass_B);
      const double halfsat=precomputed[k].c.dot(biomass_B);
      const double fact0=
	(halfsat>0
	 ?
	 1/(halfsat+s(k).handling_time_T()*s(k).attack_rate_a()*sand)
	 :
	 0 );
      
      double cf=biomass_B[k]*s(k).attack_rate_a()*fact0;
      common_factor[k]=cf;
      common_factor_thread_max=std::max(common_factor_thread_max,cf);

      time_derivative[k]=
	s(k).conversion_efficiency_eps()*s(k).attack_rate_a()*sand*fact0;
      
#ifdef NAN_DIAGNOSTIC
      if(time_derivative[k]==time_derivative[k]+1){
	time_derivative[k];
	REPORT(sand);
	REPORT(halfsat);
	REPORT(fact0);
	REPORT(biomass_B[k]);
	FATAL_ERROR("nan results");
      }
#endif

      time_derivative[k]-=s(k).turnover_rate_r();
	
#ifdef NAN_DIAGNOSTIC
      if(time_derivative[k]==time_derivative[k]+1){
	time_derivative[k];
	FATAL_ERROR("nan results");
      }
#endif

      // Compute diagnostics if requested:
      if(compute_diagnostics){
	const double u=
	  s(k).handling_time_T()*s(k).attack_rate_a()*sand/halfsat;
	s(k).the_saturation_strength=u/(1+u);

	s(k).the_GP=
	  time_derivative[k]*biomass_B(k);
      }

    }
  }else{ // no switching
    for(int k=offset;k<s.n_animals;k+=num_threads){
      const double halfsat=1;
      const double total_availability=
	precomputed[k].c.dot(biomass_B);
      const double fact0=
	1/(halfsat+s(k).handling_time_T()*s(k).attack_rate_a()
	   *total_availability);
      
      double cf=biomass_B[k]*s(k).attack_rate_a()*fact0;
      common_factor[k]=cf;
      common_factor_thread_max=std::max(common_factor_thread_max,cf);
      
      time_derivative[k]=
	s(k).conversion_efficiency_eps()*s(k).attack_rate_a()
	*total_availability*fact0;
      
#ifdef NAN_DIAGNOSTIC
      if(time_derivative[k]==time_derivative[k]+1){
	time_derivative[k];
	REPORT(sand);
	REPORT(halfsat);
	REPORT(fact0);
	REPORT(biomass_B[k]);
	FATAL_ERROR("nan results");
      }
#endif
      
      // Compute diagnostics if requested:
      if(compute_diagnostics){
	const double u=
	  s(k).handling_time_T()*s(k).attack_rate_a()*total_availability/
	  halfsat;
	s(k).the_saturation_strength=u/(1+u);
	s(k).the_GP=time_derivative[k]*biomass_B[k];
      }

      time_derivative[k]-=s(k).turnover_rate_r();

#ifdef NAN_DIAGNOSTIC
      if(time_derivative[k]==time_derivative[k]+1){
	time_derivative[k];
	FATAL_ERROR("nan results");
      }
#endif

    }
  }

#ifdef NAN_DIAGNOSTIC
  nantest(time_derivative);
#endif

  // Threads with odd offset write their max common_factor here
  // already to reduce delays in computation of overall max later.
  if((offset & 1)){
    do_dynamics_manager.common_factor_max[offset]=
      common_factor_thread_max;
  }

  
  // Plant population growth:
  for(int i=s.n_animals+offset;i<s.size();i+=num_threads){
    // See Rossberg et al. (2008) Ecology 89:567–580 for the underlying ideas.
    
    // saturation effects:
    double shadowing_sum=precomputed[i].c.dot(biomass_B);
    double light;
    if(plant_physiology_version!=3){
      light=myEXP(-shadowing_sum);
    }else{
      light=1-shadowing_sum;
    }
    
    if(compute_diagnostics){
      // !!!  These parts need re-implementation:
      // 	s[i].the_top_down_strength=being_eaten/
      // 	  (s(i).the_plant_growth_rate_sigma*s(i).the_loss_rate_over_max_production_rate_r
      // 	   *biomass_B(i) );
      //s[i].the_competition_strength=
      //         1-biomass_B(i)*c(i,i)/light_losses;
      if(plant_physiology_version==3){
	s[i].the_saturation_strength=light;
	s[i].the_light_strength=light;
	s[i].the_GP=s(i).plant_growth_rate_sigma()*light*biomass_B(i);
	s[i].the_competition_strength=
	  shadowing_sum-precomputed[i].c[i]*biomass_B[i];
      }else{
	s[i].the_saturation_strength=light;
	s[i].the_light_strength=light;
	s[i].the_GP=
	  s(i).plant_growth_rate_sigma()*
	  (1-s(i).get_loss_rate_over_max_production_rate_r())*
	  biomass_B(i);
      }
    }

    if(plant_physiology_version==3){
      time_derivative[i]=s(i).plant_growth_rate_sigma()*(light+s(i).env_effect());
    }else{
      time_derivative[i]=s(i).plant_growth_rate_sigma()*
	(light-s(i).get_loss_rate_over_max_production_rate_r());
    }
    
#ifdef NAN_DIAGNOSTIC
    if(time_derivative[i]==time_derivative[i]+1){
      time_derivative[i];
      REPORT(light);
      REPORT(biomass_B[i]);
      FATAL_ERROR("nan results");
    }
#endif

  }
  
  // Parallel computation of maximum of common_factor.  Implicitly
  // this is also a barrier for all threads.
  int other_bit=1;
  // Threads with index divisible by 2^n compute max of the
  // common_factor entries it computed itself, and those of 2^n-1
  // following thread indices.  This leads to a hierarchy of thread
  // dependencies, which is implemented by the shifting bit mask
  // other_bit.
  if((offset & 1)==0){
    // threads with odd offset have written their max already
    do{
      if( (offset | other_bit)>=num_threads )
	break;
      while(do_dynamics_manager.max_not_set(offset | other_bit)){
 	boost::this_thread::yield();
      }
      common_factor_thread_max=
	std::max(common_factor_thread_max,
		 do_dynamics_manager.common_factor_max[offset | other_bit] );
      other_bit <<= 1;
    }while((offset & other_bit)==0);
  }
  if(offset){
    // Save the max this thread is responsible for:
    do_dynamics_manager.common_factor_max[offset]=
      common_factor_thread_max;
    // Wait for overall max to be finished:
    do_dynamics_manager.barrier->wait();
  }else{
    common_factor.set_max(common_factor_thread_max);
    // Overall max is finished:
    do_dynamics_manager.barrier->wait();
  }
  //... finished computation of maximum common_factor.


#ifdef NAN_DIAGNOSTIC
  nantest(time_derivative);
#endif
  
  // Predation/consumption mortality for all species:
  for(int i=offset;i<s.size();i+=num_threads){
    double being_eaten;
    if(do_switching){
      being_eaten=
	precomputed[i].csc_being_eaten.
	sandwich_product(biomass_B,common_factor);
    }else{
      being_eaten=
	precomputed[i].cT.dot(common_factor);
    }      
    
    time_derivative[i]-=being_eaten;
    
    if(compute_diagnostics){
      if(is_plant(i)){
	if(plant_physiology_version==3){
	  s[i].the_top_down_strength=being_eaten/
	    s(i).plant_growth_rate_sigma();
	}else{
	  s[i].the_top_down_strength=being_eaten/
	    (s(i).plant_growth_rate_sigma()*
	     (1-s(i).get_loss_rate_over_max_production_rate_r()));
	}
      }else{
	s[i].the_top_down_strength=being_eaten/
	  (s[i].turnover_rate_r());
      }
    }
    //     REPORT(precomputed[i].csc_being_eaten);
    //     REPORT(log(s[i].the_mean_bodymass_M/eval("1*kilogram"))/log(10));
    //     REPORT(being_eaten*eval("1*year*meter^2")/area_per_compartment);

    /// PATCH TO AVOID EXTINCTIONS (don't use -z option, works great!):
    //time_derivative[i]+=100*s(i).turnover_rate_r()*s(i).get_threshold_biomass_B()/biomass_B(i);
  }
#ifdef NAN_DIAGNOSTIC
  nantest(time_derivative);
#endif
}


/// Main function for spawned threads:
void NewWeb::do_dynamics_manager_t::task_t::operator()(){

  NewWeb::do_dynamics_manager_t & mgr=
    dispatcher->do_dynamics_manager;


  while(true){

    {boost::mutex::scoped_lock lock(mgr.mutex);
      mgr.common_factor_max[index]=do_dynamics_manager_t::unset;  

      mgr.threads_ready++;
      mgr.condition.wait(mgr.mutex);
      mgr.threads_ready--;
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

NewWeb::do_dynamics_manager_t::do_dynamics_manager_t():
  _num_threads(1),threads_ready(0),common_factor_max(1),barrier(new boost::barrier(1))
{
}

/// Don't copy do_dynamics_manager_t when copying NewWeb:
NewWeb::do_dynamics_manager_t::
do_dynamics_manager_t(const do_dynamics_manager_t & other):
  _num_threads(1),threads_ready(0),common_factor_max(1),barrier(new boost::barrier(1)) {
}

/// Don't assign to do_dynamics_manager_t when assigning to NewWeb:
NewWeb::do_dynamics_manager_t& NewWeb::do_dynamics_manager_t::
operator=(const NewWeb::do_dynamics_manager_t& other){
}
  
void NewWeb::do_dynamics_manager_t::
initialize_threads_maybe(NewWeb * base_web){
  const int S=base_web->number_of_species();
  const int threading_threshold=
    (do_switching ? 250 : 500);
  if(S < threading_threshold){
    if(_num_threads>1){
      WARNING("Stopping multithreading.");
      stop_threads();
    }
  }else if(S > threading_threshold){
    if(_num_threads <= 1 && max_num_threads > 1){
      WARNING("Starting multithreading with " 
	      << max_num_threads << " threads.");
      start_threads(max_num_threads,base_web);
    }
  }
}

void NewWeb::do_dynamics_manager_t::start_threads(int n,NewWeb * base_web){
  if(n!=_num_threads){
    stop_threads();
    _num_threads=n;
    task.resize(_num_threads-1);
    delete barrier;
    barrier=new boost::barrier(_num_threads);
    common_factor_max.resize(_num_threads);
    stop_now=false;
    for(int i=_num_threads-1;i-->0;){
      task[i].dispatcher=base_web;
      task[i].index=i+1;
      threads.create_thread(task[i]);
    }
    while(threads_ready<_num_threads-1){
      boost::this_thread::yield();
    }
  }
}

void NewWeb::do_dynamics_manager_t::stop_threads(){
  if(_num_threads>1){
    stop_now=true;

    {boost::mutex::scoped_lock lock(mutex);
      condition.notify_all();
    }
    threads.join_all();
    delete barrier;
    barrier=new boost::barrier(1);
    _num_threads=1;
    common_factor_max.resize(_num_threads);
    task.resize(0);
  }
}

NewWeb::do_dynamics_manager_t::~do_dynamics_manager_t(){
  stop_threads();
  delete barrier;
}


int NewWeb::
dynamics_for_side_effects(ODE_vector const & state, 
			  ODE_vector & time_derivative){
  int ups=use_packed_simulation;
  packed_simulation * da=dynamics_accelerator;
  use_packed_simulation=0;
  dynamics_accelerator=0;
  int result=
   dynamics(state,time_derivative);
  use_packed_simulation=ups;
  dynamics_accelerator=da;
  return result;
}

int NewWeb::dynamics(ODE_vector const & state, 
		     ODE_vector & time_derivative)

{  
  if(dynamics_accelerator){
    return dynamics_accelerator->
      dynamics(state,time_derivative);
  }

  if(state.size()==0) return 1;

  do_dynamics_manager_t & mgr = do_dynamics_manager;

  mgr.initialize_threads_maybe(this);

  biomass_B.resize(s.size());
  common_factor.resize(s.n_animals+1);
  
  // Check if state in reasonable range:
  bool failure=false;
  for(int i=s.size();i-->0;){
    if(state[i] > log_DBL_MAX/2 /*the '2' here is the assumed
 				  switching exponent!*/  ||
       state[i] < log_DBL_MIN/2){
      WARNING("state[" << i << "] = " << state[i]);
      REPORT(assigned_column[i]);
      REPORT(s(i).turnover_rate_r());
      REPORT(s(i).get_trade_off_multiplier());
      REPORT(s(i).plant_growth_rate_sigma());
      failure=true;
      break;
    }else{
      // In principle, these exps could also be computed in parallel,
      // but the cost of meeting at a barrier afterwards is too
      // high (and it is linear in problem size!).
      biomass_B[i]=myEXP(state[i]);
    }
  }

  biomass_B.fix_max();

  if(failure){
    newton_failure=1;
    return 1; // recoverable failure calculating dynamics
  }
  
  // Pass data to other threads:
  mgr.state=&state;
  mgr.time_derivative=&time_derivative;

  mgr.common_factor_max[0] = do_dynamics_manager_t::unset;

  // Run parallel computation:
  {boost::mutex::scoped_lock lock(mgr.mutex);
    mgr.condition.notify_all();
  }
  do_dynamics(0,state,time_derivative);

  // Wait until all threads are ready to start the next iteration
  // (this serves also as exit barrier for do_dynamics).  One should
  // try to wait for all threads to become ready before, rather than
  // after doing computations.  But then an extra exit barrier would
  // be needed.
  while(mgr.threads_ready<
	mgr.get_num_threads()-1){
    boost::this_thread::yield();
  }

  return 0; // success calculating dynamics
}

// ///////////////////////////////////////////////////////////// 
// 
// The remainder of this file is for fitness computations.  The
// computations here must be adjusted if the computations in
// "do_dynamics(...)" above change.

/// Compute the linear growth rate of the species pointed to by sp.
/// This species is not yet included in the species list s.
double NewWeb::fast_linear_fitness(NewSpecies * sp){
  // ATTENTION: 
  // this function assumes common_factor and biomass_B to be set!!

  // This should be computationally cheapest for momentary linear
  // growth rate. For time averaging use version with pre-computed
  // coefficients below.
  
  const int S=s.size();
  double time_derivative;
  double matrix_truncation_epsilon=SortedMatrix().truncation_epsilon();
  
  if(sp->is_a_plant()){
    double shadowing_sum=0;
    for(int i=s.size();i-->s.n_animals;){
      shadowing_sum+=
	sp->plant_hampering_by(s(i))*
	biomass_B[i];
    }
    if(plant_physiology_version==3){
      time_derivative=
	sp->plant_growth_rate_sigma()*(1-shadowing_sum+sp->env_effect());
    }else{
      const double light=myEXP(-shadowing_sum);
      time_derivative=sp->plant_growth_rate_sigma()*
	(light - sp->get_loss_rate_over_max_production_rate_r());
    }
  }else{
    // *sp is animal:
    double functional_response;
    if(do_switching) {
      double Asum=0; // same as "total_availability" below
      double csc_sum=0;
      std::vector< int > prey(s.size()); // list of prey of a predator
      std::vector< double > A(s.size()); // Availability
      int Z=0;
      
      for(int i=s.size();i-->0;){
	double cc=sp->foraging_strength_c_on(s(i));
	if(cc > matrix_truncation_epsilon){
	  double AA=cc*biomass_B[i];
	  Asum+=AA;
	  A[Z]=AA;
	  prey[Z]=i;
	  csc_sum+=AA*AA;//contribution from diagonal term
			 //(similarity==1)
	  for(int jj=Z;jj-->0;){
	    // we make use of the symmetry here!
	    int j=prey[jj];
	    double contribution=AA*A[jj]*
	      s(i).switching_similarity_to(s(j),sp->prey_similarity_width());
	    csc_sum+=2*contribution;// also the transposed element
	  }
	  Z++;
	}
      }
      
      csc_sum*=sp->attack_rate_a();
      functional_response=
	csc_sum/
	(Asum+sp->handling_time_T()*csc_sum);
    }else{ // no switching
      double total_availability=0;
      for(int i=s.size();i-->0;){
	double cc=sp->foraging_strength_c_on(s(i));
	if(cc > matrix_truncation_epsilon){
	  total_availability+=cc*biomass_B[i];
	}
      }
      
      total_availability*=sp->attack_rate_a();
      functional_response=
	total_availability/
	(1+sp->handling_time_T()*total_availability);
    }
    
    time_derivative=
      sp->conversion_efficiency_eps()*functional_response-
      sp->turnover_rate_r();
  }
  
  // finally, the being-eaten term:
  double being_eaten=0;
  if(do_switching) {
    // compute pr.csc_being_eaten
    std::vector< double > pressure_on_sp_by(s.n_animals);
    std::vector< int > predator(s.size()); // list of predators of sp
    int Z=0;
    
    for(int j=s.n_animals;j-->0;){
      double cc=s(j).foraging_strength_c_on(*sp);
      if(cc > matrix_truncation_epsilon){
	pressure_on_sp_by[Z]=cc*common_factor[j];
	predator[Z++]=j;
      }
    }
    
    for(int i=s.size();i-->0;){
      for(int jj=Z;jj-->0;){
	int j=predator[jj];
	double contribution=precomputed[j].c[i]*
	  pressure_on_sp_by[jj]*
	  sp->switching_similarity_to(s(i),s(j).prey_similarity_width())*
	  biomass_B[i];
	being_eaten+=contribution;
      }
    }  
  } else{
    for(int j=s.n_animals;j-->0;){
      being_eaten += 
	common_factor[j]*s(j).foraging_strength_c_on(*sp);
    }
  }
  time_derivative-=being_eaten;

  return time_derivative;
}


void NewWeb::precompute(NewSpecies * sp,precomputed_entry_t & pr){
  double matrix_truncation_epsilon=SortedMatrix().truncation_epsilon();
  int S=s.size();
  pr.c.resize(S);
  pr.cT.resize(S);
  pr.csc_eating.resize(S);
  pr.csc_being_eaten.resize(S);

  if(sp->is_a_plant()){
    for(int i=s.size();i-->s.n_animals;){
      pr.c[i]=sp->plant_hampering_by(s(i));
    }
  }else{
    std::vector< int > prey(s.size()); // list of prey of a predator
    int Z=0;
    for(int i=s.size();i-->0;){
      double cc=sp->foraging_strength_c_on(s(i));
      if(cc > matrix_truncation_epsilon){
	pr.c[i]=cc;
	if(do_switching){
	  prey[Z]=i;
	  for(int jj=Z+1;jj-->0;){
	    // we make use of the symmetry here!
	    int j=prey[jj];
	    double entry=cc*pr.c[j]*
	      s(i).switching_similarity_to(s(j),sp->prey_similarity_width());
	    if(entry > matrix_truncation_epsilon){
	      pr.csc_eating[i][j]=entry;
	    }
	  }
	  Z++;
	}
      }
    }
  }

  if(do_switching){
    // compute pr.csc_being_eaten
    std::vector< double > c_on_sp_by(s.n_animals);
    std::vector< int > predator(s.size()); // list of predators of sp
    int Z=0;
    for(int i=s.n_animals;i-->0;){
      double cc=s(i).foraging_strength_c_on(*sp);
      if(cc > matrix_truncation_epsilon){
	c_on_sp_by[Z]=cc;
	predator[Z++]=i;
      }
    }
    
    for(int i=s.size();i-->0;){
      for(int jj=Z;jj-->0;){
	int j=predator[jj];
	double entry=precomputed[j].c[i]*
	  c_on_sp_by[jj]*
	  sp->switching_similarity_to(s(i),s(j).prey_similarity_width());
	if(entry > matrix_truncation_epsilon){
	  pr.csc_being_eaten[i][j]=entry;
	}
      }
    }
  }else{// not switching
    for(int i=s.n_animals;i-->0;){
      pr.cT[i]=s(i).foraging_strength_c_on(*sp);
    }
  }
}

/// High level function
double NewWeb::linear_fitness(NewSpecies * sp){
  // ATTENTION: 
  // this function assumes common_factor and biomass_B to be set!!

  // This does in principle the same as fast_linear_fitness(..)
  // above, but goes through the same approximations as used for full
  // simulations, and is therefore more reliable in predicting
  // invadibility in simulations.

  precomputed_entry_t pr;
  precompute(sp,pr);
  return linear_fitness(sp,biomass_B,common_factor,pr);
}

/// Low level function
double NewWeb::linear_fitness(NewSpecies * sp,
			      const vector_with_max &biomass_B,
			      const vector_with_max &common_factor,
			      const precomputed_entry_t & pr ){
  // Compute the linear growth rate of the species pointed to by sp.
  // This species is not yet included in the species list s.
  
  
  const int S=s.size();
  double time_derivative;

  if(sp->is_a_plant()){
    // *sp is plant:
    double shadowing_sum=pr.c.dot(biomass_B);
    if(plant_physiology_version==3){
      time_derivative=
	sp->plant_growth_rate_sigma()*(1-shadowing_sum+sp->env_effect());
    }else{
      const double light=myEXP(-shadowing_sum);
      time_derivative=sp->plant_growth_rate_sigma()*
	(light-sp->get_loss_rate_over_max_production_rate_r());
    }
  }else{
    // *sp is animal:
    double functional_response;
    if(do_switching) {
      double Asum=pr.c.dot(biomass_B);
      double csc_sum=pr.csc_eating.sandwich_product(biomass_B,biomass_B);
    
      csc_sum*=sp->attack_rate_a();
      functional_response=
	csc_sum/
	(Asum+sp->handling_time_T()*csc_sum);
    } else{
      double total_availability=pr.c.dot(biomass_B);
    
      total_availability*=sp->attack_rate_a();
      functional_response=
	total_availability/
	(1+sp->handling_time_T()*total_availability);
    }
    
    time_derivative=
      sp->conversion_efficiency_eps()*functional_response-
      sp->turnover_rate_r();
  }
  
  // finally, the being-eaten term:
  double being_eaten;
  if(do_switching) {
    being_eaten=pr.
      csc_being_eaten.sandwich_product(biomass_B,common_factor);
  }else{
    being_eaten= 
      pr.cT.dot(common_factor);
  }
  time_derivative-=being_eaten;

  return time_derivative;
}

double NewWeb::invasion_fitness(NewSpecies * sp){
  if(steady_state.size()<=1){
    double f=linear_fitness(sp);
    return f;
  }
  
  precomputed_entry_t pr;
  precompute(sp,pr);

  double growth_rate;
  sequence<double *> gr;
  gr[0]=&growth_rate;
  {
    Integrator gr_integrator(gr);

    steady_state_t::iterator start=steady_state.begin();

    if(steady_state.size() > max_steady_state_size){
      start+=steady_state.size()-max_steady_state_size;
    }
    for(steady_state_t::iterator i=start;
	i!=steady_state.end();
	++i){
      
      const saved_state & s=*i;
      growth_rate=
	linear_fitness(sp,s.biomass_B,s.common_factor,pr);
      gr_integrator.sample(i->t);
    }
    
    gr_integrator*=1/(steady_state.back().t-steady_state.front().t);
  }// write back of growth_rate at destruction of gr_integrator:
  return growth_rate;
}


// Machinery to find fit species by repeatedly creating species and
// measuring fitness.  Each cycle is called an "attempt".  If
// multithreading_only_fit==true, uses multithreading.  In this case,
// output is undetermined, because of concurrent usage of the random
// number generator.  

/// If multithreading_only_fit, all threads communicate via
/// attempts_left.  A positive value means continue.  If a thread
/// finished one attempt, it will decrement the variable.  A negative
/// value means that a fit species was found with abs(attempts_left)
/// attempts left, and its address is written to new_species.  A value
/// of zero means we have to stop trying.  The thread that sets
/// attempts_left to zero, writes the address of its current candidate
/// to new_species.  Mutex_random serves to prevent race conditions
/// for random-number generation.  We need to check thread-safety of
/// random number generators declared in random.h to ensure we can do
/// without mutex_random.
static boost::mutex mutex_attempts_left, mutex_random;
static int attempts_left;
static NewSpecies * new_species;
static NewWeb * fitness_thread_dispatcher;
static NewWeb::finder_function_t finder_function;

/// This defines the task for a thread that repeatedly lookings for a
/// fit species (usage of global variables is a bit messy at the
/// moment):
struct find_fit_task{
  double plant_fraction;
  void operator()(){
    do{
      NewSpecies * my_new_species;

      {boost::mutex::scoped_lock lock(mutex_random);
	my_new_species=
	  (fitness_thread_dispatcher->*finder_function)(plant_fraction);
      }

      double fit=fitness_thread_dispatcher->invasion_fitness(my_new_species);
      
      {boost::mutex::scoped_lock lock(mutex_attempts_left);
	if(attempts_left <= 0){
	  /* Nothing left to be done, forget what you found	*/
	  delete my_new_species;
	  break;
	}else{
	  if(fit>=0){
	    /* Success! Pass my_new_species back and make other
	       threads stop:*/
	    new_species=my_new_species;
	    attempts_left=-attempts_left;
	    break;
	  }else{
	    /* Failure! Continue searching...*/
	    attempts_left--;
	    if(!attempts_left){
	      /*... but only if attempts_left.*/
	      new_species=my_new_species;
	      break;
	    }
	    delete my_new_species;
	  }
	}
      }// end of mutex_attempts_left lock

    }while(true);
  }
};

NewSpecies * NewWeb::find_only_fit(double plant_fraction, finder_function_t f){
  const int max_attempts=max_invasion_attempts;
  new_species=0;
  attempts_left=max_attempts;
  ODE_vector state(number_of_variables()),dummy(number_of_variables());
  write_state_to(state);
  dynamics_for_side_effects(state,dummy);/* called for side effects */
  if(!multithreading_only_fit || max_num_threads <= 1){
    do{
      new_species=(this->*f)(plant_fraction);
      double f=invasion_fitness(new_species);
      TRACE(f,DYNAMICS);
      if(f>=0 ||
	 --attempts_left <=0 ) break;
      delete new_species;
    }while(true);
  }else{/*multithreading implementation:*/
    fitness_thread_dispatcher=this;
    finder_function=f;
    boost::thread_group threads;
    std::vector< find_fit_task > task(max_num_threads);
    for(int i=0; i<max_num_threads; i++){
      task[i].plant_fraction=plant_fraction;
      threads.create_thread(task[i]);
    }
    threads.join_all();
    if(attempts_left<0) attempts_left=-attempts_left;
  }
  if(attempts_left <= 0){
    WARNING("adding unfit "
	    << (new_species->is_a_plant()?"plant":"animal"));
  }
  int failed_invasion_attempts=max_attempts-attempts_left;
  REPORT(failed_invasion_attempts);
  ALWAYS_ASSERT(new_species);
  return new_species;
}

NewSpecies * NewWeb::invade_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::invade); };
NewSpecies * NewWeb::speciate_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::speciate); };
NewSpecies * NewWeb::speciate_selected_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::speciate_selected); };
NewSpecies * NewWeb::invade_or_speciate_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::invade_or_speciate); };
NewSpecies * NewWeb::invade_or_speciate_any_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::invade_or_speciate_any); };
NewSpecies * NewWeb::invade_or_speciate_any_by_biomass_only_fit(double plant_fraction){ return find_only_fit(plant_fraction,&NewWeb::invade_or_speciate_any_by_biomass); };
