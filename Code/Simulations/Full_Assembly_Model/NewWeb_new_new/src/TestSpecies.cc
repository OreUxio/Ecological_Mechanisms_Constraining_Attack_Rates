// $Id: TestSpecies.cc 1416 2009-04-28 19:19:56Z axel $
#include <string>
std::string version("$Id: TestSpecies.cc 1416 2009-04-28 19:19:56Z axel $");
extern std::string compilation_time;

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE /* glibc2 needs this for strptime */
#endif
#ifndef _BSD_SOURCE
#define _BSD_SOURCE /* so we can set TZ in struct tm */
#endif

#include <signal.h>
#include <iostream>
#include <sstream>
#include <iomanip>

//#define DO_FPE  // do floating point exceptions
#ifdef DO_FPE
#ifdef __GNUC__
#include <fenv.h>  // floating point exceptions
#endif
#ifdef SX
#include <ieeefp.h>  // floating point exceptions
#endif
#endif

#include "ODE.h"
#include "scheduler.h"
#include "parsecfg.h" // include this header file when use parsecfg lib
#include "error.h"    // my error handlers
#include "XMLStore.h"
#include "matrix_transformers.h"
#include "standardize.h"
#include "cfgPermanent.h"
#include "random.h" 
#include "Statistics.h"
#include "NewWeb.h"
#include <unistd.h>  // for getpid  and sysconf
#include <sys/times.h> //for times :)

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
#include "linpack_eigen.h"
#include "mp_sug_numthreads.h"
#include "xy_graph.h"
#include "hotspots.h"
#endif

using namespace std;

string in_file_name="NewWorld.cfg";
string web_file_name="";

static my_evaluator_t eval_first;

// adjustable parameters:
static double check_interval=eval_first("1 * year");
static double chunk_length=1.0;
static double observation_time=eval_first("100 * years");
static double time_between_printout=eval_first("1 * week");
static double time_unit_for_output=eval_first("1*year");
static int print_each_step=0;
static int random_seed=1;
static int check_data_consistency=0;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(check_interval),
  CFGDOUBLE(chunk_length),
  CFGDOUBLE(observation_time),
  CFGDOUBLE(time_between_printout),
  CFGINT(print_each_step),
  CFGDOUBLE(time_unit_for_output),
  CFGINT(random_seed),
  CFGINT(check_data_consistency),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

bool has_suffix(char * suff,char * filename){
  return 0==strcmp(suff,filename+strlen(filename)-strlen(suff));
}

int number_of_tests=10;
double plant_fraction=1;
double max_test_duration=eval_first("100*years");
bool do_print_each_step_anyway=false;
const double unit_mass=eval_first("1*kilo_*gram");
const double unit_area=eval_first("1*meter^2"); 


void read_arguments_and_data(int argc,char *argv[],NewWeb & web){
  int c;
  
  while(1){
    c=getopt(argc, argv, "Phvn:p:t:");
    if (c == -1)
      break;
    
    switch(c){
    case 'P':
      do_print_each_step_anyway=true;
      break;
    case 't':
      max_test_duration=atof(optarg);
      break;
    case 'p':
      plant_fraction=atof(optarg);
      break;
    case 'n':
      number_of_tests=atoi(optarg);
      break;
    case 'v':
      printf("%s\n",version.c_str());
      exit(0);
      break;
    case 'h':
    default:
        fprintf(stdout,"\nusage: %s [foodweb_file]\n",argv[0]);
        exit(c=='h'?0:1);
      }
  }
  argc-=optind-1;
  argv+=optind-1;

  while(argc>1){
    if(has_suffix(".cfg",argv[1])){
      in_file_name=argv[1];
      std::cout << "reading parameter file " << in_file_name << std::endl;
      cfgStruct * cfg_list=full_cfg_list();
      if (cfgParse(in_file_name.c_str(), cfg_list, CFG_SIMPLE) == -1)
	FATAL_ERROR("error reading parameter file");
      delete [] cfg_list; //this is a raw array, we need to deallocate by hand
    }else if(has_suffix(".xml",argv[1]) ||
	     has_suffix(".xml.gz",argv[1]) ||
	     has_suffix(".xml.bz2",argv[1]) ){
      web_file_name=argv[1];
      std::cout << "reading food web " << web_file_name << std::endl;
      XMLStore store(web_file_name);
      cfgPermanent parameters;
      store.get(&parameters,"adjustable_parameters");

      // make code compatibile with buggy older versions:
      string old_compilation_time;
      store.get(old_compilation_time,"compilation_time");
      REPORT(old_compilation_time);
      //typical date: May 23 2006, 14:14:16
      struct tm time_data;
      char * strptime_retval=
	strptime(old_compilation_time.c_str(),"%b %d %Y, %T",&time_data);
      if(old_compilation_time.size()>0 && !strptime_retval) 
	FATAL_ERROR("could not parse compliation time of input web");
      else if(*strptime_retval)
	REPORT(strptime_retval);
      //and now the time of the bug fix: 2006/08/01 05:45:21 UTC
      struct tm bug_fix_time;
      strptime("Aug 1 2006, 14:45:21","%b %d %Y, %T",&bug_fix_time);
      REPORT(mktime(&time_data));
      REPORT(mktime(&bug_fix_time));
      if(old_compilation_time.size()==0 || 
	 mktime(&time_data)<mktime(&bug_fix_time) ){
	WARNING("input web is from old code,");
	WARNING("emulating carrying-capacity bug!");
	emulate_carrying_capacity_bug=1;
      }

      store.get(&web,"FoodWeb");
      // this is to read older versions with separated simulatedTime
      try{
	store.get(web.current_time,"simulatedTime");
      } catch (XMLStoreException){
	WARNING("... but this variable is obsolete, anyway.");
      }
    }else
      FATAL_ERROR("unknown file type: " << argv[1]);
    argv++;
    argc--;
  }
}


void signal_handling();

int exit_now=0;
int save_now=0;
void exiter(int i){
  if(i==SIGUSR1)
    exit_now=i;
  else if(i==SIGUSR2)
    save_now=i;
  else
    exit_now=1;
  if(i==SIGFPE)
    abort();
  signal_handling();
}

void signal_handling(){
  signal(SIGUSR1,&exiter);
  signal(SIGUSR2,&exiter);
  signal(SIGXCPU,&exiter);
  signal(SIGHUP,&exiter);
  signal(SIGTERM,&exiter);
  signal(SIGFPE,&exiter);
}

static bool integrator_trouble_last_step=false; //for integrator trouble

inline double CPU_time_in_clockticks(){
  struct tms time_struct;
  times(&time_struct);
  return (time_struct.tms_utime+
	  time_struct.tms_stime+
	  time_struct.tms_cutime+
	  time_struct.tms_cstime+
	  0);
}

class frozen_web : public NewWeb{
private:
  int test_species_index;
  ODE_vector test_state;
  ODE_vector test_derivative;
  double saved_current_time;
  void prepare_test_state(){
    test_state=ODE_vector(NewWeb::number_of_variables());
    test_derivative=ODE_vector(NewWeb::number_of_variables());
    saved_current_time=current_time;
  }
public:
  frozen_web(): test_species_index(-1){};
  virtual void dynamics(ODE_vector const & state, 
			ODE_vector & time_derivative){
    ASSERT(test_species_index >= 0);
    test_state[test_species_index]=state[0];
    NewWeb::dynamics(test_state,test_derivative);
    time_derivative[0]=test_derivative[test_species_index];
  };
  virtual bool can_calculate_Jacobian(){return false;};
  virtual bool has_preconditioner(){return false;};
  virtual void write_state_to(ODE_vector & state){
    ALWAYS_ASSERT(test_species_index >= 0);
    NewWeb::write_state_to(test_state);
    state[0]=test_state[test_species_index];
  }
  virtual void read_state_from(ODE_vector & state){
    ALWAYS_ASSERT(test_species_index >= 0);
    test_state[test_species_index]=state[0];
    NewWeb::read_state_from(test_state);
  }
  virtual int number_of_variables(){
    return 1; // 
  }
#define TEST_SPECIES(fname)			\
  void fname(double plant){			\
    if(test_species_index >= 0){		\
      FATAL_ERROR("test species already set");	\
    }						\
    test_species_index=				\
      NewWeb::fname(plant);			\
    prepare_test_state();			\
  }
  TEST_SPECIES(invade_only_fit);
  TEST_SPECIES(speciate_only_fit);
  void delete_test_species(){
    ALWAYS_ASSERT(test_species_index >= 0);
    delete_species(test_species_index);
    current_time=saved_current_time;
    test_species_index=-1;
  }
  void report_test_species(){
    ALWAYS_ASSERT(test_species_index >= 0);
    double area=area_per_compartment/unit_area;
    cout << log10(s[test_species_index].the_mean_bodymass_M/unit_mass);
    cout << " " << 
      log10(s[test_species_index].the_biomass_abundance_B/unit_mass/area);
    cout << endl;
  }
};

int main(int argc,char *argv[]){
  // find out how long we can run this job
  struct tms dummy_time_struct; // we do not read it
  const clock_t wallclock_start_time=times(&dummy_time_struct);
  clock_t wallclock_when_entering_main_loop; // we get this later
#ifdef _SC_CLK_TCK
  const double clockticks_per_second=sysconf(_SC_CLK_TCK);
#else 
  const double clockticks_per_second=CLK_TCK; //obsolete
#endif
  REPORT(clockticks_per_second);
  const char * const queuename=getenv("QUEUENAME");
#ifdef SX
  if(!queuename)
    FATAL_ERROR("Environment variable QUEUENAME not set.");
#endif
  double runtime=0;
  if(queuename){
    switch(queuename[0]){
    case 'A':
      runtime=eval_first("23*h");
      break;
    case 'S':
      runtime=eval_first("25*minutes");
      break;
    case 'M':
      runtime=eval_first("47*h");
      break;
    default:
      FATAL_ERROR("Unknown QUEUENAME.");
    }
  }
  const clock_t runtime_in_clockticks=
    clock_t(runtime/eval_first("second")*clockticks_per_second);
  // finished computing runtime
	
  {// What is this good for?? Some initialization??
    frozen_web w;
  }
#ifndef ON_SX5FSV
  set_new_handler(outOfMemory); 
#endif

#ifdef DO_FPE
  // floating-point exceptions
#ifdef __GNUC__
  feenableexcept( 0 | FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); 
#endif
#ifdef SX
  fpsetmask(fpgetmask() | FP_X_DZ | FP_X_INV | FP_X_OFL );
#endif
#endif

  signal_handling();
  REPORT(version);
  REPORT(compilation_time);
  system("/bin/date");
  
  frozen_web web;

  read_arguments_and_data(argc,argv,web);
  // read parameter file:


  REPORT(random_seed);
  set_random_seed(random_seed);
  double &t=web.current_time;
  web.initialize_for_integration();
  if(check_data_consistency)
    web.check_internal_consistency();

  if(do_print_each_step_anyway)
    print_each_step=1;

  try {
    //  double time_at_load=web.time_stamp;

    for(int test_number=0;test_number<number_of_tests;test_number++){
      web.speciate_only_fit(plant_fraction);

      repeat_after_t printout(time_between_printout);
      do_after_every_t check_it(check_interval);
      scheduler_t print_scheduler;
      print_scheduler.add(printout);
      if(!print_each_step) print_scheduler.add(check_it);
      
#if defined(SX) && defined(PARALLEL)
      //start parallel processes
      if(queuename[0]=='A'){
#pragma cdir reserve=4
      }else{
#pragma cdir reserve=4
      }
#endif
      REPORT(times(&dummy_time_struct));
      wallclock_when_entering_main_loop=times(&dummy_time_struct);
      
      //main loop:
      for(int once=1;once-->0;){
	if(web.number_of_species()==0){
	  std::cout << "all species extinct, exiting..." << std::endl;
	  throw exit_now;
	}else{ // scope of ode_state:
	  web.prepare_for_integration();
	  ODE_state ode_state(&web); 
	  //web.test_Jacobian();exit(1);
	  int old_prec=std::cout.precision();
	  std::cout.precision(19);
	  std::cout << "starting new integration at t = " << t/time_unit_for_output << std::endl;
	  std::cout.precision(old_prec);
	  TRACE(ode_state,MAIN);
	  
	  //cout << ode_state.size() << " compartments"<< std::endl;
	  ODE_vector start_state(ode_state);
	
	  double t_stop=web.current_time+max_test_duration;
	
	  while(t<t_stop){
	    if(check_it.due_now(t) || true || print_each_step){
	      // nothing to check
	    }
	    if(printout.due_now(t) || 
	       print_each_step || 
	       integrator_trouble_last_step ){
	      if(t<=observation_time){
		int old_prec=std::cout.precision();
		std::cout.precision(19);
		std::cout << fixed << std::setw(20) << t/time_unit_for_output  
			  << " : " ;
		std::cout.precision(old_prec);
		std::cout << ode_state[0] << endl;
	      }
	    }
	    double target_time=min(t_stop,print_scheduler.time_of_next_call());
	    int integrator_failure;
	    TRACE(t,ODE);
	    TRACE(target_time,ODE);
	    try{
	      if(true || print_each_step || integrator_trouble_last_step)
		integrator_failure=ode_state.integrate_one_step(target_time,t);
	      else
		integrator_failure=ode_state.integrate_until(target_time,t);
	    }catch(int i){
	      integrator_failure=1;
	      REPORT(integrator_failure);
	    }
	    TRACE(t,ODE);
	    TRACE(ode_state,ODE);
	    if(integrator_failure){
	      if(integrator_trouble_last_step){
		ode_state.diagnosis();
		WARNING("Integrator trouble again");
		throw exit_now;
	      }
	      integrator_trouble_last_step=true;
	      break;
	    }else{
	      integrator_trouble_last_step=false;
	    }
	    if(t>observation_time) 
	      break;
#ifdef SX
	    if(CPU_time_in_clockticks() > runtime_in_clockticks ){
	      WARNING("CPU time is up, exiting");
	      throw exit_now;
	    }
#endif
	  }// end of ODE loop
	  std::cout << "now " 
		    << CPU_time_in_clockticks()/clockticks_per_second/3600 
		    << " h CPU time used" << std::endl; 
	  if(t>0) ode_state.diagnosis();
	  if(exit_now) {
	    std::cout.flush();
	    throw exit_now;
	  }
	} // end scope of ode_state;
	if(t>observation_time){
	  throw exit_now;
	}
	if(check_data_consistency)
	  web.check_internal_consistency();
      }
      std::cout << (times(&dummy_time_struct)-wallclock_when_entering_main_loop)/
	clockticks_per_second << " seconds in main loop" << std::endl;
      cout << "test " << test_number << " " ;
      web.report_test_species();
      web.delete_test_species();
    }//number_of_tests
  }
  catch(int){
  }
#if defined(SX) && defined(PARALLEL)
#pragma cdir release
#endif
}

