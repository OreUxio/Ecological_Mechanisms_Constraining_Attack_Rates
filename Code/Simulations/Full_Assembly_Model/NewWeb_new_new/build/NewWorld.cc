// $Id: NewWorld.cc 2460 2016-01-15 16:30:43Z axel $
#include <string>
std::string version("$Id: NewWorld.cc 2460 2016-01-15 16:30:43Z axel $");
extern std::string compilation_time;

#include <fstream>
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
#include <sys/times.h> //for times :)
#include "otherwebs.h"
#include "relax.h"

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
//#include "linpack_eigen.h"
#include "xy_graph.h" 
#include "hotspots.h"
#include "polyfit.h"
#include "SortedMatrix.h"
#include "Integrator.h"
#include "packed_simulation.h"
#ifdef TRY_ASSEMBLER_CODE
#include "packed_simulation_asm.h"
#endif
#endif

using namespace std;

string in_file_name="NewWorld.cfg";
string web_file_name="";
ofstream matrices;
ofstream list_of_inserted_deleted_species;
ofstream list_of_inserted_deleted_biomass_species;
ofstream list_of_aggressivity_values; ///< file to which all aggressivity values are written
ofstream list_of_niche_width_values; ///< file to which all niche width values are written
ofstream list_of_max_animal_size_values; ///< file to which all max size values are written
ofstream list_of_noofplantspecies_values; ///< file to which all number of plant species values are written
ofstream list_of_Di_plant_values; ///< file to which all number of plant species values are written
ofstream list_of_noofanimalspecies_values; ///< file to which all number of animal species values are written
ofstream list_of_nooffishspecies_values; ///< file to which all number of fish species values are written
ofstream list_of_plant_body_mass_values; ///< file to which mean log10 plant body mass is written
ofstream list_of_animal_body_mass_values; ///< file to which mean log10 animal body mass is written
ofstream list_of_plant_biomass_values; ///< file to which mean log10 plant body mass is written
ofstream list_of_animal_biomass_values; ///< file to which mean log10 animal body mass is written

static my_evaluator_t eval_first;

// adjustable parameters:
static double check_interval=eval_first("1 * year"); ///< obsolete 
static double chunk_length=1.0; ///< obsolete
static double observation_time=eval_first("100 * years");
static double time_between_analysis=eval_first("1 * year");
static double analysis_start=eval_first("1 * year");
static double time_between_invasions=0; ///< obsolete
static double time_between_speciations=0; ///< obsolete
static double time_between_insertions=eval_first("3 * years");
static double time_between_printout=0;///< obsolete
static double time_between_saving=eval_first("100 * years");
static double time_between_re_initialize=eval_first("0");
static double time_unit_for_output=eval_first("1*year");
static double addition_fraction=0;///< fraction of additional species inserted
static int print_each_step=0;///< obsolete
static int random_seed=1;
static int number_of_next_save=0;
static int stop_after_one_period=0;///< obsolete
static int skip_also_chaos=0;///< obsolete
static int check_data_consistency=0;
static int no_animals=0;
static int plant_abundance_printout=0;
static int invade_from_other=0;
static int min_number_of_plants=0;
static int min_number_of_animals=0;
static int max_number_of_animals=0;
static int initial_number_of_plants=0;
static int initial_number_of_animals=0;
static double plant_animal_ratio=0;///< if nonzero, fix that ratio
static double min_relaxation_time=0;///< obsolete
static double print_each_step_runlimit=0;///< obsolete
static int chaos_detection_level=0;///< obsolete
static int measure_time_in_additions=0;///< obsolete
static int n_additions=0;
static double plant_richness_saturation_criterion=0;
static int auto_extinguishing=1;
static int invade_by_biomass=0;
static int speciate_selected_species=0;///< equal spec rates for all species
static int exit_on_large_plants=1;///< exit if some plants are larger than all animals
static int do_listing_after_each_iteration=0;///< set for fine-grained info
static double x=1; ///< variable to scale parameters in .cfg file

static double cpu_time_used=0;///< not really a parameter, but permanent

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(check_interval),
  CFGDOUBLE(chunk_length),
  CFGDOUBLE(observation_time),
  CFGDOUBLE(time_between_analysis),
  CFGDOUBLE(time_between_re_initialize),
  CFGDOUBLE(analysis_start),
  CFGDOUBLE(time_between_invasions),
  CFGDOUBLE(time_between_speciations),
  CFGDOUBLE(time_between_insertions),
  CFGDOUBLE(time_between_printout),
  CFGINT(print_each_step),
  CFGDOUBLE(time_between_saving),
  CFGDOUBLE(time_unit_for_output),
  CFGINT(random_seed),
  CFGINT(exit_on_large_plants),
  CFGINT(number_of_next_save),
  CFGINT(stop_after_one_period),
  CFGINT(skip_also_chaos),
  CFGINT(check_data_consistency),
  CFGINT(no_animals),
  CFGINT(invade_from_other),
  CFGINT(plant_abundance_printout),
  CFGINT(min_number_of_plants),
  CFGINT(min_number_of_animals),
  CFGINT(max_number_of_animals),
  CFGINT(initial_number_of_plants),
  CFGINT(initial_number_of_animals),
  CFGDOUBLE(plant_animal_ratio),
  CFGDOUBLE(min_relaxation_time),
  CFGDOUBLE(print_each_step_runlimit),
  CFGDOUBLE(addition_fraction),
  CFGDOUBLE(plant_richness_saturation_criterion),
  CFGINT(chaos_detection_level),
  CFGINT(measure_time_in_additions),
  CFGINT(n_additions), // this is not actually configuration but state
  CFGINT(auto_extinguishing),
  CFGINT(invade_by_biomass),
  CFGINT(do_listing_after_each_iteration),
  CFGDOUBLE(cpu_time_used),
  CFGINT(speciate_selected_species),
  CFGDOUBLE(x),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

bool has_suffix(const char * suff,const char * filename){
  return 0==strcmp(suff,filename+strlen(filename)-strlen(suff));
}

void make_predictions(){
  NewSpecies::make_predictions();
}

static bool do_make_predictions=false;
static bool guarantee_reproducibility=true;
static double multiply_biomasses=1;
static double multiply_aggressivities=1;
static bool set_aggressivities=false;
static double new_log10_aggressivities=0;
static double add_noise=0;
static int stop_after_saving=-1;
static const char * write_parameters_to_file=0;
static double previous_cpu_time_used=0;
#ifdef _SC_CLK_TCK
const double clockticks_per_second=sysconf(_SC_CLK_TCK);
#else 
const double clockticks_per_second=CLK_TCK; ///< obsolete
#endif
static bool weed_deletion=false;
static double max_wallclock_time=0;
static double one_second=eval_first("1*second");    
static vector< string > assignments;

void do_inline_assignments(){
  do_assignments(assignments);
}

void read_arguments_and_data(int argc,char *argv[],NewWeb & web){
  
  int c; 
  const char * formatstring="k:j::hvpm:a:A:P:N:S:qwW:z";
  while(1){
    c=getopt(argc, argv, formatstring);
    if (c == -1)
      break;
    
    switch(c){ 
    case 'v':
      printf("%s\n",version.c_str());
      exit(0);
      break;
    case 'k':
      assignments.push_back(optarg);
      break;
    case 'w':
      weed_deletion=true;
      break;
    case 'p':
      do_make_predictions=true;
      break;
    case 'P':
      write_parameters_to_file=optarg;
      if(!has_suffix(".cfg",write_parameters_to_file)){
	FATAL_ERROR("Parameter output file must end in .cfg");
      }	
      break;
    case 'S':
      stop_after_saving=atoi(optarg);
      break;
    case 'j':
      max_num_threads=
	(optarg ? atoi(optarg) : 1);
      break;
    case 'm':
      multiply_biomasses=atof(optarg);
      break;
    case 'a':
      multiply_aggressivities=atof(optarg);
      break;
    case 'A':
      set_aggressivities=true;
      new_log10_aggressivities=atof(optarg);
      break;
    case 'N':
      add_noise=atof(optarg);
      break;
    case 'q':
      guarantee_reproducibility=false;
      set_cfg_parameter("multithreading_only_fit","1");
      break;
    case 'W':
      max_wallclock_time=eval(optarg);
      break;
    case 'z':
      use_packed_simulation=true;
      break;
    case 'h':
    default:
	fprintf(stdout,"usage: %s [options] [foodweb_file|config_file...]\n\n\
options:\n\
-h            : This help message\n\
-m <factor>   : Multiply all biomasses by <factor> before starting\n\
-a <factor>   : Multiply all aggressivities by <factor> before starting\n\
-A <value>    : Set all aggressivities to 10^<value> before starting\n\
-N <sigma>    : Add noise with variance <sigma> to all log-biomasses\n\
-p            : Make analytic predictions and stop\n\
-P <filename> : Generate parameter files from given food web\n\
-q            : (quick) Avoid overhead to guarantee reproducibility\n\
-j[<n>]       : Set max number of threads to n (default n=1)\n\
-k <assign>   : assign value to parameter\n\
-S <number>   : Stop after saving web<number>.xml.bz2\n\
-v            : Print version info and exit\n\
-w            : Delete plants too large to eat (weed) at startup\n\
-W <time>     : Stop after wall-clock <time> (with unit) passed\n\
-z            : Use eXperimental, optimized core code\n\
",
		argv[0]);
	exit(c=='h'?0:1);
    }
  }
  if(set_aggressivities && multiply_aggressivities != 1){
    FATAL_ERROR("Cannot use both '-a' and '-A' option");
  }
    
  argc-=optind-1;
  argv+=optind-1;

  if(argc<=1){
    read_parameters_from_file(in_file_name);
    do_inline_assignments();
    return;
  }

  while(argc>1){
    if(has_suffix(".cfg",argv[1])){
      in_file_name=argv[1];
      read_parameters_from_file(in_file_name);
      do_inline_assignments();
    }else if(has_suffix(".xml",argv[1]) ||
	     has_suffix(".xml.gz",argv[1]) ||
	     has_suffix(".xml.bz2",argv[1]) ){
      web_file_name=argv[1];
      std::cout << "reading food web " << web_file_name << std::endl;
      XMLStore store(web_file_name);
      cfgPermanent parameters;
      store.get(&parameters,"adjustable_parameters");
      previous_cpu_time_used=cpu_time_used;

      do_inline_assignments();

      if(store.older_than("Aug 1 2006, 14:45:21")){
	WARNING("input web is from old code,");
	WARNING("emulating carrying-capacity bug!");
	emulate_carrying_capacity_bug=1;
      }
      if(store.older_than("Nov 25 2009, 19:20:00")){
	WARNING("input web is from old code,");
	WARNING("emulating animal body-mass-cutoff bug!");
	emulate_animal_body_mass_cutoff_bug=1;
      }
      if(store.older_than("Feb 8 2007, 23:59:00")){
	WARNING("input web is from old code,");
	WARNING("which had been measuring time in years, not additions!");
      }
      if(store.older_than("Nov 7 2010, 13:30:00")){
	if(get_cfg_parameter("totally_random_link_strength_sigma")){
	  WARNING("input web is from old code,");
	  WARNING("emulating log-random-matching bug");
	}
	emulate_log_random_matching_bug=1;
      }
    }else
      FATAL_ERROR("unknown file type: " << argv[1]);
    argv++;
    argc--;
  }
  if(write_parameters_to_file){
    write_cfg_template(write_parameters_to_file);
    exit(0);  
  }

  // Read the food-web at the very end, after all parameters are read.
  if(!web_file_name.empty()){
    XMLStore store(web_file_name);
    store.get(&web,"FoodWeb");
    if(store.older_than("Feb 8 2007, 23:59:00")){
      // this is to read older versions with separated simulatedTime
      try{
	store.get(web.current_time,"simulatedTime");
      } catch (XMLStoreException){
	WARNING("... but this variable is obsolete, anyway.");
      }
    }else{
      web.current_time=0;
    }
  }
  if(multiply_aggressivities!=1){
    for(int i=web.number_of_species();i-->0;){
      web.s(i).
	set_aggressivity_g(web.s(i).aggressivity_g()*
			   multiply_aggressivities);
    }
  }
  if(set_aggressivities){
    double g=pow(10,new_log10_aggressivities)*eval("meter^2/kilogram");
    for(int i=web.number_of_species();i-->0;){
      web.s(i).set_aggressivity_g(g);
    }
  }    
}

static void make_sure_no_outdated_parameters_are_used(){
#define FATAL_ON_DEPRECATE(X) if((X)) WARNING("parameter " #X " has no meaning anymore")
  FATAL_ON_DEPRECATE(check_interval);
  FATAL_ON_DEPRECATE(chunk_length);
  FATAL_ON_DEPRECATE(time_between_speciations);
  FATAL_ON_DEPRECATE(time_between_invasions);
  FATAL_ON_DEPRECATE(time_between_printout);
  FATAL_ON_DEPRECATE(print_each_step);
  FATAL_ON_DEPRECATE(stop_after_one_period);
  FATAL_ON_DEPRECATE(print_each_step_runlimit);
  FATAL_ON_DEPRECATE(skip_also_chaos);
  FATAL_ON_DEPRECATE(chaos_detection_level);
  FATAL_ON_DEPRECATE(min_relaxation_time);
  FATAL_ON_DEPRECATE(measure_time_in_additions);
#undef FATAL_ON_DEPRECATE
}

static bool integrator_trouble_last_step=false; //for integrator trouble

inline double CPU_time_in_clockticks(){
  static struct tms time_struct;
  times(&time_struct);
  return (time_struct.tms_utime+
	  time_struct.tms_stime+
	  time_struct.tms_cutime+
	  time_struct.tms_cstime+
	  0);
}

static void system_level_setup(){
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
}

static void perturb_web_if_requested(NewWeb & web){
  if(multiply_biomasses!=1 || add_noise!=0){
    REPORT(multiply_biomasses);
    REPORT(add_noise);
    REPORT(gaussian(0,add_noise));
    for(int i=web.number_of_species();i-->0;){
      REPORT(web.s(i).biomass_abundance_B());
      web.s(i).
	set_biomass_abundance_B(web.s(i).biomass_abundance_B()*
				multiply_biomasses*
				exp(gaussian(0,add_noise)) );
      REPORT(web.s(i).biomass_abundance_B());
    }
  }
}

static void connect_to_other_webs_maybe(NewWeb & web){
  if(invade_from_other && guarantee_reproducibility){
    FATAL_ERROR("otherwebs currently only works with the -q option set");
    // The problem is that otherwebs are not copied if a web is
    // copied, and without -q option, we do some web copying.
  }
  web.otherwebs.activate(Otherwebs::mode(invade_from_other));
  if(number_of_next_save>=1 &&
     web.otherwebs.current_mode()==Otherwebs::by_number){
    web.otherwebs.get_webs(number_of_next_save-1,exit_now);
  }
}

static double if_nan_default_to(double def, double x){
  if(!isnormal(x)){
    return def;
  }else{
    return x;
  }
}

static void do_listing(const NewWeb & web){
  static my_evaluator_t eval_here;
  const double unit_mass=eval_here("1*kilogram"); //for output files
  const double unit_area=eval_here("1*meter^2"); //for output files
  average_meter log10_aggressivity,log10_plant_bodymass,log10_animal_bodymass,
    log10_plant_biomass,log10_animal_biomass,niche_width;
  double max_size=(web.s.size()?web.s(0).log_mean_bodymass_M()/M_LN10:0);
  double plantB=0,plantB2=0;
  const double log10_kg=log(eval("1*kilogram"))/M_LN10;
  for(int i=web.s.size();i-->0;){
    if(web.is_plant(i)){
      log10_plant_bodymass.sample(web.s(i).log_mean_bodymass_M()/M_LN10);
      log10_plant_biomass.sample(log10(web.s(i).biomass_abundance_B()));
    }else{
      log10_aggressivity.sample(log(web.s(i).aggressivity_g()*unit_mass/unit_area )/M_LN10);
      niche_width.sample(web.s(i).niche_width_wt());
      max_size=max<double>(max_size,web.s(i).log_mean_bodymass_M()/M_LN10);
      log10_animal_bodymass.sample(web.s(i).log_mean_bodymass_M()/M_LN10);
      log10_animal_biomass.sample(log10(web.s(i).biomass_abundance_B()));
    }
  }
  const double default_log10_aggressivity=
    log10(get_cfg_parameter("initial_aggressivity")*unit_mass/unit_area);
  list_of_aggressivity_values 
    << n_additions << " " 
    << if_nan_default_to(default_log10_aggressivity,log10_aggressivity)
    << " " 
    << if_nan_default_to(0,log10_aggressivity.std())
    << " " 
    << number_of_next_save << endl;
  if(get_cfg_parameter("typical_niche_width_ratio")){
    list_of_niche_width_values
      << n_additions << " " 
      << if_nan_default_to(0,niche_width)
      << " " 
      << if_nan_default_to(0,niche_width.std())
      << " " 
      << number_of_next_save 
      << endl;
    list_of_niche_width_values.flush();
    if(get_cfg_parameter("trophic_niche_width_spread") &&
       niche_width.std() < 1e-6 ){
      exit_now=true;
    }
  }
  list_of_max_animal_size_values 
    << n_additions << " " 
    << if_nan_default_to(log10(get_cfg_parameter("Mini")),max_size
			 ) - log10_kg
    << " " 
    << number_of_next_save << " " 
    << if_nan_default_to(0,max_size-log10_plant_bodymass)
    << endl;
  list_of_max_animal_size_values.flush();
  list_of_noofplantspecies_values 
    << n_additions << " " 
    << web.number_of_plants() << " " 
    << number_of_next_save  << endl;
  list_of_noofplantspecies_values.flush();
  list_of_noofanimalspecies_values 
    << n_additions << " " 
    << web.number_of_animals() << " "
    << number_of_next_save  << " "
    << (web.number_of_plants() ?
	web.number_of_animals()/double(web.number_of_plants()) :
	1 ) << endl;
  list_of_noofanimalspecies_values.flush();
  list_of_nooffishspecies_values 
    << n_additions << " " 
    << web.number_of_fish() << " "
    << number_of_next_save  << " "
    << endl;
  list_of_nooffishspecies_values.flush();
  list_of_plant_body_mass_values
    << n_additions << " " 
    << if_nan_default_to(log10(get_cfg_parameter("Mini")),
			 log10_plant_bodymass ) - log10_kg 
    << " "
    << if_nan_default_to(0,log10_plant_bodymass.std())
    << " "
    << number_of_next_save 
    << endl;
  list_of_plant_body_mass_values.flush();
  list_of_animal_body_mass_values
    << n_additions << " " 
    << if_nan_default_to(log10(get_cfg_parameter("Mini")*
			       get_cfg_parameter("Mopt"))/2,
			 log10_animal_bodymass ) - log10_kg 
    << " "
    << if_nan_default_to(0,log10_animal_bodymass.std()) 
    << " "
    << number_of_next_save 
    << endl;
  list_of_animal_body_mass_values.flush();
  list_of_plant_biomass_values
    << n_additions << " " 
    << if_nan_default_to(0,log10_plant_biomass ) - log10_kg 
    << " "
    << if_nan_default_to(0,log10_plant_biomass.std()) 
    << " "
    << number_of_next_save 
    << endl;
  list_of_plant_biomass_values.flush();
  list_of_animal_biomass_values
    << n_additions << " " 
    << if_nan_default_to(0,log10_animal_biomass ) - log10_kg 
    << " "
    << if_nan_default_to(0,log10_animal_biomass.std()) 
    << " "
    << number_of_next_save 
    << endl;
  list_of_animal_biomass_values.flush();
}

static void save_current_state(NewWeb & web){
   
  if(!do_listing_after_each_iteration){
    do_listing(web);
  }

  // original part of code
  ostringstream filename;
  random_seed=random_integer(10000);
  set_random_seed(random_seed);
  filename << "web" << number_of_next_save++ << ".xml" << ".bz2";
  std::cout << "saving as " << filename.str() << std::endl;
  cpu_time_used=
    CPU_time_in_clockticks()/(clockticks_per_second/eval("second")) + 
    previous_cpu_time_used;
  cfgPermanent parameters;
  XMLStore to_save;
  to_save
    .put(version,"version")
    .put(compilation_time,"compilation_time")
    .put(&parameters,"adjustable_parameters")
    .put(&web,"FoodWeb")
    .save(filename.str().c_str());
  if(guarantee_reproducibility){
    // Reload web so that restarting from here does the same as just
    // continuing.
    web=NewWeb();
    to_save.get(&parameters,"adjustable_parameters");
    to_save.get(&web,"FoodWeb");
    web.initialize_for_integration();
  }
  if(web.otherwebs.current_mode()==Otherwebs::by_number){
    web.otherwebs.get_webs(number_of_next_save-1,exit_now);
  }
  // this is only good for interactive use, slow, and eats up memory:
#ifndef SX
  if(web.number_of_species()<50){
    Interaction_Matrix im=threshold(in_fraction(web.intake_matrix()),0.01);
//     if(im.Number_of_Species_S())
//       im.tsort().PPrint();
//     else
//       im.PPrint();
  }else{
    cout << "web too large to print" << endl;
  }
#endif
  if(exit_now)
    throw terminal_condition("exit_now set after saving");
}

  

/// Adds a species following various protocols, returns column of added species.
static NewSpecies * find_one_species_to_add(NewWeb & web){
  NewSpecies * new_species;
  if(unified_invasion_pressure_kappa){
    new_species=web.invade_or_speciate_only_fit(0.5);
  }else if(any_invasion_pressure_q){
    bool get_a_plant=(unirand()<plant_animal_ratio/(plant_animal_ratio+1)
		      || no_animals==1);
    if(invade_by_biomass){
      new_species=
	web.invade_or_speciate_any_by_biomass_only_fit(get_a_plant);
    }else{
      new_species=
	web.invade_or_speciate_any_only_fit(get_a_plant);
    }
  }else if(plant_animal_ratio){
    if(web.number_of_plants()<min_number_of_plants || no_animals){	    
      new_species=web.invade_only_fit(NewSpecies::plant);
    }else if(web.number_of_animals()<min_number_of_animals){
      new_species=web.invade_only_fit(NewSpecies::animal);
    }else{
      if(unirand()<plant_animal_ratio/(plant_animal_ratio+1)){
	if(speciate_selected_species){
	  int parent=web.number_of_animals()+
	    random_integer(web.number_of_plants());
	  new_species=web.speciate_selected_only_fit(parent);
	}else{
	  new_species=web.speciate_only_fit(NewSpecies::plant);
	}
      }else{
	if(speciate_selected_species){
	  int parent=random_integer(web.number_of_animals());
	  new_species=web.speciate_selected_only_fit(parent);
	}else{
	  new_species=web.speciate_only_fit(NewSpecies::animal);
	}
      }
    }
  }else{// no plant_animal_ratio
    new_species=
      web.invade_or_speciate_only_fit(double(web.number_of_plants())/
				      web.number_of_species());
  }
  return new_species;
}

/// Adds some species, returns set of added species.
static NewWeb::species_set_t add_species(NewWeb & web){
  int old_prec=std::cout.precision();
  std::cout.precision(2);
  cout << "currently " 
       << web.number_of_animals() << " animals and "
       << web.number_of_plants() << " plants ("
       << web.number_of_animals()/(double(web.number_of_plants())) << "), ";
  cout << n_additions << endl;
  std::cout.precision(old_prec);

  int n_new;
  n_new=int(addition_fraction*web.number_of_species()+1);
  
  vector< NewSpecies * > species_to_add(n_new);

  // We need to find all species to add first, because inserting one
  // deletes the steady_state saved in web.
  for(int i=n_new;i-->0;){
    species_to_add[i]=find_one_species_to_add(web);
  }
  // Observe that we also have to wait with building the species set
  // until after all specie s are added, because addition may
  // reshuffle indices into s.
  
  NewWeb::species_set_t inserted;
  list_of_inserted_deleted_species << "0.1 0.1 0.1 0.1 0.1" << endl;
  int old_n_additions;
  old_n_additions=n_additions; old_n_additions++;
  for(int i=n_new;i-->0;){
    int new_species=web.insert(species_to_add[i]);
    int col = web.assigned_column(new_species);
    std::cout << web.s(new_species).taxon_name()
	      << " " << col << " inserted" << std::endl;
    inserted.insert(col);
    n_additions++;
	list_of_inserted_deleted_species << web.s(new_species).unique_id() << " " << web.s(new_species).parent_unique_id() << " " << old_n_additions << " " << web.s(new_species).aggressivity_g() << " "  << n_additions << endl;
   }
   list_of_inserted_deleted_species << "0.2 0.2 0.2 0.2 0.2" << endl;

  if(check_data_consistency)
    web.check_internal_consistency(check_data_consistency);

  return inserted;
}

// Used in main():
struct main_prey_t {
  int column;
  double initial_biomass;
  double predator_aggressivity;
  double predator_niche_width;
  double predator_body_mass;
  double predator_D;
  NewSpecies::taxon_t prey_taxon;
};

int main(int argc,char *argv[]){
  
#if 0 // turn on to verify if my_isnan works properly
  REPORT(test_my_isnan(0.0/0.0));
  REPORT(test_my_isnan(1.0/0.0));
  REPORT(test_my_isinf(0.0/0.0));
  REPORT(test_my_isinf(1.0/0.0));
  REPORT(test_my_isnan(sqrt(-1.0)));
#endif

  // Used to stop this job when there are runtime restrictions.
  struct tms dummy_time_struct; // we do not read it
  const clock_t wallclock_start_time=times(&dummy_time_struct);
	
  system_level_setup();
  signal_handling();
  REPORT(version);
  REPORT(compilation_time);

  cout << "command_line =";
  for(int i=0;i<argc;i++){
    cout << " " << argv[i];
  }
  cout << endl;
    
  if(system("/bin/date")) WARNING("/bin/date command not found.");
  NewSpecies::initialize();

  NewWeb web;
  NewSpecies::sampling_species_list=&web.s;

  #ifdef SET_CPU_AFFINITY
  max_num_threads=number_of_CPUs();
  #endif
  
  read_arguments_and_data(argc,argv,web);

  make_sure_no_outdated_parameters_are_used();

  if(web.number_of_species()==0){
    for(int i=initial_number_of_plants;i-->0;){
      web.insert(web.invade(NewSpecies::plant));
    }
    web.relax(time_between_insertions,NewWeb::species_set_t(),true);
    for(int i=initial_number_of_animals;i-->0;){
      web.insert(web.invade(NewSpecies::animal));
    }
  }

  ///////////////// patch to generate unconstrainted random webs.
  // cfgPermanent parameters;
    // XMLStore()
    //   .put(version,"version")
    //   .put(compilation_time,"compilation_time")
    //   .put(&parameters,"adjustable_parameters")
    //   .put(&web,"FoodWeb")
    //   .save("tmp.xml.bz2");
    // exit(0);


  #ifdef SET_CPU_AFFINITY 
  CPU_set_list = find_my_CPUs();//.. and restrict max_num_threads if required
  #endif

  // For keeping track of wallclock time spent computing:
  clock_t wallclock_when_entering_main_loop; // we get this later
  REPORT(clockticks_per_second);

  perturb_web_if_requested(web);

  REPORT(random_seed);
  set_random_seed(random_seed);

  double &t=web.current_time;

  connect_to_other_webs_maybe(web);

  web.delete_forbidden_species();
  if(weed_deletion) web.delete_weeds();

  web.initialize_for_integration();

  if(check_data_consistency)
    web.check_internal_consistency(check_data_consistency);

  if(do_make_predictions){
    // Try to predict some properties without simulating.
    make_predictions();
    exit(0);
  }

  // Preparation for scheduling of things to do while simulating:
  const int insertions_between_analysis=
    int(time_between_analysis/time_between_insertions+0.5);
  const int insertions_between_saving=
    int(time_between_saving/time_between_insertions+0.5);
  const int insertions_between_re_initialize=
    int(time_between_re_initialize/time_between_insertions+0.5);

  try {
    // Schedule things to be done while simulating:
    scheduler_t scheduler;
    at_multiples_of_t analysis(insertions_between_analysis,n_additions);
    at_multiples_of_t re_initialize(insertions_between_re_initialize,
				    n_additions);
    scheduler.add(analysis);
    if(insertions_between_re_initialize)
      scheduler.add(re_initialize);
    scheduler_t saving_scheduler;
    at_multiples_of_t saving(insertions_between_saving,n_additions);
    saving_scheduler.add(saving);

    //Start parallel processes if required.
#if defined(SX) && defined(PARALLEL)
    if(queuename[0]=='A'){
#pragma cdir reserve=4
    }else{
#pragma cdir reserve=4
    }
#endif
	matrices.open("matrices.txt", std::ios_base::app);//std::ios_base::app
    list_of_inserted_deleted_species.open("list_of_inserted_deleted_species.dat",ios_base::app);
    list_of_inserted_deleted_biomass_species.open("list_of_inserted_deleted_biomass_species.dat",ios_base::app);
    list_of_aggressivity_values.open ("list_of_aggressivity_values.dat",ios_base::app); //file to which aggressivity values are written
    if(get_cfg_parameter("typical_niche_width_ratio"))
      list_of_niche_width_values.open ("list_of_niche_width_values.dat",ios_base::app); //file to which aggressivity values are written
    list_of_max_animal_size_values.open ("list_of_max_animal_size_values.dat",ios_base::app); //file to which max sizes values are written
    list_of_noofplantspecies_values.open ("list_of_noofplantspecies_values.dat",ios_base::app); //file to which number of plant species values are written
    list_of_noofanimalspecies_values.open ("list_of_noofanimalspecies_values.dat",ios_base::app); //file to which number of animal species values are written
    list_of_nooffishspecies_values.open ("list_of_nooffishspecies_values.dat",ios_base::app); //file to which number of fish species values are written
    list_of_plant_body_mass_values.open ("list_of_plant_body_mass_values.dat",ios_base::app); 
    list_of_animal_body_mass_values.open ("list_of_animal_body_mass_values.dat",ios_base::app); 
    list_of_plant_biomass_values.open ("list_of_plant_biomass_values.dat",ios_base::app); 
    list_of_animal_biomass_values.open ("list_of_animal_biomass_values.dat",ios_base::app); 
    REPORT(times(&dummy_time_struct));
    wallclock_when_entering_main_loop=times(&dummy_time_struct);
    double first_round=true;
    sequence<int> inserted_taxa_counts(NewSpecies::n_taxa);
    sequence<int> speciated_taxa_counts(NewSpecies::n_taxa);

    //main loop:
    while(n_additions <  observation_time/time_between_insertions){
      clock_t wallclock_start=times(&dummy_time_struct);
      if(do_listing_after_each_iteration and
	 web.number_of_plants()>0 ){
	do_listing(web);
      }
      if(saving.due_now(n_additions)){
	if(stop_after_saving>=0){
	  if(number_of_next_save==stop_after_saving-1){
	    // Turn full debugging on!
	    set_cfg_parameter("TRACEFLAG","-1");
	  }
	}
	save_current_state(web);
	if(guarantee_reproducibility){
	  // After re-loading, no species is added in the first
	  // iteration.  So we should not do so in the next iteration,
	  // either.
	  first_round=true;
	}
	if(stop_after_saving>=0){
	  if(number_of_next_save>stop_after_saving){
	    WARNING("stopping as requested");
	    exit(0);
	  }
	}
      }
      if(re_initialize.due_now(n_additions))
	web.initialize_for_integration();
      
      if(analysis.due_now(n_additions) && t > analysis_start){
	REPORT(n_additions);
	web.fast_stats();
	link_strength_matrix ls=
	  in_fraction(web.intake_matrix(false));
	strength_distribution(ls,1e-7,"/dev/null");
      }

      NewWeb::species_set_t inserted;
      if(first_round){
	first_round=false;
      }else{
	inserted=add_species(web);
      }

      // Keep number of animal species artifically low if requested:
      if(max_number_of_animals){
	while(web.number_of_animals() > max_number_of_animals){
	  web.delete_random_animal_species();
	}
      }

      const double t_stop=t+time_between_insertions;
      if(web.number_of_species()==0){
	std::cout << "all species extinct, skipping simulation..." << std::endl;
	t=t_stop;
      }else{ 

	web.delete_all_plants_with_negative_growth_rates(inserted);

	if(exit_on_large_plants && web.large_plants()){
	  system("echo Large plants in $PWD while running on `hostname`. Exiting.|mail -s 'Large plants' $USER");
	  throw terminal_condition("large plants!");
	}

	// Remember identity and taxa of parents:
	std::set<int> parent_unique_ids;
	for(NewWeb::species_set_t::const_iterator c=inserted.begin();
	    c!=inserted.end();
	    ++c){
	  parent_unique_ids.insert(web.s[web.species_at[*c]].parent_unique_id());
	  inserted_taxa_counts[web.s[web.species_at[*c]].taxon()]++;
	}

	// Remember main prey and niche width:
	std::vector< main_prey_t > main_prey(inserted.size());
	{
	  sequence< double > diet(web.s.size());
	  int predator_index=0;
	  for(NewWeb::species_set_t::const_iterator c=inserted.begin();
	      c!=inserted.end();
	      ++c){
	    int predator=web.species_at[*c];
	    if(!web.is_plant(predator)){
	      web.get_diet(predator,diet);
	      int main_prey_of_predator=
		max_element(diet.begin(),diet.end())-diet.begin();
	      main_prey[predator_index].column=
		web.assigned_column[main_prey_of_predator];
	      main_prey[predator_index].initial_biomass=
		web.s(main_prey_of_predator).biomass_abundance_B();
	      main_prey[predator_index].predator_niche_width=
		web.s(predator).niche_width_wt();
	      main_prey[predator_index].predator_aggressivity=
		web.s(predator).aggressivity_g();
	      main_prey[predator_index].predator_body_mass=
		web.s(predator).bodymass();
	      main_prey[predator_index].predator_D=1-sum(diet*diet);
	      main_prey[predator_index].prey_taxon=
		web.s[main_prey_of_predator].taxon();
	      predator_index++;
	    }
	  }
	  main_prey.resize(predator_index);
	}

	clock_t wall_simulation_start=times(&dummy_time_struct);
	double cpu_start=CPU_time_in_clockticks()/clockticks_per_second;
	list_of_inserted_deleted_biomass_species << "1 1 1 1" << endl;
	for(int i=web.s.size();i-->0;){
		list_of_inserted_deleted_biomass_species << web.s(i).unique_id() << " " << web.s(i).aggressivity_g()  << " " << web.s(i).biomass_abundance_B()  << " " << 0 << endl;
	}
	double simulated_time_relaxing=
	  web.relax(t_stop-t,inserted,auto_extinguishing);
	int n_plants = web.number_of_plants();
	int n_animals = web.number_of_animals();
	int index_k = 0;
	int index_i = 0;
	double H1H2 = 0;
	double H = 0;
	double max_H = 0;
	double curr = 0;
	for(int i=web.s.size();i-->0;){ //-n_plants
		list_of_inserted_deleted_biomass_species << web.s(i).unique_id() << " " << web.s(i).aggressivity_g()  << " " << web.s(i).biomass_abundance_B()  << " " << n_additions << endl;
		//~ index_i =  i; //+ n_plants
			// Remember to comment out previous declaration if using this one
			//std::string nameOfFile;
			//stringstream stream;
			//stream << n_additions;
			//nameOfFile = stream.str();
			//std::string nameOfFile1 = nameOfFile + "_matrix.txt";
			//nameOfFile += ".txt";
			//std::ofstream matrices;
			//matrices.open(nameOfFile.c_str(), std::ios_base::app);//std::ios_base::app
			
			//~ matrices << n_additions << " " << web.s(i).unique_id() << " " << web.s(i).aggressivity_g() << " ";
			//~ for(int k=web.s.size();k-->0;){
				//~ index_k = k;
				//~ H1H2 = 0;
				//~ for(int j=(web.s.size());j-->0;){
					//~ H1H2+= web.precomputed[index_i].c[j]*web.precomputed[index_k].c[j];
				//~ }
			
			//~ matrices << H1H2 << " " ;
			//~ }
			
			//~ curr =0; H=0; max_H=0;
				//~ for(int j=web.s.size();j-->0;){
					//~ curr = web.precomputed[index_i].c[j];
					//~ H+= curr;
					//~ if(max_H<curr){
						//~ max_H=curr;
					//~ }
				//~ }
			//~ matrices << H << " " << max_H << " " << n_plants<< " " << n_animals << endl;
	}
	
	// For checking
				//~ std::string nameOfFile;
			//~ stringstream stream;
			//~ stream << n_additions;
			//~ nameOfFile = stream.str();
			//~ std::string nameOfFile1 = nameOfFile + "_matrix.txt";
	//~ std::ofstream during_interaction_matrix;
	//~ during_interaction_matrix.open(nameOfFile1.c_str(),ios_base::app);
	//~ for(int i=web.s.size();i-->0;){
		//~ for(int j=web.s.size();j-->0;){
			//~ if(j!=0){during_interaction_matrix << web.precomputed[i].c[j] << "," ;}else{during_interaction_matrix << web.precomputed[i].c[j] << endl;}
	  //~ }
	//~ }
	// End checking
	REPORT(simulated_time_relaxing);
	double cpu_stop=CPU_time_in_clockticks()/clockticks_per_second;
	clock_t wall_simulation_stop=times(&dummy_time_struct);
	std::cout << "now " 
		  << cpu_stop/3600 +  previous_cpu_time_used/eval("hour")
		  << " h CPU time used, " 
		  << (cpu_stop-cpu_start)/simulated_time_relaxing*eval("1000000*year")/
	  web.number_of_species()
		  << " s/MY/species" << std::endl; 
	std::cout << "walltime " 
		  << (wall_simulation_stop-wall_simulation_start)/
	  clockticks_per_second/ 
	  simulated_time_relaxing*eval("1000000*year")/
	  web.number_of_species()
		  << " s/MY/species" << std::endl;
	// 	if(t>0) ode_state.short_diagnosis();
	// 	if(exit_now) {
	// 	  std::cout.flush();
	// 	  throw terminal_condition("exit_now set after relax");
	// 	}
	
	// see which parents are still there:
	for(int i=web.s.size();i-->0;){
	  if(parent_unique_ids.count(web.s[i].unique_id()) > 0){
	    speciated_taxa_counts[web.s[i].taxon()]++;
	  }
	}
	sequence<double> speciation_fractions=speciated_taxa_counts;
	speciation_fractions/=sequence< double >(inserted_taxa_counts);
	// Speciation_fractions tells us how often a new invadere does
	// NOT replace its parent after invasion.
	REPORT(speciation_fractions);

	if(log_steady_state_stats){
	  // see how much main prey abundance changed:
	  for(int i=main_prey.size();i-->0;){
	    cout << "war";
	    cout << " " << main_prey[i].predator_aggressivity ;
	    cout << " " << main_prey[i].predator_niche_width ;
	    cout << " " << main_prey[i].predator_D;
	    cout << " ";
	    int prey=web.species_at[main_prey[i].column];
	    if(prey==NewWeb::column_unused){
	      cout << 0 ;
	    }else{
	      cout << 
		web.s(prey).biomass_abundance_B()/
		main_prey[i].initial_biomass;
	    }
	    cout << " " << int(main_prey[i].prey_taxon);
	    cout << " " << main_prey[i].predator_body_mass ;
	    cout << endl;
	  }
	} 
      }
      if(check_data_consistency)
        web.check_internal_consistency(check_data_consistency);
      clock_t wallclock_stop=times(&dummy_time_struct);
      std::cout << "time for one iteration of main loop = " 
		  << (wallclock_stop-wallclock_start)/clockticks_per_second 
		  << " s" << std::endl;
      
      // stop if wallclock limit exceeded
      if(max_wallclock_time){
	if((wallclock_stop-wallclock_start_time)/clockticks_per_second > 
	   max_wallclock_time / one_second ){
	  WARNING("wallclock_time exceeded " << max_wallclock_time/one_second << " seconds");
	  WARNING("stopping");
	  exit_now=true;
	}
      }
      if(n_additions > plant_richness_saturation_criterion){
    	WARNING("stopping by plant_richness_saturation_criterion");
    		exit_now=true;
    }
}
    std::cout << "observation time exceeded" << endl;
    std::cout << (times(&dummy_time_struct)-wallclock_when_entering_main_loop)/
    clockticks_per_second << " seconds in main loop" << std::endl;
  }
  // Save current web as SAVED.xml in case something went wrong:
  catch(terminal_condition c){
    std::cout << "caught " << c << std::endl;
    const char * outfile_name="SAVED.xml";
    std::cout 
      << (times(&dummy_time_struct)-wallclock_when_entering_main_loop)/
      clockticks_per_second << " seconds in main loop" << std::endl;
    cfgPermanent parameters;
    XMLStore()
      .put(version,"version")
      .put(compilation_time,"compilation_time")
      .put(&parameters,"adjustable_parameters")
      .put(&web,"FoodWeb")
      .save(outfile_name);
    std::cerr << "saved web as " << outfile_name << std::endl; 
    exit(0);
  }
  // Stop parallel processes if required
#if defined(SX) && defined(PARALLEL)
#pragma cdir release
#endif
  exit(0);
}
