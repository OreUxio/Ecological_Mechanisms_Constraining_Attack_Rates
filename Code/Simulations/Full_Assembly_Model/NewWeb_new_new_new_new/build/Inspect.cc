// $Id: Inspect.cc 2489 2016-12-21 13:31:36Z axel $
#include <string>
std::string version("$Revision: 2489 $");

extern std::string compilation_time;

#include <signal.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
// to_string example


//#define DO_FPE  // do floating point exceptions
#ifdef DO_FPE
#ifdef __GNUC__
#include <fenv.h>  // floating point exceptions
#endif
#ifdef SX
#include <ieeefp.h>  // floating point exceptions
#endif
#endif

#include "hotspots.h"
#include "NewSpecies.h"
#include "parsecfg.h" // include this header file when use parsecfg lib
#include "error.h"    // my error handlers
#include "XMLStore.h"
#include "cfgPermanent.h"
#include "matrix_transformers.h"
#include "standardize.h"
#include "nls_web.h"
#include "random.h"
#include "time_average.h"
#include "tlf.h"
#include "snapshot.h"
#include "perturb.h"
#include "relax.h"
#include "three_column_files.h"
#include "mass_integrator.h"
#include "Integrator.h"
#include "MSY.h"


#ifdef GET_INDIRECT_DEPENDENCIES
#include "Statistics.h"
#include "NewWeb.h"
#include "evaluate.h"
//#include "linpack_eigen.h"
#include "period_cutter.h"
#include "xy_graph.h"
#include "hotspots.h"
#include "polyfit.h"
#include "packed_simulation.h"
#ifdef TRY_ASSEMBLER_CODE
#include "packed_simulation_asm.h"
#endif
#endif

using namespace std;
ofstream fitness_only;
ofstream interaction_matrix;
ofstream agg_for_matrix;
ofstream during_agg_for_matrix;
ofstream during_interaction_matrix;
ofstream during_agg_for_matrix1;
ofstream during_interaction_matrix1;
ofstream before_agg_for_matrix;
ofstream before_interaction_matrix;
ofstream after_agg_for_matrix;
ofstream after_interaction_matrix;
ofstream after_agg_for_matrix1;
ofstream after_interaction_matrix1;
ofstream biomass_interaction_matrix;
ofstream list_of_inserted_deleted_biomass_species_w;
ofstream list_of_inserted_deleted_biomass_species_ww;
ofstream list_of_inserted_deleted_biomass_species_www;

static int random_seed=1;
static int number_of_next_save=0;
static double abundance_multiplier=1;
static bool needs_flows=false;
static int coldelS=0;

 
// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGINT(random_seed),
  CFGINT(number_of_next_save),
  CFGDOUBLE(abundance_multiplier),
  CFGINT(coldelS),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

enum file_type_t {xml,nls,tlf,dot_web};
static file_type_t file_type=xml;

bool has_suffix(const char * suff,const char * filename){
  return 0==strcmp(suff,filename+strlen(filename)-strlen(suff));
}

static my_evaluator_t eval_now;

static string in_file_name="";
static string web_file_name="";
static string dot_graph_file="";
static string pajek_graph_file="";
static bool use_c=false;
static const double default_threshold_value=0.01;
static double threshold_value=default_threshold_value;
static bool print_properties=false;
static string graph_file_name="graph.eps";
static bool standardized=false;
static bool lump=false;
static bool print=false;
static bool print_flows=false;
static bool lump_lowest_level=false;
static bool get_msy=false;
static char* rank_abundance_plot=0;
static char* size_spectrum_file=0;
static char* strength_spectrum_file=0;
static char* species_table=0;
static char* species_table2=0;
static char* species_table3=0;
static char* link_table=0;
static char* link_table2=0;
static char* pgm_file=0;
static char* eigenvalue_file=0;
static char* generality_file=0;
static char* biomass_spectrum_file=0;
static char* vulnerability_file=0;
static char* distance_file=0;
static vector< string > assignments;
static bool print_insertion_count;
static bool choose_other=false; // for debugging
static double time_averaging_time=0;
static bool fast_stats=false;
static bool report_pressures=false;
static int iterations_for_pressures=100;
static double nls_area=0;
static char* log_availability_structure_file=0;
static char* rectification_base_name=0;
static char* perturbation_type=0;
static char* waiting_time=0;
static char* writeback_file=0;
static char* invasion_fitness_curve=0;
static int fitness_samples=0;
static char* link_coefficient_matrix_file=0;
static char* competition_file=0;
static bool print_indices=false;
static double min_size_of_things_we_care_about=0;
static bool jolly=false;
static double mass_integral_start=0;
static double mass_integral_start_LSI=0;
static bool biodiv_indices_start=false;
static double all_indices_start=0;
static bool scale_flows=false;
static bool rate_of_change_indices=false;
static bool snap_to_fixed_point=false;
static bool forbid_fishing=false;
static char* fished_web_name=0;
static bool enforce_equilibrium=false;
static char* aggressivenes_file=0;
//~ static bool aggressivenes_file=false;


void do_inline_assignments(){
  do_assignments(assignments);
}


void print_usage(char ** argv,const char * formatstring){
  fprintf(stdout,"usage: %s [options] foodweb_file\n\n\
options (some work only with special input file types):\n\
-a [<time>]  : time average over non-constant states for at most <time>\n\
-A <area>    : set area for SCOR files, e.g. \"10*ha\"\n\
-b           : print connection matrix\n\
-B <filename>: generate a biomass spectrum\n\
-c           : use potential link strength rather than realized link strength\n\
-C <filename>: output log availability structure to file\n\
-d <filename>: generate a dot file for graphs\n\
-D <basename>: write neighbors-within-distance curves to files ending in basename\n\
-e           : output invasion fitness curve\n\
-E <filename>: write eigenvalues of Jacobian ( * -1) to file\n\
-f           : print strengths of connections\n\
-F <n>       : strengths of evolutionary pressures using <n> iterations\n\
-g [<name>]  : generate web graph (dysfunctional)\n\
-G <filename>: write cumulative generality distribution to file\n\
-H <filename>: generate a dat file to draw a food-web graph using Pajek\n\
-i           : print the number of insertions into the web\n\
-I           : print macroecological ecological indices\n\
-I <size>    : restrict some indices to size larger than this\n\
-j           : jolly: don't do things that are computationally heavy\n\
-J           : write dynamics of variables related to rate of change of LFI, LSI and biodiv indices\n\
-k <assign>  : assign value to parameter\n\
-K <filename>: generate PGM image of competition matrix\n\
-l           : lump lowest level\n\
-L <filename>: generate a table of link properties\n\
-m           : print some fast statistics (e.g. on bodymasses)\n\
-M <size>    : write LFI dynamics with size as the threshold\n\
-n <filename>: write link-strength matrix to file\n\
-N <size>    : write LSI dynamics with size as the threshold\n\
-o           : analyze state at fix point\n\
-O           : write dynamics for biodiversity indices\n\
-p           : compute and print topological properties\n\
-P <filename>: generate a PGM file of the connection matrix\n\
-q <filename>: generate link properties \n\
-Q <size>    : write dynamics for LFI, LSI and biodiversity indices\n\
-r <filename>: generate a rank-abundance plot (affected by -l)\n\
-R <filename>: rectify the interaction matrix, write into various files\n\
-s           : standardize topology\n\
-S <filename>: write diet partitioning function (DPF) to file\n\
-t <value>   : set threshold value for link assignment (default %g)\n\
-T <filename>: generate a table of species properties (affected by -l)\n\
-u <filename>: generate a table of properties for fish species, together with supporting tables (affected by -l)\n\
-U <filename>: generate a table of properties for fish species, together with supporting tables (affected by -l)\n\
-v           : print version information\n\
-V <filename>: write cumulative vulnerability distribution to file\n\
-w <time>    : wait <time> to relax before averaging\n\
-W <filename>  write relaxed state to <filename>\n\
-x           : lump trophic species for topology\n\
-X <type>    : apply perturbation <type> before waiting and averaging\n\
               (use type \"help\" for help on types)\n\
-y           : Do MSY computation (experimental)\n\
-Y           : divide flows by resource turnover rate before using\n\
-z [<ncpus>] : user experimental core code (as -z in NewWorld)\n\
-Z <filename>: generate a size spectrum\n\
-0           : set fishing mortalities to zero\n\
-3 <n>       : generate |n| samples of animal inv. fitness (n<0 for plants)\n\
-4           : enforce equilibrium state by adjusting turnover rates\n\
-5 <webfile> : copy fishing mortalities from web file\n\
",argv[0],double(default_threshold_value));
}

void read_arguments(int argc,char *argv[]){
  int c;
  const char * formatstring= 
    "k:K:Q:N:M:n:ma::B:E:L:q:G:V:P:T:u:U:D:ihvg::sxbfdH:ct:plr:Z:S:12yYA:C:R:X:w:W:I::jF:e:OoJz::3:05:4";

  while(1){
    c=getopt(argc, argv, formatstring);
    if (c == -1)
      break;
    
    switch(c){
    case 'v':
      printf("Main code version: %s\n",version.c_str());
      printf("Compilation time:  %s\n",compilation_time.c_str());
      printf("Copyright (C) 2007,2008 Axel G. Rossberg\n");
      printf("License GPL 3\n");
      exit(0);
      break;
    case 'R':
      rectification_base_name=optarg;
      needs_flows=true;
      break;
    case 's':
      standardized=true;
      break;
    case 'x':
      lump=true;
      break;
    case 'm':
      fast_stats=true;
      break;
    case 'M':
      mass_integral_start=eval(optarg);
      break;
    case 'N':
      mass_integral_start_LSI=eval(optarg);
      break;
    case 'o':
      snap_to_fixed_point=true;
      break;
    case 'O':
      biodiv_indices_start=true;
      break;
    case 'Q':
      all_indices_start=eval(optarg);
      break;
    case 'J':
      rate_of_change_indices=true;
      break;
    case 'l':
      lump_lowest_level=true;
      break; 
    case 'b':
      print=true;
      needs_flows=true;
      break;
    case 'y':
      get_msy=true;
      break;
    case 'Y':
      scale_flows=true;
      needs_flows=true;
      break;
    case 'i':
      print_insertion_count=true;
      break;
    case 'I':
      print_indices=true;
      needs_flows=true;
      if(optarg)
	min_size_of_things_we_care_about=eval(optarg);
      break;
    case 'j':
      jolly=true;
      break;
    case 'f':
      print_flows=true;
      needs_flows=true;
      break;
    case 'F':
      report_pressures=true;
      iterations_for_pressures=eval(optarg);
      break;
    case 'P':
      pgm_file=optarg;
      needs_flows=true;
      break;
    case 'C':
      log_availability_structure_file=optarg;
      break;
    case 'E':
      eigenvalue_file=optarg;
      break;
    case 'e':
      invasion_fitness_curve=optarg;//(char *)"*";
      break;
    case 'G':
      generality_file=optarg;
      needs_flows=true;
      break;
    case 'V':
      vulnerability_file=optarg;
      needs_flows=true;
      break;
    case 'd':
      dot_graph_file=optarg;
      needs_flows=true;
      break;
    case 'r':
      rank_abundance_plot=optarg;
      break;
    case 'Z':
      size_spectrum_file=optarg;
      break;
    case 'z':
      use_packed_simulation=true;
      if(optarg){
	max_num_threads=atoi(optarg);
      }
      break;
    case 'B':
      biomass_spectrum_file=optarg;
      break;
    case 'D':
      distance_file=optarg;
      break;
    case 'S':
      strength_spectrum_file=optarg;
      needs_flows=true;
      break;
    case 'K':
      competition_file=optarg;
      break;
    case 'k':
      assignments.push_back(optarg);
      break;
    case 'T':
      species_table=optarg;
      needs_flows=true;
      break;
    case 'u':
      species_table2=optarg;
      needs_flows=true;
      break;
    case 'U':
      species_table3=optarg;
      needs_flows=true;
      break;
    case 'X':
      perturbation_type=optarg;
      break;
    case 'w':
      waiting_time=optarg;
      break;
    case 'W':
      //~ writeback_file=optarg;
      aggressivenes_file=optarg;
      needs_flows=true;
      break;
    case 'L':
      link_table=optarg;
      needs_flows=true;
      break;
    case 'q':
      link_table2=optarg;
      needs_flows=true;
      break;
    case 'n': 
      link_coefficient_matrix_file=optarg;
      break;
    case 'a':
      if(optarg)
	time_averaging_time=eval(optarg);
      else
	time_averaging_time=eval("1*year");
      REPORT(time_averaging_time);
      break;
    case 'c':
      use_c=true;
      break;
    case 'p':
      print_properties=true;
      needs_flows=true;
      break;
    case '5':
      fished_web_name=optarg;
      break;
    case '4':
      enforce_equilibrium=true;
      break;
    case '3':
      fitness_samples=atof(optarg);
      break;
    case '2':
      choose_other=true;
      break;
    case '1':
      choose_other=false;
      break;
    case '0':
      forbid_fishing=true;
      break;
    case 't':
      threshold_value=atof(optarg);
      break;
    case 'A':
      nls_area=eval(optarg);
      break;
    case 'H':
      pajek_graph_file=optarg;
      needs_flows=true;
      break;
    case 'h':
    default:
      print_usage(argv,formatstring);
      exit(c=='h'?0:1);
    }
  }
  argc-=optind-1;
  argv+=optind-1;
  
  while(argc>1){
    if(has_suffix(".cfg",argv[1])){
      in_file_name=argv[1];
    }else if(has_suffix(".dat",argv[1]) ||
	     has_suffix(".DAT",argv[1]) ||
	     has_suffix(".tlf",argv[1]) ||
	     has_suffix(".web",argv[1]) ||
	     has_suffix(".xml",argv[1]) ||
	     has_suffix(".xml.gz",argv[1]) ||
	     has_suffix(".xml.bz2",argv[1]) ){
      web_file_name=argv[1];
      if(has_suffix(".dat",argv[1])||has_suffix(".DAT",argv[1])){
	// SCOR files
	file_type=nls;
      }else if(has_suffix(".tlf",argv[1])){
	// custom "Tuesday Lake" format 
	file_type=tlf;
      }else if(has_suffix(".web",argv[1])){
	file_type=dot_web;
      }
    }else
      FATAL_ERROR("unknown file type: " << argv[1]);
    argv++;
    argc--;
  }
  if(web_file_name.size()==0){
    print_usage(argv,formatstring);
    exit(1);
  }
}


string order_of_magnitude(double x){
  int y=int(-log(x)/log(10.0));
  ostringstream s;
  if(y>=6 || x==0){
    s << ".";
  }else if(y<=-10){
    s << "X";
  }else{
    s << y;
  }
  return s.str();
}    

extern "C" {
  void openblas_set_num_threads(int num_threads);
}

int main(int argc,char *argv[]){
  std::set_new_handler(outOfMemory); 
  signal_handling();

#ifdef DO_FPE
  // floating-point exceptions
#ifdef __GNUC__
  feenableexcept( 0 | FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW ); 
#endif
#ifdef SX
  fpsetmask(fpgetmask() | FP_X_DZ | FP_X_INV | FP_X_OFL );
#endif
#endif

#ifdef SET_CPU_AFFINITY 
  max_num_threads=number_of_CPUs();
#endif
  
  read_arguments(argc,argv);

  do_inline_assignments();

#ifdef SET_CPU_AFFINITY 
  CPU_set_list = find_my_CPUs();//.. and restrict max_num_threads if required
#endif

  openblas_set_num_threads(max_num_threads);
  
  // read parameter file for the first time:
  set_cfg_parameter("DEFAULT_RELATIVE_TOLERANCE","0.001");
  if(in_file_name[0])
    if (cfgParse(in_file_name.c_str(), full_cfg_list(), CFG_SIMPLE) == -1)
      FATAL_ERROR("error reading parameter file");
  REPORT(random_seed);
  set_random_seed(random_seed);
  NewSpecies::initialize();
  
  cfgPermanent parameters;
  NewWeb web,baseline_web;
  link_strength_matrix ls;
  sequence<double> biomass_B;

  nls_web nlsweb;
  Snapshot web_data;
  
  switch(file_type){
  case xml:
    {// NewWeb simulation result
      XMLStore store(web_file_name);
      store.get(&parameters,"adjustable_parameters");
      do_inline_assignments();
      // read parameter file again to overwrite if desired!
      if(in_file_name[0])
	if (cfgParse(in_file_name.c_str(), full_cfg_list(), CFG_SIMPLE) == -1)
	  FATAL_ERROR("error reading parameter file");
      REPORT(random_seed);
      set_random_seed(random_seed);

      if(store.older_than("Aug 1 2006, 14:45:21")){
	WARNING("input web is from old code,");
	WARNING("emulating carrying-capacity bug!");
	emulate_carrying_capacity_bug=1;
      }
      if(store.older_than("Nov 7 2010, 13:30:00")){
	if(get_cfg_parameter("totally_random_link_strength_sigma")){
	  WARNING("input web is from old code,");
	  WARNING("emulating log-random-matching bug");
	}
	emulate_log_random_matching_bug=1;
      }

      if(get_cfg_parameter("totally_random_link_strength_sigma")){
	double as_ratio=
	  -log(get_cfg_parameter("conversion_efficiency_epsilon")*
	       get_cfg_parameter("attack_rate_allometric_prefactor")/
	       get_cfg_parameter("respiration_rate_allometric_prefactor"))/
	  get_cfg_parameter("totally_random_link_strength_sigma");
	double Sr_guess=
	  exp(pow(as_ratio,2)/2)*as_ratio*sqrt(2*M_PI)+1;
	REPORT(Sr_guess);
	double ic_guess=sqrt(2*log(Sr_guess))/
	  get_cfg_parameter("totally_random_link_strength_sigma");
	REPORT(ic_guess);
      }

      store.get(&web,"FoodWeb");
      
      for(int i=web.s.size();i-->0;){
	web.s(i).
	  set_biomass_abundance_B(web.s(i).biomass_abundance_B()*
				  abundance_multiplier );
      }

      web.initialize_for_integration();
      web.compute_diagnostics=true;

      if(enforce_equilibrium){
	web.enforce_equilibrium();
      }
      if(fished_web_name){
	NewWeb fished_web;
	XMLStore(fished_web_name).get(&fished_web,"FoodWeb");
	web.copy_fishing_mortalities_from(fished_web);
      }
      
      if(print_indices){
	baseline_web=web;
	if(time_averaging_time){
	  time_average avr(baseline_web,time_averaging_time); 
	}else{
	  baseline_web.get_top_down_strength();
	}
	if(snap_to_fixed_point){
	  FATAL_ERROR("fixed point snapping and index calculation cannot be combined, yet");
	}
      }

      if(forbid_fishing){
	web.forbid_fishing();
      }

      if(perturbation_type){
	perturb(web,perturbation_type);
      }

      //~ if(writeback_file){
	//~ cfgPermanent parameters;
	//~ XMLStore()
	  //~ .put(version,"version")
	  //~ .put(compilation_time,"compilation_time")
	  //~ .put(&parameters,"adjustable_parameters")
	  //~ .put(&web,"FoodWeb")
	  //~ .save(writeback_file);
      //~ }
      
      
      do_inline_assignments();
      

  //~ FATAL_ERROR("EXPERIMENTAL CODE HERE");

      if(fast_stats){
	web.fast_stats();
      }

      //web.initialize_for_integration();
      web.compute_diagnostics=true;
      if(time_averaging_time and not snap_to_fixed_point){
	time_average avr(web,time_averaging_time); // here the
	// averages are
	// computed!!
	if(needs_flows) ls=web.fx;
	biomass_B.resize(web.biomass_B.size());
	std::copy(web.biomass_B.begin(),web.biomass_B.end(),biomass_B.begin());
      }else{
	if(snap_to_fixed_point){
	  if(time_averaging_time){
	    time_average avr(web,time_averaging_time);
	  }
	  // Snap web to fixed point:
	  {	
	    WARNING("snapping");
	    fixed_point_analyzer snapper(&web,-10);
	    snapper.set_time_scale_to(1e-5);
	    snapper.snap_to_fixed_point();
	  }
	  cfgPermanent parameters;
	  XMLStore()
	    .put(version,"version")
	    .put(compilation_time,"compilation_time")
	    .put(&parameters,"adjustable_parameters")
	    .put(&web,"FoodWeb")
	    .save("snapped.xml.bz2");
	}
	if(needs_flows){
	  web.get_top_down_strength();// called for side effects.
	  if(scale_flows)
	    ls=web.scaled_intake_matrix(not jolly);//don't precompute if jolly
	  else
	    ls=web.intake_matrix(not jolly);//don't precompute if jolly
	}
	biomass_B=web.get_biomass_B();
	web.get_top_down_strength();// called for side effects.
      }
      if(use_c && needs_flows){
	ls=web.matching_matrix();
      }
      
   //~ for(int i=10;i-->0;){
    //~ NewSpecies * s = web.speciate(NewSpecies::animal);
    //~ s->set_aggressivity_g(0.06); // you need to adjust this!
    //~ REPORT(web.invasion_fitness(s));
  //~ }
  //~ FATAL_ERROR("EXPERIMENTAL CODE HERE");
  
        if(waiting_time){
	//~ web.compute_diagnostics=false;
	//~ cout << "relaxing ... ";
	//~ cout.flush();
	//~ double T=eval(waiting_time);
	//~ if(T>0){
	  //~ double simulated_time_relaxing=
	    //~ web.relax(T);
	  //~ REPORT(simulated_time_relaxing);
	//~ }else{
	  //~ T=-T;
//~ #if 1
	  //~ double simulated_time_relaxing=
	    //~ web.relax(T,NewWeb::species_set_t(),false);	  
	  //~ REPORT(simulated_time_relaxing);
//~ #else // experimental code to make movies
	  //~ int frame_number=0;
	  //~ double t=0;
	  //~ for(double next_t=0.001;next_t<T;next_t*=1.1 ){
	    //~ biomass_B=web.get_biomass_B();
	    //~ ostringstream filename;
	    //~ filename  << "frame" << setfill('0') << setw(5) 
		      //~ << frame_number << ".dat";
	    //~ Mass_Integrator().spectrum(1e-12*web.get_biomass_B(),
				       //~ web.get_bodymass_M(),
				       //~ filename.str().c_str());
	    //~ REPORT(frame_number);
	    //~ double simulated_time_relaxing=
	      //~ web.relax(next_t-t,
			//~ NewWeb::species_set_t(),false);	  	    
	    //~ REPORT(simulated_time_relaxing);
	    //~ t=next_t;
	    //~ REPORT(t);
	    //~ frame_number++;
	  //~ }
//~ #endif // experimental code to make movies
	//~ }
	//~ cout << "finished" << endl;
		//~ ifstream file(aggressivenes_file);
		//~ string line;
		//~ getline(file, line);
		//~ stringstream geek(line);
		//~ double agg = 0; geek >> agg;
		//~ std::fstream myfile(waiting_time, std::ios_base::in);
		//~ float a;
		//~ while (myfile >> a){
			//~ printf("%f ", a);
		//~ }
		//~ getchar();
		std::string line;
		ifstream f (waiting_time);
		fitness_only.open("fitness_only.dat",ios_base::app);
		double agg = 0;
		while (std::getline(f, line)) {
			stringstream ss(line);
			ss >> agg;
			cout << agg << endl;
			
			
			//~ NewSpecies * new_species = web.speciate_only_fit(NewSpecies::animal);
			NewSpecies * new_species = web.speciate(NewSpecies::animal);
			new_species->set_aggressivity_g(agg);
			web.initialize_for_integration();
			double fitness = web.invasion_fitness(new_species);
			fitness_only << agg << " " << number_of_next_save << " " << random_seed << " "<< fitness << endl;
			//~ cout << agg << " " << number_of_next_save << " " << random_seed << " "<< fitness << endl;
		}
      }


      
      web_data.set_number_of_compartments(web.number_of_species());
      web_data.set_biomasses(biomass_B);
      web_data.set_bodymasses(web.get_bodymass_M());
      web_data.set_bottom(web.get_is_plant());
      if(needs_flows){
	web_data.set_flows(ls);
	web_data.adjust_links_given_flows(threshold_value);
      }
      web_data.set_area(area_per_compartment);

      cout << "web_number " << number_of_next_save-1 << endl;
      int old_prec=std::cout.precision();
      std::cout.precision(10);
      cout << "current_time " << web.current_time << endl;
      std::cout.precision(old_prec);
    }
    break;
  case nls:
    {// nls web in SCORE format:
      nlsweb=nls_web(web_file_name.c_str());
      biomass_B=nlsweb.biomass;
      biomass_B.resize(nlsweb.number_of_living_compartments);
      ls=nlsweb.trophic_intake_matrix();
 
      if(nls_area)
	web_data.set_area(nls_area);
      else{
	WARNING("area not set, assuming area of 1 meter^2");
	web_data.set_area(eval_now("1*meter^2"));
      }

      web_data.
	set_number_of_compartments(nlsweb.number_of_living_compartments);
      web_data.set_biomasses(biomass_B);
      web_data.set_flows(ls);
      web_data.adjust_links_given_flows(threshold_value);
      web_data.adjust_bottom_given_links();
    }
    break;
  case tlf:
    {// custom "Tuesday Lake" format
      read_tlf_file(web_file_name.c_str(),web_data);
    }
    break;
  case dot_web:
    web_data.set_links(load_three_column_file(web_file_name));
  }

  if(fast_stats and file_type!=xml){
    web_data.fast_stats();
  }
  
        	//~ for(int i=10;i-->0;){
		//~ fitness_only2.open("fitness_only2.dat",ios_base::app);
		//~ NewSpecies * s = web.speciate(NewSpecies::animal);
		//~ s->set_aggressivity_g(0.06); // you need to adjust this!
		//~ fitness_only2 << web.invasion_fitness(s) << endl;
		//~ REPORT(web.invasion_fitness(s));
	//~ }
	

  
  if(mass_integral_start){
    Mass_Integrator integrator;
    double mass_integral=
      integrator.integrate(biomass_B,web_data.get_bodymasses(),
			   mass_integral_start);
    REPORT(mass_integral);
  }
  
  if(print_insertion_count){
    cout << "attempted insertions " 
	 << web.get_invasion_counter() + web.get_speciation_counter() 
	 << endl;
  }
  
  if(strength_spectrum_file){
    if(choose_other || !file_type==nls)
      strength_distribution_new(web_data.get_ifrac(),
			    threshold_value/100.0,
			    strength_spectrum_file,
			    web_data.get_bodymasses());
    else
      nlsweb.
	strength_distribution(threshold_value/100.0,strength_spectrum_file);
  }
  if(distance_file && file_type==xml){
    web.distance_lists(distance_file);
  }
  
  Interaction_Matrix im=
    (needs_flows ? web_data.get_im() : Interaction_Matrix());
  if(lump_lowest_level){
    im=im.lump_lowest_level();
  }
  if(standardized){
    im=standardize(im);
  }
  if(lump){
    im=im.trophic();
  }
  if(print || pgm_file){
    if(true || standardized || lump)
      im=forgiving_tsort(im);
    if(pgm_file){
      im.pgm_write(pgm_file);
    }
    if(print){
      im.Print();
    }
  }
    if(aggressivenes_file){
    
	stringstream ss(aggressivenes_file);
	double agg = 0; ss >> agg;
	//string str = ss.str();
	//string random;
	
	//~ srand((unsigned)time(0));
	//~ int random_integer = rand();
	//~ int Number = rand();     // number to be converted to a string
	//~ ostringstream convert;   // stream used for the conversion
	//~ convert << random_seed;       
	//~ random = convert.str();

	//~ str = "_"+ str + "_"+ random;
	string out1 = "before_interaction_matrix.dat";
	//~ out1 = out1 + str + ".dat";
	string out2 = "before_agg_for_matrix.dat";
	//~ out2 = out2 + str + ".dat";
	before_interaction_matrix.open(out1.c_str(),ios_base::app);
	before_agg_for_matrix.open(out2.c_str(),ios_base::app);
	
	
		web.initialize_for_integration();

		//~ Print out original community matrix
		
		for(int i=web.s.size();i-->0;){
			before_agg_for_matrix << web.s(i).unique_id() << " "<<web.s(i).aggressivity_g() << " " << web.s(i).attack_rate_a() << " " << web.s(i).biomass_abundance_B() << endl;
			for(int j=web.s.size();j-->0;){
				if(j!=0){before_interaction_matrix << web.precomputed[i].c[j] << ",";}else{before_interaction_matrix << web.precomputed[i].c[j] << endl;}
			}
		}
}  
  if(print_flows){
	interaction_matrix.open("interaction_matrix.dat",ios_base::app);
	biomass_interaction_matrix.open("biomass_interaction_matrix.dat",ios_base::app);
	agg_for_matrix.open("agg_for_matrix.dat",ios_base::app);
    const link_strength_matrix & m=ls;
    std::cout << "   ";
    for(unsigned int j=0;j<m.size();j++)
      std::cout << std::setw(2) << (j/10)%10;
    std::cout << endl;
    std::cout << "   ";
    for(unsigned int j=0;j<m.size();j++)
      std::cout << std::setw(2) << j%10;
    std::cout << endl;
    std::cout << "   ";
    for(unsigned int j=0;j<m.size();j++)
      std::cout << "--";
    std::cout << std::endl;
    for(unsigned int i=0;i<m.size();i++){
      agg_for_matrix << web.s(i).aggressivity_g() << endl;
      std::cout << std::setw(2 ) << i%100 << "|";
      for(unsigned int j=0;j<m.size();j++){
	std::cout << std::setw(2) << order_of_magnitude(m[j][i]);
	interaction_matrix << m[i][j]/web.s(i).aggressivity_g() << endl;
	biomass_interaction_matrix << m[i][j]/web.s(i).aggressivity_g()*web.s(j).biomass_abundance_B() << endl;
      }
      std::cout << endl;
    }
  }
  if(dot_graph_file.size()){
    im.dot_graph(dot_graph_file);
  }
  if(pajek_graph_file.size()){
    web_data.pajek_graph(pajek_graph_file,threshold_value);
  }
  if(print_properties){
    sequence<const char *> names=Interaction_Matrix::prop_names();
    Interaction_Matrix::distribution pr=im.props();

    for(unsigned int j=0;j<pr.size();j++){
      std::cout << names[j] << " " 
	   << pr[j]
	   << endl;
    }
  }
  if(eigenvalue_file){
    web.eigenvalue_list(eigenvalue_file);
  }
  if(generality_file){
    ofstream os(generality_file);
    im.cumulative_prey_hist(os);
  }
  if(vulnerability_file){
    ofstream os(vulnerability_file);
    im.cumulative_predator_hist(os);
  }
  if(rank_abundance_plot){
    web.rank_abundance_plot(rank_abundance_plot,!lump_lowest_level);
  }    
  if(size_spectrum_file){
    if(file_type==nls)
      FATAL_ERROR("cannot compute size spectrum for SCOR files");
    ostringstream p,a,t;
    p << "p" << size_spectrum_file;
    a << "a" << size_spectrum_file;
    t << "t" << size_spectrum_file;
    web_data.size_spectrum(p.str().c_str(),1,0);
    web_data.size_spectrum(a.str().c_str(),0,1);
    web_data.size_spectrum(t.str().c_str(),1,1);
    Mass_Integrator().spectrum(web_data.get_biomasses(),
			       web_data.get_bodymasses(),
			       size_spectrum_file,/* take log 10?*/true);
  }    
  if(biomass_spectrum_file){
    if(file_type==nls)
      FATAL_ERROR("cannot compute biomass spectrum for SCOR files");
    ostringstream p,a,t;
    p << "p" << biomass_spectrum_file;
    a << "a" << biomass_spectrum_file;
    t << "t" << biomass_spectrum_file;
    web_data.biomass_spectrum(p.str().c_str(),1,0);
    web_data.biomass_spectrum(a.str().c_str(),0,1);
    web_data.biomass_spectrum(t.str().c_str(),1,1);
  }    
  if(rectification_base_name){
    ostringstream i,t;
    i << rectification_base_name << ".pgm";
    t << rectification_base_name << ".dat";
    permutation p=web_data.rectification();
    threshold(web_data.ifrac(),threshold_value).
      permute(p.inverse()).
      pgm_write(i.str().c_str());
    cout << "rectified_niche_width " << sqrt(web_data.Omega(p,2.2)) 
	 << "                           " << endl;
  }
  if(file_type==xml && competition_file){
    web.write_competition_matrix(competition_file);
  }
  if(species_table){
    if(file_type==xml){
      web.species_table(species_table,true,
			threshold_value,web_data,jolly);
    }else{
      web_data.species_table(species_table,true,
			     threshold_value);
    }
  }
  if(species_table2){
    if(file_type==xml){
      web.species_table2(species_table2,true,
			threshold_value,web_data,ls,coldelS,jolly);
    }else{
      web_data.species_table2(species_table2,true,
			     threshold_value);
    }
  }
  if(species_table3){
    if(file_type==xml){
      web.species_table3(species_table3,true,
			 threshold_value,web_data,ls,jolly);
    }else{
      web_data.species_table3(species_table3,true,
			     threshold_value);
    }
  }
  if(link_coefficient_matrix_file){
    if(file_type==xml){
      web.write_link_coefficient_matrix(link_coefficient_matrix_file);
    }else{
      FATAL_ERROR("Option -n not implemented for this filetype.");
    }
  }
  if(link_table){
    if(file_type==xml){
      web.link_table(link_table,threshold_value,ls,biomass_B);
    }else{
      web_data.link_table(link_table, threshold_value);
    }
  }
  if(link_table2){
    if(file_type==xml){
      web.link_table2(link_table2,threshold_value,ls,biomass_B);
    }else{
      web_data.link_table2(link_table2, threshold_value);
    }
  }
  if(file_type==xml && log_availability_structure_file){
    multinormal_distribution log_avs=
      web.log_availability_structure();
    log_avs.save(log_availability_structure_file);
    std::cout << "var_log_availability " 
	      << log_avs.var_of_sum() 
	      << std::endl;
  }
  if(print_indices){
    bool with_plants=!lump_lowest_level;

    double S=(with_plants?
	      web_data.number_of_species():
	      web_data.number_of_animals() );
      
    double H=web_data.Shanon_diversity_H(with_plants,true);

    std::cout << S << "(S)" << std::endl;
    std::cout << 1/web_data.Simpson_s_diversity_D(with_plants,true) << "(1/D)" << std::endl;
    std::cout << 1-web_data.Simpson_s_diversity_D(with_plants,true) << "(1-D)" << std::endl;
    std::cout << H << "(H)" << std::endl;
    std::cout << exp(H) << "(expH)" << std::endl;
    std::cout << H/log(S) << "(E)" << std::endl; // Pielou's evenness index
    std::cout << web_data.max_level() << "(hmax)" << std::endl;
    std::cout << web_data.mean_level() << "(hmean)" << std::endl;
    std::cout << web_data.biomass_weighted_mean_level() << "(hBmean)" << std::endl;
    if(file_type==xml){
      std::cout << 
	web.biomass_larger(min_size_of_things_we_care_about)/
	baseline_web.biomass_larger(min_size_of_things_we_care_about) << 
	"(FishI)" << std::endl;
      std::cout << web.GPP()/baseline_web.GPP() << "(GPPI)" << std::endl;
      std::cout << 
	web.Living_Planet_Index(baseline_web,with_plants,true) << "(LPI)" << std::endl;
      std::cout << 
	web.Biodiversity_Intactness_Index(baseline_web,with_plants,true) 
		<< "(BII)" << std::endl;
    }
  }
  if(report_pressures || invasion_fitness_curve || fitness_samples){
    if(!waiting_time){
      REPORT(waiting_time);
      WARNING("Cannot report evolutionary pressures");
      WARNING("Must set -w option");
    }else if(file_type!=xml){
      WARNING("Cannot report evolutionary pressures");
      WARNING("Works only for xml files");
    }else{
      if(report_pressures){
	web.report_evolutionary_pressures(iterations_for_pressures);
      }
      if(invasion_fitness_curve){
	web.invasion_fitness_curve(invasion_fitness_curve);
      }
      if(fitness_samples){
	web.invasion_fitness_samples(fitness_samples);
      }
    }
  }
  if(get_msy){
#if 0
    vector<string> catchability_function(1);
#define MMIN "(1*gram)"
    catchability_function[0]=
      "(1-1/(1+(Mmat/" MMIN ")^4))" "*(Mmat/" MMIN ")^(-1/4)";
    //    catchability_function[1]="exp(-log(Mmat/(100*gram))^2/(2*log(3)^2))";
   //    catchability_function[0]="(1-1/(1+(Mmat/(1000*gram))^4))*sin(1000*log(Mmat))^2";
   //    catchability_function[1]="(1-1/(1+(Mmat/(1000*gram))^4))*sin(7327*log(Mmat))^2";
    //catchability_function[1]="10*exp(-log(Mmat/(100*gram))^2/(2*log(3)^2))";
    //catchability_function[0]="exp(-log(Mmat/(10*gram))^2/(2*log(10)^2))";
    MSY_fleets msy(web,catchability_function);
#else 
    vector<int> fished;
    double cutoff=get_cfg_parameter("catchability_cutoff");
    for(int i=web.number_of_species();i-->0;){
      if(web.s[i].mean_bodymass_M() >= cutoff){
	fished.push_back(web.assigned_column(i));
	// if(!cin.eof()){
	//   cout << "Pipe in list of mortalities, or hit [return]" << endl;
	//   // Read fishing mortalities from standard input if
	//   // available.
	//   double F;
	//   cin >> F;
	//   web.s[i].set_fishing_mortality(F);
	// }
      }
    }
    MSY_species msy(web,fished);
    if(1){
      //individual_transposed_interaction_controller_t HCR(msy);
      //individual_target_pressure_controller_t HCR(msy);
      individual_productive_state_controller_t HCR(msy);
      //soft_individual_productive_state_controller_t HCR(msy);
      //individual_B_transposed_interaction_controller_t HCR(msy);
      //individual_B_target_pressure_controller_t HCR(msy);
      //individual_B_target_pressure_controller_t HCR(msy);
      //individual_B_productive_state_controller_t HCR(msy);
      //productive_state_controller_t HCR(msy);
      //transposed_interaction_controller_t HCR(msy);
      //target_pressure_controller_t HCR(msy);
      //growth_rate_controller_t HCR(msy);
      //CFP_controller_t HCR(msy);
      //DPF_controller_t HCR(msy);
      HCR.recompute_strategy();
      for(int rep=31;rep-->0;){
	HCR.relax(eval("50*year"),NewWeb::species_set_t(),false);
	HCR.recompute_strategy();
      }
      HCR.relax(eval("50*year"),NewWeb::species_set_t(),false);
      NewWeb test=HCR;
      cfgPermanent parameters;
      XMLStore()
	.put(version,"version")
	.put(compilation_time,"compilation_time")
	.put(&parameters,"adjustable_parameters")
	.put(&HCR,"FoodWeb")
	.save("HCR-web.xml.bz2");
      exit(0);
    }
#endif
    msy.find_MSY();
    REPORT(msy.yield());
    //    exit_now=false; // we are exiting soon anyway!
    msy.optimally_exploited_web(web);

    {
      cfgPermanent parameters;
      XMLStore()
	.put(version,"version")
	.put(compilation_time,"compilation_time")
	.put(&parameters,"adjustable_parameters")
	.put(&web,"FoodWeb")
	.save("MSY-web.xml.bz2");
    }
  }
}

