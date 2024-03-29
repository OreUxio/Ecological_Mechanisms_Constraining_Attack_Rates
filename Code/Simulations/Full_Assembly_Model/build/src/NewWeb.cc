// $Id: NewWeb.cc 2503 2017-02-27 17:37:46Z axel $
#include <string.h>
#include <stdlib.h>
#include <unistd.h>  // for getpid, sysconf and number of CPUs
#include <sys/types.h> // for number of CPUs on MAC
#include <sys/sysctl.h> // for number of CPUs on MAC
#include <math.h> // this should include NAN
#ifndef NAN
//#warning "defining my own NAN"
static const double my_nan=sqrt(-1.0);
#define NAN my_nan
#endif

#include <fstream>
#include <iostream>

#include "NewWeb.h"
#include "random.h"
//#include "TOP.h"
#include "evaluate.h"
#include "random_pick_field.h"
#include "matrix_transformers.h"
//#include "linpack_eigen.h"
#include "xy_graph.h"
#include "link_strength_matrix.h"
#include "CompMatrix.h"
#include "packed_simulation.h"
//#include <values.h>


using namespace std;

static my_evaluator_t eval_here;
static double animal_invasion_pressure_kappa_a=1.4;
static double plant_invasion_pressure_kappa_p=1.4;
static int number_of_plant_species=0; // UNUSED!! 
double any_invasion_pressure_q=0;
double unified_invasion_pressure_kappa=0; //extern
static int mean_field_threshold=0; // At what abundance to use mean
				   // field invasions. 0=NEVER
static double speciation_rate=1;  // UNUSED!!
static double invasion_rate=1;    // UNUSED!!
static double n_neutral=1;
static double conservation_tolerance=10000;///< prevents inserted species from being delete right away
int max_steady_state_size=0; ///< 0 mean don't use
static double fishall_intercept=0*eval_here("1/year");
static double fishall_slope=0*eval_here("1/kilogram/year");
static double fishall_residual_se=0*eval_here("1/year");
static double fishMrangemax=10000000000.0;
static double fishTLrangemax=10000000000.0;
static double const_fishing_mortality_rate=0.1;
static double Fmax=1*eval_here("1/year");
static double fishMrangeCSFfactor=1.0;
static int addrandF=0;
static double randFstd=0.1*eval_here("1/year");
static double log_small_value_tolerance=log(1.001);
static int dumb_preconditioner=true;
static double fishMlowerthreshold=eval_here("1*gram");
static double fishMupperthreshold=eval_here("345000*gram");
int use_packed_simulation=0; ///< Use experimental, optimized code.
#ifdef SET_CPU_AFFINITY
CPU_set_list_t CPU_set_list=CPU_set_list_t(0);
#endif
int log_steady_state_stats=0;


// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  {
    CFGDOUBLE(speciation_rate),
    CFGDOUBLE(invasion_rate),
    CFGDOUBLE(animal_invasion_pressure_kappa_a),
    CFGDOUBLE(plant_invasion_pressure_kappa_p),
    CFGDOUBLE(unified_invasion_pressure_kappa),
    CFGDOUBLE(any_invasion_pressure_q),
    CFGDOUBLE(n_neutral),
    CFGINT(number_of_plant_species),
    CFGINT(mean_field_threshold),
    CFGINT(max_steady_state_size),
    CFGDOUBLE(fishall_intercept),
    CFGDOUBLE(fishall_slope),
    CFGDOUBLE(fishall_residual_se),
    CFGDOUBLE(fishMrangemax),
    CFGDOUBLE(fishTLrangemax),
    CFGDOUBLE(const_fishing_mortality_rate),
    CFGDOUBLE(Fmax),
    CFGDOUBLE(log_small_value_tolerance),
    CFGINT(addrandF),
    CFGDOUBLE(randFstd),
    CFGDOUBLE(fishMrangeCSFfactor),
    CFGINT(dumb_preconditioner),
    CFGINT(log_steady_state_stats),
    {0, CFG_END, 0}
  };
static cfg_add dummy(cfg);


const double log_DBL_MAX=5*log(DBL_MAX)/7.0; // we added some safety margin
const double log_DBL_MIN=5*log(DBL_MIN)/7.0; // we added some safety margin

NewWeb::NewWeb():
  invasion_counter(0),
  speciation_counter(0),
  max_ever_plant_richness(0),
  dynamics_accelerator(0),
  compute_diagnostics(false)
{
  // nothing to do
}

NewWeb::~NewWeb(){
  while(precomputed.begin()!=precomputed.end()){
    delete precomputed.detach(precomputed.end()-1);
  }
};

// This copy constructure is not essential, but save us from copying 
// some parts which are not essential.
NewWeb::NewWeb(const NewWeb& other):
  biomass_B(other.biomass_B),
  common_factor(other.common_factor),
  dynamics_accelerator(0),
  fx(other.fx),
  compute_diagnostics(false),
  speciation_counter(other.speciation_counter),
  invasion_counter(other.invasion_counter),
  max_ever_plant_richness(other.max_ever_plant_richness),
  s(other.s),
  species_at(other.species_at),
  assigned_column(other.assigned_column),
  otherwebs(other.otherwebs) {
}
	
NewWeb& NewWeb::operator=(const NewWeb& other){
  if(this==&other) return *this;

  if(dynamics_accelerator) {// Perhaps not necessary but safer.
    delete dynamics_accelerator;
    dynamics_accelerator = 0;
  }
  while(precomputed.begin()!=precomputed.end()){
    delete precomputed.detach(precomputed.end()-1);
  }

  biomass_B=other.biomass_B;
  common_factor=other.common_factor;
  fx=other.fx;
  compute_diagnostics=false;
  current_time=other.current_time;
  speciation_counter=other.speciation_counter;
  invasion_counter=other.invasion_counter;
  max_ever_plant_richness=other.max_ever_plant_richness;
  s=other.s;
  species_at=other.species_at;
  assigned_column=other.assigned_column;
  otherwebs=other.otherwebs;
  return *this;
}
  

void NewWeb::data_mapping(){
  PERMANENT(current_time); //from ODE_dynamical_object
  PERMANENT(invasion_counter);
  PERMANENT(speciation_counter);
  PERMANENT(max_ever_plant_richness);
  if(the_remember->read_or_write()==rememberWrite){
    int new_foraging_vectors=1;
    PERMANENT(new_foraging_vectors);
    int used_species=0;
    for(unsigned int i=0;i<s.size();i++){
      NewSpecies & Species = s[i];
      PERMANENTP(Species);
      used_species++;
    }
    std::cout << 
      "**** WRITING " << used_species << " SPECIES ****" << std::endl;
    assigned_column.resize(s.size());
    for(Ancestry::iterator i=ancestors.begin();
	i!=ancestors.end();
	++i){
      if(i->second.column==parent_not_extant){
	NewSpecies & Ancestor=i->second.extinct_ancestor;
	PERMANENTP(Ancestor);
      }
    }
    PERMANENT(assigned_column);
  }else if(the_remember->read_or_write()==rememberRead){
    int new_foraging_vectors=0;
    PERMANENT(new_foraging_vectors);
    NewSpecies Species;
    if(PERMANENT(assigned_column)!=0){
      // there were no assigned columns, so we do it by hand:
      for(int i=s.size();i-->0;)
	assigned_column[i]=i;
    }
    // species will be reordered, so we should remember
    // assigned_column elsewhere
    sequence<int> assigned_column_hold(assigned_column);
    assigned_column.resize(0);
    int i_before_saving=0;
    while(PERMANENTP(Species)==0){
      // This is a bit time consuming right now, because we have to
      // reallocate memory with each species inserted:
      int i_now;
      s[i_now=get_unused_index_and_allocate_memory(Species.is_a_plant())]=
	Species;
      assigned_column[i_now]=assigned_column_hold(i_before_saving);
      species_at[assigned_column_hold(i_before_saving)];//resize species_at
      i_before_saving++;
    }
    NewSpecies Ancestor;
    while(PERMANENTP(Ancestor)==0){
      ancestors[Ancestor.unique_id()]=ancestral_record(Ancestor);
    }    
    // the species_at array is messed up, we have to recompute it:
    species_at.resize(0);
    for(int i=s.size();i-->0;){
      species_at[assigned_column(i)]=i+1;
    }
    ASSERT(column_unused==-1);
    species_at-=1;
    // fix foraging vectors if they are old
    if(!new_foraging_vectors){
      for(int i=s.size();i-->0;){
	s[i].fix_foraging_vector();
      }
    }
  }
};

void NewWeb::write_state_to(ODE_vector & state) const{
  int S=s.size();
  for(int i=S;i-->0;){
    state[i]=log(s(i).biomass_abundance_B());
  }
}

void NewWeb::read_state_from(const ODE_vector & state){
  int S=s.size();
  for(int i=S;i-->0;){
    s(i).set_biomass_abundance_B(exp(state[i]));
  }
}


int NewWeb::
number_of_variables() const{
  return s.size();
}

const int NewWeb::column_unused=-1;
const int NewWeb::parent_not_extant=-1;

int NewWeb::column_allocate(){
  int i=0;
  while(true){
    if(i==species_at.size()){
      species_at.resize(i+1);
      return i;
    }
    if(species_at(i)!=column_unused)
      i++;
    else
      return i;
  }
}

void NewWeb::column_free(int i){
  species_at[i]=column_unused;
}

int NewWeb::get_unused_index_and_allocate_memory(int plant){
  int si=s.size();
  int i;
  // enlarge arrays (for a simple_vector this is a cheap operation if
  // going back and forth):
  precomputed.resize(si+1,s);
  fx.resize(si+1);
  TRACE(fx.size(),FLOWS);
  biomass_B.resize(si+1);
  s.resize(si+1);
  assigned_column.resize(si+1);
  if(plant){
    i=si;
  }else{//animal
    move_species(s.n_animals,si);
    i=s.n_animals;
    s.n_animals++;
  }
  species_at[assigned_column(i)=column_allocate()]=i;
  return i;
}

void NewWeb::move_species(int old_i,int new_i){
  if(old_i==new_i)
    return;

  species_at(assigned_column(new_i)=assigned_column(old_i))=new_i;
  // Note that this also moved c(old_i,old_i) to c(new_i,new_i)!
  s(new_i)=s(old_i);
  precomputed.move(old_i,new_i,s);
}
  
int NewWeb::delete_species_and_report(int i){
	ofstream list_of_inserted_deleted_species;
	list_of_inserted_deleted_species.open("list_of_inserted_deleted_species.dat",ios_base::app);
	list_of_inserted_deleted_species << s(i).unique_id() << " " << s(i).parent_unique_id() << " " << 0 << " " << s(i).aggressivity_g()  << " " << 0 << endl;
	return delete_species(i);
}

int NewWeb::delete_species(int i){
  ASSERT(0<=i);
  ASSERT(i<s.size());
  TRACE(s[i].unique_id(),SPECIES);
  TRACE(s[i].clade_id(),SPECIES);
  int end=s.size()-1;

  if(log_steady_state_stats && !is_plant(i)){
    cout << "note " << s(i).aggressivity_g();
    cout << " " << s(i).niche_width_wt();
    cout << " " << s(i).bodymass();
    cout << " " << s(i).number_of_invading_offspring();
    cout << " " << s(i).age();
    sequence< double > diet(s.size());
    get_diet(i,diet);
    int main_prey_of_predator=
      max_element(diet.begin(),diet.end())-diet.begin();
    cout << " " << is_plant(main_prey_of_predator);
    cout << endl;
  }

  column_free(assigned_column(i));

  if(is_plant(i)){
    move_species(end,i);
  }else{
    move_species(--s.n_animals,i);
    move_species(end,s.n_animals);
  }
  s.resize(end);
  fx.resize(end);
  precomputed.resize(end,s);
  steady_state.clear();

  return i < end;
}

void NewWeb::delete_forbidden_species(){
  for(int i=s.size();i-->0;){ // Loop must go backwards!
    if(s(i).outside_allowed_trait_space()){
      delete_species(i);
    }
  }
}

// Note that insert() deletes the object the argument points to!
int NewWeb::insert(const NewSpecies * sp){
  int i=get_unused_index_and_allocate_memory(sp->is_a_plant());
  s[i]=*sp;
  delete sp;
  s[i].set_sequential_id();
  precomputed.update(i,s);
  steady_state.clear();
  max_ever_plant_richness=
    max(max_ever_plant_richness,number_of_plants());

  int parent_i=find_extant_parent(* sp);
  if(parent_i != parent_not_extant){
    s(parent_i).notify_speciation(* sp);
  }
  return i;
}

// Helper to find a parent species:
int NewWeb::find_extant_parent(const NewSpecies & sp){
  int puid=sp.parent_unique_id();
  for(int i=s.size();i-->0;){
    if(s[i].unique_id()==puid){
      return i;
    }
  }
  return parent_not_extant;
}

void NewWeb::Ancestry::note_speciation(const NewSpecies & offspring,
				       const NewWeb & web){
  int puid=offspring.parent_unique_id();
  Ancestry::iterator parent_p=find(puid);
  if(parent_p==this->end()){
    // Parent not recorded in ancestry, yet
    for(int i=web.s.size();i-->0;){
      if(web.s[i].unique_id()==puid){
	// Found the parent
	parent_p=
	  Ancestry::
	  insert(pair< int, ancestral_record >
		 (puid,ancestral_record(web.assigned_column[i])) ).first;
	break;
      }else{
	WARNING("Speciating species has no parent!");
	return;
      }
    }
  }
  parent_p->second.number_of_known_offspring++;
  return;
}

void NewWeb::Ancestry::
note_impending_extinction(const NewSpecies & going_extinct){
  int uid=going_extinct.unique_id();
  
  Ancestry::iterator record_p=find(uid);
  if(record_p->second.number_of_known_offspring == 0){
    int puid=going_extinct.parent_unique_id();
    // Remove entire branch from tree:
    record_p=Ancestry::find(puid);
    FATAL_ERROR("REMOVAL OF ANCESTRAL BRANCHES NOT IMPLEMENTED YET");
    if(record_p != Ancestry::end()){
      record_p->second.number_of_known_offspring--;
      //... stuff missing here, probaby need some recursive function or so.
    }
    return;
  }else{
    // explicitly record among ancestors, overwriting implicit record:
    FATAL_ERROR("number_of_known_offspring WOULD GET LOST IN NEXT LINE!");
    record_p->second=ancestral_record(going_extinct);
  }
  return;
}

void NewWeb::clear(){
  for(int i=s.size();i-->0;){
    delete_species(i);
  }
}

void NewWeb::record_for_steady_state(){
  if(not dynamics_accelerator){
    steady_state.
      push_back(new saved_state(current_time,biomass_B,common_factor));
  }else{
    steady_state.
      push_back(new saved_state
		(current_time,
		 vector_with_max(dynamics_accelerator->biomass_B.begin(),
				 dynamics_accelerator->biomass_B.end() ),
		 vector_with_max(dynamics_accelerator->common_factor.begin(),
				 dynamics_accelerator->common_factor.end() )
		 ));
  }    
  if(max_steady_state_size && 
     steady_state.size()==2*max_steady_state_size){
    steady_state_t::iterator start=
      steady_state.begin();
    steady_state_t::iterator rest=
      steady_state.begin()+max_steady_state_size;
    for(;rest!=steady_state.end();++start,++rest){
      ptr_swap(start,rest);
    }
    while(steady_state.size()>max_steady_state_size){
      delete steady_state.detach(steady_state.end()-1);
    }
  }
}

void NewWeb::delete_weeds(){
  double max_animal_size=0;
  for(int i=s.n_animals;i-->0;){
    max_animal_size=
      max<double>(max_animal_size,s(i).mean_bodymass_M());
  }
  const double weed_threshold_mass=
    max_animal_size/
    pow(10,get_cfg_parameter("log10_bodymass_ratio_window_center"));

  for(int i=s.size();i-->s.n_animals;){
    if(s(i).mean_bodymass_M()>weed_threshold_mass){
      delete_species(i);
    }
  }
}

	
void NewWeb::steady_state_is_fixed_point(){
  WARNING("deleting steady_state up to final point");
  steady_state.clear();
  steady_state.push_back(new saved_state(current_time,biomass_B,common_factor));
}
	
double NewWeb::biomass_larger(double threshold) const{
  double sum=0;
  for(int i=s.size();i-->0;){
    if(s[i].bodymass()>threshold)
      sum+=s[i].biomass_abundance_B();
  }
  return sum;
}

double NewWeb::biomass_yield(){
  double Y=0;
  for(int i=s.size();i-->0;){
    Y+=s(i).fishing_yield();
  }
  return Y;
}


double NewWeb::plant_biomass() const{
  double sum=0;
  for(int i=s.size();i-->s.n_animals;){
    sum+=s[i].biomass_abundance_B();
  }
  return sum;
}

double NewWeb::animal_biomass() const{
  double sum=0;
  for(int i=s.n_animals;i-->0;){
    sum+=s[i].biomass_abundance_B();
  }
  return sum;
}

void NewWeb::check_internal_consistency(int check_data_consistency){
  if(check_data_consistency>1){
    REPORT(s.n_animals);
    REPORT(assigned_column);
    REPORT(species_at);
  }
  for(int i=species_at.size();i-->0;){
    if(species_at(i)!=column_unused){
      ALWAYS_ASSERT(assigned_column(species_at(i))==i);
    }
  }
  ALWAYS_ASSERT(s.n_animals<=s.size());
  for(int i=s.size();i-->0;){
    ALWAYS_ASSERT(s(i).is_a_plant()==is_plant(i));
    ALWAYS_ASSERT(is_plant(i)==(i>=s.n_animals));
  }
  
  
  precomputed_t pre2=precomputed;
  initialize_precomputed();
  
  for(int i=0;i<s.size();i++){
    pre2.update(i,s);
  }
  if(do_switching){
    for(int i=0;i<s.size();i++){
      if(check_data_consistency>1){
	REPORT(i);
	REPORT(pre2[i].csc_eating);
	REPORT(precomputed[i].csc_eating);
      }
      ALWAYS_ASSERT(pre2[i].csc_eating==precomputed[i].csc_eating);
      if(check_data_consistency>1){
	REPORT(pre2[i].csc_being_eaten);
	REPORT(precomputed[i].csc_being_eaten);
      }
      ALWAYS_ASSERT(pre2[i].csc_being_eaten==precomputed[i].csc_being_eaten);
      if(check_data_consistency>1){
	REPORT(pre2[i].c);
	REPORT(precomputed[i].c);
      }
      ALWAYS_ASSERT(pre2[i].c==precomputed[i].c);
    }
  }else{// not switching:
    // .. not implemented ..
  }
}

NewSpecies * NewWeb::speciate_selected(double parent){
  speciation_counter++;
  return s[parent].speciate();
}  

// Find a random species and speciate from it:
NewSpecies * NewWeb::speciate(double plant_fraction){
  bool speciate_plant;
  int parent;
  if(unified_invasion_pressure_kappa){
    parent=random_integer(s.size());
    speciate_plant=s[parent].is_a_plant();
  }else{
    if(plant_fraction <=0 && number_of_animals()==0){
      FATAL_ERROR("No animals to speciate");
      return 0;
    }
    if(plant_fraction >=1 && number_of_plants()==0){
      FATAL_ERROR("No plants to speciate");
      return 0;
    }
    if(number_of_animals()==0){
      if(number_of_plants()>0)
	speciate_plant=1;
      else{
	FATAL_ERROR("No species to speciate");
	return 0;
      }
    }
    else if(number_of_plants()==0)
      speciate_plant=0;
    else speciate_plant=unirand()<plant_fraction;
    do{
      parent=random_integer(s.size());
    }while(s[parent].is_a_plant()!=speciate_plant);
  }
  return speciate_selected(parent);
}
  

NewSpecies * NewWeb::speciate_by_biomass(int get_a_plant){
  int parent;
  if(number_of_species()==0){
    FATAL_ERROR("no species to speciate");
  }
  TRACE(get_a_plant,SPECIES);
  double biomass_sum=(get_a_plant?plant_biomass():animal_biomass());
  TRACE(biomass_sum,SPECIES);
  double specindex=biomass_sum*unirand();
  double cumulative_mass_sum=0;
  parent=0;
  if(get_a_plant){
    for(int j=s.size();j-->s.n_animals;){
      cumulative_mass_sum+=s[j].biomass_abundance_B();
      if(cumulative_mass_sum>specindex){
	parent=j;
	break;
      }
    }
  }else{// get an animal
    for(int j=s.n_animals;j-->0;){
      cumulative_mass_sum+=s[j].biomass_abundance_B();
      if(cumulative_mass_sum>specindex){
	parent=j;
	break;
      }
    }
  }
  TRACE(parent,SPECIES);
  TRACE(s[parent].unique_id(),SPECIES);
  return speciate_selected(parent);
}
  

NewSpecies * NewWeb::invade(double plant_fraction){
  bool get_a_plant=(unirand()<plant_fraction);
  NewSpecies * species_p;
  if(mean_field_threshold){
      int sample_size=(get_a_plant ? s.size()-s.n_animals : s.n_animals)-1;
      // (subtract 1 because the last is the empty slot 1)
      if(sample_size < mean_field_threshold){
	species_p=new NewSpecies(get_a_plant);
      }else{
	int offset=(get_a_plant ? s.n_animals : 0);
	species_p=NewSpecies(s,offset,sample_size).speciate();
      }
  }else{
    // invade from other webs?
    if(otherwebs.current_mode()==Otherwebs::newest)
      otherwebs.get_newest_webs();
    int n_elsewhere=
      (get_a_plant?
       otherwebs.total_number_of_plants():
       otherwebs.total_number_of_animals() );
    if(otherwebs.current_mode()!=Otherwebs::off &&  
       n_elsewhere>0 && unirand()<n_elsewhere/double(n_neutral+n_elsewhere) ){
      species_p=
	otherwebs.
	get_random_species((get_a_plant ? 
			    NewSpecies::plant : NewSpecies::animal)).
	speciate_as_new();
    }else{
      species_p=new NewSpecies(get_a_plant);
    }
  }
  invasion_counter++;
  return species_p;
}

NewSpecies * NewWeb::invade_by_biomass(double plant_fraction){
  bool get_a_plant=(unirand()<plant_fraction);
  NewSpecies * species_p;
  // invade from other webs?
  if(otherwebs.current_mode()==Otherwebs::newest)
    otherwebs.get_newest_webs();
  if(get_a_plant ?
     otherwebs.total_number_of_plants()>0 :
     otherwebs.total_number_of_animals()>0 ){
    species_p=
      otherwebs.
      get_random_species_by_biomass(get_a_plant ? 
				    NewSpecies::plant : NewSpecies::animal).
      speciate_as_new();
  }else{
    species_p=new NewSpecies(get_a_plant);
  }
  invasion_counter++;
  return species_p;
}
  
NewSpecies * NewWeb::invade_or_speciate(double plant_fraction){
  if(plant_fraction!=0 && plant_fraction!=1){
    plant_fraction=(unirand()<plant_fraction);
  }
  int number_of_relevant_species;
  if(unified_invasion_pressure_kappa){
    number_of_relevant_species=number_of_species();
  }else{
    number_of_relevant_species=
      ( plant_fraction==0 ? number_of_animals() : number_of_plants() );
  }
  double kappa;
  if(unified_invasion_pressure_kappa){
    kappa=unified_invasion_pressure_kappa;
  }else{
    kappa=
      ( plant_fraction==0 ? animal_invasion_pressure_kappa_a : plant_invasion_pressure_kappa_p );
  }
  TRACE(kappa,SPECIES);
  if(!number_of_relevant_species ||
     kappa/(kappa + number_of_relevant_species) >= 
     unirand() ){
    //    cout << "invading " << plant_fraction << endl;
    return invade(plant_fraction);
  }else{
    //    cout << "speciating " << plant_fraction << endl;
    ALWAYS_ASSERT(number_of_relevant_species);
    return speciate(plant_fraction);
  }
}

NewSpecies * NewWeb::invade_or_speciate_any(double get_a_plant){
  if(int(get_a_plant)!=get_a_plant){
    get_a_plant=(unirand()<get_a_plant);
  }
  if(!any_invasion_pressure_q){
    FATAL_ERROR("invade_or_speciate_any needs any_invasion_pressure_q");
  }
  int plants_to_invade=otherwebs.total_number_of_plants();
  int animals_to_invade=otherwebs.total_number_of_animals();
  int species_to_invade=(get_a_plant?plants_to_invade:animals_to_invade);
  int species_here=(get_a_plant?number_of_plants():number_of_animals());
  int total_number_of_species=
    species_to_invade+species_here;

  if(total_number_of_species==0){
    return new NewSpecies(get_a_plant);
  }
    
  
  double probability_of_invasion=
    (any_invasion_pressure_q*species_to_invade)/total_number_of_species;

  if(unirand()<probability_of_invasion || !species_here){
    // invade:
    return invade(get_a_plant);
  }else{
    // speciate:
    return speciate(get_a_plant);
  }
}

NewSpecies * NewWeb::invade_or_speciate_any_by_biomass(double get_a_plant){
  if(int(get_a_plant)!=get_a_plant){
    get_a_plant=(unirand()<get_a_plant);
  }
  if(!any_invasion_pressure_q){
    FATAL_ERROR
      ("invade_or_speciate_any_by_biomass needs any_invasion_pressure_q");
  }
  double plants_to_invade=otherwebs.total_plant_biomass();
  double animals_to_invade=otherwebs.total_animal_biomass();
  double species_to_invade=
    (get_a_plant?
     otherwebs.total_plant_biomass() :
     otherwebs.total_animal_biomass() );
  double species_here=
    (get_a_plant?
     plant_biomass() :
     animal_biomass() );

  int total_number_of_species=
    species_to_invade+species_here;
  
  if(total_number_of_species==0){
    return new NewSpecies(get_a_plant);
  }
  
  double probability_of_invasion=
    (any_invasion_pressure_q*species_to_invade)/total_number_of_species;

  if(unirand()<probability_of_invasion || !species_here){
    // invade:
    return invade_by_biomass(get_a_plant);
  }else{
    // speciate:
    return speciate_by_biomass(get_a_plant);
  }
}

// Some dead code (still here, because annealing was really cool!):
 
// #define ONLY_FIT_AFTER_ANNEALING_DEFINE(X)              \
// int NewWeb::X##_only_fit_after_annealing(double plant_fraction){	\
//   species_list_t save_s=s;				\
//   int new_species;					\
//   int attempts_left=999;				\
//   do{							\
//     new_species=insert(X(plant_fraction));		\
//     if(fitness_after_annealing(new_species)>=0 ||	\
//        attempts_left-- <=0 ) break;			\
//     delete_species(new_species);			\
//     s=save_s;						\
//   }while(true);						\
//   if(attempts_left <=0){				\
//     WARNING(#X " unfit " << (plant_fraction?"plant":(plant_fraction==0?"animal":""))); \
//   }							\
//   return new_species;					\
// }

// ONLY_FIT_AFTER_ANNEALING_DEFINE(invade);
// ONLY_FIT_AFTER_ANNEALING_DEFINE(speciate);
// ONLY_FIT_AFTER_ANNEALING_DEFINE(invade_or_speciate);
 
// int NewWeb::invade_only_fit_or_speciate_only_fit_after_annealing(double plant_fraction){
//   if(plant_fraction!=0 && plant_fraction!=1){
//     plant_fraction=(unirand()<plant_fraction);
//   }
//   int number_of_relevant_species;
//   if(unified_invasion_pressure_kappa){
//     number_of_relevant_species=number_of_species();
//   }else{
//     number_of_relevant_species=
//       ( plant_fraction==0 ? number_of_animals() : number_of_plants() );
//   }
//   double kappa;
//   if(unified_invasion_pressure_kappa){
//     kappa=unified_invasion_pressure_kappa;
//   }else{
//     kappa=
//       ( plant_fraction==0 ? animal_invasion_pressure_kappa_a : plant_invasion_pressure_kappa_p );
//   }
//   TRACE(kappa,SPECIES);
//   REPORT(kappa);
//   if(kappa/(kappa + number_of_relevant_species) >= 
//      unirand() ){
//     cout << "invading " << plant_fraction << endl;
//     return invade_only_fit_after_annealing(plant_fraction);
//   }else{
//     cout << "speciating " << plant_fraction << endl;
//     return speciate_only_fit_after_annealing(plant_fraction);
//   }
// }

// int NewWeb::invade_fit_or_speciate_fit (double plant_fraction){
//   // In this function FIRST the species to speciate is decided, and
//   // THEN a fit descendant is sought.
//  if(plant_fraction!=0 && plant_fraction!=1){
//     plant_fraction=(unirand()<plant_fraction);
//   }
//   // this assumes unified_invasion_pressure_kappa to be set
//   ALWAYS_ASSERT(unified_invasion_pressure_kappa);
//   double kappa=unified_invasion_pressure_kappa;
//   if(kappa/(kappa + s.size() )  >=   unirand()){
//     cout << "invading " << int(plant_fraction) << endl;
//     return invade_only_fit(plant_fraction);
//   }else{
//     int parent=random_integer(s.size());
//     cout << "speciating " << is_plant(parent) << endl;
//     int child=get_unused_index_and_allocate_memory(is_plant(parent));
//     int attempts_left = 999;
//     do{
//       TRACE(parent,SPECIES);
//       TRACE(child,SPECIES);
//       s[child]=s[parent].speciate();
//       precomputed.update(child,s);
//       speciation_counter++;
//       if(fitness_after_annealing(child)>=0) break;
//     }while(attempts_left-- > 0);
//     if(!attempts_left)
//       WARNING("speciating unfit " << (plant_fraction?"plant":(plant_fraction==0?"animal":""))); 
//     return child;
//   }
// }

// double NewWeb::fitness_after_annealing(int test){
//   double best_fitness=fitness(test);
//   NewSpecies best_species=s[test];
//   const double annealing_slowdown=0.1;
//   const int annealing_attempts=niche_space_dimensions_D*10;
//   int attempts_left=annealing_attempts;
//   int annealing_steps=0;
//   do{
//     s[test]=best_species.speciate(annealing_slowdown);
//     precomputed.update(test,s);
//     double new_fitness=fitness(test);
//     if(new_fitness>best_fitness){
//       best_fitness=new_fitness;
//       best_species=s[test];
//       attempts_left=annealing_attempts;
//       annealing_steps++;
//     }
//   }while(--attempts_left);
//   REPORT(annealing_steps);
//   s[test]=best_species;
//   return best_fitness;
// }


// void NewWeb::update_link_strength(int i){
//   // i on other:
//   if(!s(i).is_plant){
//     for(int j=s.size();j-->0;){
//       log_c(j,i)=s(i).log_foraging_strength_c_on(s(j));
//       c(j,i)=(log_c(j,i)<log_DBL_MIN ? 0 : exp(log_c(j,i)));
//     }
//     TRACE(i,LINKS);
//     //    TRACE(log_c(i),LINKS);
//   }else{
//     for(int j=s.size();j-->0;){
//       if(s(j).is_plant)
// 	c(j,i)=s(i).plant_hampering_by(s(j));
//       else 
// 	c(j,i)=0;
//     }
//     TRACE(i,LINKS);
//     //    TRACE(c(i),LINKS);
//   }
//   // other on i
//   for(int j=s.size();j-->0;){
//     if(!s(j).is_plant){
//       log_c(i,j)=s(j).log_foraging_strength_c_on(s(i));
//       c(i,j)=(log_c(i,j)<log_DBL_MIN ? 0 : exp(log_c(i,j)));
//     }
//     else if(s(i).is_plant)
//       c(i,j)=s(j).plant_hampering_by(s(i));
//     else 
//       c(i,j)=0;
//   }
// }

void NewWeb::initialize_for_integration(){
  NewSpecies::initialize();
  for(int i=s.size();i-->0;){
    s(i).adjust_trait_vector_dimensions();
  }

  initialize_precomputed();
}

void NewWeb::prepare_for_integration(){
  if(dynamics_accelerator) {
    delete dynamics_accelerator;
    dynamics_accelerator = 0;
  }
  if(use_packed_simulation and
     number_of_animals() >= max_num_threads and
     number_of_plants()  >= max_num_threads and
     do_switching
     ){
    dynamics_accelerator=new packed_simulation(*this);
  }
}
void NewWeb::cleanup_after_integration(){
  if(dynamics_accelerator) {
    delete dynamics_accelerator;
    dynamics_accelerator = 0;
  }
}

void NewWeb::precondition(ODE_vector const & state,
		  ODE_vector const & in,
		  ODE_vector & out,
		  realtype gamma,
		  bool left_rather_than_right){
//   The preconditioner computes
//   out = P_i^(-1) in,
//   where i=1,2 and P=P_1*P_2 is an approximation of the Newton Matrix N
//   N = I - gamma J
//   where J is the Jacobian of the dynamics f.

//   A dumb preconditioner does not do much harm,
//   so we check that we are at least better than this:
  
  if(dumb_preconditioner || !left_rather_than_right){
    memcpy(&out[0],&in[0],in.size()*sizeof(in[0]));
  }else{
    if(dynamics_accelerator) {
      dynamics_accelerator->precondition(state,in,out,gamma);
    }else{
      for(int i=s.size();i-->0;){
	out[i]=in[i]/(1+gamma*s(i).turnover_rate_r());
      }
    }
  }
}

void NewWeb::get_inherent_rates(ODE_vector & rates){
  for(int i=number_of_variables();i-->0;){
    rates[i]=s(i).turnover_rate_r();
  }
}


// the next is dead code??
link_strength_matrix select(const link_strength_matrix & m,
			    const sequence<short int> & u){
  link_strength_matrix m0;
  int i0=0;
  for(unsigned int i=0;i<m.size();i++){
    if(u[i]){
      int j0=0;
      for(unsigned int j=0;j<m.size();j++){
	if(u[j]){
	  m0(j0,i0)=m(j,i);
	  j0++;
	}
      }
      i0++;
    }
  }
  return m0;
}

link_strength_matrix NewWeb::matching_matrix(){
  prepare_for_integration();
  link_strength_matrix m;
  m.resize(s.size());
  for(int i=s.size();i-->0;){
    if(s[i].is_a_plant())
      for(int j=s.size();j-->0;){
	m(j,i)=0;
      }
    else // j is the forager:
      for(int j=s.size();j-->0;){
	m(j,i)=precomputed[j].c[i];
      }
  }
  return m;
}
		   
const link_strength_matrix & NewWeb::get_intake_matrix(){
  return fx;
}

link_strength_matrix NewWeb::intake_matrix(bool initialize){
  prepare_for_integration();
  ODE_vector state(number_of_variables());
  NewWeb::write_state_to(state);
  compute_flows(state,initialize); // sets fx as a side effect
  return fx;
}
		   
link_strength_matrix NewWeb::scaled_intake_matrix(bool initialize){
  link_strength_matrix f=
    intake_matrix(initialize);
  for(int i=s.size();i-->0;){
    for(int j=s.size();j-->0;){
      f(i,j)/=s[i].turnover_rate_r();
    }
  }
  return f;
}
	
#if 1
void NewWeb::get_top_down_strength(){
  ODE_vector ddt(number_of_variables());
  ODE_state ode_state(this); 
  dynamics(ode_state,ddt);
  return;
}
#else
void NewWeb::get_top_down_strength(){
  int hold=get_cfg_parameter("record_relaxation_dynamics");
  set_cfg_parameter("record_relaxation_dynamics","0");
  relax(0);
  set_cfg_parameter("record_relaxation_dynamics",
		    (hold==0?"0":"1") );
  return;
}
#endif
		   
void NewWeb::enforce_equilibrium(){
  // ... and delete species for which it cannot be enforced
  int bad_species;
  do{
    bad_species=-1; // -1 means no species is bad
    {
      ODE_vector state(number_of_variables()),
	time_derivative(number_of_variables());
      NewWeb::write_state_to(state);
      NewWeb::dynamics(state,time_derivative);
      for(int i=s.size();i-->0;){
	if(!s[i].increase_loss_rate_by(time_derivative[i])){
	  bad_species=i;
	  break;
	}
      }
    }
    if(bad_species>=0){
      WARNING("deleting species " << assigned_column[bad_species]);
      delete_species(bad_species);
    }else{
      return;
    }
  }while(true);
}

double NewWeb::fitness(int test){
  ODE_vector state(number_of_variables()),
    time_derivative(number_of_variables());
  NewWeb::write_state_to(state);
  NewWeb::dynamics(state,time_derivative);
  return time_derivative[test];
}

double NewWeb::fitness(NewSpecies & test_species){
  // species should have number of individuals set!
  int test=get_unused_index_and_allocate_memory(test_species.is_a_plant());
  s[test]=test_species;
  precomputed.update(test,s);
  double f=fitness(test);
  delete_species(test);
  return f;
}

void NewWeb::delete_random_animal_species(){
  delete_species(unirand()*number_of_animals());
}

// delete x randomly chosen fish species
void NewWeb::delete_fish_species_at_random(double x){
  int nooffishspecies=0;
  double dummyindex=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      // Choose a random number representing the fish species to be deleted
      dummyindex = nooffishspecies*unirand();
      // int() always rounds down to nearest integer, regardless of size of fraction
      deletedfishspeciesindex = int(dummyindex);
      //cout << "dummyindex = " << dummyindex << endl;
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}

// delete x randomly chosen fish species
void NewWeb::delete_fish_species_by_size(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double largestM=0;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      largestM=0;
      // Choose the largest existing fish species and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the largest existing fish species
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(s[i].bodymass()>largestM){
	    largestM=s[i].bodymass();
	    deletedfishspeciesindex=dummyindex;
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}

// delete x randomly chosen fish species
void NewWeb::delete_fish_species_by_size_reverse(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double smallestM=0;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      smallestM=fishMupperthreshold;
      // Choose the smallest existing fish species and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the smallest existing fish species
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(s[i].bodymass()<smallestM){
	    smallestM=s[i].bodymass();
	    deletedfishspeciesindex=dummyindex;
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}

// delete x fish species chosen in order of decreasing biomass
void NewWeb::delete_fish_species_by_biomass(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double largestB=0;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      largestB=0;
      // Choose the existing fish species with largest biomass and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with largest biomass
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(s(i).biomass_abundance_B()>largestB){
	    largestB=s(i).biomass_abundance_B();
	    deletedfishspeciesindex=dummyindex;
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}


// delete x fish species chosen in order of increasing biomass
void NewWeb::delete_fish_species_by_biomass_reverse(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double smallestB=0;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      smallestB=0;
      // Choose the existing fish species with smallest biomass and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with smallest biomass
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(dummyindex==0){
	    smallestB=s(i).biomass_abundance_B();
	    deletedfishspeciesindex=dummyindex;
	  }
	  if(dummyindex>0){
	    if(s(i).biomass_abundance_B()<smallestB){
	      smallestB=s(i).biomass_abundance_B();
	      deletedfishspeciesindex=dummyindex;
	    }
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}


// delete x fish species chosen in order of decreasing TL
void NewWeb::delete_fish_species_by_TL(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double largestTL=0;
  // Work out TLs
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewVector level;
  NewMatrix dietF;
  dietF=ifrac;
  NewVector ones(S);
  for(int i=0;i<S;i++){
    ones[i]=1;
  }
  level=ones;
  NewVector level2;
  NewVector leveldiff;
  double dummycounter = 0;
  double dummydiff = 0;
  do{
    level2=level;
    level=dietF*level+ones;
    leveldiff=level-level2;
    dummydiff=NewVectorAbs(leveldiff);
    dummycounter=dummycounter+1;
    //cout << "dummycounter = " << dummycounter << endl;
    //cout << "dummydiff = " << dummydiff << endl;
  }while(dummydiff>=0.0001);
  cout << "dummycounter = " << dummycounter << endl;
  cout << "dummydiff = " << dummydiff << endl;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      largestTL=0;
      // Choose the existing fish species with largest TL and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with largest TL
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(level[i]>largestTL){
	    largestTL=level[i];
	    deletedfishspeciesindex=dummyindex;
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}


// delete x fish species chosen in order of increasing TL
void NewWeb::delete_fish_species_by_TL_reverse(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double smallestTL=0;
  // Work out TLs
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewVector level;
  NewMatrix dietF;
  dietF=ifrac;
  NewVector ones(S);
  for(int i=0;i<S;i++){
    ones[i]=1;
  }
  level=ones;
  NewVector level2;
  NewVector leveldiff;
  double dummycounter = 0;
  double dummydiff = 0;
  do{
    level2=level;
    level=dietF*level+ones;
    leveldiff=level-level2;
    dummydiff=NewVectorAbs(leveldiff);
    dummycounter=dummycounter+1;
    //cout << "dummycounter = " << dummycounter << endl;
    //cout << "dummydiff = " << dummydiff << endl;
  }while(dummydiff>=0.0001);
  cout << "dummycounter = " << dummycounter << endl;
  cout << "dummydiff = " << dummydiff << endl;
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      smallestTL=0;
      // Choose the existing fish species with smallest TL and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with smallest TL
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(dummyindex==0){
	    smallestTL=level[i];
	    deletedfishspeciesindex=dummyindex;
	  }
	  if(dummyindex>0){
	    if(level[i]<smallestTL){
	      smallestTL=level[i];
	      deletedfishspeciesindex=dummyindex;
	    }
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}


// delete x fish species chosen in order of decreasing Conn
void NewWeb::delete_fish_species_by_Conn(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double largestConn=0;
  // Work out Conns
  const double gram=eval_here("1*gram"); //get a dimensionless gram
  const double unit_mass=eval_here("1*kilogram"); //for output files
  const double meter2=eval_here("1*meter^2");
  const double unit_area=eval_here("1*meter^2"); //for output files
  const double year=eval_here("1*year");
  double area=area_per_compartment/unit_area;
  //cout << "area_per_compartment = " << area_per_compartment << endl;
  //cout << "meter2 = " << meter2 << endl;
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  sequence<double> b_idotI;
  sequence<double> b_dotjI;
  double b_dotdotI = 0.0;
  double dummy3 = 0.0;
  double dummy4 = 0.0;
  double resolog10M = 0.0;
  double conslog10M = 0.0;
  sequence<double> Gen_I;
  sequence<double> Vul_I;
  sequence<double> Conn_I;
  double dietfractionthreshold = 0.01; 
  // First populate b_idotI, b_dotjI.
  // Also calcuate b_dotdotI.
  for(int i=0;i<S;i++){
    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    dummy3 = 0;
    dummy4 = 0;
    for(int j=0;j<S;j++){
      NewSpecies & cons=s[j];
      conslog10M = log10(cons.bodymass()/unit_mass);
      dummy3 += intake(i,j)/unit_mass/area/(1/year);
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      // End loop for j.
    }
    b_idotI[i] = dummy3;
    b_dotjI[i] = dummy4;
    // End loop for i.
  }
  // Now populate Gen_I, Vul_I and 
  // Conn_I (= Gen + Vul)
  // for each species.
  double dummy100 = 0.0;
  double dummy101 = 0.0;
  double sumGen = 0.0;
  double sumVul = 0.0;
  double sumConn = 0.0;
  for(int i=0;i<S;i++){
    sumGen = 0;
    sumVul = 0;
    sumConn = 0;
    for(int j=0;j<S;j++){
      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);
      // Work out Gen, Vul and Conn.
      if(b_dotjI[i]>0){
	if((dummy101/b_dotjI[i])>dietfractionthreshold){
	  sumGen = sumGen + 1;
	  sumConn = sumConn + 1;
	}
      }
      if(b_dotjI[j]>0){
	if((dummy100/b_dotjI[j])>dietfractionthreshold){
	  sumVul = sumVul + 1;
	  sumConn = sumConn + 1;
	}
      }
      // End loop for j.
    }
    Gen_I[i] = sumGen;
    Vul_I[i] = sumVul;
    Conn_I[i] = sumConn;
    // End loop for i.
  }
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      largestConn=0;
      // Choose the existing fish species with largest Conn and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with largest Conn
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(Conn_I[i]>largestConn){
	    largestConn=Conn_I[i];
	    deletedfishspeciesindex=dummyindex;
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}


// delete x fish species chosen in order of increasing Conn
void NewWeb::delete_fish_species_by_Conn_reverse(double x){
  int dummyindex=0;
  int nooffishspecies=0;
  int deletedfishspeciesindex=0;
  int indextracker=0;
  double smallestConn=0;
  // Work out Conns
  const double gram=eval_here("1*gram"); //get a dimensionless gram
  const double unit_mass=eval_here("1*kilogram"); //for output files
  const double meter2=eval_here("1*meter^2");
  const double unit_area=eval_here("1*meter^2"); //for output files
  const double year=eval_here("1*year");
  double area=area_per_compartment/unit_area;
  //cout << "area_per_compartment = " << area_per_compartment << endl;
  //cout << "meter2 = " << meter2 << endl;
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  sequence<double> b_idotI;
  sequence<double> b_dotjI;
  double b_dotdotI = 0.0;
  double dummy3 = 0.0;
  double dummy4 = 0.0;
  double resolog10M = 0.0;
  double conslog10M = 0.0;
  sequence<double> Gen_I;
  sequence<double> Vul_I;
  sequence<double> Conn_I;
  double dietfractionthreshold = 0.01; 
  // First populate b_idotI, b_dotjI.
  // Also calcuate b_dotdotI.
  for(int i=0;i<S;i++){
    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    dummy3 = 0;
    dummy4 = 0;
    for(int j=0;j<S;j++){
      NewSpecies & cons=s[j];
      conslog10M = log10(cons.bodymass()/unit_mass);
      dummy3 += intake(i,j)/unit_mass/area/(1/year);
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      // End loop for j.
    }
    b_idotI[i] = dummy3;
    b_dotjI[i] = dummy4;
    // End loop for i.
  }
  // Now populate Gen_I, Vul_I and 
  // Conn_I (= Gen + Vul)
  // for each species.
  double dummy100 = 0.0;
  double dummy101 = 0.0;
  double sumGen = 0.0;
  double sumVul = 0.0;
  double sumConn = 0.0;
  for(int i=0;i<S;i++){
    sumGen = 0;
    sumVul = 0;
    sumConn = 0;
    for(int j=0;j<S;j++){
      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);
      // Work out Gen, Vul and Conn.
      if(b_dotjI[i]>0){
	if((dummy101/b_dotjI[i])>dietfractionthreshold){
	  sumGen = sumGen + 1;
	  sumConn = sumConn + 1;
	}
      }
      if(b_dotjI[j]>0){
	if((dummy100/b_dotjI[j])>dietfractionthreshold){
	  sumVul = sumVul + 1;
	  sumConn = sumConn + 1;
	}
      }
      // End loop for j.
    }
    Gen_I[i] = sumGen;
    Vul_I[i] = sumVul;
    Conn_I[i] = sumConn;
    // End loop for i.
  }
  // Loop to delete fish species - one fish species deleted per iteration
  for(int i=0;i<x;i++){
    // First work out the total number of fish species. If this is zero, then a message is printed out
    nooffishspecies=0;
    for(int i=s.size();i-->0;){
      if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	nooffishspecies+=1;
      }
    }
    //cout << "nooffishspecies = " << nooffishspecies << endl;
    if(nooffishspecies==0){
      cout << "No fish species left" << endl;
    }else{
      dummyindex=0;
      smallestConn=0;
      // Choose the existing fish species with smallest Conn and delete it
      // First find the index (deletedfishspeciesindex) corresponding to the existing fish species
      // with smallest Conn
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  if(dummyindex==0){
	    smallestConn=Conn_I[i];
	    deletedfishspeciesindex=dummyindex;
	  }
	  if(dummyindex>0){
	    if(Conn_I[i]<smallestConn){
	      smallestConn=Conn_I[i];
	      deletedfishspeciesindex=dummyindex;
	    }
	  }
	  dummyindex+=1;
	}
      }
      //cout << "deletedfishspeciesindex = " << deletedfishspeciesindex << endl;
      indextracker=0;
      for(int i=s.size();i-->0;){
	if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
	  // If indextracker is equal to deletedfishspeciesindex, then current fish species 
	  // is deleted and indextracker is incremented by 1; else, indextracker is 
	  // incremented by 1 without species deletion
	  if(indextracker==deletedfishspeciesindex){
	    //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
	    cout << "deleting " 
		 << s[i].taxon_name()
		 << " " << assigned_column[i] << endl;
	    delete_species(i);
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }else{
	    indextracker+=1;
	    //cout << "indextracker = " << indextracker << endl;
	  }
	}
      }
    }
  }
}

// delete last fish species making up some % of total fish biomass when arranged in descending 
// order of biomass
void NewWeb::delete_rare_fish_species(double percentage){
  sequence< double > fishbiomasses;
  double totalfishbiomass=0;
  double fractionthreshold=(100-percentage)/100;
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
      fishbiomasses.push_back(s(i).biomass_abundance_B());
      totalfishbiomass+=s(i).biomass_abundance_B();
    }    
  }
  //cout << "totalfishbiomass = " << totalfishbiomass << endl;
  double totalfishbiomassthreshold=fractionthreshold*totalfishbiomass;
  double fishbiomassthreshold=0;
  int fishbiomassthresholdindex=0;
  //cout << "totalfishbiomassthreshold = " << totalfishbiomassthreshold << endl;
  // sort fish biomasses in descending order (because we are working backwards in i)
  sort(fishbiomasses.begin(),fishbiomasses.end());
  totalfishbiomass=0;
  for(int i=fishbiomasses.size();i-->0;){
    //cout << "fishbiomasses.at(" << i << ") = " << fishbiomasses.at(i) << endl;   
    if((totalfishbiomass>totalfishbiomassthreshold)&&(fishbiomassthresholdindex==0)){
      fishbiomassthreshold=fishbiomasses.at(i);
      fishbiomassthresholdindex=1;
    }
    totalfishbiomass+=fishbiomasses.at(i);
    //cout << "totalfishbiomass = " << totalfishbiomass << endl;   
  }
  //cout << "fishbiomassthreshold = " << fishbiomassthreshold << endl;
  //cout << "fishbiomassthresholdindex = " << fishbiomassthresholdindex << endl;
  // now delete rare fish species
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)&&(s(i).biomass_abundance_B()<=fishbiomassthreshold)){
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

// delete last fish species making up 5% of total fish biomass when arranged in descending 
// order of biomass, but using a positive sloping cut-off rather than a horizontal cut-off
void NewWeb::delete_rare_fish_species_5per_2(double slope){
  // Rotated fish biomasses are fish biomasses rotated clockwise through angle defined by slope 
  // We only care about the rotated biomass, not the rotated body mass,
  // but the rotated biomass depends on original body mass.
  // Since line that does cutting-off is in log10 kg and log10 kg m-2, 
  // need to work in these units.
  // Rotated fish biomasses vector also has original fish biomass attached to each value.
  // There is also a vector of rotated biomasses for all species, required at the end.
  sequence< double > fishbiomasses;
  double totalfishbiomass=0;
  double percentage = 5;
  double fractionthreshold=(100-percentage)/100;
  sequence< pair<double,double> > rotatedlog10fishbiomassesunitmassunitarea;
  sequence< double > rotatedlog10biomassesunitmassunitarea;
  double theta=atan(slope);
  double massunit = 1000;
  double areaunit = 1;
  //cout << "theta = " << theta << endl;
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
      fishbiomasses.push_back(s(i).biomass_abundance_B());
      rotatedlog10fishbiomassesunitmassunitarea.push_back(make_pair(
								    -(log10((s(i).bodymass())/massunit))*sin(theta)
								    +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta)
								    ,s(i).biomass_abundance_B()));
      totalfishbiomass+=s(i).biomass_abundance_B();
    }
    rotatedlog10biomassesunitmassunitarea.push_back(
						    -(log10((s(i).bodymass())/massunit))*sin(theta)
						    +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta)
						    );
  }
  // Re-populate vector rotatedlog10biomassesunitmassunitarea to get elements in right order - this time, at(i) works because
  // there are existing elements
  for(int i=s.size();i-->0;){
    rotatedlog10biomassesunitmassunitarea.at(i)=-(log10((s(i).bodymass())/massunit))*sin(theta)
      +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta);
  }
  // For rotated fish biomasses vector,
  // sort in descending order (because we are working backwards in i)
  // then work out cut-off point above which there is (100-percentage)% of fish biomass
  // this is done in a loop where we go from fish species with highest biomass
  // down to cut-off point (as opposed to other way round), so always get
  // at least (100-percentage)% of total fish biomass
  double totalfishbiomassthreshold=fractionthreshold*totalfishbiomass;
  double rotatedlog10fishbiomassthresholdunitmassunitarea=0;
  int rotatedlog10fishbiomassthresholdunitmassunitareaindex=0;
  //cout << "totalfishbiomass = " << totalfishbiomass << endl;
  //cout << "totalfishbiomassthreshold = " << totalfishbiomassthreshold << endl;
  sort(rotatedlog10fishbiomassesunitmassunitarea.begin(),rotatedlog10fishbiomassesunitmassunitarea.end());
  totalfishbiomass=0;
  for(int i=rotatedlog10fishbiomassesunitmassunitarea.size();i-->0;){
    //cout << "rotatedlog10fishbiomassesunitmassunitarea.at(" << i << ") = " << rotatedlog10fishbiomassesunitmassunitarea.at(i).first << endl;   
    if((totalfishbiomass>totalfishbiomassthreshold)&&(rotatedlog10fishbiomassthresholdunitmassunitareaindex==0)){
      rotatedlog10fishbiomassthresholdunitmassunitarea=rotatedlog10fishbiomassesunitmassunitarea.at(i).first;
      rotatedlog10fishbiomassthresholdunitmassunitareaindex=1;
    }
    totalfishbiomass+=rotatedlog10fishbiomassesunitmassunitarea.at(i).second;
    //cout << "totalfishbiomass = " << totalfishbiomass << endl;   
  }
  //cout << "rotatedlog10fishbiomassthresholdunitmassunitarea = " << rotatedlog10fishbiomassthresholdunitmassunitarea << endl;
  //cout << "rotatedlog10fishbiomassthresholdunitmassunitareaindex = " << rotatedlog10fishbiomassthresholdunitmassunitareaindex << endl;
  // now delete rare fish species
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)&&(rotatedlog10biomassesunitmassunitarea.at(i)<=rotatedlog10fishbiomassthresholdunitmassunitarea)){
      //cout << "rotatedlog10biomassesunitmassunitarea.at(" << i << ") = " << rotatedlog10biomassesunitmassunitarea.at(i) << endl;
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

// delete last fish species making up 10% of total fish biomass when arranged in descending 
// order of biomass, but using a positive sloping cut-off rather than a horizontal cut-off
void NewWeb::delete_rare_fish_species_10per_2(double slope){
  // Rotated fish biomasses are fish biomasses rotated clockwise through angle defined by slope 
  // We only care about the rotated biomass, not the rotated body mass,
  // but the rotated biomass depends on original body mass.
  // Since line that does cutting-off is in log10 kg and log10 kg m-2, 
  // need to work in these units.
  // Rotated fish biomasses vector also has original fish biomass attached to each value.
  // There is also a vector of rotated biomasses for all species, required at the end.
  sequence< double > fishbiomasses;
  double totalfishbiomass=0;
  double percentage = 10;
  double fractionthreshold=(100-percentage)/100;
  sequence< pair<double,double> > rotatedlog10fishbiomassesunitmassunitarea;
  sequence< double > rotatedlog10biomassesunitmassunitarea;
  double theta=atan(slope);
  double massunit = 1000;
  double areaunit = 1;
  //cout << "theta = " << theta << endl;
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
      fishbiomasses.push_back(s(i).biomass_abundance_B());
      rotatedlog10fishbiomassesunitmassunitarea.push_back(make_pair(
								    -(log10((s(i).bodymass())/massunit))*sin(theta)
								    +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta)
								    ,s(i).biomass_abundance_B()));
      totalfishbiomass+=s(i).biomass_abundance_B();
    }
    rotatedlog10biomassesunitmassunitarea.push_back(
						    -(log10((s(i).bodymass())/massunit))*sin(theta)
						    +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta)
						    );
  }
  // Re-populate vector rotatedlog10biomassesunitmassunitarea to get elements in right order - this time, at(i) works because
  // there are existing elements
  for(int i=s.size();i-->0;){
    rotatedlog10biomassesunitmassunitarea.at(i)=-(log10((s(i).bodymass())/massunit))*sin(theta)
      +(log10(((s(i).biomass_abundance_B())/massunit)/areaunit))*cos(theta);
  }
  // For rotated fish biomasses vector,
  // sort in descending order (because we are working backwards in i)
  // then work out cut-off point above which there is (100-percentage)% of fish biomass
  // this is done in a loop where we go from fish species with highest biomass
  // down to cut-off point (as opposed to other way round), so always get
  // at least (100-percentage)% of total fish biomass
  double totalfishbiomassthreshold=fractionthreshold*totalfishbiomass;
  double rotatedlog10fishbiomassthresholdunitmassunitarea=0;
  int rotatedlog10fishbiomassthresholdunitmassunitareaindex=0;
  //cout << "totalfishbiomass = " << totalfishbiomass << endl;
  //cout << "totalfishbiomassthreshold = " << totalfishbiomassthreshold << endl;
  sort(rotatedlog10fishbiomassesunitmassunitarea.begin(),rotatedlog10fishbiomassesunitmassunitarea.end());
  totalfishbiomass=0;
  for(int i=rotatedlog10fishbiomassesunitmassunitarea.size();i-->0;){
    //cout << "rotatedlog10fishbiomassesunitmassunitarea.at(" << i << ") = " << rotatedlog10fishbiomassesunitmassunitarea.at(i).first << endl;   
    if((totalfishbiomass>totalfishbiomassthreshold)&&(rotatedlog10fishbiomassthresholdunitmassunitareaindex==0)){
      rotatedlog10fishbiomassthresholdunitmassunitarea=rotatedlog10fishbiomassesunitmassunitarea.at(i).first;
      rotatedlog10fishbiomassthresholdunitmassunitareaindex=1;
    }
    totalfishbiomass+=rotatedlog10fishbiomassesunitmassunitarea.at(i).second;
    //cout << "totalfishbiomass = " << totalfishbiomass << endl;   
  }
  //cout << "rotatedlog10fishbiomassthresholdunitmassunitarea = " << rotatedlog10fishbiomassthresholdunitmassunitarea << endl;
  //cout << "rotatedlog10fishbiomassthresholdunitmassunitareaindex = " << rotatedlog10fishbiomassthresholdunitmassunitareaindex << endl;
  // now delete rare fish species
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)&&(rotatedlog10biomassesunitmassunitarea.at(i)<=rotatedlog10fishbiomassthresholdunitmassunitarea)){
      //cout << "rotatedlog10biomassesunitmassunitarea.at(" << i << ") = " << rotatedlog10biomassesunitmassunitarea.at(i) << endl;
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

void NewWeb::delete_all_species_with_biomass_less_than(double m){
  for(int i=s.size();i-->0;){
    if(!isnormal(s(i).biomass_abundance_B()) || s(i).biomass_abundance_B() < m ){
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i]
	   << " with body mass " << s[i].bodymass() << endl;
      delete_species(i);
    }
  }
}

// delete all fish species except the one with assigned column number x
void NewWeb::delete_fish_species_except_1(double x){
  //cout << "x = " << x << endl;
  for(int i=s.size();i-->0;){
    //cout << "i = " << i << endl;
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)&&(assigned_column[i]!=x)){
      //cout << "deleted species s[i].bodymass() = " << s[i].bodymass() << endl;
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

void NewWeb::
delete_all_plants_with_negative_growth_rates(species_set_t& inserted){
  for(int i=s.size();i-->s.n_animals;){
    if(s(i).plant_growth_rate_sigma()<0){
      inserted.erase(assigned_column[i]);
      delete_species(i);
      species_set_t corrected;
    }
  }
}

NewWeb::species_set_t
NewWeb::delete_species_larger_than_exp(const sequence<double> & si,
				       const species_set_t&  conserved){
  species_set_t deleted;
  for(int i=s.size();i-->0;){
    // must be downwards loop to work!!
    if(!isnormal(s(i).biomass_abundance_B()) || 
       si(i)<(log(s[i].get_threshold_biomass_B()/
		  (conserved.count(assigned_column[i])?conservation_tolerance:1) 
		  ))){
      std::cout << s(i).taxon_name() << " "
		<< assigned_column(i) << " went extinct (2)" << std::endl;
      deleted.insert(assigned_column[i]);
      delete_species(i);
    }
  }
  return deleted;
}
  

relaxing_dynamical_object::species_set_t NewWeb::
delete_all_species_with_less_than_one_individual(const species_set_t
						 & conserved){
  const double tolerance = 1e-5;
  bool OK=true;

  // not OK, do cleanup
  species_set_t deleted;
  for(int i=s.size();i-->0;){
    if(!isnormal(s(i).biomass_abundance_B()) 
       || 
       (1+tolerance)*s(i).biomass_abundance_B()*
       (conserved.count(assigned_column[i])?conservation_tolerance:1) < 
       s(i).get_threshold_biomass_B() ){
      std::cout << s(i).taxon_name() << " "
	   << assigned_column(i) << " went extinct (1)" << std::endl;
      deleted.insert(assigned_column[i]);
      delete_species(i);
    }
  }
  return deleted;
}

void NewWeb::delete_all_animals(){
  for(int i=s.size();i-->0;){
    // must be downwards loop to work!!
    if(!is_plant(i)){
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

void NewWeb::delete_all_species_with_level_higher_than(double h){

  // compute the trophic level vector:
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewMatrix F=ifrac;
  NewVector ones(S);
  int err;
  for(int i=S;i-->0;)
    ones[i]=1;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=solve(D,ones);
  if(err) FATAL_ERROR("could not compute trophic levels");

  for(int i=s.size();i-->0;){
    // must be downwards loop to work!!
    if(level[i]>h){
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

void NewWeb::scale_primary_productivity(double x){
  for(int i=s.size();i-->s.n_animals;){
    s[i].set_plant_growth_rate_sigma(x*s[i].plant_growth_rate_sigma());
  }
}


// See if there are plants larger than the largest animals::
bool NewWeb::
large_plants(){
  if(number_of_animals()==0) return false;

  double max_animal_body_mass=0;
  for(int i=s.n_animals;i-->0;){
    max_animal_body_mass=
      max(s[i].bodymass(),max_animal_body_mass);
  }
  for(int i=s.size();i-->s.n_animals;){
    if(s[i].bodymass()>max_animal_body_mass)
      return true;
  }
  return false;
}

// quick and dirty implementation!!
void NewWeb::fish(double n){
  // increase mortality of an abundand large fish species by n
  double max_fishing_value=0;
  int max_species=0;
  for(int i=s.n_animals;i-->0;){
    double fishing_value=s[i].biomass_abundance_B()*
      s[i].bodymass();
    if(fishing_value>max_fishing_value){
      max_fishing_value=fishing_value;
      max_species=i;
    }
  }
  cout << "fishing target is " << assigned_column[max_species] 
       << " (internal " << max_species << ")"
       << endl;
  //double m=s[max_species].turnover_rate_r()*(n);
  s[max_species].set_fishing_mortality(n);
}

void NewWeb::forbid_fishing(){
  for(int i=s.size();i-->0;){
    s(i).set_fishing_mortality(0);
  }
}

void NewWeb::copy_fishing_mortalities_from(const NewWeb & other){
  for(int i=other.s.size();i-->0;){
    const int c=other.assigned_column[i];
    if(c >= species_at.size() ||
       species_at[c]==column_unused ){
      WARNING("no correspondence for fished species " << c);
    }else{
      s(species_at[c]).set_fishing_mortality(other.s(i).fishing_mortality());
    }
  }
}


void NewWeb::fishF(double F){
  // Fish all animals > 1g with F
  for(int i=s.n_animals;i-->0;){
    if(s[i].bodymass() >= eval_here("1*gram")){
      s[i].set_fishing_mortality(F);
    }
  } 
}

void NewWeb::fishall(double n){
  // increase mortality of all fish species with bodymass >n kg
  double fishing_value_intercept=0;
  double fishing_value=0;
  for(int i=s.n_animals;i-->0;){
    if(s[i].bodymass() > n*eval_here("kilogram")){
      fishing_value_intercept=fishall_intercept+gaussian(0,fishall_residual_se);
      fishing_value=pow(10,fishing_value_intercept)*pow((s[i].bodymass())/eval_here("kilogram"),fishall_slope);
      //cout << "fishing_value is " << fishing_value << endl;
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

void NewWeb::fishMrange(double n){
  // increase mortality of all fish species with bodymass >n kg
  // and <10^fishMrangemax kg
  double fishing_value_intercept=0;
  double fishing_value=0;
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      fishing_value_intercept=fishall_intercept+gaussian(0,fishall_residual_se);
      fishing_value=pow(10,fishing_value_intercept)*pow((s[i].bodymass())/eval_here("kilogram"),fishall_slope);
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

void NewWeb::custom_perturbation(double x){
  // increase mortality of all fish species with bodymass >n kg
  double fishing_value_intercept=0;
  double fishing_value=0;
  for(int i=s.n_animals;i-->0;){
    if(s[i].bodymass() > eval_here("10*gram")){
      s[i].set_fishing_mortality(eval_here("0.5/year"));
     }
  } 
}


void NewWeb::fishTLrange(double n){
  // increase mortality of all fish species with TL >n and < fishTLrangemax
  double fishing_value_intercept=0;
  double fishing_value=0;
  // compute the trophic level vector:
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewMatrix F=ifrac;
  NewVector ones(S);
  for(int i=S;i-->0;)
    ones[i]=1;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=solve(D,ones);
  // apply fishing
  for(int i=s.n_animals;i-->0;){
    if((level[i]>n)&&(level[i]<fishTLrangemax)){
      fishing_value_intercept=fishall_intercept+gaussian(0,fishall_residual_se);
      fishing_value=pow(10,fishing_value_intercept)*pow((s[i].bodymass())/eval_here("kilogram"),fishall_slope);
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

void NewWeb::fishMrangeconst(double n){
  // increase mortality of all fish species with bodymass >n kg
  // and <fishMrangemax kg
  double fishing_value=0;
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      fishing_value=const_fishing_mortality_rate;
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

// fishMrangeconst2 same as fishMrangeconst except that F can be chosen
void NewWeb::fishMrangeconst2(double n, double p){
  // increase mortality of all fish species with bodymass >n kg
  // and <fishMrangemax kg
  double fishing_value=0;
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      fishing_value=p;
      //cout << "fishing_value = " << p << endl;
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

void NewWeb::fishTLrangeconst(double n){
  // increase mortality of all fish species with TL >n and < fishTLrangemax
  double fishing_value=0;
  // compute the trophic level vector:
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewMatrix F=ifrac;
  NewVector ones(S);
  for(int i=S;i-->0;)
    ones[i]=1;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=solve(D,ones);
  // apply fishing
  for(int i=s.n_animals;i-->0;){
    if((level[i]>n)&&(level[i]<fishTLrangemax)){
      fishing_value=const_fishing_mortality_rate;
      s[i].set_fishing_mortality(fishing_value);
     }
  } 
}

void NewWeb::fishMrangeline(double n){
  // increase mortality of all fish species with bodymass >n kg
  // and <fishMrangemax kg
  // An upper limit Fmax is set to the fishing mortality
  double fishing_value=0;
  double fishing_value2=0;
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      fishing_value=pow(10,fishall_intercept)*pow((s[i].bodymass())/eval_here("kilogram"),fishall_slope);
      fishing_value2=min(fishing_value,Fmax);
      s[i].set_fishing_mortality(fishing_value2);
     }
  } 
}

void NewWeb::fishTLrangeline(double n){
  // increase mortality of all fish species with TL >n and < fishTLrangemax
  // An upper limit Fmax is set to the fishing mortality
  double fishing_value=0;
  double fishing_value2=0;
  // compute the trophic level vector:
  link_strength_matrix intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  NewMatrix F=ifrac;
  NewVector ones(S);
  for(int i=S;i-->0;)
    ones[i]=1;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=solve(D,ones);
  // apply fishing
  for(int i=s.n_animals;i-->0;){
    if((level[i]>n)&&(level[i]<fishTLrangemax)){
      fishing_value=pow(10,fishall_intercept)*pow((s[i].bodymass())/eval_here("kilogram"),fishall_slope);
      fishing_value2=min(fishing_value,Fmax);
      s[i].set_fishing_mortality(fishing_value2);
     }
  } 
}


void NewWeb::fishMrangeCS(double n){
  // species1Mdiff is the list of body mass differences between species 1 and all other fish species.
  // species1Mdiffsorted is species1Mdiff sorted in ascending order.
  // species1M is the body mass of species 1.
  // species1F is F for species 1.
  // indices is the list of indices for the fish species. 
  // newindices1 is the list of indices for species 1 once they have been sorted.
  sequence<double> species1Mdiff;
  sequence<double> species2Mdiff;
  sequence<double> species3Mdiff;
  sequence<double> species4Mdiff;
  sequence<double> species5Mdiff;
  sequence<double> species6Mdiff;
  sequence<double> species1Mdiffsorted;
  sequence<double> species2Mdiffsorted;
  sequence<double> species3Mdiffsorted;
  sequence<double> species4Mdiffsorted;
  sequence<double> species5Mdiffsorted;
  sequence<double> species6Mdiffsorted;
  double species1M = 0.0691*eval_here("kilogram");
  double species2M = 0.144*eval_here("kilogram");
  double species3M = 0.535*eval_here("kilogram");
  double species4M = 0.626*eval_here("kilogram");
  double species5M = 1.01*eval_here("kilogram");
  double species6M = 2.10*eval_here("kilogram");
  double species1F = fishMrangeCSFfactor*0.437*eval_here("1/year");
  double species2F = fishMrangeCSFfactor*0.277*eval_here("1/year");
  double species3F = fishMrangeCSFfactor*0.566*eval_here("1/year");
  double species4F = fishMrangeCSFfactor*0.255*eval_here("1/year");
  double species5F = fishMrangeCSFfactor*0.714*eval_here("1/year");
  double species6F = fishMrangeCSFfactor*0.234*eval_here("1/year");
  sequence<int> indices;
  sequence<int> newindices1;
  sequence<int> newindices2;
  sequence<int> newindices3;
  sequence<int> newindices4;
  sequence<int> newindices5;
  sequence<int> newindices6;
  int count1 = 0;

  // For each model fish species, compute differences between body mass and body masses of 6 species
  // Note that the body mass of each species is in g. *eval_here("kilogram") above 
  // converts all the body masses into g automatically by multiplying by 1000. 
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      species1Mdiff[count1]=fabs(s[i].bodymass()-species1M);
      species2Mdiff[count1]=fabs(s[i].bodymass()-species2M);
      species3Mdiff[count1]=fabs(s[i].bodymass()-species3M);
      species4Mdiff[count1]=fabs(s[i].bodymass()-species4M);
      species5Mdiff[count1]=fabs(s[i].bodymass()-species5M);
      species6Mdiff[count1]=fabs(s[i].bodymass()-species6M);
      species1Mdiffsorted[count1]=fabs(s[i].bodymass()-species1M);
      species2Mdiffsorted[count1]=fabs(s[i].bodymass()-species2M);
      species3Mdiffsorted[count1]=fabs(s[i].bodymass()-species3M);
      species4Mdiffsorted[count1]=fabs(s[i].bodymass()-species4M);
      species5Mdiffsorted[count1]=fabs(s[i].bodymass()-species5M);
      species6Mdiffsorted[count1]=fabs(s[i].bodymass()-species6M);
      indices[count1]=i;
      count1++;      
     }
  } 

  // Now sort each vector of body mass differences (ascending order). 
  sort(species1Mdiffsorted.begin(),species1Mdiffsorted.end());
  sort(species2Mdiffsorted.begin(),species2Mdiffsorted.end());
  sort(species3Mdiffsorted.begin(),species3Mdiffsorted.end());
  sort(species4Mdiffsorted.begin(),species4Mdiffsorted.end());
  sort(species5Mdiffsorted.begin(),species5Mdiffsorted.end());
  sort(species6Mdiffsorted.begin(),species6Mdiffsorted.end());

  // assign new indices by comparing sorted 
  for(int i=0;i<count1;i++){
    for(int j=0;j<count1;j++){
      if(species1Mdiff[i]==species1Mdiffsorted[j]) {
	newindices1[j]=indices[i];
      }
      if(species2Mdiff[i]==species2Mdiffsorted[j]) {
	newindices2[j]=indices[i];
      }
      if(species3Mdiff[i]==species3Mdiffsorted[j]) {
	newindices3[j]=indices[i];
      }
      if(species4Mdiff[i]==species4Mdiffsorted[j]) {
	newindices4[j]=indices[i];
      }
      if(species5Mdiff[i]==species5Mdiffsorted[j]) {
	newindices5[j]=indices[i];
      }
      if(species6Mdiff[i]==species6Mdiffsorted[j]) {
	newindices6[j]=indices[i];
      }
    }
  }

  // Now we calculate the sum of the (absolute) differences for each combination of the first 6 elements of each 
  // difference list (stored in sumdiff). Only combinations that contain no duplicate species are considered.
  // count2 is used to check if the combinations have duplicates - if it is >0, then there is at least one duplicate.
  // combination is the vector of values for a combination. combinationindices is vector of corresponding species indices.
  // lowestdiff gives the body mass differences giving lowestsumdiff
  // bestindices is the list of indices that gives the smallest sum of differences.
  double sumdiff = 0.0; 
  double lowestsumdiff = 1e100;
  int count2 = 0;
  sequence<double> combination;
  sequence<int> combinationindices;
  sequence<double> lowestdiff;
  sequence<int> bestindices;
  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      for(int k=0;k<6;k++){
	for(int l=0;l<6;l++){
	  for(int m=0;m<6;m++){
	    for(int n=0;n<6;n++){
	      combination[0]=species1Mdiffsorted[i];
	      combination[1]=species2Mdiffsorted[j];
	      combination[2]=species3Mdiffsorted[k];
	      combination[3]=species4Mdiffsorted[l];
	      combination[4]=species5Mdiffsorted[m];
	      combination[5]=species6Mdiffsorted[n];
	      combinationindices[0]=newindices1[i];
	      combinationindices[1]=newindices2[j];
	      combinationindices[2]=newindices3[k];
	      combinationindices[3]=newindices4[l];
	      combinationindices[4]=newindices5[m];
	      combinationindices[5]=newindices6[n];
	      // Check for duplicates for current combination.
	      for(int o=0;o<6;o++){
		for(int p=o+1;p<6;p++){
		  if(combinationindices[o]==combinationindices[p]){count2++;}
		}  
	      }
	      // If there are no duplicates, then work out the sum of differences.
	      // If there are duplicates, then sumdiff is set to a very large number.
	      if(count2==0){
		for(int q=0;q<6;q++){
		  sumdiff += combination[q]; 
		}
	      } else {
		sumdiff = 1e100;
	      }
	      // Update lowestsumdiff and bestindices if necessary.
	      if(sumdiff<lowestsumdiff){
		lowestsumdiff=sumdiff;
		lowestdiff[0]=combination[0];
		lowestdiff[1]=combination[1];
		lowestdiff[2]=combination[2];
		lowestdiff[3]=combination[3];
		lowestdiff[4]=combination[4];
		lowestdiff[5]=combination[5];
		bestindices[0]=combinationindices[0];
		bestindices[1]=combinationindices[1];
		bestindices[2]=combinationindices[2];
		bestindices[3]=combinationindices[3];
		bestindices[4]=combinationindices[4];
		bestindices[5]=combinationindices[5];
	      }
	      sumdiff = 0.0;
	      count2 = 0;
	    }
	  }
	}
      }
    }
  }

  cout << "lowestsumdiff = " << lowestsumdiff << endl;

  int bestindex = 0;
  for(int i=0;i<6;i++){
    bestindex = bestindices[i];
    cout << "body mass of species giving lowest difference with " << i << "th species = " << s[bestindex].bodymass() << endl;
    cout << i << "th element of lowestdiff = " << lowestdiff[i] << endl;
  }

  // Now we apply F to all model fish species. 
  double fishing_value=0;
  // First, F is applied to all species according to straight lines between the 6 known F values. Equation of a straight line used. 
  // For model fish species smaller than the smallest species with known F, the line derived from the 2 smallest species with known F is used. 
  // For model fish species larger than the largest species with known F, the line derived from the 2 largest species with known F is used.
  // A random gaussian variable is also added if addrandF is 1, with mean 0 and std randFstd.
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass() > n*eval_here("kilogram")) 
       && (s[i].bodymass() < fishMrangemax*eval_here("kilogram"))){
      if(s[i].bodymass() <= species2M){ fishing_value = ((species2F-species1F)/(species2M-species1M))*(s[i].bodymass()-species1M)+species1F; }
      if((s[i].bodymass() > species2M)&&(s[i].bodymass() <= species3M)){ fishing_value = ((species3F-species2F)/(species3M-species2M))*(s[i].bodymass()-species2M)+species2F; }
      if((s[i].bodymass() > species3M)&&(s[i].bodymass() <= species4M)){ fishing_value = ((species4F-species3F)/(species4M-species3M))*(s[i].bodymass()-species3M)+species3F; }
      if((s[i].bodymass() > species4M)&&(s[i].bodymass() <= species5M)){ fishing_value = ((species5F-species4F)/(species5M-species4M))*(s[i].bodymass()-species4M)+species4F; }
      if(s[i].bodymass() > species5M){ fishing_value = ((species6F-species5F)/(species6M-species5M))*(s[i].bodymass()-species5M)+species5F; }
      if(addrandF==1){ fishing_value = fishing_value + gaussian(0,randFstd); }
      // The final line has a negative slope, so we need to ensure that the fishing_value does not become negative. 
      // Also, the addition of a negative Guassian random number may give a negative fishing value, so again, we need to ensure it does not become negative.
      if(fishing_value<0){ fishing_value=0; }; 
      // cout << "fishing_value = " << fishing_value << endl;
      s[i].set_fishing_mortality(fishing_value);
    }
  }
  // Now we apply the known F values to the 6 model species with closest body mass.
  int index1 = 0;
  index1 = bestindices[0];
  s[index1].set_fishing_mortality(species1F);
  index1 = bestindices[1];
  s[index1].set_fishing_mortality(species2F);
  index1 = bestindices[2];
  s[index1].set_fishing_mortality(species3F);
  index1 = bestindices[3];
  s[index1].set_fishing_mortality(species4F);
  index1 = bestindices[4];
  s[index1].set_fishing_mortality(species5F);
  index1 = bestindices[5];
  s[index1].set_fishing_mortality(species6F);
  
}


void NewWeb::delete_plant_fraction(double r){

  ALWAYS_ASSERT(r<1);
  int plants_to_delete=int(r*number_of_plants());
  
  while(plants_to_delete--){
    int n=s.n_animals+random_integer(number_of_plants());
    cout << "deleting plant " << assigned_column[n] << endl;
    delete_species(n);
  }
}

void NewWeb::
delete_lower_fraction_of_property(const sequence<double> & property,
				  double r){
  //quick and dirty implementation
  ALWAYS_ASSERT(property.size()==number_of_species());
  ALWAYS_ASSERT(0<=r and r<1);
  sequence<double> prop_copy=property;
  sort(prop_copy.begin(),prop_copy.end());
  const double threshold=prop_copy[number_of_species()*r];
  REPORT(threshold);
  // This must be a backwards loop!
  for(int i=number_of_species();i-->0;){
    if(property[i]<threshold){
      cout << "deleting " 
	   << s[i].taxon_name()
	   << " " << assigned_column[i] << endl;
      delete_species(i);
    }
  }
}

void NewWeb::delete_plants_by_vulnerability_trait(double r){
  sequence< double > V1(number_of_species());
  double V1_max=0;// it's no problem if all plant V1 are < 0 ...
  for(int i=number_of_species();i-->s.n_animals;){
    V1[i]=s[i].vulnerability_V()[1];
    if(V1[i]>V1_max) V1_max=V1[i];
  }
  for(int i=s.n_animals;i-->0;){
    V1[i]=V1_max;
  }
  double fraction=
    (r*number_of_plants())/
    number_of_species();
  delete_lower_fraction_of_property(V1,fraction);
}

void NewWeb::delete_rare_species_fraction(double r){

  // compute cutoff body mass:
  sequence<double> Biomass(s.size());
  for(int i=s.size();i-->0;){
    Biomass(i)=s(i).biomass_abundance_B();
  }
  delete_lower_fraction_of_property(Biomass,r);
}

void NewWeb::delete_large_species_fraction(double r){

  // compute cutoff body mass:
  sequence<double> neg_bodymass=get_bodymass_M()*-1;
  delete_lower_fraction_of_property(neg_bodymass,r);
}

void NewWeb::delete_10_species_near_size(double size){
  // compute cutoff body mass:
  sequence<double> diff_bodymass=get_bodymass_M()-size;
  diff_bodymass*=diff_bodymass;
  delete_lower_fraction_of_property(diff_bodymass,10.5/number_of_species());
}


Interaction_Matrix NewWeb::trivial_loop_removal(const Interaction_Matrix & im){
  ALWAYS_ASSERT(im.size()==s.size());
  const int S=s.size();
  Interaction_Matrix im2(im);
  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      if((!s[i].is_a_plant()) && (!s[j].is_a_plant()) && 
	 s[i].bodymass() <= s[j].bodymass() ){
	im2[i][j]=NetworkAnalysis::none;
      }
    }
  }
  return im2;
}

bool NewWeb::small_values_in(ODE_vector & state,
			     const species_set_t&  conserved){
  const double log_tolerance = log_small_value_tolerance;
  double log_conservation_tolerance=log(conservation_tolerance);
  for(int j=0;j<s.size();j++){
    //cout << state[j] << ", " << s(j).the_log_mean_bodymass_M << std::endl;
    if(state[j] < log(s(j).get_threshold_biomass_B()) - log_tolerance -
       (conserved.count(assigned_column[j])?log_conservation_tolerance:0))
      return true;
  }
  return false;
}

int NewWeb::number_of_species() const{
  return s.size();
}

  
int NewWeb::number_of_animals() const{
  return s.n_animals;
}

int NewWeb::number_of_plants() const{
  return s.size()-s.n_animals;
}

int NewWeb::number_of_fish() const{
  int fishcount=0;  
  for(int i=s.size();i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
      fishcount++;
    }
  }
  return fishcount;
}

const double gram=eval_here("1*gram"); //get a dimensionless gram
const double unit_mass=eval_here("1*kilogram"); //for output files
const double meter2=eval_here("1*meter^2");
const double unit_area=eval_here("1*meter^2"); //for output files
const double year=eval_here("1*year");

void NewWeb::rank_abundance_plot(const char* filename,bool with_plants){
  std::ofstream os(filename);
  sequence<double> abund;
  for(int i=s.size();i-->0;){
    if(with_plants || ! s[i].is_a_plant()){
      abund[abund.size()]=s[i].biomass_abundance_B()/
	s[i].bodymass();
    }
  }
  sort(abund.begin(),abund.end());
  for(int i=abund.size();i-->0;){
    os << abund.size()-i << " " << abund[i] << std::endl;
  }

  abund.resize(0);
  for(int i=s.n_animals;i-->0;){
    abund[abund.size()]=s[i].biomass_abundance_B();
  }
  sort(abund.begin(),abund.end());
  std::ostringstream afname;
  afname << "a" << filename;
  std::ofstream as(afname.str().c_str());
  for(int i=abund.size();i-->0;){
    as << abund.size()-i << " " << abund[i] << std::endl;
  }

  if(with_plants){
    abund.resize(0);
    for(int i=s.size();i-->s.n_animals;){
      abund[abund.size()]=s[i].biomass_abundance_B();
    }
    sort(abund.begin(),abund.end());
    std::ostringstream pfname;
    pfname << "p" << filename;
    std::ofstream ps(pfname.str().c_str());
    for(int i=abund.size();i-->0;){
      ps << abund.size()-i << " " << abund[i] << std::endl;
    }
  }
}

void NewWeb::size_spectrum(const char* filename,bool plants,bool animals){
  sequence<double> abundance;
  sequence<double> bodymass;
  double area=area_per_compartment/unit_area;

  int j=0;
  for(int i=s.size();i-->0;){
    if((s[i].is_a_plant() && plants) || (!s[i].is_a_plant() && animals)){
      abundance[j]=s[i].biomass_abundance_B()/s[i].bodymass()/
	area;
      bodymass[j]=s[i].bodymass()/
	unit_mass;
      j++;
    }
  }
  std::ofstream os(filename);
  os << log_spectrum(xy_graph(bodymass,abundance),10).log_xy();
}

void NewWeb::biomass_spectrum(const char* filename,bool plants,bool animals){
  sequence<double> biomass;
  sequence<double> bodymass;
  double area=area_per_compartment/unit_area;

  int j=0;
  for(int i=s.size();i-->0;){
    if((s[i].is_a_plant() && plants) || (!s[i].is_a_plant() && animals)){
      biomass[j]=s[i].biomass_abundance_B()/area/unit_mass;
      bodymass[j]=s[i].bodymass()/unit_mass;
      j++;
    }
  }
  std::ofstream os(filename);
  os << log_spectrum(xy_graph(bodymass,biomass),10).log_xy();
}

link_strength_matrix NewWeb::trophic_distance_matrix(){
  int S=number_of_species();
  link_strength_matrix d;
  d.resize(S);

  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      d(i,j)=s[j].foraging_distance_to(s[i]);
    }
  }
  return d;
}

void NewWeb::species_table(const char* filename,
			   bool with_plants,double th,Snapshot & snap,
			   bool jolly){

  bool links_are_random=
    (get_cfg_parameter("totally_random_link_strength_sigma")!=0);

  sequence<std::string> additional_headings;
  sequence<std::string> additional_data;

  int column=0;

  int S=s.size();
  double area=area_per_compartment/unit_area;
  
  // compute niche widths and mean interaction strength
  link_strength_matrix tdm=trophic_distance_matrix();
  sequence<double> niche_width(S);
  sequence<double> meanC(S);
  for(int j=s.n_animals;j-->0;){
    for(int i=S;i-->0;){
      niche_width(j)+=tdm(i,j)*snap.ifrac()(i,j);
      meanC(j)+=s[j].aggressivity_g()*
	s[j].foraging_strength_c_on(s[i])*snap.ifrac()(i,j);
    }
  }
  meanC*=unit_mass/unit_area;

  // compute GPP
  double GPP=0;
  for(int i=s.size();i-->s.n_animals;){
    GPP+=s[i].the_GP;
  }
  cout << "GPP " << GPP/area*year/gram << endl; 
  cout << "(in units of g m^-2 yr^-1)" << endl;

  // compute remaining niche space main axis:
  simple_vector<bool> selection(niche_space_dimensions_D,true);
  selection[0]=false;
  simple_vector<double> VV1(niche_space_dimensions_D,0);
  simple_vector<double> FFa1(niche_space_dimensions_D,0);
  simple_vector<double> FFp1(niche_space_dimensions_D,0);
  sequence<double> pa_separation(niche_space_dimensions_D,0);
  if(!links_are_random){
    chi_square_meter VV(selection);
    chi_square_meter VVp(selection);
    chi_square_meter VVa(selection);
    chi_square_meter FFp(selection);
    chi_square_meter FFa(selection);
    for(int n=0;n<S;n++){
      VV.sample(s[n].vulnerability_V());
      if(is_plant(n)){
	VVp.sample(s[n].vulnerability_V());
	FFp.sample(s[n].niche_G());
      }else{
	VVa.sample(s[n].vulnerability_V());
	FFa.sample(s[n].foraging_F());
      }
    }
    VV1=VV.estimate().main_axis();
    FFa1=FFa.estimate().main_axis();
    FFp1=FFp.estimate().main_axis();
    pa_separation=sequence<double>(VVp.mean())-VVa.mean();
    pa_separation[0]=0;
    pa_separation/=abs(pa_separation);
    VV1.prepend(0);
    FFa1.prepend(0);
    FFp1.prepend(0);
  }


  // compute some switching exponents
  sequence<double> switching_b(S);
  if(!jolly){
    // By the way: this algorithm is not optimized. We need only one
    // column of the intake matrix in each iteration, not the full
    // thing.
    if(S>=2){    
      prepare_for_integration();
      ODE_vector state(number_of_variables()),
	time_derivative(number_of_variables());
      NewWeb::write_state_to(state);
      compute_flows(state,true); // initialize flow computation
      for(int j=s.n_animals;j-->0;){
	// first get two most important diet items:
	int first=0,second=1;
	if(snap.ifrac()(0,j)>snap.ifrac()(1,j)){
	  second=0;
	  first=1;
	}
	for(int i=S;i-->2;){
	  if(snap.ifrac()(i,j)>snap.ifrac()(first,j)){
	    // new first
	    second=first;
	    first=i;
	  }else if(snap.ifrac()(i,j)>snap.ifrac()(second,j)){
	    // new second
	    second=i;
	  }
	}
	// compute switching b between first and second
	double delta=0.01;
	NewWeb::write_state_to(state);
	state[first]-=delta/2;
	state[second]+=delta/2;
	compute_flows(state,false);
	double diet_ratio=fx(first,j)/fx(second,j);
	state[first]+=delta;
	state[second]-=delta;
	compute_flows(state,false);
	double diet_ratio1=fx(first,j)/fx(second,j);
	switching_b(j)=log(diet_ratio1/diet_ratio)/(2*delta);
      }
    }
  }

  simple_vector<weighted_average_meter> mean_link_strength(number_of_animals());
  // compute mean attack rate
  for(int i=number_of_animals();i-->0;){
    for(int j=number_of_species();j-->0;){
      mean_link_strength[i].
	sample(s[i].foraging_strength_c_on(s[j]),
	       exp(s[i].log_mass_matching(s[j])) 
	       );
    }
  }

  additional_headings[column++]=
    "logAg (log10 aggressivity)";
  additional_headings[column++]=
    "V0    (first component of vulnerability vector)";
  additional_headings[column++]=
    "F0    (first component of foraging vector)";
  additional_headings[column++]=
    "id    (sequential ID for each species)";
  additional_headings[column++]=
    "cid   (unique ID for each clade)";
  additional_headings[column++]=
    "pn    (passive(!) ingestion efficiency)";
  additional_headings[column++]=
    "f     (feeding level)";
  additional_headings[column++]=
    "cs    (competition strength)";
  additional_headings[column++]=
    "ls    (light strength)";
  additional_headings[column++]=
    (not jolly ? 
     "b     (switching exponent)" :  
     "b     (NOT COMPUTED! switching exponent)" );
  additional_headings[column++]=
    "nw    (niche width)";
  additional_headings[column++]=
    "col   (assigned column)";
  additional_headings[column++]=
    "V1    (second component of vulnerability)";
  additional_headings[column++]=
    "F1    (second component of foraging)";
  additional_headings[column++]=
    "Vm    (main component of remaining vulnerability)";
  additional_headings[column++]=
    "Fm    (corresponding component of foraging vector)";
  additional_headings[column++]=
    "Vr    (length of vulnerability vector)";
  additional_headings[column++]=
    "Fr    (length of foraging vector)";
  additional_headings[column++]=
    "GP    (yearly gross production density)";
  additional_headings[column++]=
    "tm    (trade-off multiplier)";
  additional_headings[column++]=
    "Vpa   (animal-plant separating component)";
  additional_headings[column++]=
    "Fpa   (corresponding foraging component)";
  additional_headings[column++]=
    "logC  (log10 diet-weighted interaction strength)";
  additional_headings[column++]=
    "wt    (trophic niche width)";
  additional_headings[column++]=
    "mc    (log10(mean link strength))";
  additional_headings[column++]=
    "F     (fishing mortality)";


  for(int n=0;n<S;n++){
    std::ostringstream os;
    //os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(7);
    NewSpecies & sp=s[n];
    if(with_plants || !sp.is_a_plant()){
      os.width(5);
      os << (sp.is_a_plant() ? 0 : 
	     log10(sp.aggressivity_g()*unit_mass/unit_area)) << " ";
      os.width(5);
      os << sp.vulnerability_V()[0] << " ";
      os.width(5);
      os << sp.foraging_F()[0] << " ";
      os << sp.sequential_id() << " ";
      os << sp.clade_id() << " ";
      os.width(9);
      os << sp.the_top_down_strength/(1+sp.the_top_down_strength) 
	 << " "; //= ingestion efficiency
      os.width(9);
      os << sp.the_saturation_strength << " "; //= feeding level
      os.width(9);
      os << sp.the_competition_strength << " ";
      os.width(9);
      os << sp.the_light_strength << " ";
      os.width(9);
      os << (not jolly ? switching_b(n) : 0) << " ";
      os.width(9);
      os << niche_width(n) << " ";
      os.width(9);
      os << assigned_column(n) << " ";
      os.width(5);
      os << (niche_space_dimensions_D<2?0:sp.vulnerability_V()[1]) << " ";
      os.width(5);
      os << (niche_space_dimensions_D<2?0:sp.foraging_F()[1]) << " ";
      os.width(5);
      os << dot(sp.vulnerability_V(),VV1) << " ";
      os.width(5);
      os << dot(sp.foraging_F(),sp.is_a_plant() ? FFp1 : VV1) << " ";
      os.width(5);
      os << sqrt(dot(sp.vulnerability_V(),sp.vulnerability_V())) << " ";
      os.width(5);
      os << sqrt(dot(sp.foraging_F(),sp.foraging_F())) << " ";
      os.width(9);
      os << sp.the_GP/area*year << " ";
      os.width(9);
      os << sp.get_trade_off_multiplier() << " ";
      os.width(5);
      os << dot(sp.vulnerability_V(),pa_separation) << " ";
      os.width(5);
      os << dot(sp.foraging_F(),pa_separation) << " ";
      os.width(5);
      os << ( sp.is_a_plant() ? 0 : log10(meanC(n)) ) << " ";
      os.width(5);
      os << ( sp.is_a_plant() ? 0 : sp.niche_width_wt() ) << " ";
      os.width(5);
      os << ( sp.is_a_plant() ? 0 : log10(mean_link_strength[n]) ) << " ";
      os.width(5);
      os << ( sp.fishing_mortality()*year ) << " ";
      additional_data[n]=os.str();
    }
  }
  snap.species_table(filename,with_plants,th,
		     additional_headings,additional_data);
}

void NewWeb::species_table2(const char* filename,
			   bool with_plants,double th,Snapshot & snap,
			   link_strength_matrix & intake, int colno,
			   bool jolly){
  double area=area_per_compartment/unit_area;

  if(intake.size()==0)
    intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  // Work out number of fish species
  double Sfish = 0.0;
  for(int n=0;n<S;n++){
    if(log10(s[n].bodymass()/unit_mass)>-3){
      Sfish += 1;
    }
  }
  //cout << "Sfish = " << Sfish << endl;

  sequence<std::string> additional_headings;
  sequence<std::string> additional_data;

  int column=0;

  //int S=s.size();
  //double area=area_per_compartment/unit_area;

  additional_headings[column++]=
    "col   (assigned column)";
  additional_headings[column++]=
    "prod   (productivity)";
  additional_headings[column++]=
    "pred   (predation mortality)";
  additional_headings[column++]=
    "comp   (mean dietary overlap for fish species)";
  additional_headings[column++]=
    "nin   (nin for fish species, considering all species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "nout   (nout for fish species, considering all species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninw   (nin for fish species, considering all species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutw   (nout for fish species, considering all species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninf   (nin for fish species, considering just fish species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutf   (nout for fish species, considering just fish species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninwf   (nin for fish species, considering just fish species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutwf   (nout for fish species, considering just fish species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "gen   (generality for fish species, using diet fraction threshold of default 0.01)";
  additional_headings[column++]=
    "vul   (generality for fish species, using diet fraction threshold of default 0.01)";
  additional_headings[column++]=
    "conn   (connectivity for fish species, using diet fraction threshold of default 0.01)";
  additional_headings[column++]=
    "Mxsim   (max trophic similarity for fish species, using diet fraction threshold of default 0.01)";
  additional_headings[column++]=
    "MPL   (min path length to deleted fish species for fish species, using diet fraction threshold of default 0.01)";
  
  // Calculate indicators for each fish species
  sequence<double> b_idotI;
  sequence<double> b_dotjI;
  sequence<double> b_idotIf;
  sequence<double> b_dotjIf;
  double b_dotdotI = 0.0;
  double b_dotdotIf = 0.0;
  sequence<double> H_ink;
  sequence<double> H_outk;
  sequence<double> n_ink;
  sequence<double> n_outk;
  sequence<double> H_inkf;
  sequence<double> H_outkf;
  sequence<double> n_inkf;
  sequence<double> n_outkf;
  double dummy3 = 0.0;
  double dummy4 = 0.0;
  double dummy5 = 0.0;
  double dummy6 = 0.0;
  double dummy7 = 0.0;
  double resolog10M = 0.0;
  double conslog10M = 0.0;
  sequence<double> Gen_I;
  sequence<double> Vul_I;
  sequence<double> Conn_I;
  double dietfractionthreshold = 0.01; 
  // These are matrices where ijth entry is number of prey shared by species i and j,
  // number of predators shared by species i and j,
  // number of prey and predator species shared by i and j divided by the total number
  // of prey and predators for both species (trophic similarity), all using diet threshold of 0.01.
  link_strength_matrix Npreyshared=intake;
  link_strength_matrix Npredshared=intake;
  link_strength_matrix Trophsim=intake;
  // Vector of max TrophSim for each species
  sequence<double> MaxTrophsim;
  // This is matrix indicating whether there are two links, one link or no links between two species,
  // using diet threshold of 0.01 (2, 1 or 0).
  // Diagonal entries are set to 0 as we are interested in path lengths between species, and setting
  // diagonal entries to non-zero means adding an extra link that should not be counted as a link in a path
  // between species.
  link_strength_matrix Linkornot=intake;
  // First populate b_idotI, b_dotjI, b_idotIf and b_dotjIf.
  // Also calcuate b_dotdotI and b_dotdotIf.
  for(int i=0;i<S;i++){
    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    dummy3 = 0;
    dummy4 = 0;
    dummy5 = 0;
    dummy6 = 0;
    for(int j=0;j<S;j++){
      NewSpecies & cons=s[j];
      conslog10M = log10(cons.bodymass()/unit_mass);
      dummy3 += intake(i,j)/unit_mass/area/(1/year);
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      if(resolog10M > -3){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M > -3)){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M <= -3)){
	dummy7 += intake(i,j)/unit_mass/area/(1/year);
      }
      // End loop for j.
    }
    b_idotI[i] = dummy3;
    b_dotjI[i] = dummy4;
    b_idotIf[i] = dummy5;
    b_dotjIf[i] = dummy6;
    // End loop for i.
  }
  b_dotdotIf = b_dotdotI - dummy7;
  // Now populate H_ink, H_outk, n_ink, n_outk.
  // Also populate H_inkf, H_outkf, n_inkf, n_outkf.
  // Also populate Gen_I, Vul_I and 
  // Conn_I (= Gen + Vul)
  // for each species.
  // Also calculate Conn_tot_I. To avoid double-counting, only count links for each species as a consumer.
  double dummy100 = 0.0;
  double dummy101 = 0.0;
  double dummy102 = 0.0;
  double dummy103 = 0.0;
  double dummy104 = 0.0;
  double dummy105 = 0.0;
  double sumGen = 0.0;
  double sumVul = 0.0;
  double sumConn = 0.0;
  double sumConntot = 0.0;
  double Conn_tot_I = 0.0;
  double sumpreyshared = 0.0;
  double sumpredshared = 0.0;
  double sumpreypred = 0.0;
  double dummy106 = 0.0;
  double dummy107 = 0.0;
  double dummy108 = 0.0;
  double dummy109 = 0.0;
  double sumGentot = 0.0;
  double sumGensquaredtot = 0.0;
  double sumVultot = 0.0;
  double sumVulsquaredtot = 0.0; 
  //double sumMaxTrophsim = 0.0;
  for(int i=0;i<S;i++){
    dummy102 = 0.0;
    dummy103 = 0.0;
    dummy104 = 0.0;
    dummy105 = 0.0;
    sumGen = 0;
    sumVul = 0;
    sumConn = 0;
    dummy109 = 0.0;
    for(int j=0;j<S;j++){
      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);
      if((dummy101>0)&&(b_dotjI[i]>0)){
	dummy102 = dummy102 - ((dummy101/b_dotjI[i])*
			       (log10(dummy101/b_dotjI[i])/log10(2)));
      }
      if((dummy100>0)&&(b_idotI[i]>0)){
	dummy103 = dummy103 - ((dummy100/b_idotI[i])*
			       (log10(dummy100/b_idotI[i])/log10(2)));
      }
      if((dummy101>0)&&(b_dotjIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)>-3)){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
      }
      if((dummy100>0)&&(b_idotIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)<=-3)){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
      }  
      // Work out Gen, Vul and Conn.
      if(b_dotjI[i]>0){
	if((dummy101/b_dotjI[i])>dietfractionthreshold){
	  sumGen = sumGen + 1;
	  sumConn = sumConn + 1;
	}
      }
      if(b_dotjI[j]>0){
	if((dummy100/b_dotjI[j])>dietfractionthreshold){
	  sumVul = sumVul + 1;
	  sumConn = sumConn + 1;
	}
      }
      // For species i and j, work out Npreyshared, Npredshared and Trophsim.
      // Only do this for pairs of fish species.
      if((log10(s[i].bodymass()/unit_mass)>-3)&&(log10(s[j].bodymass()/unit_mass)>-3)){
	sumpreyshared = 0.0;
	sumpredshared = 0.0;
	sumpreypred = 0.0;
	for(int k=0;k<S;k++){
	  dummy106 = intake(k,i)/unit_mass/area/(1/year);
	  dummy107 = intake(k,j)/unit_mass/area/(1/year);
	  dummy108 = 0.0;
	  if((b_dotjI[i]>0)&&((dummy106/b_dotjI[i])>dietfractionthreshold)&&(b_dotjI[j]>0)&&((dummy107/b_dotjI[j])>dietfractionthreshold)){
	    sumpreyshared = sumpreyshared + 1;
	  }
	  if((b_dotjI[i]>0)&&((dummy106/b_dotjI[i])>dietfractionthreshold)){
	    dummy108 = 1;
	  }
	  if((b_dotjI[j]>0)&&((dummy107/b_dotjI[j])>dietfractionthreshold)){
	    dummy108 = 1;
	  }
	  if(dummy108==1){
	    sumpreypred = sumpreypred + 1;
	  }
	  dummy106 = intake(i,k)/unit_mass/area/(1/year);
	  dummy107 = intake(j,k)/unit_mass/area/(1/year);
	  dummy108 = 0.0;
	  if((b_dotjI[k]>0)&&((dummy106/b_dotjI[k])>dietfractionthreshold)&&((dummy107/b_dotjI[k])>dietfractionthreshold)){
	    sumpredshared = sumpredshared + 1;
	  }
	  if((b_dotjI[k]>0)&&((dummy106/b_dotjI[k])>dietfractionthreshold)){
	    dummy108 = 1;
	  }
	  if((b_dotjI[k]>0)&&((dummy107/b_dotjI[k])>dietfractionthreshold)){
	    dummy108 = 1;
	  }
	  if(dummy108==1){
	    sumpreypred = sumpreypred + 1;
	  }
	}
	Npreyshared(i,j)=sumpreyshared;
	Npredshared(i,j)=sumpredshared;
	Trophsim(i,j)=(sumpreyshared+sumpredshared)/sumpreypred;
	// Update max value of TrophSim if necessary.
	// Remember constraint that i is not equal to j (Williams and Martinez, 2000).
	if((Trophsim(i,j)>dummy109)&&(i!=j)){
	  dummy109=Trophsim(i,j);
	}
      }
      // Update Linkornot matrix.
      if((b_dotjI[j]>0)&&((dummy100/b_dotjI[j])>dietfractionthreshold)&&(b_dotjI[i]>0)&&((dummy101/b_dotjI[i])>dietfractionthreshold)){
	Linkornot(i,j)=2;
      }else{
	if((b_dotjI[j]>0)&&((dummy100/b_dotjI[j])>dietfractionthreshold)){
	  Linkornot(i,j)=1;
	}else{
	  if((b_dotjI[i]>0)&&((dummy101/b_dotjI[i])>dietfractionthreshold)){
	    Linkornot(i,j)=1;
	  }else{
	    Linkornot(i,j)=0;
	  }
	}
      }
      // If diagonal entry, set Linkornot to 0.
      if(i==j){
	Linkornot(i,j)=0;
      }
      // End loop for j.
    }
    H_ink[i] = dummy102;
    H_outk[i] = dummy103;
    if(b_idotI[i]>0){
      n_outk[i] = pow(2,H_outk[i]);
    }else{
      n_outk[i] = 0;
    }
    if(b_dotjI[i]>0){
      n_ink[i] = pow(2,H_ink[i]);
    }else{
      n_ink[i] = 0;
    }
    H_inkf[i] = dummy104;
    H_outkf[i] = dummy105;
    if(b_idotIf[i]>0){
      n_outkf[i] = pow(2,H_outkf[i]);
    }else{
      n_outkf[i] = 0;
    }
    if(b_dotjIf[i]>0){
      n_inkf[i] = pow(2,H_inkf[i]);
    }else{
      n_inkf[i] = 0;
    }
    Gen_I[i] = sumGen;
    Vul_I[i] = sumVul;
    Conn_I[i] = sumConn;
    Conn_tot_I = Conn_tot_I + sumGen;
    //MaxTrophsim[i] = dummy109;
    sumGentot = sumGentot + sumGen;
    sumGensquaredtot = sumGensquaredtot + (sumGen*sumGen);
    sumVultot = sumVultot + sumVul;
    sumVulsquaredtot = sumVulsquaredtot + (sumVul*sumVul);
    //sumMaxTrophsim = sumMaxTrophsim + dummy109;
    // End loop for i.
  }
  // Work out Conn and topological LD and C for whole web.
  // Also work out topological mean standardised Gen, s.d. standardised Gen, 
  // mean standardised Vul and s.d. standardised Vul for whole web.
  // Print these values out.
  cout << "Conn_tot_I = " << Conn_tot_I << endl;
  cout << "LD_tot_I = " << Conn_tot_I/S << endl;
  cout << "C_tot_I = " << Conn_tot_I/(S*S) << endl;
  cout << "Mean_stand_Gen_I = " << sumGentot/Conn_tot_I << endl;
  cout << "sd_stand_Gen_I = " << ((sumGensquaredtot)*(S/(Conn_tot_I*Conn_tot_I)))-((sumGentot/Conn_tot_I)*(sumGentot/Conn_tot_I)) << endl;
  cout << "Mean_stand_Vul_I = " << sumVultot/Conn_tot_I << endl;
  cout << "sd_stand_Vul_I = " << ((sumVulsquaredtot)*(S/(Conn_tot_I*Conn_tot_I)))-((sumVultot/Conn_tot_I)*(sumVultot/Conn_tot_I)) << endl;
  //cout << "Mean_MaxTrophsim = " << sumMaxTrophsim/S << endl;
  // Work out entries for matrix for whether there is a path between two species of length 2,3,4,...,10, using diet threshold of 0.01.
  // Only calculate for entries in column corresponding to assigned_column for deleted species.
  link_strength_matrix Linkornot2=intake;
  link_strength_matrix Linkornot3=intake;
  link_strength_matrix Linkornot4=intake;
  link_strength_matrix Linkornot5=intake;
  link_strength_matrix Linkornot6=intake;
  link_strength_matrix Linkornot7=intake;
  link_strength_matrix Linkornot8=intake;
  link_strength_matrix Linkornot9=intake;
  link_strength_matrix Linkornot10=intake;
  // First find index corresponding to assigned_column of deleted species.
  int indexdelS = -1;
  //cout << "colno = " << colno << endl;
  for(int i=0;i<S;i++){
    if(assigned_column(i)==colno){
      indexdelS=i;
    }
  }
  //cout << "indexdelS = " << indexdelS << endl;

  double sumMatentry=0.0;
  for(int i=0;i<S;i++){
    //cout << i << endl;
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot(k,indexdelS));
    }
    Linkornot2(i,indexdelS)=sumMatentry;
    //cout << "Linkornot2(i,indexdelS) = " << Linkornot2(i,indexdelS) << endl;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot2(k,indexdelS));
    }
    Linkornot3(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot3(k,indexdelS));
    }
    Linkornot4(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot4(k,indexdelS));
    }
    Linkornot5(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot5(k,indexdelS));
    }
    Linkornot6(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot6(k,indexdelS));
    }
    Linkornot7(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot7(k,indexdelS));
    }
    Linkornot8(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot8(k,indexdelS));
    }
    Linkornot9(i,indexdelS)=sumMatentry;
  }
  for(int i=0;i<S;i++){
    sumMatentry=0.0;
    for(int k=0;k<S;k++){
      sumMatentry=sumMatentry+(Linkornot(i,k)*Linkornot9(k,indexdelS));
    }
    Linkornot10(i,indexdelS)=sumMatentry;
  }
  // Work out minimum path length between two species, using diet threshold of 0.01. 11 represents >=11.
  // Only calculate for entries in column corresponding to assigned_column for deleted species.
  // For diagonal entries, Minpathlength set to 0.
  link_strength_matrix Minpathlength=intake;
  int Minpathlenghentry=0;
  for(int i=0;i<S;i++){
    Minpathlenghentry=0;
    if(Linkornot(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=1;
      }
    }
    if(Linkornot2(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=2;
      }
    }
    if(Linkornot3(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=3;
      }
    }
    if(Linkornot4(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=4;
      }
    }
    if(Linkornot5(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=5;
      }
    }
    if(Linkornot6(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=6;
      }
    }
    if(Linkornot7(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=7;
      }
    }
    if(Linkornot8(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=8;
      }
    }
    if(Linkornot9(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=9;
      }
    }
    if(Linkornot10(i,indexdelS)>0){
      if(Minpathlenghentry==0){
	Minpathlenghentry=10;
      }
    }   
    if((Minpathlenghentry==0)&&(i!=indexdelS)){
      Minpathlenghentry=11;
    }
    if(i==indexdelS){
      Minpathlenghentry=0;
    }
    Minpathlength(i,indexdelS)=Minpathlenghentry;
  }
  //cout << "Minpathlength(indexdelS,indexdelS) = " << Minpathlength(indexdelS,indexdelS) << endl;

  // Work out TLs
  NewVector level;
  NewMatrix dietF;
  dietF=ifrac;
  NewVector ones(S);
  for(int i=0;i<S;i++){
    ones[i]=1;
  }
  level=ones;
  NewVector level2;
  NewVector leveldiff;
  double dummycounter = 0;
  double dummydiff = 0;
  do{
    level2=level;
    level=prod(dietF,level)+ones;
    leveldiff=level-level2;
    dummydiff=NewVectorAbs(leveldiff);
    dummycounter=dummycounter+1;
    //cout << "dummycounter = " << dummycounter << endl;
    //cout << "dummydiff = " << dummydiff << endl;
  }while(dummydiff>=0.0001);
  cout << "dummycounter = " << dummycounter << endl;
  cout << "dummydiff = " << dummydiff << endl;

  // Calculate indicators involving TL.
  double dummy1 = 0.0;
  int dummy2 = 0;
  sequence<double> rho_TL;
  for(int i=0;i<45;i++){
    rho_TL[i]=0;
  }
  for(int i=0;i<S;i++){
    dummy1 = (level[i] - 1)/0.1;
    //cout << "dummy1 = " << dummy1 << endl;
    dummy2 = floor(dummy1);
    //cout << "dummy2 = " << dummy2 << endl;
    rho_TL[dummy2] += biomass_B(i)/unit_mass/area;
  }
  //for(int i=0;i<45;i++){
  //cout << "rho_TL[" << i << "] = " << rho_TL[i] << endl;
  //}
  cout << "rho_1_1point1_TL = " << rho_TL[0] << endl;
  cout << "rho_1point1_1point2_TL = " << rho_TL[1] << endl;
  cout << "rho_1point2_1point3_TL = " << rho_TL[2] << endl;
  cout << "rho_1point3_1point4_TL = " << rho_TL[3] << endl;
  cout << "rho_1point4_1point5_TL = " << rho_TL[4] << endl;
  cout << "rho_1point5_1point6_TL = " << rho_TL[5] << endl;
  cout << "rho_1point6_1point7_TL = " << rho_TL[6] << endl;
  cout << "rho_1point7_1point8_TL = " << rho_TL[7] << endl;
  cout << "rho_1point8_1point9_TL = " << rho_TL[8] << endl;
  cout << "rho_1point9_2_TL = " << rho_TL[9] << endl;
  cout << "rho_2_2point1_TL = " << rho_TL[10] << endl;
  cout << "rho_2point1_2point2_TL = " << rho_TL[11] << endl;
  cout << "rho_2point2_2point3_TL = " << rho_TL[12] << endl;
  cout << "rho_2point3_2point4_TL = " << rho_TL[13] << endl;
  cout << "rho_2point4_2point5_TL = " << rho_TL[14] << endl;
  cout << "rho_2point5_2point6_TL = " << rho_TL[15] << endl;
  cout << "rho_2point6_2point7_TL = " << rho_TL[16] << endl;
  cout << "rho_2point7_2point8_TL = " << rho_TL[17] << endl;
  cout << "rho_2point8_2point9_TL = " << rho_TL[18] << endl;
  cout << "rho_2point9_3_TL = " << rho_TL[19] << endl;
  cout << "rho_3_3point1_TL = " << rho_TL[20] << endl;
  cout << "rho_3point1_3point2_TL = " << rho_TL[21] << endl;
  cout << "rho_3point2_3point3_TL = " << rho_TL[22] << endl;
  cout << "rho_3point3_3point4_TL = " << rho_TL[23] << endl;
  cout << "rho_3point4_3point5_TL = " << rho_TL[24] << endl;
  cout << "rho_3point5_3point6_TL = " << rho_TL[25] << endl;
  cout << "rho_3point6_3point7_TL = " << rho_TL[26] << endl;
  cout << "rho_3point7_3point8_TL = " << rho_TL[27] << endl;
  cout << "rho_3point8_3point9_TL = " << rho_TL[28] << endl;
  cout << "rho_3point9_4_TL = " << rho_TL[29] << endl;
  cout << "rho_4_4point1_TL = " << rho_TL[30] << endl;
  cout << "rho_4point1_4point2_TL = " << rho_TL[31] << endl;
  cout << "rho_4point2_4point3_TL = " << rho_TL[32] << endl;
  cout << "rho_4point3_4point4_TL = " << rho_TL[33] << endl;
  cout << "rho_4point4_4point5_TL = " << rho_TL[34] << endl;
  cout << "rho_4point5_4point6_TL = " << rho_TL[35] << endl;
  cout << "rho_4point6_4point7_TL = " << rho_TL[36] << endl;
  cout << "rho_4point7_4point8_TL = " << rho_TL[37] << endl;
  cout << "rho_4point8_4point9_TL = " << rho_TL[38] << endl;
  cout << "rho_4point9_5_TL = " << rho_TL[39] << endl;
  cout << "rho_5_5point1_TL = " << rho_TL[40] << endl;
  cout << "rho_5point1_5point2_TL = " << rho_TL[41] << endl;
  cout << "rho_5point2_5point3_TL = " << rho_TL[42] << endl;
  cout << "rho_5point3_5point4_TL = " << rho_TL[43] << endl;
  cout << "rho_5point4_5point5_TL = " << rho_TL[44] << endl;

  // Flows from organisms of lower or equal TL to those of higher TL
  // Flows from organisms of higher TL to those of lower or equal TL
  // Then work out b_xTL.
  double S_c = number_of_animals();
  double S_f = 0.0;
  sequence<double> b_lowtohigh_TL;
  sequence<double> b_lowtohigh_TL_Mcorrected;
  sequence<double> b_hightolow_TL;
  sequence<double> b_hightolow_TL_Mcorrected;
  for(int i=0;i<36;i++){
    b_lowtohigh_TL[i]=0;
    b_lowtohigh_TL_Mcorrected[i]=0;
    b_hightolow_TL[i]=0;
    b_hightolow_TL_Mcorrected[i]=0;
  }
  resolog10M = 0.0;
  conslog10M = 0.0;
  //sequence<double> b_dotjI;
  //sequence<double> b_dotjIf;
  b_dotdotI = 0.0;
  b_dotdotIf = 0.0;
  double resoTL = 0.0;
  double consTL = 0.0;
  dummy4 = 0.0;
  dummy6 = 0.0;
  dummy7 = 0.0;
  
  for(int i=0;i<S;i++){
    
    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    resoTL = level[i];
    dummy4 = 0;
    dummy6 = 0;     
    
    for(int j=0;j<S;j++){
      
      NewSpecies & cons=s[j];	
      conslog10M = log10(cons.bodymass()/unit_mass);
      consTL = level[j];
      // Update sums for calculation of flows for species
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      if(resolog10M > -3){
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M > -3)){
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M <= -3)){
	dummy7 += intake(i,j)/unit_mass/area/(1/year);
      }
      
      if((resoTL <= 1.9)&&(consTL > 1.9)){
	b_lowtohigh_TL[0] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2)&&(consTL > 2)){
	b_lowtohigh_TL[1] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.1)&&(consTL > 2.1)){
	b_lowtohigh_TL[2] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.2)&&(consTL > 2.2)){
	b_lowtohigh_TL[3] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.3)&&(consTL > 2.3)){
	b_lowtohigh_TL[4] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.4)&&(consTL > 2.4)){
	b_lowtohigh_TL[5] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.5)&&(consTL > 2.5)){
	b_lowtohigh_TL[6] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.6)&&(consTL > 2.6)){
	b_lowtohigh_TL[7] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.7)&&(consTL > 2.7)){
	b_lowtohigh_TL[8] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.8)&&(consTL > 2.8)){
	b_lowtohigh_TL[9] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 2.9)&&(consTL > 2.9)){
	b_lowtohigh_TL[10] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3)&&(consTL > 3)){
	b_lowtohigh_TL[11] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.1)&&(consTL > 3.1)){
	b_lowtohigh_TL[12] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.2)&&(consTL > 3.2)){
	b_lowtohigh_TL[13] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.3)&&(consTL > 3.3)){
	b_lowtohigh_TL[14] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.4)&&(consTL > 3.4)){
	b_lowtohigh_TL[15] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.5)&&(consTL > 3.5)){
	b_lowtohigh_TL[16] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.6)&&(consTL > 3.6)){
	b_lowtohigh_TL[17] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.7)&&(consTL > 3.7)){
	b_lowtohigh_TL[18] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.8)&&(consTL > 3.8)){
	b_lowtohigh_TL[19] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 3.9)&&(consTL > 3.9)){
	b_lowtohigh_TL[20] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4)&&(consTL > 4)){
	b_lowtohigh_TL[21] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.1)&&(consTL > 4.1)){
	b_lowtohigh_TL[22] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.2)&&(consTL > 4.2)){
	b_lowtohigh_TL[23] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.3)&&(consTL > 4.3)){
	b_lowtohigh_TL[24] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.4)&&(consTL > 4.4)){
	b_lowtohigh_TL[25] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.5)&&(consTL > 4.5)){
	b_lowtohigh_TL[26] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.6)&&(consTL > 4.6)){
	b_lowtohigh_TL[27] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.7)&&(consTL > 4.7)){
	b_lowtohigh_TL[28] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.8)&&(consTL > 4.8)){
	b_lowtohigh_TL[29] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 4.9)&&(consTL > 4.9)){
	b_lowtohigh_TL[30] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 5)&&(consTL > 5)){
	b_lowtohigh_TL[31] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 5.1)&&(consTL > 5.1)){
	b_lowtohigh_TL[32] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 5.2)&&(consTL > 5.2)){
	b_lowtohigh_TL[33] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 5.3)&&(consTL > 5.3)){
	b_lowtohigh_TL[34] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL <= 5.4)&&(consTL > 5.4)){
	b_lowtohigh_TL[35] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_TL_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      
      if((resoTL > 1.9)&&(consTL <= 1.9)){
	b_hightolow_TL[0] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2)&&(consTL <= 2)){
	b_hightolow_TL[1] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.1)&&(consTL <= 2.1)){
	b_hightolow_TL[2] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.2)&&(consTL <= 2.2)){
	b_hightolow_TL[3] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.3)&&(consTL <= 2.3)){
	b_hightolow_TL[4] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.4)&&(consTL <= 2.4)){
	b_hightolow_TL[5] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.5)&&(consTL <= 2.5)){
	b_hightolow_TL[6] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.6)&&(consTL <= 2.6)){
	b_hightolow_TL[7] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.7)&&(consTL <= 2.7)){
	b_hightolow_TL[8] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.8)&&(consTL <= 2.8)){
	b_hightolow_TL[9] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 2.9)&&(consTL <= 2.9)){
	b_hightolow_TL[10] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3)&&(consTL <= 3)){
	b_hightolow_TL[11] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.1)&&(consTL <= 3.1)){
	b_hightolow_TL[12] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.2)&&(consTL <= 3.2)){
	b_hightolow_TL[13] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.3)&&(consTL <= 3.3)){
	b_hightolow_TL[14] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.4)&&(consTL <= 3.4)){
	b_hightolow_TL[15] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.5)&&(consTL <= 3.5)){
	b_hightolow_TL[16] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.6)&&(consTL <= 3.6)){
	b_hightolow_TL[17] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.7)&&(consTL <= 3.7)){
	b_hightolow_TL[18] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.8)&&(consTL <= 3.8)){
	b_hightolow_TL[19] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 3.9)&&(consTL <= 3.9)){
	b_hightolow_TL[20] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4)&&(consTL <= 4)){
	b_hightolow_TL[21] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.1)&&(consTL <= 4.1)){
	b_hightolow_TL[22] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.2)&&(consTL <= 4.2)){
	b_hightolow_TL[23] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.3)&&(consTL <= 4.3)){
	b_hightolow_TL[24] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.4)&&(consTL <= 4.4)){
	b_hightolow_TL[25] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.5)&&(consTL <= 4.5)){
	b_hightolow_TL[26] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.6)&&(consTL <= 4.6)){
	b_hightolow_TL[27] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.7)&&(consTL <= 4.7)){
	b_hightolow_TL[28] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.8)&&(consTL <= 4.8)){
	b_hightolow_TL[29] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 4.9)&&(consTL <= 4.9)){
	b_hightolow_TL[30] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 5)&&(consTL <= 5)){
	b_hightolow_TL[31] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 5.1)&&(consTL <= 5.1)){
	b_hightolow_TL[32] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 5.2)&&(consTL <= 5.2)){
	b_hightolow_TL[33] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 5.3)&&(consTL <= 5.3)){
	b_hightolow_TL[34] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resoTL > 5.4)&&(consTL <= 5.4)){
	b_hightolow_TL[35] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_TL_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      
      // close loop for j
    }
    
    // Update vectors of flows for species
    b_dotjI[i] = dummy4;
    b_dotjIf[i] = dummy6;
    
    // Update species richness
    if(resolog10M > -3){
      S_f += 1;
    }
      
    // close loop for i 
  }

  b_dotdotIf = b_dotdotI - dummy7;

  cout << "b_1point9_TL = " << b_lowtohigh_TL[0]-b_hightolow_TL[0] << endl;
  cout << "b_2_TL = " << b_lowtohigh_TL[1]-b_hightolow_TL[1] << endl;
  cout << "b_2point1_TL = " << b_lowtohigh_TL[2]-b_hightolow_TL[2] << endl;
  cout << "b_2point2_TL = " << b_lowtohigh_TL[3]-b_hightolow_TL[3] << endl;
  cout << "b_2point3_TL = " << b_lowtohigh_TL[4]-b_hightolow_TL[4] << endl;
  cout << "b_2point4_TL = " << b_lowtohigh_TL[5]-b_hightolow_TL[5] << endl;
  cout << "b_2point5_TL = " << b_lowtohigh_TL[6]-b_hightolow_TL[6] << endl;
  cout << "b_2point6_TL = " << b_lowtohigh_TL[7]-b_hightolow_TL[7] << endl;
  cout << "b_2point7_TL = " << b_lowtohigh_TL[8]-b_hightolow_TL[8] << endl;
  cout << "b_2point8_TL = " << b_lowtohigh_TL[9]-b_hightolow_TL[9] << endl;
  cout << "b_2point9_TL = " << b_lowtohigh_TL[10]-b_hightolow_TL[10] << endl;
  cout << "b_3_TL = " << b_lowtohigh_TL[11]-b_hightolow_TL[11] << endl;
  cout << "b_3point1_TL = " << b_lowtohigh_TL[12]-b_hightolow_TL[12] << endl;
  cout << "b_3point2_TL = " << b_lowtohigh_TL[13]-b_hightolow_TL[13] << endl;
  cout << "b_3point3_TL = " << b_lowtohigh_TL[14]-b_hightolow_TL[14] << endl;
  cout << "b_3point4_TL = " << b_lowtohigh_TL[15]-b_hightolow_TL[15] << endl;
  cout << "b_3point5_TL = " << b_lowtohigh_TL[16]-b_hightolow_TL[16] << endl;
  cout << "b_3point6_TL = " << b_lowtohigh_TL[17]-b_hightolow_TL[17] << endl;
  cout << "b_3point7_TL = " << b_lowtohigh_TL[18]-b_hightolow_TL[18] << endl;
  cout << "b_3point8_TL = " << b_lowtohigh_TL[19]-b_hightolow_TL[19] << endl;
  cout << "b_3point9_TL = " << b_lowtohigh_TL[20]-b_hightolow_TL[20] << endl;
  cout << "b_4_TL = " << b_lowtohigh_TL[21]-b_hightolow_TL[21] << endl;
  cout << "b_4point1_TL = " << b_lowtohigh_TL[22]-b_hightolow_TL[22] << endl;
  cout << "b_4point2_TL = " << b_lowtohigh_TL[23]-b_hightolow_TL[23] << endl;
  cout << "b_4point3_TL = " << b_lowtohigh_TL[24]-b_hightolow_TL[24] << endl;
  cout << "b_4point4_TL = " << b_lowtohigh_TL[25]-b_hightolow_TL[25] << endl;
  cout << "b_4point5_TL = " << b_lowtohigh_TL[26]-b_hightolow_TL[26] << endl;
  cout << "b_4point6_TL = " << b_lowtohigh_TL[27]-b_hightolow_TL[27] << endl;
  cout << "b_4point7_TL = " << b_lowtohigh_TL[28]-b_hightolow_TL[28] << endl;
  cout << "b_4point8_TL = " << b_lowtohigh_TL[29]-b_hightolow_TL[29] << endl;
  cout << "b_4point9_TL = " << b_lowtohigh_TL[30]-b_hightolow_TL[30] << endl;
  cout << "b_5_TL = " << b_lowtohigh_TL[31]-b_hightolow_TL[31] << endl;
  cout << "b_5point1_TL = " << b_lowtohigh_TL[32]-b_hightolow_TL[32] << endl;
  cout << "b_5point2_TL = " << b_lowtohigh_TL[33]-b_hightolow_TL[33] << endl;
  cout << "b_5point3_TL = " << b_lowtohigh_TL[34]-b_hightolow_TL[34] << endl;
  cout << "b_5point4_TL = " << b_lowtohigh_TL[35]-b_hightolow_TL[35] << endl;
  cout << "b_1point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[0]-b_hightolow_TL_Mcorrected[0] << endl;
  cout << "b_2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[1]-b_hightolow_TL_Mcorrected[1] << endl;
  cout << "b_2point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[2]-b_hightolow_TL_Mcorrected[2] << endl;
  cout << "b_2point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[3]-b_hightolow_TL_Mcorrected[3] << endl;
  cout << "b_2point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[4]-b_hightolow_TL_Mcorrected[4] << endl;
  cout << "b_2point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[5]-b_hightolow_TL_Mcorrected[5] << endl;
  cout << "b_2point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[6]-b_hightolow_TL_Mcorrected[6] << endl;
  cout << "b_2point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[7]-b_hightolow_TL_Mcorrected[7] << endl;
  cout << "b_2point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[8]-b_hightolow_TL_Mcorrected[8] << endl;
  cout << "b_2point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[9]-b_hightolow_TL_Mcorrected[9] << endl;
  cout << "b_2point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[10]-b_hightolow_TL_Mcorrected[10] << endl;
  cout << "b_3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[11]-b_hightolow_TL_Mcorrected[11] << endl;
  cout << "b_3point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[12]-b_hightolow_TL_Mcorrected[12] << endl;
  cout << "b_3point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[13]-b_hightolow_TL_Mcorrected[13] << endl;
  cout << "b_3point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[14]-b_hightolow_TL_Mcorrected[14] << endl;
  cout << "b_3point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[15]-b_hightolow_TL_Mcorrected[15] << endl;
  cout << "b_3point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[16]-b_hightolow_TL_Mcorrected[16] << endl;
  cout << "b_3point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[17]-b_hightolow_TL_Mcorrected[17] << endl;
  cout << "b_3point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[18]-b_hightolow_TL_Mcorrected[18] << endl;
  cout << "b_3point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[19]-b_hightolow_TL_Mcorrected[19] << endl;
  cout << "b_3point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[20]-b_hightolow_TL_Mcorrected[20] << endl;
  cout << "b_4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[21]-b_hightolow_TL_Mcorrected[21] << endl;
  cout << "b_4point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[22]-b_hightolow_TL_Mcorrected[22] << endl;
  cout << "b_4point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[23]-b_hightolow_TL_Mcorrected[23] << endl;
  cout << "b_4point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[24]-b_hightolow_TL_Mcorrected[24] << endl;
  cout << "b_4point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[25]-b_hightolow_TL_Mcorrected[25] << endl;
  cout << "b_4point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[26]-b_hightolow_TL_Mcorrected[26] << endl;
  cout << "b_4point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[27]-b_hightolow_TL_Mcorrected[27] << endl;
  cout << "b_4point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[28]-b_hightolow_TL_Mcorrected[28] << endl;
  cout << "b_4point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[29]-b_hightolow_TL_Mcorrected[29] << endl;
  cout << "b_4point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[30]-b_hightolow_TL_Mcorrected[30] << endl;
  cout << "b_5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[31]-b_hightolow_TL_Mcorrected[31] << endl;
  cout << "b_5point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[32]-b_hightolow_TL_Mcorrected[32] << endl;
  cout << "b_5point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[33]-b_hightolow_TL_Mcorrected[33] << endl;
  cout << "b_5point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[34]-b_hightolow_TL_Mcorrected[34] << endl;
  cout << "b_5point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[35]-b_hightolow_TL_Mcorrected[35] << endl;
  
  sequence<double> omn_j;
  dummy100 = 0.0;
  dummy101 = 0.0;
  double dummy110 = 0.0;
  double dummy111 = 0.0;
  double dummy112 = 0.0;
  double dummy113 = 0.0;
  double dummy114 = 0.0;
  double dummy115 = 0.0;

  for(int i=0;i<S;i++){
    
    // Reset sums for vectors
    dummy110 = 0;
    
    for(int j=0;j<S;j++){
      
      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);
      
      if(b_dotjI[i]>0){
	dummy110 += (dummy101/b_dotjI[i])*pow(level[j]-level[i],2);
      }
      
      if((level[i]>=1.5)&&(level[i]<2.5)&&(level[j]>=2.5)&&(level[j]<3.5)){
	dummy113 += dummy100;
      }
      if((level[i]>=2.5)&&(level[i]<3.5)&&(level[j]>=3.5)&&(level[j]<4.5)){
	dummy114 += dummy100;
      }
      if((level[i]>=3.5)&&(level[i]<4.5)&&(level[j]>=4.5)&&(level[j]<5.5)){
	dummy115 += dummy100;
      }
      
      // end loop for j
    }
    
    omn_j[i] = sqrt(dummy110);
    
    if(level[i]==1){ dummy112 += (s[i].the_GP)/unit_mass/area/(1/year); } 
    
    // end loop for i
  }

  double omn_prime = 0.0;
  double omn = 0.0;
  double omn_f_prime = 0.0;
  double omn_f = 0.0;
  double epsilon_3 = 0.0;
  double epsilon_4 = 0.0;
  double epsilon_5 = 0.0;
  double dummy204 = 0.0;
  double dummy205 = 0.0;
  double dummy206 = 0.0;
  double dummy207 = 0.0;
  double dummy208 = 0.0;    
  
  for(int i=0;i<S;i++){
    if(log10(s[i].bodymass()/unit_mass)>-3){
      dummy208 += b_dotjI[i];
    }
  }
  for(int i=0;i<S;i++){
    dummy204 += omn_j[i];
    dummy205 += (b_dotjI[i]/b_dotdotI)*omn_j[i];
    if(log10(s[i].bodymass()/unit_mass)>-3){
      dummy206 += omn_j[i];
      dummy207 += (b_dotjI[i]/dummy208)*omn_j[i];
    }
  }
  // It has been checked that dividing by S gives decimal
  omn_prime = dummy204/S_c;
  omn = dummy205;
  omn_f_prime = dummy206/S_f;
  omn_f = dummy207;
  epsilon_3 = dummy113/dummy112;
  epsilon_4 = dummy114/dummy113;
  epsilon_5 = dummy115/dummy114;
  
  cout << "omn_prime = " << omn_prime << endl;
  cout << "omn = " << omn << endl;
  cout << "omn_f_prime = " << omn_f_prime << endl;
  cout << "omn_f = " << omn_f << endl;

  
  for(int n=0;n<S;n++){
    std::ostringstream os;
    double prod = 0.0;
    double pred = 0.0;
    double comp = 0.0;
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;
    int counter4 = 0;
    //os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(7);
    NewSpecies & sp=s[n];     
    if(log10(sp.bodymass()/unit_mass)>-3){
      if(with_plants || !sp.is_a_plant()){
	os.width(1);
	os << assigned_column(n) << " ";
	for(int i=0;i<S;i++){
	  NewSpecies & sp2=s[i];
	  prod += 0.6*(intake(i,n)/unit_mass/area/(1/year));
	  pred += intake(n,i)/unit_mass/area/(1/year);
	  if(log10(sp2.bodymass()/unit_mass)>-3){
	    double comp1 = 0.0;
	    double comp2 = 0.0;
	    double comp3 = 0.0;
	    for(int j=0;j<S;j++){
	      comp1 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	      comp2 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,n)/unit_mass/area/(1/year));
	      comp3 += (intake(j,i)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	      // End loop for j.
	    }
	    comp += comp1/(sqrt(comp2*comp3));
	  }
	  // End loop for i.
   	}
	// Work out mean of competition coefficients
	comp = comp/Sfish;
	os.width(5);
	os << prod << " ";
	os.width(5);
	os << pred << " ";
	os.width(5);
	os << comp << " ";
	os.width(5);
	os << n_ink[n] << " ";
	os.width(5);
	os << n_outk[n] << " ";
	os.width(5);
	os << (b_dotjI[n]/b_dotdotI)*n_ink[n] << " ";
	os.width(5);
	os << (b_idotI[n]/b_dotdotI)*n_outk[n] << " ";
	os.width(5);
	os << n_inkf[n] << " ";
	os.width(5);
	os << n_outkf[n] << " ";
	os.width(5);
	os << (b_dotjIf[n]/b_dotdotIf)*n_inkf[n] << " ";
	os.width(5);
	os << (b_idotIf[n]/b_dotdotIf)*n_outkf[n] << " ";
	os.width(5);
	os << Gen_I[n] << " ";
	os.width(5);
	os << Vul_I[n] << " ";
	os.width(5);
	os << Conn_I[n] << " ";
	os.width(5);
	os << MaxTrophsim[n] << " ";
	os.width(5);
	os << Minpathlength(n,indexdelS) << " ";
	additional_data[n]=os.str();
      }
    }
  }
  snap.species_table2(filename,with_plants,th,
		     additional_headings,additional_data);

  // Work out matrix of consumption terms for fish species.
  // Each row gives, for each fish species, consumption on
  // other fish species in order of increasing body mass.
  std::ofstream os("consum_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_consum;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  bodymass_consum.push_back(make_pair(sp2.bodymass(),intake(i,n)/unit_mass/area/(1/year)));
	}
      }
      sort(bodymass_consum.begin(),bodymass_consum.end());
      os << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_consum.size();i++){
	//cout << bodymass_consum.at(i).first << endl;
	os << " " << bodymass_consum.at(i).second;
      }
      os << endl;
    }
  }

  // Work out matrix of diet overlaps for fish species.
  // Each row gives, for each fish species, diet overlap
  // with other fish species in order of increasing body mass.
  std::ofstream os2("dietoverlap_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os2.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_dietoverlap;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  double comp1 = 0.0;
	  double comp2 = 0.0;
	  double comp3 = 0.0;
	  double comp = 0.0;
	  for(int j=0;j<S;j++){
	    comp1 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	    comp2 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,n)/unit_mass/area/(1/year));
	    comp3 += (intake(j,i)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	  }
	  comp = comp1/(sqrt(comp2*comp3));
	  bodymass_dietoverlap.push_back(make_pair(sp2.bodymass(),comp));
	}
      }
      sort(bodymass_dietoverlap.begin(),bodymass_dietoverlap.end());
      os2 << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_dietoverlap.size();i++){
	//cout << bodymass_dietoverlap.at(i).first << endl;
	os2 << " " << bodymass_dietoverlap.at(i).second;
      }
      os2 << endl;
    }
  }

  // Work out matrix of c coefficient overlaps for fish species.
  // Each row gives, for each fish species, c coefficient overlap
  // with other fish species in order of increasing body mass.
  std::ofstream os3("ccoeffoverlap_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os3.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_ccoeffoverlap;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  double ccoeff1 = 0.0;
	  double ccoeff2 = 0.0;
	  double ccoeff3 = 0.0;
	  double ccoeff = 0.0;
	  for(int j=0;j<S;j++){
	    ccoeff1 += (s[n].foraging_strength_c_on(s[j]))*(s[i].foraging_strength_c_on(s[j]));
	    ccoeff2 += (s[n].foraging_strength_c_on(s[j]))*(s[n].foraging_strength_c_on(s[j]));
	    ccoeff3 += (s[i].foraging_strength_c_on(s[j]))*(s[i].foraging_strength_c_on(s[j]));
	  }
	  ccoeff = ccoeff1/(sqrt(ccoeff2*ccoeff3));
	  bodymass_ccoeffoverlap.push_back(make_pair(sp2.bodymass(),ccoeff));
	}
      }
      sort(bodymass_ccoeffoverlap.begin(),bodymass_ccoeffoverlap.end());
      os3 << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_ccoeffoverlap.size();i++){
	//cout << bodymass_ccoeffoverlap.at(i).first << endl;
	os3 << " " << bodymass_ccoeffoverlap.at(i).second;
      }
      os3 << endl;
    }
  }

  // Work out matrix of Npreyshared for fish species.
  // Each row gives, for each fish species, Npreyshared
  // with other fish species in order of increasing body mass.
  std::ofstream os4("Npreyshared_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os4.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_Npreyshared;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  bodymass_Npreyshared.push_back(make_pair(sp2.bodymass(),Npreyshared(n,i)));
	}
      }
      sort(bodymass_Npreyshared.begin(),bodymass_Npreyshared.end());
      os4 << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_Npreyshared.size();i++){
	//cout << bodymass_Npreyshared.at(i).first << endl;
	os4 << " " << bodymass_Npreyshared.at(i).second;
      }
      os4 << endl;
    }
  }

  // Work out matrix of Npredshared for fish species.
  // Each row gives, for each fish species, Npredshared
  // with other fish species in order of increasing body mass.
  std::ofstream os5("Npredshared_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os5.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_Npredshared;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  bodymass_Npredshared.push_back(make_pair(sp2.bodymass(),Npredshared(n,i)));
	}
      }
      sort(bodymass_Npredshared.begin(),bodymass_Npredshared.end());
      os5 << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_Npredshared.size();i++){
	//cout << bodymass_Npredshared.at(i).first << endl;
	os5 << " " << bodymass_Npredshared.at(i).second;
      }
      os5 << endl;
    }
  }

  // Work out matrix of Trophsim for fish species.
  // Each row gives, for each fish species, Trophsim
  // with other fish species in order of increasing body mass.
  std::ofstream os6("Trophsim_matrix.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os6.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    if(log10(sp.bodymass()/unit_mass)>-3){
      sequence< pair<double,double> > bodymass_Trophsim;
      for(int i=0;i<S;i++){
	NewSpecies & sp2=s[i];
	if(log10(sp2.bodymass()/unit_mass)>-3){
	  bodymass_Trophsim.push_back(make_pair(sp2.bodymass(),Trophsim(n,i)));
	}
      }
      sort(bodymass_Trophsim.begin(),bodymass_Trophsim.end());
      os6 << log10(sp.bodymass()/unit_mass);
      for(int i=0;i<bodymass_Trophsim.size();i++){
	//cout << bodymass_Trophsim.at(i).first << endl;
	os6 << " " << bodymass_Trophsim.at(i).second;
      }
      os6 << endl;
    }
  }

  // Work out matrix oflog10(bodymass), log10(biomass density) and TL  
  // for each species.
  std::ofstream os7("log10M_log10B_TL.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os7.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    os7 << log10(sp.bodymass()/unit_mass);
    os7 << " " << log10(biomass_B(n)/unit_mass/area);
    os7 << " " << level[n];
    os7 << endl;
  }

  // Work out matrix oflog10(bodymass), log10(biomass density), TL, 
  // Gen, Vul and Conn for each species.
  std::ofstream os8("log10M_log10B_TL_Gen_Vul_Conn.dat");
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os8.precision(7);
  for(int n=0;n<S;n++){
    NewSpecies & sp=s[n];
    os8 << log10(sp.bodymass()/unit_mass);
    os8 << " " << log10(biomass_B(n)/unit_mass/area);
    os8 << " " << level[n];
    os8 << " " << Gen_I[n];
    os8 << " " << Vul_I[n];
    os8 << " " << Conn_I[n];
    os8 << endl;
  }

}

void NewWeb::species_table3(const char* filename,
			   bool with_plants,double th,Snapshot & snap,
			   link_strength_matrix & intake,
			   bool jolly){
  double area=area_per_compartment/unit_area;

  if(intake.size()==0)
    intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  // Work out number of fish species
  double Sfish = 0.0;
  for(int n=0;n<S;n++){
    if(log10(s[n].bodymass()/unit_mass)>-3){
      Sfish += 1;
    }
  }
  //cout << "Sfish = " << Sfish << endl;

  sequence<std::string> additional_headings;
  sequence<std::string> additional_data;

  int column=0;

  //int S=s.size();
  //double area=area_per_compartment/unit_area;

  additional_headings[column++]=
    "col   (assigned column)";
  additional_headings[column++]=
    "prod   (productivity)";
  additional_headings[column++]=
    "pred   (predation mortality)";
  additional_headings[column++]=
    "comp   (mean dietary overlap for fish species)";
  additional_headings[column++]=
    "nin   (nin for fish species, considering all species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "nout   (nout for fish species, considering all species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninw   (nin for fish species, considering all species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutw   (nout for fish species, considering all species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninf   (nin for fish species, considering just fish species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutf   (nout for fish species, considering just fish species; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "ninwf   (nin for fish species, considering just fish species and weighted; see Bersier et al. (2002) for definition)";
  additional_headings[column++]=
    "noutwf   (nout for fish species, considering just fish species and weighted; see Bersier et al. (2002) for definition)";
  
  // Calculate indicators for each fish species
  sequence<double> b_idotI;
  sequence<double> b_dotjI;
  sequence<double> b_idotIf;
  sequence<double> b_dotjIf;
  double b_dotdotI = 0.0;
  double b_dotdotIf = 0.0;
  sequence<double> H_ink;
  sequence<double> H_outk;
  sequence<double> n_ink;
  sequence<double> n_outk;
  sequence<double> H_inkf;
  sequence<double> H_outkf;
  sequence<double> n_inkf;
  sequence<double> n_outkf;
  double dummy3 = 0.0;
  double dummy4 = 0.0;
  double dummy5 = 0.0;
  double dummy6 = 0.0;
  double dummy7 = 0.0;
  double resolog10M = 0.0;
  double conslog10M = 0.0;
  // First populate b_idotI, b_dotjI, b_idotIf and b_dotjIf.
  // Also calcuate b_dotdotI and b_dotdotIf.
  for(int i=0;i<S;i++){
    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    dummy3 = 0;
    dummy4 = 0;
    dummy5 = 0;
    dummy6 = 0;
    for(int j=0;j<S;j++){
      NewSpecies & cons=s[j];
      conslog10M = log10(cons.bodymass()/unit_mass);
      dummy3 += intake(i,j)/unit_mass/area/(1/year);
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      if(resolog10M > -3){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M > -3)){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M <= -3)){
	dummy7 += intake(i,j)/unit_mass/area/(1/year);
      }
      // End loop for j.
    }
    b_idotI[i] = dummy3;
    b_dotjI[i] = dummy4;
    b_idotIf[i] = dummy5;
    b_dotjIf[i] = dummy6;
    // End loop for i.
  }
  b_dotdotIf = b_dotdotI - dummy7;
  // Now populate H_ink, H_outk, n_ink, n_outk.
  // Also populate H_inkf, H_outkf, n_inkf, n_outkf.
  double dummy100 = 0.0;
  double dummy101 = 0.0;
  double dummy102 = 0.0;
  double dummy103 = 0.0;
  double dummy104 = 0.0;
  double dummy105 = 0.0;
  for(int i=0;i<S;i++){
    dummy102 = 0;
    dummy103 = 0;
    dummy104 = 0;
    dummy105 = 0;
    for(int j=0;j<S;j++){
      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);
      if((dummy101>0)&&(b_dotjI[i]>0)){
	dummy102 = dummy102 - ((dummy101/b_dotjI[i])*
			       (log10(dummy101/b_dotjI[i])/log10(2)));
      }
      if((dummy100>0)&&(b_idotI[i]>0)){
	dummy103 = dummy103 - ((dummy100/b_idotI[i])*
			       (log10(dummy100/b_idotI[i])/log10(2)));
      }
      if((dummy101>0)&&(b_dotjIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)>-3)){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
      }
      if((dummy100>0)&&(b_idotIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)<=-3)){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
      }  
      // End loop for j.
    }
    H_ink[i] = dummy102;
    H_outk[i] = dummy103;
    if(b_idotI[i]>0){
      n_outk[i] = pow(2,H_outk[i]);
    }else{
      n_outk[i] = 0;
    }
    if(b_dotjI[i]>0){
      n_ink[i] = pow(2,H_ink[i]);
    }else{
      n_ink[i] = 0;
    }
    H_inkf[i] = dummy104;
    H_outkf[i] = dummy105;
    if(b_idotIf[i]>0){
      n_outkf[i] = pow(2,H_outkf[i]);
    }else{
      n_outkf[i] = 0;
    }
    if(b_dotjIf[i]>0){
      n_inkf[i] = pow(2,H_inkf[i]);
    }else{
      n_inkf[i] = 0;
    }
    // End loop for i.
  }
  
  for(int n=0;n<S;n++){
    std::ostringstream os;
    double prod = 0.0;
    double pred = 0.0;
    double comp = 0.0;
    int counter1 = 0;
    int counter2 = 0;
    int counter3 = 0;
    int counter4 = 0;
    //os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(7);
    NewSpecies & sp=s[n];     
    if(log10(sp.bodymass()/unit_mass)>-3){
      if(with_plants || !sp.is_a_plant()){
	os.width(1);
	os << assigned_column(n) << " ";
	for(int i=0;i<S;i++){
	  NewSpecies & sp2=s[i];
	  prod += 0.6*(intake(i,n)/unit_mass/area/(1/year));
	  pred += intake(n,i)/unit_mass/area/(1/year);
	  if(log10(sp2.bodymass()/unit_mass)>-3){
	    double comp1 = 0.0;
	    double comp2 = 0.0;
	    double comp3 = 0.0;
	    for(int j=0;j<S;j++){
	      comp1 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	      comp2 += (intake(j,n)/unit_mass/area/(1/year))*(intake(j,n)/unit_mass/area/(1/year));
	      comp3 += (intake(j,i)/unit_mass/area/(1/year))*(intake(j,i)/unit_mass/area/(1/year));
	      // End loop for j.
	    }
	    comp += comp1/(sqrt(comp2*comp3));
	  }
	  // End loop for i.
   	}
	// Work out mean of competition coefficients
	comp = comp/Sfish;
	os.width(5);
	os << prod << " ";
	os.width(5);
	os << pred << " ";
	os.width(5);
	os << comp << " ";
	os.width(5);
	os << n_ink[n] << " ";
	os.width(5);
	os << n_outk[n] << " ";
	os.width(5);
	os << (b_dotjI[n]/b_dotdotI)*n_ink[n] << " ";
	os.width(5);
	os << (b_idotI[n]/b_dotdotI)*n_outk[n] << " ";
	os.width(5);
	os << n_inkf[n] << " ";
	os.width(5);
	os << n_outkf[n] << " ";
	os.width(5);
	os << (b_dotjIf[n]/b_dotdotIf)*n_inkf[n] << " ";
	os.width(5);
	os << (b_idotIf[n]/b_dotdotIf)*n_outkf[n] << " ";
	additional_data[n]=os.str();
      }
    }
  }
  snap.species_table3(filename,with_plants,th,
		     additional_headings,additional_data);

  // Work out mean, standard deviation and skewness of interaction strengths (measured as flows)
  // Exclude interaction strengths of magnitude 0
  // Formula for skewness has been verified using a set of data and R code
  double mean_interaction_strengths = 0.0;
  double sd_interaction_strengths = 0.0;
  double skewness_interaction_strengths = 0.0;
  double sumnoofentries = 0.0;
  double sumX = 0.0;
  double sumX2 = 0.0;
  double sumX3 = 0.0;
  double interaction_strength_dummy = 0.0;
  for(int i=0;i<S;i++){
    for(int j=0;j<S;j++){
      if((intake(i,j)/area/(1/year))>0){
	interaction_strength_dummy = intake(i,j)/area/(1/year);
	sumnoofentries = sumnoofentries+1;
	sumX = sumX+interaction_strength_dummy;
	sumX2 = sumX2+(interaction_strength_dummy*interaction_strength_dummy);
	sumX3 = sumX3+(interaction_strength_dummy*interaction_strength_dummy*interaction_strength_dummy);
      }
    }
  }
  mean_interaction_strengths = sumX/sumnoofentries;
  sd_interaction_strengths = sqrt((sumX2/sumnoofentries)-(mean_interaction_strengths*mean_interaction_strengths));
  skewness_interaction_strengths = ((sumX3/sumnoofentries)-(3*mean_interaction_strengths*sd_interaction_strengths*sd_interaction_strengths)
				    -(mean_interaction_strengths*mean_interaction_strengths*mean_interaction_strengths))/
    (sd_interaction_strengths*sd_interaction_strengths*sd_interaction_strengths);
  cout << "Ltot = " << sumnoofentries << endl;
  cout << "LDtot = " << sumnoofentries/S << endl;
  cout << "Ctot = " << sumnoofentries/(S*S) << endl;
  cout << "mean_interaction_strengths = " << mean_interaction_strengths << endl;
  cout << "sd_interaction_strengths = " << sd_interaction_strengths << endl;
  cout << "skewness_interaction_strengths = " << skewness_interaction_strengths << endl;
  
  // Work out mean, standard deviation and skewness of interaction strengths (measured as flows)
  // divided by biomass densities of predator and prey
  // Exclude interaction strengths of magnitude 0
  // Formula for skewness has been verified using a set of data and R code
  double mean_interaction_strengths2 = 0.0;
  double sd_interaction_strengths2 = 0.0;
  double skewness_interaction_strengths2 = 0.0;
  sumnoofentries = 0.0;
  sumX = 0.0;
  sumX2 = 0.0;
  sumX3 = 0.0;
  interaction_strength_dummy = 0.0;
  for(int i=0;i<S;i++){
    for(int j=0;j<S;j++){
      if((intake(i,j)/area/(1/year))>0){
	interaction_strength_dummy = (intake(i,j)/area/(1/year))/((biomass_B(i)/unit_mass/area)*(biomass_B(j)/unit_mass/area));
	sumnoofentries = sumnoofentries+1;
	sumX = sumX+interaction_strength_dummy;
	sumX2 = sumX2+(interaction_strength_dummy*interaction_strength_dummy);
	sumX3 = sumX3+(interaction_strength_dummy*interaction_strength_dummy*interaction_strength_dummy);
      }
    }
  }
  mean_interaction_strengths2 = sumX/sumnoofentries;
  sd_interaction_strengths2 = sqrt((sumX2/sumnoofentries)-(mean_interaction_strengths2*mean_interaction_strengths2));
  skewness_interaction_strengths2 = ((sumX3/sumnoofentries)-(3*mean_interaction_strengths2*sd_interaction_strengths2*sd_interaction_strengths2)
				    -(mean_interaction_strengths2*mean_interaction_strengths2*mean_interaction_strengths2))/
    (sd_interaction_strengths2*sd_interaction_strengths2*sd_interaction_strengths2);
  cout << "mean_interaction_strengths2 = " << mean_interaction_strengths2 << endl;
  cout << "sd_interaction_strengths2 = " << sd_interaction_strengths2 << endl;
  cout << "skewness_interaction_strengths2 = " << skewness_interaction_strengths2 << endl;

  // Work out mean, standard deviation and skewness of interaction strengths (measured as flows)
  // divided by biomass density of predator
  // Exclude interaction strengths of magnitude 0
  // Formula for skewness has been verified using a set of data and R code
  double mean_interaction_strengths3 = 0.0;
  double sd_interaction_strengths3 = 0.0;
  double skewness_interaction_strengths3 = 0.0;
  sumnoofentries = 0.0;
  sumX = 0.0;
  sumX2 = 0.0;
  sumX3 = 0.0;
  interaction_strength_dummy = 0.0;
  for(int i=0;i<S;i++){
    for(int j=0;j<S;j++){
      if((intake(i,j)/area/(1/year))>0){
	interaction_strength_dummy = (intake(i,j)/area/(1/year))/(biomass_B(j)/unit_mass/area);
	sumnoofentries = sumnoofentries+1;
	sumX = sumX+interaction_strength_dummy;
	sumX2 = sumX2+(interaction_strength_dummy*interaction_strength_dummy);
	sumX3 = sumX3+(interaction_strength_dummy*interaction_strength_dummy*interaction_strength_dummy);
      }
    }
  }
  mean_interaction_strengths3 = sumX/sumnoofentries;
  sd_interaction_strengths3 = sqrt((sumX2/sumnoofentries)-(mean_interaction_strengths3*mean_interaction_strengths3));
  skewness_interaction_strengths3 = ((sumX3/sumnoofentries)-(3*mean_interaction_strengths3*sd_interaction_strengths3*sd_interaction_strengths3)
				    -(mean_interaction_strengths3*mean_interaction_strengths3*mean_interaction_strengths3))/
    (sd_interaction_strengths3*sd_interaction_strengths3*sd_interaction_strengths3);
  cout << "mean_interaction_strengths3 = " << mean_interaction_strengths3 << endl;
  cout << "sd_interaction_strengths3 = " << sd_interaction_strengths3 << endl;
  cout << "skewness_interaction_strengths3 = " << skewness_interaction_strengths3 << endl;
  
}

void NewWeb::link_table(const char* filename,double th,
			link_strength_matrix & intake,
			const sequence<double> & biomass_B){
  double area=area_per_compartment/unit_area;

  if(intake.size()==0)
    intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();
  
  // compute the trophic level vector:
  NewMatrix F=ifrac;
  NewVector ones(S);
  for(int i=S;i-->0;)
    ones[i]=1;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=solve(D,ones);
  
  std::cout << "writing a link table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram;"
	    << " area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "the meaning of the columns is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": hc (consumer trophic level)" << std::endl;
  std::cout << "column " << column++ << ": hr (resource trophic level)" << std::endl;
  std::cout << "column " << column++ << ": Mc (log10 consumer body mass)" << std::endl;
  std::cout << "column " << column++ << ": Mr (log10 resource body mass)" << std::endl;
  std::cout << "column " << column++ << ": Bc (log10 consumer biomass density)" << std::endl;
  std::cout << "column " << column++ << ": Br (log10 resource biomass density)" << std::endl;
  std::cout << "column " << column++ << ": cc (consumer column)" << std::endl;
  std::cout << "column " << column++ << ": cr (resource column)" << std::endl;
  std::cout << "column " << column++ << ": F  (biomass flow density [per year])" << std::endl;
  std::cout << "column " << column++ << ": f  (consumer diet fraction)" << std::endl;
  std::cout << "column " << column++ << ": c  (log10 coupling coefficient)" << std::endl;

  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */


  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(4);
  
  for(int n=0;n<S;n++){
    for(int m=0;m<S;m++){
      if(ifrac(m,n)>th){
	NewSpecies & cons=s[n];
	NewSpecies & reso=s[m];
	os.width(9);
	os << level[n] << " ";
	os.width(9);
	os << level[m] << " ";
	os.width(9);
	os << log10(cons.bodymass()/unit_mass) << " ";
	os.width(9);
	os << log10(reso.bodymass()/unit_mass) << " ";
	os.width(9);
	os << log10(biomass_B(n)/unit_mass/area) << " ";
	os.width(9);
	os << log10(biomass_B(m)/unit_mass/area) << " ";
	os.width(9);
	os << assigned_column[n] << " ";
	os.width(9);
	os << assigned_column[m] << " ";
	os.width(9);
	os << intake(m,n)/unit_mass/area/(1/year) << " ";
	os.width(9);
	os << ifrac(m,n) << " ";
 	os.width(9);
 	os << log10(precomputed[n].c[m]) << " ";
	os << std::endl;
      }
    }
  }
}

void NewWeb::link_table2(const char* filename,double th,
			link_strength_matrix & intake,
			const sequence<double> & biomass_B){
  double area=area_per_compartment/unit_area;

  if(intake.size()==0)
    intake=intake_matrix();
  link_strength_matrix ifrac=in_fraction(intake);
  int S=ifrac.size();

  // Determine whether to calculate TL indicators or not
  int compute_TL_indicators = 0;
  
  // compute the trophic level vector:
  if(compute_TL_indicators==1){

    NewMatrix F=ifrac;
    NewVector ones(S);
    for(int i=S;i-->0;)
      ones[i]=1;
    NewMatrix D=NewIdentityMatrix(S,S)-F;
    NewVector level=solve(D,ones);

    double dummy1 = 0.0;
    int dummy2 = 0;
    sequence<double> rho_TL;
    for(int i=0;i<45;i++){
      rho_TL[i]=0;
    }
    for(int i=0;i<S;i++){
      dummy1 = (level[i] - 1)/0.1;
      //cout << "dummy1 = " << dummy1 << endl;
      dummy2 = floor(dummy1);
      //cout << "dummy2 = " << dummy2 << endl;
      rho_TL[dummy2] += biomass_B(i)/unit_mass/area;
    }
    //for(int i=0;i<45;i++){
    //cout << "rho_TL[" << i << "] = " << rho_TL[i] << endl;
    //}
    cout << "rho_1_1point1_TL = " << rho_TL[0] << endl;
    cout << "rho_1point1_1point2_TL = " << rho_TL[1] << endl;
    cout << "rho_1point2_1point3_TL = " << rho_TL[2] << endl;
    cout << "rho_1point3_1point4_TL = " << rho_TL[3] << endl;
    cout << "rho_1point4_1point5_TL = " << rho_TL[4] << endl;
    cout << "rho_1point5_1point6_TL = " << rho_TL[5] << endl;
    cout << "rho_1point6_1point7_TL = " << rho_TL[6] << endl;
    cout << "rho_1point7_1point8_TL = " << rho_TL[7] << endl;
    cout << "rho_1point8_1point9_TL = " << rho_TL[8] << endl;
    cout << "rho_1point9_2_TL = " << rho_TL[9] << endl;
    cout << "rho_2_2point1_TL = " << rho_TL[10] << endl;
    cout << "rho_2point1_2point2_TL = " << rho_TL[11] << endl;
    cout << "rho_2point2_2point3_TL = " << rho_TL[12] << endl;
    cout << "rho_2point3_2point4_TL = " << rho_TL[13] << endl;
    cout << "rho_2point4_2point5_TL = " << rho_TL[14] << endl;
    cout << "rho_2point5_2point6_TL = " << rho_TL[15] << endl;
    cout << "rho_2point6_2point7_TL = " << rho_TL[16] << endl;
    cout << "rho_2point7_2point8_TL = " << rho_TL[17] << endl;
    cout << "rho_2point8_2point9_TL = " << rho_TL[18] << endl;
    cout << "rho_2point9_3_TL = " << rho_TL[19] << endl;
    cout << "rho_3_3point1_TL = " << rho_TL[20] << endl;
    cout << "rho_3point1_3point2_TL = " << rho_TL[21] << endl;
    cout << "rho_3point2_3point3_TL = " << rho_TL[22] << endl;
    cout << "rho_3point3_3point4_TL = " << rho_TL[23] << endl;
    cout << "rho_3point4_3point5_TL = " << rho_TL[24] << endl;
    cout << "rho_3point5_3point6_TL = " << rho_TL[25] << endl;
    cout << "rho_3point6_3point7_TL = " << rho_TL[26] << endl;
    cout << "rho_3point7_3point8_TL = " << rho_TL[27] << endl;
    cout << "rho_3point8_3point9_TL = " << rho_TL[28] << endl;
    cout << "rho_3point9_4_TL = " << rho_TL[29] << endl;
    cout << "rho_4_4point1_TL = " << rho_TL[30] << endl;
    cout << "rho_4point1_4point2_TL = " << rho_TL[31] << endl;
    cout << "rho_4point2_4point3_TL = " << rho_TL[32] << endl;
    cout << "rho_4point3_4point4_TL = " << rho_TL[33] << endl;
    cout << "rho_4point4_4point5_TL = " << rho_TL[34] << endl;
    cout << "rho_4point5_4point6_TL = " << rho_TL[35] << endl;
    cout << "rho_4point6_4point7_TL = " << rho_TL[36] << endl;
    cout << "rho_4point7_4point8_TL = " << rho_TL[37] << endl;
    cout << "rho_4point8_4point9_TL = " << rho_TL[38] << endl;
    cout << "rho_4point9_5_TL = " << rho_TL[39] << endl;
    cout << "rho_5_5point1_TL = " << rho_TL[40] << endl;
    cout << "rho_5point1_5point2_TL = " << rho_TL[41] << endl;
    cout << "rho_5point2_5point3_TL = " << rho_TL[42] << endl;
    cout << "rho_5point3_5point4_TL = " << rho_TL[43] << endl;
    cout << "rho_5point4_5point5_TL = " << rho_TL[44] << endl;

    // Flows from organisms of lower or equal TL to those of higher TL
    // Flows from organisms of higher TL to those of lower or equal TL
    // Then work out b_xTL.
    double S_c = number_of_animals();
    double S_f = 0.0;
    sequence<double> b_lowtohigh_TL;
    sequence<double> b_lowtohigh_TL_Mcorrected;
    sequence<double> b_hightolow_TL;
    sequence<double> b_hightolow_TL_Mcorrected;
    for(int i=0;i<36;i++){
      b_lowtohigh_TL[i]=0;
      b_lowtohigh_TL_Mcorrected[i]=0;
      b_hightolow_TL[i]=0;
      b_hightolow_TL_Mcorrected[i]=0;
    }
    double resolog10M = 0.0;
    double conslog10M = 0.0;
    sequence<double> b_dotjI;
    sequence<double> b_dotjIf;
    double b_dotdotI = 0.0;
    double b_dotdotIf = 0.0;
    double resoTL = 0.0;
    double consTL = 0.0;
    double dummy4 = 0.0;
    double dummy6 = 0.0;
    double dummy7 = 0.0;

    for(int i=0;i<S;i++){

      NewSpecies & reso=s[i];
      resolog10M = log10(reso.bodymass()/unit_mass);
      resoTL = level[i];
      dummy4 = 0;
      dummy6 = 0;     

      for(int j=0;j<S;j++){
	
	NewSpecies & cons=s[j];	
	conslog10M = log10(cons.bodymass()/unit_mass);
	consTL = level[j];
	// Update sums for calculation of flows for species
	dummy4 += intake(j,i)/unit_mass/area/(1/year);
	b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
	if(resolog10M > -3){
	  dummy6 += intake(j,i)/unit_mass/area/(1/year);
	}
	if((resolog10M <= -3)&&(conslog10M > -3)){
	  dummy6 += intake(j,i)/unit_mass/area/(1/year);
	}
	if((resolog10M <= -3)&&(conslog10M <= -3)){
	dummy7 += intake(i,j)/unit_mass/area/(1/year);
	}
         
	if((resoTL <= 1.9)&&(consTL > 1.9)){
	  b_lowtohigh_TL[0] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2)&&(consTL > 2)){
	  b_lowtohigh_TL[1] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.1)&&(consTL > 2.1)){
	  b_lowtohigh_TL[2] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.2)&&(consTL > 2.2)){
	  b_lowtohigh_TL[3] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.3)&&(consTL > 2.3)){
	  b_lowtohigh_TL[4] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.4)&&(consTL > 2.4)){
	  b_lowtohigh_TL[5] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.5)&&(consTL > 2.5)){
	  b_lowtohigh_TL[6] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.6)&&(consTL > 2.6)){
	  b_lowtohigh_TL[7] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.7)&&(consTL > 2.7)){
	  b_lowtohigh_TL[8] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.8)&&(consTL > 2.8)){
	  b_lowtohigh_TL[9] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 2.9)&&(consTL > 2.9)){
	  b_lowtohigh_TL[10] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3)&&(consTL > 3)){
	  b_lowtohigh_TL[11] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.1)&&(consTL > 3.1)){
	  b_lowtohigh_TL[12] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.2)&&(consTL > 3.2)){
	  b_lowtohigh_TL[13] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.3)&&(consTL > 3.3)){
	  b_lowtohigh_TL[14] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.4)&&(consTL > 3.4)){
	  b_lowtohigh_TL[15] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.5)&&(consTL > 3.5)){
	  b_lowtohigh_TL[16] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.6)&&(consTL > 3.6)){
	  b_lowtohigh_TL[17] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.7)&&(consTL > 3.7)){
	  b_lowtohigh_TL[18] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.8)&&(consTL > 3.8)){
	  b_lowtohigh_TL[19] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 3.9)&&(consTL > 3.9)){
	  b_lowtohigh_TL[20] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4)&&(consTL > 4)){
	  b_lowtohigh_TL[21] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.1)&&(consTL > 4.1)){
	  b_lowtohigh_TL[22] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.2)&&(consTL > 4.2)){
	  b_lowtohigh_TL[23] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.3)&&(consTL > 4.3)){
	  b_lowtohigh_TL[24] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.4)&&(consTL > 4.4)){
	  b_lowtohigh_TL[25] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.5)&&(consTL > 4.5)){
	  b_lowtohigh_TL[26] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.6)&&(consTL > 4.6)){
	  b_lowtohigh_TL[27] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.7)&&(consTL > 4.7)){
	  b_lowtohigh_TL[28] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.8)&&(consTL > 4.8)){
	  b_lowtohigh_TL[29] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 4.9)&&(consTL > 4.9)){
	  b_lowtohigh_TL[30] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 5)&&(consTL > 5)){
	  b_lowtohigh_TL[31] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 5.1)&&(consTL > 5.1)){
	  b_lowtohigh_TL[32] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 5.2)&&(consTL > 5.2)){
	  b_lowtohigh_TL[33] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 5.3)&&(consTL > 5.3)){
	  b_lowtohigh_TL[34] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL <= 5.4)&&(consTL > 5.4)){
	  b_lowtohigh_TL[35] += intake(i,j)/unit_mass/area/(1/year);
	  b_lowtohigh_TL_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	
	if((resoTL > 1.9)&&(consTL <= 1.9)){
	  b_hightolow_TL[0] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2)&&(consTL <= 2)){
	  b_hightolow_TL[1] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.1)&&(consTL <= 2.1)){
	  b_hightolow_TL[2] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.2)&&(consTL <= 2.2)){
	  b_hightolow_TL[3] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.3)&&(consTL <= 2.3)){
	  b_hightolow_TL[4] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.4)&&(consTL <= 2.4)){
	  b_hightolow_TL[5] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.5)&&(consTL <= 2.5)){
	  b_hightolow_TL[6] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.6)&&(consTL <= 2.6)){
	  b_hightolow_TL[7] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.7)&&(consTL <= 2.7)){
	  b_hightolow_TL[8] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.8)&&(consTL <= 2.8)){
	  b_hightolow_TL[9] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 2.9)&&(consTL <= 2.9)){
	  b_hightolow_TL[10] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3)&&(consTL <= 3)){
	  b_hightolow_TL[11] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.1)&&(consTL <= 3.1)){
	  b_hightolow_TL[12] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.2)&&(consTL <= 3.2)){
	  b_hightolow_TL[13] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.3)&&(consTL <= 3.3)){
	  b_hightolow_TL[14] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.4)&&(consTL <= 3.4)){
	  b_hightolow_TL[15] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.5)&&(consTL <= 3.5)){
	  b_hightolow_TL[16] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.6)&&(consTL <= 3.6)){
	  b_hightolow_TL[17] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.7)&&(consTL <= 3.7)){
	  b_hightolow_TL[18] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.8)&&(consTL <= 3.8)){
	  b_hightolow_TL[19] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 3.9)&&(consTL <= 3.9)){
	  b_hightolow_TL[20] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4)&&(consTL <= 4)){
	  b_hightolow_TL[21] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.1)&&(consTL <= 4.1)){
	  b_hightolow_TL[22] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.2)&&(consTL <= 4.2)){
	  b_hightolow_TL[23] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.3)&&(consTL <= 4.3)){
	  b_hightolow_TL[24] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.4)&&(consTL <= 4.4)){
	  b_hightolow_TL[25] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.5)&&(consTL <= 4.5)){
	  b_hightolow_TL[26] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.6)&&(consTL <= 4.6)){
	  b_hightolow_TL[27] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.7)&&(consTL <= 4.7)){
	  b_hightolow_TL[28] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.8)&&(consTL <= 4.8)){
	  b_hightolow_TL[29] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 4.9)&&(consTL <= 4.9)){
	  b_hightolow_TL[30] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 5)&&(consTL <= 5)){
	  b_hightolow_TL[31] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 5.1)&&(consTL <= 5.1)){
	  b_hightolow_TL[32] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 5.2)&&(consTL <= 5.2)){
	  b_hightolow_TL[33] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 5.3)&&(consTL <= 5.3)){
	  b_hightolow_TL[34] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}
	if((resoTL > 5.4)&&(consTL <= 5.4)){
	  b_hightolow_TL[35] += intake(i,j)/unit_mass/area/(1/year);
	  b_hightolow_TL_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
	}

	// close loop for j
      }

      // Update vectors of flows for species
      b_dotjI[i] = dummy4;
      b_dotjIf[i] = dummy6;
      
      // Update species richness
      if(resolog10M > -3){
	S_f += 1;
      }
      
      // close loop for i 
    }

    b_dotdotIf = b_dotdotI - dummy7;

    cout << "b_1point9_TL = " << b_lowtohigh_TL[0]-b_hightolow_TL[0] << endl;
    cout << "b_2_TL = " << b_lowtohigh_TL[1]-b_hightolow_TL[1] << endl;
    cout << "b_2point1_TL = " << b_lowtohigh_TL[2]-b_hightolow_TL[2] << endl;
    cout << "b_2point2_TL = " << b_lowtohigh_TL[3]-b_hightolow_TL[3] << endl;
    cout << "b_2point3_TL = " << b_lowtohigh_TL[4]-b_hightolow_TL[4] << endl;
    cout << "b_2point4_TL = " << b_lowtohigh_TL[5]-b_hightolow_TL[5] << endl;
    cout << "b_2point5_TL = " << b_lowtohigh_TL[6]-b_hightolow_TL[6] << endl;
    cout << "b_2point6_TL = " << b_lowtohigh_TL[7]-b_hightolow_TL[7] << endl;
    cout << "b_2point7_TL = " << b_lowtohigh_TL[8]-b_hightolow_TL[8] << endl;
    cout << "b_2point8_TL = " << b_lowtohigh_TL[9]-b_hightolow_TL[9] << endl;
    cout << "b_2point9_TL = " << b_lowtohigh_TL[10]-b_hightolow_TL[10] << endl;
    cout << "b_3_TL = " << b_lowtohigh_TL[11]-b_hightolow_TL[11] << endl;
    cout << "b_3point1_TL = " << b_lowtohigh_TL[12]-b_hightolow_TL[12] << endl;
    cout << "b_3point2_TL = " << b_lowtohigh_TL[13]-b_hightolow_TL[13] << endl;
    cout << "b_3point3_TL = " << b_lowtohigh_TL[14]-b_hightolow_TL[14] << endl;
    cout << "b_3point4_TL = " << b_lowtohigh_TL[15]-b_hightolow_TL[15] << endl;
    cout << "b_3point5_TL = " << b_lowtohigh_TL[16]-b_hightolow_TL[16] << endl;
    cout << "b_3point6_TL = " << b_lowtohigh_TL[17]-b_hightolow_TL[17] << endl;
    cout << "b_3point7_TL = " << b_lowtohigh_TL[18]-b_hightolow_TL[18] << endl;
    cout << "b_3point8_TL = " << b_lowtohigh_TL[19]-b_hightolow_TL[19] << endl;
    cout << "b_3point9_TL = " << b_lowtohigh_TL[20]-b_hightolow_TL[20] << endl;
    cout << "b_4_TL = " << b_lowtohigh_TL[21]-b_hightolow_TL[21] << endl;
    cout << "b_4point1_TL = " << b_lowtohigh_TL[22]-b_hightolow_TL[22] << endl;
    cout << "b_4point2_TL = " << b_lowtohigh_TL[23]-b_hightolow_TL[23] << endl;
    cout << "b_4point3_TL = " << b_lowtohigh_TL[24]-b_hightolow_TL[24] << endl;
    cout << "b_4point4_TL = " << b_lowtohigh_TL[25]-b_hightolow_TL[25] << endl;
    cout << "b_4point5_TL = " << b_lowtohigh_TL[26]-b_hightolow_TL[26] << endl;
    cout << "b_4point6_TL = " << b_lowtohigh_TL[27]-b_hightolow_TL[27] << endl;
    cout << "b_4point7_TL = " << b_lowtohigh_TL[28]-b_hightolow_TL[28] << endl;
    cout << "b_4point8_TL = " << b_lowtohigh_TL[29]-b_hightolow_TL[29] << endl;
    cout << "b_4point9_TL = " << b_lowtohigh_TL[30]-b_hightolow_TL[30] << endl;
    cout << "b_5_TL = " << b_lowtohigh_TL[31]-b_hightolow_TL[31] << endl;
    cout << "b_5point1_TL = " << b_lowtohigh_TL[32]-b_hightolow_TL[32] << endl;
    cout << "b_5point2_TL = " << b_lowtohigh_TL[33]-b_hightolow_TL[33] << endl;
    cout << "b_5point3_TL = " << b_lowtohigh_TL[34]-b_hightolow_TL[34] << endl;
    cout << "b_5point4_TL = " << b_lowtohigh_TL[35]-b_hightolow_TL[35] << endl;
    cout << "b_1point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[0]-b_hightolow_TL_Mcorrected[0] << endl;
    cout << "b_2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[1]-b_hightolow_TL_Mcorrected[1] << endl;
    cout << "b_2point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[2]-b_hightolow_TL_Mcorrected[2] << endl;
    cout << "b_2point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[3]-b_hightolow_TL_Mcorrected[3] << endl;
    cout << "b_2point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[4]-b_hightolow_TL_Mcorrected[4] << endl;
    cout << "b_2point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[5]-b_hightolow_TL_Mcorrected[5] << endl;
    cout << "b_2point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[6]-b_hightolow_TL_Mcorrected[6] << endl;
    cout << "b_2point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[7]-b_hightolow_TL_Mcorrected[7] << endl;
    cout << "b_2point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[8]-b_hightolow_TL_Mcorrected[8] << endl;
    cout << "b_2point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[9]-b_hightolow_TL_Mcorrected[9] << endl;
    cout << "b_2point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[10]-b_hightolow_TL_Mcorrected[10] << endl;
    cout << "b_3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[11]-b_hightolow_TL_Mcorrected[11] << endl;
    cout << "b_3point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[12]-b_hightolow_TL_Mcorrected[12] << endl;
    cout << "b_3point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[13]-b_hightolow_TL_Mcorrected[13] << endl;
    cout << "b_3point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[14]-b_hightolow_TL_Mcorrected[14] << endl;
    cout << "b_3point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[15]-b_hightolow_TL_Mcorrected[15] << endl;
    cout << "b_3point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[16]-b_hightolow_TL_Mcorrected[16] << endl;
    cout << "b_3point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[17]-b_hightolow_TL_Mcorrected[17] << endl;
    cout << "b_3point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[18]-b_hightolow_TL_Mcorrected[18] << endl;
    cout << "b_3point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[19]-b_hightolow_TL_Mcorrected[19] << endl;
    cout << "b_3point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[20]-b_hightolow_TL_Mcorrected[20] << endl;
    cout << "b_4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[21]-b_hightolow_TL_Mcorrected[21] << endl;
    cout << "b_4point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[22]-b_hightolow_TL_Mcorrected[22] << endl;
    cout << "b_4point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[23]-b_hightolow_TL_Mcorrected[23] << endl;
    cout << "b_4point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[24]-b_hightolow_TL_Mcorrected[24] << endl;
    cout << "b_4point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[25]-b_hightolow_TL_Mcorrected[25] << endl;
    cout << "b_4point5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[26]-b_hightolow_TL_Mcorrected[26] << endl;
    cout << "b_4point6_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[27]-b_hightolow_TL_Mcorrected[27] << endl;
    cout << "b_4point7_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[28]-b_hightolow_TL_Mcorrected[28] << endl;
    cout << "b_4point8_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[29]-b_hightolow_TL_Mcorrected[29] << endl;
    cout << "b_4point9_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[30]-b_hightolow_TL_Mcorrected[30] << endl;
    cout << "b_5_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[31]-b_hightolow_TL_Mcorrected[31] << endl;
    cout << "b_5point1_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[32]-b_hightolow_TL_Mcorrected[32] << endl;
    cout << "b_5point2_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[33]-b_hightolow_TL_Mcorrected[33] << endl;
    cout << "b_5point3_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[34]-b_hightolow_TL_Mcorrected[34] << endl;
    cout << "b_5point4_TL_Mcorrected = " << b_lowtohigh_TL_Mcorrected[35]-b_hightolow_TL_Mcorrected[35] << endl;

    sequence<double> omn_j;
    double dummy100 = 0.0;
    double dummy101 = 0.0;
    double dummy110 = 0.0;
    double dummy111 = 0.0;
    double dummy112 = 0.0;
    double dummy113 = 0.0;
    double dummy114 = 0.0;
    double dummy115 = 0.0;

    for(int i=0;i<S;i++){

      // Reset sums for vectors
      dummy110 = 0;

      for(int j=0;j<S;j++){

	dummy100 = intake(i,j)/unit_mass/area/(1/year);
	dummy101 = intake(j,i)/unit_mass/area/(1/year);
      
	if(b_dotjI[i]>0){
	  dummy110 += (dummy101/b_dotjI[i])*pow(level[j]-level[i],2);
	}

	if((level[i]>=1.5)&&(level[i]<2.5)&&(level[j]>=2.5)&&(level[j]<3.5)){
	  dummy113 += dummy100;
	}
	if((level[i]>=2.5)&&(level[i]<3.5)&&(level[j]>=3.5)&&(level[j]<4.5)){
	  dummy114 += dummy100;
	}
	if((level[i]>=3.5)&&(level[i]<4.5)&&(level[j]>=4.5)&&(level[j]<5.5)){
	  dummy115 += dummy100;
	}
      
	// end loop for j
      }

      omn_j[i] = sqrt(dummy110);
      
      if(level[i]==1){ dummy112 += (s[i].the_GP)/unit_mass/area/(1/year); } 
    
      // end loop for i
    }

    double omn_prime = 0.0;
    double omn = 0.0;
    double omn_f_prime = 0.0;
    double omn_f = 0.0;
    double epsilon_3 = 0.0;
    double epsilon_4 = 0.0;
    double epsilon_5 = 0.0;
    double dummy204 = 0.0;
    double dummy205 = 0.0;
    double dummy206 = 0.0;
    double dummy207 = 0.0;
    double dummy208 = 0.0;    

    for(int i=0;i<S;i++){
      if(log10(s[i].bodymass()/unit_mass)>-3){
	dummy208 += b_dotjI[i];
      }
    }
    for(int i=0;i<S;i++){
      dummy204 += omn_j[i];
      dummy205 += (b_dotjI[i]/b_dotdotI)*omn_j[i];
      if(log10(s[i].bodymass()/unit_mass)>-3){
	dummy206 += omn_j[i];
	dummy207 += (b_dotjI[i]/dummy208)*omn_j[i];
      }
    }
    // It has been checked that dividing by S gives decimal
    omn_prime = dummy204/S_c;
    omn = dummy205;
    omn_f_prime = dummy206/S_f;
    omn_f = dummy207;
    epsilon_3 = dummy113/dummy112;
    epsilon_4 = dummy114/dummy113;
    epsilon_5 = dummy115/dummy114;

    cout << "omn_prime = " << omn_prime << endl;
    cout << "omn = " << omn << endl;
    cout << "omn_f_prime = " << omn_f_prime << endl;
    cout << "omn_f = " << omn_f << endl;
  

    // end if statement for TL indicators
  }
  

  /*std::cout << "writing a link table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram;"
	    << " area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "the meaning of the columns is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": hc (consumer trophic level)" << std::endl;
  std::cout << "column " << column++ << ": hr (resource trophic level)" << std::endl;
  std::cout << "column " << column++ << ": Mc (log10 consumer body mass)" << std::endl;
  std::cout << "column " << column++ << ": Mr (log10 resource body mass)" << std::endl;
  std::cout << "column " << column++ << ": Bc (log10 consumer biomass density)" << std::endl;
  std::cout << "column " << column++ << ": Br (log10 resource biomass density)" << std::endl;
  std::cout << "column " << column++ << ": cc (consumer column)" << std::endl;
  std::cout << "column " << column++ << ": cr (resource column)" << std::endl;
  std::cout << "column " << column++ << ": F  (biomass flow density [per year])" << std::endl;
  std::cout << "column " << column++ << ": f  (consumer diet fraction)" << std::endl;
  std::cout << "column " << column++ << ": c  (log10 coupling coefficient)" << std::endl;*/

  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */


  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(4);

  // Print out species richness
  cout << "S = " << S << endl;
  double S_p = number_of_plants();
  double S_c = number_of_animals();
  cout << "S_p = " << S_p << endl;
  cout << "S_c = " << S_c << endl;

  // The code below calculates and prints out the stock indicators

  sequence<double> rho_M;
  double dummy1 = 0.0;
  int dummy2 = 0;
  for(int i=0;i<71;i++){
    rho_M[i]=0;
  }
  for(int i=0;i<S;i++){
    NewSpecies & speciesM=s[i];
    //cout << "log10(speciesM.bodymass()/unit_mass)[" << i << "] = " << log10(speciesM.bodymass()/unit_mass) << endl;
    dummy1 = (log10(speciesM.bodymass()/unit_mass) + 15)/0.25;
    //cout << "dummy1 = " << dummy1 << endl;
    dummy2 = floor(dummy1);
    //cout << "dummy2 = " << dummy2 << endl;
    rho_M[dummy2] += biomass_B(i)/unit_mass/area;
  }
  //for(int i=0;i<71;i++){
  //cout << "rho_M[" << i << "] = " << rho_M[i] << endl;
  //}
  cout << "rho_minus15_minus14point75_M = " << rho_M[0] << endl;
  cout << "rho_minus14point75_minus14point5_M = " << rho_M[1] << endl;
  cout << "rho_minus14point5_minus14point25_M = " << rho_M[2] << endl;
  cout << "rho_minus14point25_minus14_M = " << rho_M[3] << endl;
  cout << "rho_minus14_minus13point75_M = " << rho_M[4] << endl;
  cout << "rho_minus13point75_minus13point5_M = " << rho_M[5] << endl;
  cout << "rho_minus13point5_minus13point25_M = " << rho_M[6] << endl;
  cout << "rho_minus13point25_minus13_M = " << rho_M[7] << endl;
  cout << "rho_minus13_minus12point75_M = " << rho_M[8] << endl;
  cout << "rho_minus12point75_minus12point5_M = " << rho_M[9] << endl;
  cout << "rho_minus12point5_minus12point25_M = " << rho_M[10] << endl;
  cout << "rho_minus12point25_minus12_M = " << rho_M[11] << endl;
  cout << "rho_minus12_minus11point75_M = " << rho_M[12] << endl;
  cout << "rho_minus11point75_minus11point5_M = " << rho_M[13] << endl;
  cout << "rho_minus11point5_minus11point25_M = " << rho_M[14] << endl;
  cout << "rho_minus11point25_minus11_M = " << rho_M[15] << endl;
  cout << "rho_minus11_minus10point75_M = " << rho_M[16] << endl;
  cout << "rho_minus10point75_minus10point5_M = " << rho_M[17] << endl;
  cout << "rho_minus10point5_minus10point25_M = " << rho_M[18] << endl;
  cout << "rho_minus10point25_minus10_M = " << rho_M[19] << endl;
  cout << "rho_minus10_minus9point75_M = " << rho_M[20] << endl;
  cout << "rho_minus9point75_minus9point5_M = " << rho_M[21] << endl;
  cout << "rho_minus9point5_minus9point25_M = " << rho_M[22] << endl;
  cout << "rho_minus9point25_minus9_M = " << rho_M[23] << endl;
  cout << "rho_minus9_minus8point75_M = " << rho_M[24] << endl;
  cout << "rho_minus8point75_minus8point5_M = " << rho_M[25] << endl;
  cout << "rho_minus8point5_minus8point25_M = " << rho_M[26] << endl;
  cout << "rho_minus8point25_minus8_M = " << rho_M[27] << endl;
  cout << "rho_minus8_minus7point75_M = " << rho_M[28] << endl;
  cout << "rho_minus7point75_minus7point5_M = " << rho_M[29] << endl;
  cout << "rho_minus7point5_minus7point25_M = " << rho_M[30] << endl;
  cout << "rho_minus7point25_minus7_M = " << rho_M[31] << endl;
  cout << "rho_minus7_minus6point75_M = " << rho_M[32] << endl;
  cout << "rho_minus6point75_minus6point5_M = " << rho_M[33] << endl;
  cout << "rho_minus6point5_minus6point25_M = " << rho_M[34] << endl;
  cout << "rho_minus6point25_minus6_M = " << rho_M[35] << endl;
  cout << "rho_minus6_minus5point75_M = " << rho_M[36] << endl;
  cout << "rho_minus5point75_minus5point5_M = " << rho_M[37] << endl;
  cout << "rho_minus5point5_minus5point25_M = " << rho_M[38] << endl;
  cout << "rho_minus5point25_minus5_M = " << rho_M[39] << endl;
  cout << "rho_minus5_minus4point75_M = " << rho_M[40] << endl;
  cout << "rho_minus4point75_minus4point5_M = " << rho_M[41] << endl;
  cout << "rho_minus4point5_minus4point25_M = " << rho_M[42] << endl;
  cout << "rho_minus4point25_minus4_M = " << rho_M[43] << endl;
  cout << "rho_minus4_minus3point75_M = " << rho_M[44] << endl;
  cout << "rho_minus3point75_minus3point5_M = " << rho_M[45] << endl;
  cout << "rho_minus3point5_minus3point25_M = " << rho_M[46] << endl;
  cout << "rho_minus3point25_minus3_M = " << rho_M[47] << endl;
  cout << "rho_minus3_minus2point75_M = " << rho_M[48] << endl;
  cout << "rho_minus2point75_minus2point5_M = " << rho_M[49] << endl;
  cout << "rho_minus2point5_minus2point25_M = " << rho_M[50] << endl;
  cout << "rho_minus2point25_minus2_M = " << rho_M[51] << endl;
  cout << "rho_minus2_minus1point75_M = " << rho_M[52] << endl;
  cout << "rho_minus1point75_minus1point5_M = " << rho_M[53] << endl;
  cout << "rho_minus1point5_minus1point25_M = " << rho_M[54] << endl;
  cout << "rho_minus1point25_minus1_M = " << rho_M[55] << endl;
  cout << "rho_minus1_minus0point75_M = " << rho_M[56] << endl;
  cout << "rho_minus0point75_minus0point5_M = " << rho_M[57] << endl;
  cout << "rho_minus0point5_minus0point25_M = " << rho_M[58] << endl;
  cout << "rho_minus0point25_0_M = " << rho_M[59] << endl;
  cout << "rho_0_0point25_M = " << rho_M[60] << endl;
  cout << "rho_0point25_0point5_M = " << rho_M[61] << endl;
  cout << "rho_0point5_0point75_M = " << rho_M[62] << endl;
  cout << "rho_0point75_1_M = " << rho_M[63] << endl;
  cout << "rho_1_1point25_M = " << rho_M[64] << endl;
  cout << "rho_1point25_1point5_M = " << rho_M[65] << endl;
  cout << "rho_1point5_1point75_M = " << rho_M[66] << endl;
  cout << "rho_1point75_2_M = " << rho_M[67] << endl;
  cout << "rho_2_2point25_M = " << rho_M[68] << endl;
  cout << "rho_2point25_2point5_M = " << rho_M[69] << endl;
  cout << "rho_2point5_2point75_M = " << rho_M[70] << endl;

  double rho_f = 0.0;
  for(int i=48;i<71;i++){
    rho_f += rho_M[i];
  }
  cout << "rho_f = " << rho_f << endl;

  
  // The code below calculates and prints out the flow indicators, 
  // as well as some species richness indicators
  // S_t is the number of top predators
  // S_i is number of intermediate species
  // There are also vectors showing flows for a species
  double S_f = 0.0;
  double S_t = 0.0;
  double S_i = 0.0;
  //sequence<double> b_idot;
  //sequence<double> b_dotj;
  sequence<double> b_idotI;
  sequence<double> b_dotjI;
  //sequence<double> b_idotf;
  //sequence<double> b_dotjf;
  sequence<double> b_idotIf;
  sequence<double> b_dotjIf;
  double dummy3 = 0.0;
  double dummy4 = 0.0;
  double dummy5 = 0.0;
  double dummy6 = 0.0;
  double dummy7 = 0.0;
  double dummy8 = 0.0;
  
  // Flows from organisms of lower or equal M to those of higher M
  // Flows from organisms of higher M to those of lower or equal M
  // Then work out b_xM.
  sequence<double> b_lowtohigh;
  sequence<double> b_lowtohigh_Mcorrected;
  sequence<double> b_hightolow;
  sequence<double> b_hightolow_Mcorrected;
  for(int i=0;i<70;i++){
    b_lowtohigh[i]=0;
    b_lowtohigh_Mcorrected[i]=0;
    b_hightolow[i]=0;
    b_hightolow_Mcorrected[i]=0;
  }
  double resolog10M = 0.0;
  double conslog10M = 0.0;
  // Total consumption by fish species
  double Q_f = 0.0;
  // TST values
  //double b_dotdot = 0.0;
  double b_dotdotI = 0.0;
  //double b_dotdotf = 0.0;
  double b_dotdotIf = 0.0;

  for(int i=0;i<S;i++){

    NewSpecies & reso=s[i];
    resolog10M = log10(reso.bodymass()/unit_mass);
    dummy3 = 0;
    dummy4 = 0;
    dummy5 = 0;
    dummy6 = 0;

    for(int j=0;j<S;j++){

      NewSpecies & cons=s[j];
      conslog10M = log10(cons.bodymass()/unit_mass);
      // Update sums for calculation of flows for species
      // Update TST and TST minus respiration
      dummy3 += intake(i,j)/unit_mass/area/(1/year);
      dummy4 += intake(j,i)/unit_mass/area/(1/year);
      //b_dotdot += intake(i,j)/unit_mass/area/(1/year);
      b_dotdotI += intake(i,j)/unit_mass/area/(1/year);
      if(resolog10M > -3){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M > -3)){
	dummy5 += intake(i,j)/unit_mass/area/(1/year);
	dummy6 += intake(j,i)/unit_mass/area/(1/year);
      }
      if((resolog10M <= -3)&&(conslog10M <= -3)){
	dummy7 += intake(i,j)/unit_mass/area/(1/year);
      }

      if((resolog10M <= -14.75)&&(conslog10M > -14.75)){
	b_lowtohigh[0] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -14.5)&&(conslog10M > -14.5)){
	b_lowtohigh[1] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -14.25)&&(conslog10M > -14.25)){
	b_lowtohigh[2] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -14)&&(conslog10M > -14)){
	b_lowtohigh[3] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -13.75)&&(conslog10M > -13.75)){
	b_lowtohigh[4] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -13.5)&&(conslog10M > -13.5)){
	b_lowtohigh[5] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -13.25)&&(conslog10M > -13.25)){
	b_lowtohigh[6] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -13)&&(conslog10M > -13)){
	b_lowtohigh[7] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -12.75)&&(conslog10M > -12.75)){
	b_lowtohigh[8] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -12.5)&&(conslog10M > -12.5)){
	b_lowtohigh[9] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -12.25)&&(conslog10M > -12.25)){
	b_lowtohigh[10] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -12)&&(conslog10M > -12)){
	b_lowtohigh[11] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -11.75)&&(conslog10M > -11.75)){
	b_lowtohigh[12] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -11.5)&&(conslog10M > -11.5)){
	b_lowtohigh[13] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -11.25)&&(conslog10M > -11.25)){
	b_lowtohigh[14] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -11)&&(conslog10M > -11)){
	b_lowtohigh[15] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -10.75)&&(conslog10M > -10.75)){
	b_lowtohigh[16] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -10.5)&&(conslog10M > -10.5)){
	b_lowtohigh[17] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -10.25)&&(conslog10M > -10.25)){
	b_lowtohigh[18] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -10)&&(conslog10M > -10)){
	b_lowtohigh[19] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -9.75)&&(conslog10M > -9.75)){
	b_lowtohigh[20] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -9.5)&&(conslog10M > -9.5)){
	b_lowtohigh[21] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -9.25)&&(conslog10M > -9.25)){
	b_lowtohigh[22] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -9)&&(conslog10M > -9)){
	b_lowtohigh[23] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -8.75)&&(conslog10M > -8.75)){
	b_lowtohigh[24] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -8.5)&&(conslog10M > -8.5)){
	b_lowtohigh[25] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -8.25)&&(conslog10M > -8.25)){
	b_lowtohigh[26] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -8)&&(conslog10M > -8)){
	b_lowtohigh[27] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -7.75)&&(conslog10M > -7.75)){
	b_lowtohigh[28] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -7.5)&&(conslog10M > -7.5)){
	b_lowtohigh[29] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -7.25)&&(conslog10M > -7.25)){
	b_lowtohigh[30] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -7)&&(conslog10M > -7)){
	b_lowtohigh[31] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -6.75)&&(conslog10M > -6.75)){
	b_lowtohigh[32] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -6.5)&&(conslog10M > -6.5)){
	b_lowtohigh[33] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -6.25)&&(conslog10M > -6.25)){
	b_lowtohigh[34] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -6)&&(conslog10M > -6)){
	b_lowtohigh[35] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -5.75)&&(conslog10M > -5.75)){
	b_lowtohigh[36] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[36] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -5.5)&&(conslog10M > -5.5)){
	b_lowtohigh[37] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[37] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -5.25)&&(conslog10M > -5.25)){
	b_lowtohigh[38] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[38] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -5)&&(conslog10M > -5)){
	b_lowtohigh[39] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[39] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -4.75)&&(conslog10M > -4.75)){
	b_lowtohigh[40] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[40] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -4.5)&&(conslog10M > -4.5)){
	b_lowtohigh[41] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[41] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -4.25)&&(conslog10M > -4.25)){
	b_lowtohigh[42] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[42] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -4)&&(conslog10M > -4)){
	b_lowtohigh[43] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[43] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -3.75)&&(conslog10M > -3.75)){
	b_lowtohigh[44] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[44] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -3.5)&&(conslog10M > -3.5)){
	b_lowtohigh[45] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[45] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -3.25)&&(conslog10M > -3.25)){
	b_lowtohigh[46] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[46] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -3)&&(conslog10M > -3)){
	b_lowtohigh[47] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[47] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -2.75)&&(conslog10M > -2.75)){
	b_lowtohigh[48] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[48] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -2.5)&&(conslog10M > -2.5)){
	b_lowtohigh[49] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[49] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -2.25)&&(conslog10M > -2.25)){
	b_lowtohigh[50] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[50] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -2)&&(conslog10M > -2)){
	b_lowtohigh[51] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[51] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -1.75)&&(conslog10M > -1.75)){
	b_lowtohigh[52] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[52] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -1.5)&&(conslog10M > -1.5)){
	b_lowtohigh[53] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[53] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -1.25)&&(conslog10M > -1.25)){
	b_lowtohigh[54] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[54] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -1)&&(conslog10M > -1)){
	b_lowtohigh[55] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[55] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -0.75)&&(conslog10M > -0.75)){
	b_lowtohigh[56] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[56] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -0.5)&&(conslog10M > -0.5)){
	b_lowtohigh[57] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[57] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= -0.25)&&(conslog10M > -0.25)){
	b_lowtohigh[58] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[58] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 0)&&(conslog10M > 0)){
	b_lowtohigh[59] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[59] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 0.25)&&(conslog10M > 0.25)){
	b_lowtohigh[60] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[60] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 0.5)&&(conslog10M > 0.5)){
	b_lowtohigh[61] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[61] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 0.75)&&(conslog10M > 0.75)){
	b_lowtohigh[62] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[62] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 1)&&(conslog10M > 1)){
	b_lowtohigh[63] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[63] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 1.25)&&(conslog10M > 1.25)){
	b_lowtohigh[64] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[64] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 1.5)&&(conslog10M > 1.5)){
	b_lowtohigh[65] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[65] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 1.75)&&(conslog10M > 1.75)){
	b_lowtohigh[66] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[66] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 2)&&(conslog10M > 2)){
	b_lowtohigh[67] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[67] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 2.25)&&(conslog10M > 2.25)){
	b_lowtohigh[68] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[68] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M <= 2.5)&&(conslog10M > 2.5)){
	b_lowtohigh[69] += intake(i,j)/unit_mass/area/(1/year);
	b_lowtohigh_Mcorrected[69] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }

      if((resolog10M > -14.75)&&(conslog10M <= -14.75)){
	b_hightolow[0] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[0] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -14.5)&&(conslog10M <= -14.5)){
	b_hightolow[1] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[1] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -14.25)&&(conslog10M <= -14.25)){
	b_hightolow[2] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[2] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -14)&&(conslog10M <= -14)){
	b_hightolow[3] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[3] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -13.75)&&(conslog10M <= -13.75)){
	b_hightolow[4] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[4] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -13.5)&&(conslog10M <= -13.5)){
	b_hightolow[5] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[5] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -13.25)&&(conslog10M <= -13.25)){
	b_hightolow[6] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[6] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -13)&&(conslog10M <= -13)){
	b_hightolow[7] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[7] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -12.75)&&(conslog10M <= -12.75)){
	b_hightolow[8] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[8] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -12.5)&&(conslog10M <= -12.5)){
	b_hightolow[9] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[9] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -12.25)&&(conslog10M <= -12.25)){
	b_hightolow[10] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[10] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -12)&&(conslog10M <= -12)){
	b_hightolow[11] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[11] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -11.75)&&(conslog10M <= -11.75)){
	b_hightolow[12] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[12] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -11.5)&&(conslog10M <= -11.5)){
	b_hightolow[13] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[13] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -11.25)&&(conslog10M <= -11.25)){
	b_hightolow[14] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[14] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -11)&&(conslog10M <= -11)){
	b_hightolow[15] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[15] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -10.75)&&(conslog10M <= -10.75)){
	b_hightolow[16] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[16] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -10.5)&&(conslog10M <= -10.5)){
	b_hightolow[17] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[17] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -10.25)&&(conslog10M <= -10.25)){
	b_hightolow[18] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[18] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -10)&&(conslog10M <= -10)){
	b_hightolow[19] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[19] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -9.75)&&(conslog10M <= -9.75)){
	b_hightolow[20] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[20] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -9.5)&&(conslog10M <= -9.5)){
	b_hightolow[21] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[21] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -9.25)&&(conslog10M <= -9.25)){
	b_hightolow[22] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[22] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -9)&&(conslog10M <= -9)){
	b_hightolow[23] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[23] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -8.75)&&(conslog10M <= -8.75)){
	b_hightolow[24] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[24] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -8.5)&&(conslog10M <= -8.5)){
	b_hightolow[25] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[25] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -8.25)&&(conslog10M <= -8.25)){
	b_hightolow[26] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[26] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -8)&&(conslog10M <= -8)){
	b_hightolow[27] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[27] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -7.75)&&(conslog10M <= -7.75)){
	b_hightolow[28] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[28] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -7.5)&&(conslog10M <= -7.5)){
	b_hightolow[29] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[29] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -7.25)&&(conslog10M <= -7.25)){
	b_hightolow[30] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[30] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -7)&&(conslog10M <= -7)){
	b_hightolow[31] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[31] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -6.75)&&(conslog10M <= -6.75)){
	b_hightolow[32] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[32] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -6.5)&&(conslog10M <= -6.5)){
	b_hightolow[33] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[33] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -6.25)&&(conslog10M <= -6.25)){
	b_hightolow[34] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[34] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -6)&&(conslog10M <= -6)){
	b_hightolow[35] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[35] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -5.75)&&(conslog10M <= -5.75)){
	b_hightolow[36] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[36] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -5.5)&&(conslog10M <= -5.5)){
	b_hightolow[37] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[37] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -5.25)&&(conslog10M <= -5.25)){
	b_hightolow[38] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[38] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -5)&&(conslog10M <= -5)){
	b_hightolow[39] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[39] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -4.75)&&(conslog10M <= -4.75)){
	b_hightolow[40] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[40] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -4.5)&&(conslog10M <= -4.5)){
	b_hightolow[41] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[41] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -4.25)&&(conslog10M <= -4.25)){
	b_hightolow[42] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[42] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -4)&&(conslog10M <= -4)){
	b_hightolow[43] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[43] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -3.75)&&(conslog10M <= -3.75)){
	b_hightolow[44] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[44] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -3.5)&&(conslog10M <= -3.5)){
	b_hightolow[45] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[45] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -3.25)&&(conslog10M <= -3.25)){
	b_hightolow[46] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[46] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -3)&&(conslog10M <= -3)){
	b_hightolow[47] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[47] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -2.75)&&(conslog10M <= -2.75)){
	b_hightolow[48] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[48] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -2.5)&&(conslog10M <= -2.5)){
	b_hightolow[49] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[49] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -2.25)&&(conslog10M <= -2.25)){
	b_hightolow[50] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[50] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -2)&&(conslog10M <= -2)){
	b_hightolow[51] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[51] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -1.75)&&(conslog10M <= -1.75)){
	b_hightolow[52] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[52] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -1.5)&&(conslog10M <= -1.5)){
	b_hightolow[53] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[53] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -1.25)&&(conslog10M <= -1.25)){
	b_hightolow[54] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[54] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -1)&&(conslog10M <= -1)){
	b_hightolow[55] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[55] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -0.75)&&(conslog10M <= -0.75)){
	b_hightolow[56] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[56] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -0.5)&&(conslog10M <= -0.5)){
	b_hightolow[57] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[57] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > -0.25)&&(conslog10M <= -0.25)){
	b_hightolow[58] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[58] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 0)&&(conslog10M <= 0)){
	b_hightolow[59] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[59] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 0.25)&&(conslog10M <= 0.25)){
	b_hightolow[60] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[60] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 0.5)&&(conslog10M <= 0.5)){
	b_hightolow[61] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[61] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 0.75)&&(conslog10M <= 0.75)){
	b_hightolow[62] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[62] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 1)&&(conslog10M <= 1)){
	b_hightolow[63] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[63] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 1.25)&&(conslog10M <= 1.25)){
	b_hightolow[64] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[64] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 1.5)&&(conslog10M <= 1.5)){
	b_hightolow[65] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[65] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 1.75)&&(conslog10M <= 1.75)){
	b_hightolow[66] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[66] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 2)&&(conslog10M <= 2)){
	b_hightolow[67] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[67] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 2.25)&&(conslog10M <= 2.25)){
	b_hightolow[68] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[68] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }
      if((resolog10M > 2.5)&&(conslog10M <= 2.5)){
	b_hightolow[69] += intake(i,j)/unit_mass/area/(1/year);
	b_hightolow_Mcorrected[69] += ((intake(i,j)/unit_mass/area/(1/year))/pow(cons.bodymass()/unit_mass,-0.25));
      }

      if(conslog10M > -3){
	Q_f += intake(i,j)/unit_mass/area/(1/year);
      }

      // close loop for j  
    }
    
    // Update vectors of flows for species
    b_idotI[i] = dummy3;
    b_dotjI[i] = dummy4;
    //b_dotj[i] = dummy4;
    b_idotIf[i] = dummy5;
    b_dotjIf[i] = dummy6;
    //b_dotjf[i] = dummy6;

    // Update species richness
    if(resolog10M > -3){
      S_f += 1;
    }
    if(b_idotI[i]==0){
      S_t += 1;
    }

    // Update TST and flows with respiration
    /*if(b_dotj[i]==0){
      // is a producer
      b_dotdot += (s[i].the_GP)/unit_mass/area/(1/year);
      b_idot[i] = dummy3 + (s[i].the_GP)/unit_mass/area/(1/year);
    }else{
      // is a consumer
      b_dotdot += 0.259*pow(reso.bodymass()/unit_mass,-0.25)*(biomass_B(i)/unit_mass/area);
      b_idot[i] = dummy3 + 0.259*pow(reso.bodymass()/unit_mass,-0.25)*(biomass_B(i)/unit_mass/area);
      if(resolog10M > -3){ 
	dummy8 += 0.259*pow(reso.bodymass()/unit_mass,-0.25)*(biomass_B(i)/unit_mass/area);
	b_idotf[i] = dummy5 + 0.259*pow(reso.bodymass()/unit_mass,-0.25)*(biomass_B(i)/unit_mass/area);
      }
      } */

    // close loop for i
  }

  b_dotdotIf = b_dotdotI - dummy7;
  //b_dotdotf = b_dotdotIf + dummy8;
  // Update flow vectors for compartment S+1
  /*b_idot[S] = 0;
  b_dotj[S] = b_dotdot - b_dotdotI;
  b_idotf[S] = 0;
  b_dotjf[S] = b_dotdotf - b_dotdotIf;*/

  /*for(int i=0;i<S;i++){
    cout << "b_idot[" << i << "] = " << b_idot[i] << endl;
  }
  for(int i=0;i<S;i++){
    cout << "b_dotj[" << i << "] = " << b_dotj[i] << endl;
    }*/

  cout << "S_f = " << S_f << endl;
  cout << "S_t = " << S_t << endl;
  S_i = S_c-S_t;
  cout << "S_i = " << S_i << endl;

  cout << "b_minus14point75_M = " << b_lowtohigh[0]-b_hightolow[0] << endl;
  cout << "b_minus14point5_M = " << b_lowtohigh[1]-b_hightolow[1] << endl;
  cout << "b_minus14point25_M = " << b_lowtohigh[2]-b_hightolow[2] << endl;
  cout << "b_minus14_M = " << b_lowtohigh[3]-b_hightolow[3] << endl;
  cout << "b_minus13point75_M = " << b_lowtohigh[4]-b_hightolow[4] << endl;
  cout << "b_minus13point5_M = " << b_lowtohigh[5]-b_hightolow[5] << endl;
  cout << "b_minus13point25_M = " << b_lowtohigh[6]-b_hightolow[6] << endl;
  cout << "b_minus13_M = " << b_lowtohigh[7]-b_hightolow[7] << endl;
  cout << "b_minus12point75_M = " << b_lowtohigh[8]-b_hightolow[8] << endl;
  cout << "b_minus12point5_M = " << b_lowtohigh[9]-b_hightolow[9] << endl;
  cout << "b_minus12point25_M = " << b_lowtohigh[10]-b_hightolow[10] << endl;
  cout << "b_minus12_M = " << b_lowtohigh[11]-b_hightolow[11] << endl;
  cout << "b_minus11point75_M = " << b_lowtohigh[12]-b_hightolow[12] << endl;
  cout << "b_minus11point5_M = " << b_lowtohigh[13]-b_hightolow[13] << endl;
  cout << "b_minus11point25_M = " << b_lowtohigh[14]-b_hightolow[14] << endl;
  cout << "b_minus11_M = " << b_lowtohigh[15]-b_hightolow[15] << endl;
  cout << "b_minus10point75_M = " << b_lowtohigh[16]-b_hightolow[16] << endl;
  cout << "b_minus10point5_M = " << b_lowtohigh[17]-b_hightolow[17] << endl;
  cout << "b_minus10point25_M = " << b_lowtohigh[18]-b_hightolow[18] << endl;
  cout << "b_minus10_M = " << b_lowtohigh[19]-b_hightolow[19] << endl;
  cout << "b_minus9point75_M = " << b_lowtohigh[20]-b_hightolow[20] << endl;
  cout << "b_minus9point5_M = " << b_lowtohigh[21]-b_hightolow[21] << endl;
  cout << "b_minus9point25_M = " << b_lowtohigh[22]-b_hightolow[22] << endl;
  cout << "b_minus9_M = " << b_lowtohigh[23]-b_hightolow[23] << endl;
  cout << "b_minus8point75_M = " << b_lowtohigh[24]-b_hightolow[24] << endl;
  cout << "b_minus8point5_M = " << b_lowtohigh[25]-b_hightolow[25] << endl;
  cout << "b_minus8point25_M = " << b_lowtohigh[26]-b_hightolow[26] << endl;
  cout << "b_minus8_M = " << b_lowtohigh[27]-b_hightolow[27] << endl;
  cout << "b_minus7point75_M = " << b_lowtohigh[28]-b_hightolow[28] << endl;
  cout << "b_minus7point5_M = " << b_lowtohigh[29]-b_hightolow[29] << endl;
  cout << "b_minus7point25_M = " << b_lowtohigh[30]-b_hightolow[30] << endl;
  cout << "b_minus7_M = " << b_lowtohigh[31]-b_hightolow[31] << endl;
  cout << "b_minus6point75_M = " << b_lowtohigh[32]-b_hightolow[32] << endl;
  cout << "b_minus6point5_M = " << b_lowtohigh[33]-b_hightolow[33] << endl;
  cout << "b_minus6point25_M = " << b_lowtohigh[34]-b_hightolow[34] << endl;
  cout << "b_minus6_M = " << b_lowtohigh[35]-b_hightolow[35] << endl;
  cout << "b_minus5point75_M = " << b_lowtohigh[36]-b_hightolow[36] << endl;
  cout << "b_minus5point5_M = " << b_lowtohigh[37]-b_hightolow[37] << endl;
  cout << "b_minus5point25_M = " << b_lowtohigh[38]-b_hightolow[38] << endl;
  cout << "b_minus5_M = " << b_lowtohigh[39]-b_hightolow[39] << endl;
  cout << "b_minus4point75_M = " << b_lowtohigh[40]-b_hightolow[40] << endl;
  cout << "b_minus4point5_M = " << b_lowtohigh[41]-b_hightolow[41] << endl;
  cout << "b_minus4point25_M = " << b_lowtohigh[42]-b_hightolow[42] << endl;
  cout << "b_minus4_M = " << b_lowtohigh[43]-b_hightolow[43] << endl;
  cout << "b_minus3point75_M = " << b_lowtohigh[44]-b_hightolow[44] << endl;
  cout << "b_minus3point5_M = " << b_lowtohigh[45]-b_hightolow[45] << endl;
  cout << "b_minus3point25_M = " << b_lowtohigh[46]-b_hightolow[46] << endl;
  cout << "b_minus3_M = " << b_lowtohigh[47]-b_hightolow[47] << endl;
  cout << "b_minus2point75_M = " << b_lowtohigh[48]-b_hightolow[48] << endl;
  cout << "b_minus2point5_M = " << b_lowtohigh[49]-b_hightolow[49] << endl;
  cout << "b_minus2point25_M = " << b_lowtohigh[50]-b_hightolow[50] << endl;
  cout << "b_minus2_M = " << b_lowtohigh[51]-b_hightolow[51] << endl;
  cout << "b_minus1point75_M = " << b_lowtohigh[52]-b_hightolow[52] << endl;
  cout << "b_minus1point5_M = " << b_lowtohigh[53]-b_hightolow[53] << endl;
  cout << "b_minus1point25_M = " << b_lowtohigh[54]-b_hightolow[54] << endl;
  cout << "b_minus1_M = " << b_lowtohigh[55]-b_hightolow[55] << endl;
  cout << "b_minus0point75_M = " << b_lowtohigh[56]-b_hightolow[56] << endl;
  cout << "b_minus0point5_M = " << b_lowtohigh[57]-b_hightolow[57] << endl;
  cout << "b_minus0point25_M = " << b_lowtohigh[58]-b_hightolow[58] << endl;
  cout << "b_0_M = " << b_lowtohigh[59]-b_hightolow[59] << endl;
  cout << "b_0point25_M = " << b_lowtohigh[60]-b_hightolow[60] << endl;
  cout << "b_0point5_M = " << b_lowtohigh[61]-b_hightolow[61] << endl;
  cout << "b_0point75_M = " << b_lowtohigh[62]-b_hightolow[62] << endl;
  cout << "b_1_M = " << b_lowtohigh[63]-b_hightolow[63] << endl;
  cout << "b_1point25_M = " << b_lowtohigh[64]-b_hightolow[64] << endl;
  cout << "b_1point5_M = " << b_lowtohigh[65]-b_hightolow[65] << endl;
  cout << "b_1point75_M = " << b_lowtohigh[66]-b_hightolow[66] << endl;
  cout << "b_2_M = " << b_lowtohigh[67]-b_hightolow[67] << endl;
  cout << "b_2point25_M = " << b_lowtohigh[68]-b_hightolow[68] << endl;
  cout << "b_2point5_M = " << b_lowtohigh[69]-b_hightolow[69] << endl;
  cout << "b_minus14point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[0]-b_hightolow_Mcorrected[0] << endl;
  cout << "b_minus14point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[1]-b_hightolow_Mcorrected[1] << endl;
  cout << "b_minus14point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[2]-b_hightolow_Mcorrected[2] << endl;
  cout << "b_minus14_M_Mcorrected = " << b_lowtohigh_Mcorrected[3]-b_hightolow_Mcorrected[3] << endl;
  cout << "b_minus13point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[4]-b_hightolow_Mcorrected[4] << endl;
  cout << "b_minus13point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[5]-b_hightolow_Mcorrected[5] << endl;
  cout << "b_minus13point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[6]-b_hightolow_Mcorrected[6] << endl;
  cout << "b_minus13_M_Mcorrected = " << b_lowtohigh_Mcorrected[7]-b_hightolow_Mcorrected[7] << endl;
  cout << "b_minus12point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[8]-b_hightolow_Mcorrected[8] << endl;
  cout << "b_minus12point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[9]-b_hightolow_Mcorrected[9] << endl;
  cout << "b_minus12point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[10]-b_hightolow_Mcorrected[10] << endl;
  cout << "b_minus12_M_Mcorrected = " << b_lowtohigh_Mcorrected[11]-b_hightolow_Mcorrected[11] << endl;
  cout << "b_minus11point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[12]-b_hightolow_Mcorrected[12] << endl;
  cout << "b_minus11point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[13]-b_hightolow_Mcorrected[13] << endl;
  cout << "b_minus11point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[14]-b_hightolow_Mcorrected[14] << endl;
  cout << "b_minus11_M_Mcorrected = " << b_lowtohigh_Mcorrected[15]-b_hightolow_Mcorrected[15] << endl;
  cout << "b_minus10point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[16]-b_hightolow_Mcorrected[16] << endl;
  cout << "b_minus10point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[17]-b_hightolow_Mcorrected[17] << endl;
  cout << "b_minus10point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[18]-b_hightolow_Mcorrected[18] << endl;
  cout << "b_minus10_M_Mcorrected = " << b_lowtohigh_Mcorrected[19]-b_hightolow_Mcorrected[19] << endl;
  cout << "b_minus9point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[20]-b_hightolow_Mcorrected[20] << endl;
  cout << "b_minus9point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[21]-b_hightolow_Mcorrected[21] << endl;
  cout << "b_minus9point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[22]-b_hightolow_Mcorrected[22] << endl;
  cout << "b_minus9_M_Mcorrected = " << b_lowtohigh_Mcorrected[23]-b_hightolow_Mcorrected[23] << endl;
  cout << "b_minus8point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[24]-b_hightolow_Mcorrected[24] << endl;
  cout << "b_minus8point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[25]-b_hightolow_Mcorrected[25] << endl;
  cout << "b_minus8point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[26]-b_hightolow_Mcorrected[26] << endl;
  cout << "b_minus8_M_Mcorrected = " << b_lowtohigh_Mcorrected[27]-b_hightolow_Mcorrected[27] << endl;
  cout << "b_minus7point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[28]-b_hightolow_Mcorrected[28] << endl;
  cout << "b_minus7point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[29]-b_hightolow_Mcorrected[29] << endl;
  cout << "b_minus7point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[30]-b_hightolow_Mcorrected[30] << endl;
  cout << "b_minus7_M_Mcorrected = " << b_lowtohigh_Mcorrected[31]-b_hightolow_Mcorrected[31] << endl;
  cout << "b_minus6point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[32]-b_hightolow_Mcorrected[32] << endl;
  cout << "b_minus6point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[33]-b_hightolow_Mcorrected[33] << endl;
  cout << "b_minus6point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[34]-b_hightolow_Mcorrected[34] << endl;
  cout << "b_minus6_M_Mcorrected = " << b_lowtohigh_Mcorrected[35]-b_hightolow_Mcorrected[35] << endl;
  cout << "b_minus5point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[36]-b_hightolow_Mcorrected[36] << endl;
  cout << "b_minus5point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[37]-b_hightolow_Mcorrected[37] << endl;
  cout << "b_minus5point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[38]-b_hightolow_Mcorrected[38] << endl;
  cout << "b_minus5_M_Mcorrected = " << b_lowtohigh_Mcorrected[39]-b_hightolow_Mcorrected[39] << endl;
  cout << "b_minus4point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[40]-b_hightolow_Mcorrected[40] << endl;
  cout << "b_minus4point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[41]-b_hightolow_Mcorrected[41] << endl;
  cout << "b_minus4point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[42]-b_hightolow_Mcorrected[42] << endl;
  cout << "b_minus4_M_Mcorrected = " << b_lowtohigh_Mcorrected[43]-b_hightolow_Mcorrected[43] << endl;
  cout << "b_minus3point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[44]-b_hightolow_Mcorrected[44] << endl;
  cout << "b_minus3point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[45]-b_hightolow_Mcorrected[45] << endl;
  cout << "b_minus3point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[46]-b_hightolow_Mcorrected[46] << endl;
  cout << "b_minus3_M_Mcorrected = " << b_lowtohigh_Mcorrected[47]-b_hightolow_Mcorrected[47] << endl;
  cout << "b_minus2point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[48]-b_hightolow_Mcorrected[48] << endl;
  cout << "b_minus2point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[49]-b_hightolow_Mcorrected[49] << endl;
  cout << "b_minus2point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[50]-b_hightolow_Mcorrected[50] << endl;
  cout << "b_minus2_M_Mcorrected = " << b_lowtohigh_Mcorrected[51]-b_hightolow_Mcorrected[51] << endl;
  cout << "b_minus1point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[52]-b_hightolow_Mcorrected[52] << endl;
  cout << "b_minus1point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[53]-b_hightolow_Mcorrected[53] << endl;
  cout << "b_minus1point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[54]-b_hightolow_Mcorrected[54] << endl;
  cout << "b_minus1_M_Mcorrected = " << b_lowtohigh_Mcorrected[55]-b_hightolow_Mcorrected[55] << endl;
  cout << "b_minus0point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[56]-b_hightolow_Mcorrected[56] << endl;
  cout << "b_minus0point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[57]-b_hightolow_Mcorrected[57] << endl;
  cout << "b_minus0point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[58]-b_hightolow_Mcorrected[58] << endl;
  cout << "b_0_M_Mcorrected = " << b_lowtohigh_Mcorrected[59]-b_hightolow_Mcorrected[59] << endl;
  cout << "b_0point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[60]-b_hightolow_Mcorrected[60] << endl;
  cout << "b_0point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[61]-b_hightolow_Mcorrected[61] << endl;
  cout << "b_0point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[62]-b_hightolow_Mcorrected[62] << endl;
  cout << "b_1_M_Mcorrected = " << b_lowtohigh_Mcorrected[63]-b_hightolow_Mcorrected[63] << endl;
  cout << "b_1point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[64]-b_hightolow_Mcorrected[64] << endl;
  cout << "b_1point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[65]-b_hightolow_Mcorrected[65] << endl;
  cout << "b_1point75_M_Mcorrected = " << b_lowtohigh_Mcorrected[66]-b_hightolow_Mcorrected[66] << endl;
  cout << "b_2_M_Mcorrected = " << b_lowtohigh_Mcorrected[67]-b_hightolow_Mcorrected[67] << endl;
  cout << "b_2point25_M_Mcorrected = " << b_lowtohigh_Mcorrected[68]-b_hightolow_Mcorrected[68] << endl;
  cout << "b_2point5_M_Mcorrected = " << b_lowtohigh_Mcorrected[69]-b_hightolow_Mcorrected[69] << endl;

  cout << "Q_f = " << Q_f << endl;
  //cout << "b_dotdot = " << b_dotdot << endl;
  cout << "b_dotdotI = " << b_dotdotI << endl;
  //cout << "b_dotdotf = " << b_dotdotf << endl;
  cout << "b_dotdotIf = " << b_dotdotIf << endl;


  // More flow indicators
  // There are some vectors too for H_ink, H_outk, n_ink, n_outk

  //double H_btot = 0.0;
  //double Asc = 0.0;
  //double AMI = 0.0;
  //double Cd = 0.0;
  double Asc_I = 0.0;
  //double Asc_E = 0.0;
  double Cd_I = 0.0;
  //double Cd_E = 0.0;
  double H_btotI = 0.0;
  //double H_btotE = 0.0;
  double theta_I = 0.0;
  //double theta_E = 0.0;
  double AMI_I = 0.0;
  double effC = 0.0;
  double effF = 1.0;
  double effN = 0.0;
  double effR = 0.0;
  double H_bout = 0.0;
  double H_bin = 0.0;
  sequence<double> H_ink;
  sequence<double> H_outk;
  sequence<double> n_ink;
  sequence<double> n_outk;
  double LD_q_prime = 0.0;
  double LD_q = 0.0;
  double C_q_prime = 0.0;
  double C_q = 0.0;
  //double H_btotf = 0.0;
  //double Asc_f = 0.0;
  //double Cd_f = 0.0;
  //double theta_f = 0.0;
  double Asc_If = 0.0;
  //double Asc_Ef = 0.0;
  double Cd_If = 0.0;
  //double Cd_Ef = 0.0;
  double H_btotIf = 0.0;
  //double H_btotEf = 0.0;
  double theta_If = 0.0;
  //double theta_Ef = 0.0;
  double AMI_If = 0.0;
  double effC_f = 0.0;
  double effF_f = 1.0;
  double effN_f = 0.0;
  double effR_f = 0.0;
  double H_boutf = 0.0;
  double H_binf = 0.0;
  sequence<double> H_inkf;
  sequence<double> H_outkf;
  sequence<double> n_inkf;
  sequence<double> n_outkf;
  double LD_qf_prime = 0.0;
  double LD_qf = 0.0;
  double C_qf_prime = 0.0;
  double C_qf = 0.0;
  double G_q_prime = 0.0;
  double G_q = 0.0;
  double V_q_prime = 0.0;
  double V_q = 0.0;

  double dummy100 = 0.0;
  double dummy101 = 0.0;
  double dummy102 = 0.0;
  double dummy103 = 0.0;
  double dummy104 = 0.0;
  double dummy105 = 0.0;
  double dummy106 = 0.0;
  double dummy107 = 0.0;
  double dummy108 = 0.0;
  double dummy109 = 0.0;
  double dummy110 = 0.0;
  double dummy111 = 0.0;
  double dummy112 = 0.0;
  double dummy113 = 0.0;
  double dummy114 = 0.0;
  double dummy115 = 0.0;

  for(int i=0;i<S;i++){

    // Reset sums for vectors
    dummy102 = 0;
    dummy103 = 0;
    dummy104 = 0;
    dummy105 = 0;

    for(int j=0;j<S;j++){

      dummy100 = intake(i,j)/unit_mass/area/(1/year);
      dummy101 = intake(j,i)/unit_mass/area/(1/year);

      if(((b_idotI[i]*b_dotjI[j])>0)&&(dummy100>0)){
	Asc_I += dummy100*(log10((dummy100*b_dotdotI)/(b_idotI[i]*b_dotjI[j]))/log10(2));
      }
      if(dummy100>0){
	H_btotI = H_btotI - ((dummy100/b_dotdotI)*
			   (log10(dummy100/b_dotdotI)/log10(2)));
	effF = effF*pow(dummy100/b_dotdotI,-dummy100/b_dotdotI);
      }
      // Updating sums for vectors
      if((dummy101>0)&&(b_dotjI[i]>0)){
	dummy102 = dummy102 - ((dummy101/b_dotjI[i])*
			       (log10(dummy101/b_dotjI[i])/log10(2)));
      }
      if((dummy100>0)&&(b_idotI[i]>0)){
	dummy103 = dummy103 - ((dummy100/b_idotI[i])*
			       (log10(dummy100/b_idotI[i])/log10(2)));
      }
      // Fish indicators
      if(((b_idotIf[i]*b_dotjIf[j])>0)&&(dummy100>0)){
	if((log10(s[i].bodymass()/unit_mass)>-3)&&(log10(s[j].bodymass()/unit_mass)<=-3)){
	  Asc_If += dummy100*(log10((dummy100*b_dotdotIf)/(b_idotIf[i]*b_dotjIf[j]))/log10(2));
	}
	if(log10(s[j].bodymass()/unit_mass)>-3){
	  Asc_If += dummy100*(log10((dummy100*b_dotdotIf)/(b_idotIf[i]*b_dotjIf[j]))/log10(2));
	}
      }
      if(dummy100>0){
	if((log10(s[i].bodymass()/unit_mass)>-3)&&(log10(s[j].bodymass()/unit_mass)<=-3)){
	  H_btotIf = H_btotIf - ((dummy100/b_dotdotIf)*
				 (log10(dummy100/b_dotdotIf)/log10(2)));
	  effF_f = effF_f*pow(dummy100/b_dotdotIf,-dummy100/b_dotdotIf);
	}
	if(log10(s[j].bodymass()/unit_mass)>-3){
	  H_btotIf = H_btotIf - ((dummy100/b_dotdotIf)*
				 (log10(dummy100/b_dotdotIf)/log10(2)));
	  effF_f = effF_f*pow(dummy100/b_dotdotIf,-dummy100/b_dotdotIf);
	}
      }
      // Updating sums for vectors
      if((dummy101>0)&&(b_dotjIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)>-3)){
	  dummy104 = dummy104 - ((dummy101/b_dotjIf[i])*
				 (log10(dummy101/b_dotjIf[i])/log10(2)));
	}
      }
      if((dummy100>0)&&(b_idotIf[i]>0)){
	if(log10(s[i].bodymass()/unit_mass)>-3){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
	if((log10(s[i].bodymass()/unit_mass)<=-3)&&(log10(s[j].bodymass()/unit_mass)<=-3)){
	  dummy105 = dummy105 - ((dummy100/b_idotIf[i])*
				 (log10(dummy100/b_idotIf[i])/log10(2)));
	}
      }
       
      // end loop for j
    }

    H_ink[i] = dummy102;
    H_outk[i] = dummy103;
    if(b_idotI[i]>0){
      H_bout = H_bout - ((b_idotI[i]/b_dotdotI)*
			 (log10(b_idotI[i]/b_dotdotI)/log10(2)));
      n_outk[i] = pow(2,H_outk[i]);
    }else{
      n_outk[i] = 0;
    }
    if(b_dotjI[i]>0){
      H_bin = H_bin - ((b_dotjI[i]/b_dotdotI)*
		       (log10(b_dotjI[i]/b_dotdotI)/log10(2)));
      n_ink[i] = pow(2,H_ink[i]);
    }else{
      n_ink[i] = 0;
    }
    // It has been checked that 1.0/(2.0*S) is a decimal
    LD_q_prime += (1.0/(2.0*S))*(n_outk[i]+n_ink[i]);
    LD_q += (0.5)*(((b_idotI[i]/b_dotdotI)*n_outk[i])+((b_dotjI[i]/b_dotdotI)*n_ink[i]));
    // It has been checked that e.g. 1.0/(S_t+S_i) is a decimal
    G_q_prime += (1.0/(S_t+S_i))*n_ink[i];
    G_q += (b_dotjI[i]/b_dotdotI)*n_ink[i];
    V_q_prime += (1.0/(S_i+S_p))*n_outk[i];
    V_q += (b_idotI[i]/b_dotdotI)*n_outk[i];

    H_inkf[i] = dummy104;
    H_outkf[i] = dummy105;
    if(b_idotIf[i]>0){
      H_boutf = H_boutf - ((b_idotIf[i]/b_dotdotIf)*
			 (log10(b_idotIf[i]/b_dotdotIf)/log10(2)));
      n_outkf[i] = pow(2,H_outkf[i]);
    }else{
      n_outkf[i] = 0;
    }
    if(b_dotjIf[i]>0){
      H_binf = H_binf - ((b_dotjIf[i]/b_dotdotIf)*
		       (log10(b_dotjIf[i]/b_dotdotIf)/log10(2)));
      n_inkf[i] = pow(2,H_inkf[i]);
    }else{
      n_inkf[i] = 0;
    }
    // It has been checked that 1.0/(2.0*S) is a decimal
    LD_qf_prime += (1.0/(2.0*S))*(n_outkf[i]+n_inkf[i]);
    LD_qf += (0.5)*(((b_idotIf[i]/b_dotdotIf)*n_outkf[i])+((b_dotjIf[i]/b_dotdotIf)*n_inkf[i]));

    dummy106 += n_ink[i];
    dummy107 += b_dotjI[i]*n_ink[i];
    dummy108 += n_outk[i];
    dummy109 += b_idotI[i]*n_outk[i];

    // end loop for i
  }

  /*for(int i=0;i<S;i++){
    cout << "H_ink[" << i << "] = " << H_ink[i] << endl;
    cout << "H_outk[" << i << "] = " << H_outk[i] << endl;
    cout << "n_ink[" << i << "] = " << n_ink[i] << endl;
    cout << "n_outk[" << i << "] = " << n_outk[i] << endl;
    }
  for(int i=0;i<S_c;i++){
    cout << "omn_j[" << i << "] = " << omn_j[i] << endl;
  }
  for(int i=0;i<S_f;i++){
    cout << "omn_jf[" << i << "] = " << omn_jf[i] << endl;
    }*/

  Cd_I = b_dotdotI*H_btotI;
  theta_I = Cd_I - Asc_I;
  AMI_I = Asc_I/b_dotdotI;
  effC = pow(2,theta_I/(2.0*b_dotdotI));
  effN = effF/effC;
  effR = pow(2,AMI_I);
  C_q_prime = LD_q_prime/S;
  C_q = LD_q/S;
  
  Cd_If = b_dotdotIf*H_btotIf;
  theta_If = Cd_If - Asc_If;
  AMI_If = Asc_If/b_dotdotIf;
  effC_f = pow(2,theta_If/(2.0*b_dotdotIf));
  effN_f = effF_f/effC_f;
  effR_f = pow(2,AMI_If);
  C_qf_prime = LD_qf_prime/S;
  C_qf = LD_qf/S;

  cout << "Asc_I = " << Asc_I << endl;
  cout << "Cd_I = " << Cd_I << endl;
  cout << "H_btotI = " << H_btotI << endl;
  cout << "theta_I = " << theta_I << endl;
  cout << "AMI_I = " << AMI_I << endl;
  cout << "effC = " << effC << endl;
  cout << "effF = " << effF << endl;
  cout << "effN = " << effN << endl;
  cout << "effR = " << effR << endl;
  cout << "H_bout = " << H_bout << endl;
  cout << "H_bin = " << H_bin << endl;
  cout << "LD_q_prime = " << LD_q_prime << endl;
  cout << "LD_q = " << LD_q << endl;
  cout << "C_q_prime = " << C_q_prime << endl;
  cout << "C_q = " << C_q << endl;

  cout << "Asc_If = " << Asc_If << endl;
  cout << "Cd_If = " << Cd_If << endl;
  cout << "H_btotIf = " << H_btotIf << endl;
  cout << "theta_If = " << theta_If << endl;
  cout << "AMI_If = " << AMI_If << endl;
  cout << "effC_f = " << effC_f << endl;
  cout << "effF_f = " << effF_f << endl;
  cout << "effN_f = " << effN_f << endl;
  cout << "effR_f = " << effR_f << endl;
  cout << "H_boutf = " << H_boutf << endl;
  cout << "H_binf = " << H_binf << endl;
  cout << "LD_qf_prime = " << LD_qf_prime << endl;
  cout << "LD_qf = " << LD_qf << endl;
  cout << "C_qf_prime = " << C_qf_prime << endl;
  cout << "C_qf = " << C_qf << endl;
  
  cout << "G_q_prime = " << G_q_prime << endl;
  cout << "G_q = " << G_q << endl;
  cout << "V_q_prime = " << V_q_prime << endl;
  cout << "V_q = " << V_q << endl;

  // Final set of flow indicators
  // There are some vectors too for g_k_prime, g_k, v_k_prime, v_k, o_j and o_jf

  sequence<double> g_k_prime;
  sequence<double> g_k;
  sequence<double> v_k_prime;
  sequence<double> v_k;
  double sdG_q_prime = 0.0;
  double sdG_q = 0.0;
  double sdV_q_prime = 0.0;
  double sdV_q = 0.0;

  double dummy200 = 0.0;
  double dummy201 = 0.0;
  double dummy202 = 0.0;
  double dummy203 = 0.0;

  for(int i=0;i<S;i++){
    g_k_prime[i] = (S*n_ink[i])/dummy106;
    g_k[i] = (S*b_dotjI[i]*n_ink[i])/dummy107;
    v_k_prime[i] = (S*n_outk[i])/dummy108;
    v_k[i] = (S*b_idotI[i]*n_outk[i])/dummy109;
  }
  // Note that from formulae, the means of above are always 1
  for(int i=0;i<S;i++){
    dummy200 += pow(g_k_prime[i]-1,2.0);
    dummy201 += pow(g_k[i]-1,2.0);
    dummy202 += pow(v_k_prime[i]-1,2.0);
    dummy203 += pow(v_k[i]-1,2.0);  
  }
  // It has been checked that dividing dummy by S gives decimal
  sdG_q_prime = sqrt(dummy200/S);
  sdG_q = sqrt(dummy201/S);
  sdV_q_prime = sqrt(dummy202/S);
  sdV_q = sqrt(dummy203/S);

  cout << "sdG_q_prime = " << sdG_q_prime << endl;
  cout << "sdG_q = " << sdG_q << endl;
  cout << "sdV_q_prime = " << sdV_q_prime << endl;
  cout << "sdV_q = " << sdV_q << endl; 


  // Old link table file no longer produced
  /* for(int n=0;n<S;n++){
    for(int m=0;m<S;m++){
      if(ifrac(m,n)>th){
	NewSpecies & cons=s[n];
	NewSpecies & reso=s[m];
	os.width(9);
	os << level[n] << " ";
	os.width(9);
	os << level[m] << " ";
	os.width(9);
	os << log10(cons.bodymass()/unit_mass) << " ";
	os.width(9);
	os << log10(reso.bodymass()/unit_mass) << " ";
	os.width(9);
	os << log10(biomass_B(n)/unit_mass/area) << " ";
	os.width(9);
	os << log10(biomass_B(m)/unit_mass/area) << " ";
	os.width(9);
	os << assigned_column[n] << " ";
	os.width(9);
	os << assigned_column[m] << " ";
	os.width(9);
	os << intake(m,n)/unit_mass/area/(1/year) << " ";
	os.width(9);
	os << ifrac(m,n) << " ";
 	os.width(9);
 	os << log10(precomputed[n].c[m]) << " ";
	os << std::endl;
      }
    }
    } */

}

int NewWeb::get_invasion_counter(){
  return invasion_counter;
}
int NewWeb::get_speciation_counter(){
  return speciation_counter;
}
int NewWeb::get_max_ever_plant_richness(){
  return max_ever_plant_richness;
}

link_strength_matrix NewWeb::numerical_Jacobian(double dx){
  link_strength_matrix J;
  J.resize(s.size());
  ODE_dynamical_object::numerical_Jacobian(J,dx);
  return J;
}

void NewWeb::eigenvalue_list(const char* filename,double dx){
#if 1
  link_strength_matrix J;
  if(dx){
    J=numerical_Jacobian(dx);
  }else{
    J=numerical_Jacobian();
  }
#else // do eigenvalues of the plant competition matrix instead.
  link_strength_matrix J;
  J.resize(number_of_plants());
  for(int i=number_of_plants();i-->0;){
    for(int j=number_of_plants();j-->0;){
      J[i][j]=-s(i+s.n_animals).plant_hampering_by(s(j+s.n_animals))/
	sqrt(s(i+s.n_animals).plant_hampering_by(s(i+s.n_animals)))/
	sqrt(s(j+s.n_animals).plant_hampering_by(s(j+s.n_animals)));
    }
  }
#endif
  
  
  NewMatrix mat(J);
#if 0 // compute time-scaled community matrix instead
  ublas::diagonal_matrix<double> iB(s.size()),ir(s.size());
  for(int i=s.size();i-->0;)
    iB[i][i]=1/s(i).biomass_abundance_B();
  ublas::diagonal_matrix<double> iB(s.size()),ir(s.size());
  for(int i=s.size();i-->0;){
    iB[i][i]=1/s(i).biomass_abundance_B();
    ir[i][i]=1/(is_plant(i)?
		s(i).plant_growth_rate_sigma() : 
		s(i).turnover_rate_r() );
  }
  mat=ir*mat*iB/year*(animal_biomass()+plant_biomass())/s.size();
#endif

  mat*=1/year; // make sure we use units of 1/year.
  
  eigen_report(mat, filename);

  {
    ostringstream F,FS;
    F << "F" << filename;
    FS << "FS" << filename;

    int S=number_of_species();
    
    vector<int> fish(0);
    for(int i=S;i-->0;){
      if(!is_plant(i) && s(i).bodymass() > fishMlowerthreshold &&
	 s(i).bodymass() < fishMupperthreshold ){
	fish.push_back(i);
      }
    }
    arma::uvec sel(fish.size());
    if(fish.size()){
      copy(fish.begin(),fish.end(),sel.begin());
      NewMatrix J_fish=Schur_complement(mat,sel);
      eigen_report(J_fish,FS.str().c_str());
      eigen_report(mat(sel,sel),F.str().c_str());
    }
  }

  return;
}

sequence<double> NewWeb::get_biomass_B() const{
  sequence<double> b;
  for(int i=s.size();i-->0;){
    b[i]=s[i].biomass_abundance_B();
  }
  return b;
}

sequence<double> NewWeb::get_bodymass_M() const{
  sequence<double> m;
  for(int i=s.size();i-->0;){
    m[i]=s[i].bodymass();
  }
  return m;
}

sequence<double> NewWeb::get_log_bodymass_M() const{
  sequence<double> m;
  for(int i=s.size();i-->0;){
    m[i]=s[i].log_mean_bodymass_M();
  }
  return m;
}

sequence<int> NewWeb::get_is_plant(){
  sequence<int> m;
  for(int i=s.size();i-->0;){
    m[i]=is_plant(i);
  }
  return m;
}

static void 
sample_squares_not_0_dimension(average_meter & am,
			       const NewSpecies::trait_t & t){
  for(int i=t.size();i-->1;){
    double val=t(i);
    am.sample(val*val);
  }
}

void NewWeb::fast_stats(){
  // print some stats that are easily computed
  cout << "n_animals " << number_of_animals() << endl;
  cout << "n_plants " << number_of_plants() << endl;

  double max_log_bodymass=(s.size()?s(0).log_mean_bodymass_M():-9999);
  average_meter log10_plant_bodymass;
  average_meter log10_animal_bodymass;
  average_meter log10_aggressivity;
  double log10=log(10);
  double log10unitmass=log(unit_mass)/log10;
  for(int i=s.size();i-->0;){
    if(is_plant(i)){
      log10_plant_bodymass.sample(s(i).log_mean_bodymass_M()/log10);
    }else{
      log10_animal_bodymass.sample(s(i).log_mean_bodymass_M()/log10);
      log10_aggressivity.sample(log(s(i).aggressivity_g()*unit_mass/unit_area )/log10);
    }
    if(s(i).log_mean_bodymass_M() > max_log_bodymass){
      max_log_bodymass=s(i).log_mean_bodymass_M();
    }
  }
  cout << "max_log10_bodymass " << max_log_bodymass/log10-log10unitmass<< endl;
  cout << "log10_animal_bodymass " << log10_animal_bodymass.readout() -log10unitmass
       << " +/- " << log10_animal_bodymass.std() << endl;
  cout << "log10_plant_bodymass " << log10_plant_bodymass.readout() -log10unitmass
       << " +/- " << log10_plant_bodymass.std() << endl;
  weighted_average_meter log10_bodymass;
  log10_bodymass+=log10_plant_bodymass;
  log10_bodymass+=log10_animal_bodymass;
  cout << "log10_bodymass " << log10_bodymass.readout() -log10unitmass
       << " +/- " << log10_bodymass.sample_std() << endl;
  cout << "log10_aggressivity " 
       << log10_aggressivity.readout()
       << " +/- " << log10_aggressivity.error() 
       << " std " << log10_aggressivity.sample_std() << endl;

  cout << "yield " << biomass_yield()*year << endl;

  if(get_cfg_parameter("totally_random_link_strength_sigma"))
    return;

  average_meter plant_absVrel,animal_absVrel,absV;
  for(int n=0;n<s.size();n++){
    if(is_plant(n)){
      plant_absVrel.
	sample(abs(s[n].relative_vulnerability()));
    }else{
      animal_absVrel.
	sample(abs(s[n].relative_vulnerability()));
    }
    absV.sample(abs(s[n].vulnerability_V()));
  }
  cout << "absV " << absV.readout() << endl;;
  cout << "plant_absVrel " << plant_absVrel.readout() << endl;
  cout << "animal_absVrel " << animal_absVrel.readout() << endl;

  simple_vector<bool> selection(niche_space_dimensions_D,true);
  chi_square_meter VV(selection);
  chi_square_meter VVp(selection);
  chi_square_meter VVa(selection);
  chi_square_meter FFp(selection);
  chi_square_meter FFa(selection);
  average_meter V2,Vp2,Va2,F2,G2;

  for(int n=0;n<s.size();n++){
    VV.sample(s[n].vulnerability_V());
    sample_squares_not_0_dimension(V2,s[n].vulnerability_V());
    if(is_plant(n)){
      VVp.sample(s[n].vulnerability_V());
      sample_squares_not_0_dimension(Vp2,s[n].vulnerability_V());
      FFp.sample(s[n].niche_G());
      sample_squares_not_0_dimension(G2,s[n].niche_G());
    }else{
      VVa.sample(s[n].vulnerability_V());
      sample_squares_not_0_dimension(Va2,s[n].vulnerability_V());
      FFa.sample(s[n].foraging_F());
      sample_squares_not_0_dimension(F2,s[n].foraging_F());
    }
  }
  REPORT(V2);
  REPORT(Vp2);
  REPORT(Va2);
  REPORT(F2);
  REPORT(G2);

  if(get_cfg_parameter("niche_space_dimensions_D") >= number_of_plants())
    return;
  if(get_cfg_parameter("niche_space_dimensions_D") >= number_of_animals())
    return;
  return;
  
  double log_det,log_det_a,log_det_p;
  FFa.chi_square(&log_det,FFa.mean());
  cout << "foraging_trait_volume " << exp(log_det/2) << endl;
  VV.chi_square(&log_det,VV.mean());
  cout << "vulnerability_trait_volume " << exp(log_det/2) << endl;
  VVa.chi_square(&log_det_a,VVa.mean());
  cout << "animal vulnerability_trait_volume " << exp(log_det_a/2) << endl;
  VVp.chi_square(&log_det_p,VVp.mean());
  cout << "plant vulnerability_trait_volume " << exp(log_det_p/2) << endl;
  cout << endl;
  cout << "animal/plant volume ratio " << exp((log_det_a-log_det_p)/2) << endl;
  FFp.chi_square(&log_det,FFp.mean());
  cout << "plant_trait_volume " << exp(log_det/2) << endl;

  return;
}

void NewWeb::distance_histogram(const char * filename){
  ostringstream p,a;
  p << "p" << filename;
  a << "a" << filename;
  Histogram_Estimator hp(0.0,0.02);
  for(int i=s.size();i-->s.n_animals;){
    for(int j=i;j-->s.n_animals;){
      sequence<double> diff=s[i].niche_G()-s[j].niche_G();
      hp.sample(abs(diff));
    }
  }
  Histogram_Estimator ha(0.0,0.2);
  for(int i=s.n_animals;i-->0;){
    for(int j=i;j-->0;){
      sequence<double> diff=s[i].foraging_F()-s[j].foraging_F();
      ha.sample(abs(diff));
    }
  }
  hp.save(p.str().c_str());
  ha.save(a.str().c_str());
}

class distanced_list_helper{
  ofstream file;
  vector<double> values;
  double normalizer;
public:
  distanced_list_helper(const char * filename,double n):
    file(filename),normalizer(n){
  }
  ~distanced_list_helper(){
    sort(values.begin(),values.end());
    for(int i=0;i<values.size();i++){
      if(!(values[i]+1 == values[i])){//< poor man's test for nan/inf
	file  << values[i] << " " << (i+1)*normalizer << endl;
      }
    }
  }
  void operator<<(double x){
    values.push_back(x);
  }
};
    
    
  

void NewWeb::distance_lists(const char * filename){
  ostringstream aV,pV,mV,F,G,FV,aVfish,Ffish;
  aV << "aV" << filename; 
  pV << "pV" << filename; 
  mV << "mV" << filename; 
  F << "F" << filename; 
  G << "G" << filename; 
  FV << "FV" << filename; 
  aVfish << "aVfish" << filename; 
  Ffish << "Ffish" << filename; 

  distanced_list_helper aVs(aV.str().c_str(),2.0/(number_of_animals()));
  distanced_list_helper pVs(pV.str().c_str(),2.0/(number_of_plants()));
  distanced_list_helper mVs(mV.str().c_str(),1);
  distanced_list_helper Fs(F.str().c_str(),2.0/(number_of_animals()));
  distanced_list_helper Gs(G.str().c_str(),2.0/(number_of_plants()));
  distanced_list_helper FVs(FV.str().c_str(),1.0/number_of_animals());
  distanced_list_helper aVfishs(aVfish.str().c_str(),2.0/(number_of_fish()));
  distanced_list_helper Ffishs(Ffish.str().c_str(),2.0/(number_of_fish()));

  sequence<double> diff(niche_space_dimensions_D);
  
  for(int i=s.size();i-->s.n_animals;){
    for(int j=i;j-->s.n_animals;){
      diff=s[i].niche_G()-s[j].niche_G();
      Gs << abs(diff);
      diff=s[i].vulnerability_V()-s[j].vulnerability_V();
      pVs << abs(diff);
    }
    for(int j=s.n_animals;j-->0;){
      diff=s[i].vulnerability_V()-s[j].vulnerability_V();
      mVs << abs(diff);
    }
  }
  for(int i=s.n_animals;i-->0;){
    for(int j=i;j-->0;){
      diff=s[i].foraging_F()-s[j].foraging_F();
      Fs << abs(diff);
      diff=s[i].vulnerability_V()-s[j].vulnerability_V();
      aVs << abs(diff);
    }
    int FV_lower_bound=
      (get_cfg_parameter("no_carnivores") ?
       s.n_animals : 0);
    for(int j=s.size();j-->FV_lower_bound;){
      diff=s[i].foraging_F()-s[j].vulnerability_V();
      FVs << abs(diff);
    }
  }
  for(int i=s.n_animals;i-->0;){
    if((s[i].bodymass()>fishMlowerthreshold)&&(s[i].bodymass()<fishMupperthreshold)){
      for(int j=i;j-->0;){
	if((s[j].bodymass()>fishMlowerthreshold)&&(s[j].bodymass()<fishMupperthreshold)){
	  diff=s[i].foraging_F()-s[j].foraging_F();
	  Ffishs << abs(diff);
	  diff=s[i].vulnerability_V()-s[j].vulnerability_V();
	  aVfishs << abs(diff);
	}
      }
    }
  }
}

     
void NewWeb::compute_biomass_action_products(){
  if(!biomass_action_products.size()){
    biomass_action_products.resize(s.size());
  }
  for(int i=s.n_animals;i-->0;){
    for(int j=s.size();j-->0;){
      // NOT CURRENTLY WORKING:
//       biomass_action_products[i][j]=
//  	biomass_B(i)*pow(biomass_B(j),s[i].switching_exponent_b());
    }
  }
}

multinormal_distribution NewWeb::log_availability_structure(){
  // Compute correlation structure of terms contributing to the
  // logarithm of the availability c_ik B_i.  Include a "temporal
  // niche overlap" factor correcting for temporal correlations
  // between B_i^b_k and B_k as well as the interchange of temporal
  // averaging and potentiation for B_i.

  FATAL_ERROR("log_availability_structure is not adjusted to new switching, yet");
  
  if(!biomass_action_products.size()){
    compute_biomass_action_products();
  }
  
  simple_vector<bool> selection(log_av_n_components,true);
  sequence<double> component(log_av_n_components);

  chi_square_meter term_structure(selection);
  
  for(int i=s.n_animals;i-->0;){
    for(int j=s.size();j-->0;){
      if(s(i).bodymass()>s(j).bodymass()){
	component(log_av_niche)=
	  s(i).log_niche_matching(s(j));
	cerr << component(log_av_niche) << endl;
	component(log_av_small_prey)=
	  s(i).log_small_prey_matching(s(j));
	component(log_av_large_prey)=
	  s(i).log_large_prey_matching(s(j));
	component(log_av_other)=
	  s(i).log_other_matching(s(j));
	component(log_av_biomass)=log(biomass_B(j));
 	// note that this has any effect only after computing time_average:
	// NOT CURRENTLY WORKING:
 	component(log_av_temporal)=0;
//  	  log(biomass_action_products(j,i)/
//  	      (biomass_B(i)*pow(biomass_B(j),s(i).switching_exponent_b())) )/
//  	  s(i).switching_exponent_b();
	term_structure.sample(component);
      }
    }
  }
  return term_structure.estimate();
}
  
void NewWeb::model_environmenatal_change(double strength){
  for(int i=s.size();i-->0;){
    s[i].perturb_physiological_parameters(exp(gaussian(0,strength)));
  }
}
  
double NewWeb::Living_Planet_Index(NewWeb & baseline,
				   bool with_plants, bool with_animals){
  const double regularizer=0.01;// assume 1% of population unseen
  double delta_log_B_sum=0;
  int count=0;
  for(int i=baseline.species_at.size();i-->0;){
    if((baseline.species_at[i]!=column_unused)&&
       (baseline.is_plant(baseline.species_at[i])?with_plants:with_animals)){
      double B2=(species_at[i]==column_unused?
		 0:
		 biomass_B(species_at[i]));
      double B1=baseline.biomass_B(baseline.species_at[i]);
      double mean=(B1+B2)/2;
      delta_log_B_sum+=log(B2+regularizer*mean)-log(B1+regularizer*mean);
      count++;
    }
  }
  return exp(delta_log_B_sum/count);
}

double NewWeb::Biodiversity_Intactness_Index(NewWeb & baseline,
					     bool with_plants, 
					     bool with_animals){
  double I_sum=0;
  int count=0;
  for(int i=baseline.species_at.size();i-->0;){
    if((baseline.species_at[i]!=column_unused)&&
       (baseline.is_plant(baseline.species_at[i])?with_plants:with_animals)){
      if(species_at[i]==column_unused){
	count++;
      }else{
	double B2=biomass_B(species_at[i]);
	double B1=baseline.biomass_B(baseline.species_at[i]);
	I_sum+=B2/B1;
	//cout << B2/B1 << endl;
	count++;
      }
    }
  }
  return I_sum/count;
}

double NewWeb::GPP(){
  double GPP=0;
  for(int i=s.size();i-->s.n_animals;){
    GPP+=s[i].the_GP;
  }
  return GPP;
}

void NewWeb::report_evolutionary_pressures(const int iterations){

  if(steady_state.begin()==steady_state.end()){
    FATAL_ERROR("Must compute steady_state first: do relax()");
  }
  
  typedef NewSpecies::taxon_t taxon_t;
  typedef NewSpecies::trait_selection_t trait_selection_t;
  const int n_taxa=NewSpecies::n_taxa;
  const int S=number_of_species();
  vector< int > total_mutations(n_taxa,0),total_invasions(n_taxa,0);
  bool use_fast_fitness=
    (++steady_state.begin()==steady_state.end());
  if(use_fast_fitness){
    WARNING("Using fast_linear_fitness()");
  }

  // Get an array of size n_taxa X n_traits:
  vector< vector< weighted_average_meter > > pressure(n_taxa);
  vector< average_meter > baseline(n_taxa);
  for(int x=0;x<n_taxa;x++){
    pressure[x].resize(NewSpecies::number_of_kinds_of_traits(taxon_t(x)));
  }

  // Measure pressures
  int nan_found=0;//how many #nan we found
  for(int ii=iterations;ii-->0;){
    cout << ii << " \r"; cout.flush();

    int i=random_integer(S);
    double f_baseline;
    if(use_fast_fitness){
      f_baseline=fast_linear_fitness(&s[i]);
    }else{
      f_baseline=invasion_fitness(&s[i]);
    }
    if(f_baseline+1==f_baseline){
      nan_found++;
    }else{
      double fixed=f_baseline/s[i].turnover_rate_r();
      baseline[s[i].taxon()].sample(fixed*fixed);
    }

    int n_kinds=
      NewSpecies::number_of_kinds_of_traits(s(i).taxon());

    vector< NewSpecies::trait_t > mutating(n_kinds);
    vector< average_meter > mutating_av0(n_kinds);
    vector< NewSpecies::trait_t > invading(n_kinds);
    int mutations=0,invasions=0;

    // Sample traits after mutation/invasion for this species:
    int n=iterations;// invasions to be found.
    while(n>0){
      NewSpecies test_species=s[i];
      double f;
      test_species.mutate_selected_traits(trait_selection_t(-1));
      if(use_fast_fitness){
	f=fast_linear_fitness(&test_species);
      }else{
	f=invasion_fitness(&test_species);
      }
      f-=f_baseline;
      if(f+1==f){
	nan_found++;
      }else{
	mutations++;
	if(f>0){
	  invasions++;
	  n--;
	}
	for(int trait=n_kinds;trait-->0;){
	  NewSpecies::trait_t t=test_species.get_trait_values(trait);
	  mutating[trait]+=t;
	  mutating_av0[trait].sample(t[0]);
	  if(f>0){
	    invading[trait]+=t;
	  }
	}
      }
    }

    total_mutations[s(i).taxon()]+=mutations;
    total_invasions[s(i).taxon()]+=invasions;

    // Compute squared differences between sucessful and mutated trait:
    if(invasions>3){
      for(int trait=n_kinds;trait-->0;){
	mutating[trait]/=mutations;
	invading[trait]/=invasions;
	double estimated_errors_squared=
	  mutating_av0[trait].var()/s[i].trait_dimensions(trait)
	  // The following line takes correlations between the means
	  // "mutating" and "invading" into account.  It assumes
	  // invasion independent of the trait value(!) and takes a
	  // few lines to verify. 
	  *(mutations-invasions)/(double(mutations)*invasions);
	double estimated_relative_distance_square_of_means_per_dimension=
	  (mutating[trait].dist_square(invading[trait])
	   /s[i].trait_dimensions(trait)
	   - estimated_errors_squared
	   ) / mutating_av0[trait].var();
	pressure[s[i].taxon()][trait].
	  sample(estimated_relative_distance_square_of_means_per_dimension,
		 invasions/double(mutations));
      }
    }
  }    
  
  // Report pressures:
  for(int x=0;x<n_taxa;x++){
    NewSpecies dummy(x);
    
    cout << "Taxon: " 
	 << dummy.taxon_name()
	 << endl;
    
    double invasion_success_rate=double(total_invasions[x])/total_mutations[x];
    REPORT(invasion_success_rate);

    cout << "baseline fitness^2: " << baseline[x] << endl;
    
    double sum=0;
    for(int trait=0;trait<NewSpecies::number_of_kinds_of_traits(taxon_t(x));
	++trait){
      cout << dummy.trait_name(trait) << ": "
	   << pressure[x][trait].readout() << " (SD " 
	   << pressure[x][trait].sample_std() << ")" 
	   << " x " << dummy.trait_dimensions(trait) 
	   << endl;
      sum+=pressure[x][trait].readout()*dummy.trait_dimensions(trait);
    }
    REPORT(sum);
  }
  if(nan_found){
    WARNING(nan_found << " NAN found");
  }
}
  
void NewWeb::invasion_fitness_curve(const char * details){
  if(steady_state.begin()==steady_state.end()){
    FATAL_ERROR("Must compute steady_state first: do relax()");
  }
  
  //int focal_species=s.n_animals+random_integer(number_of_plants());
  int focal_species=random_integer(s.n_animals);
  NewSpecies::trait_sweep_t sweep=
    s[focal_species].trait_sweep(details);

  ofstream os("landscape.dat");
  for(NewSpecies::trait_sweep_t::iterator i=sweep.begin();
      i!=sweep.end();
      ++i){
    os << i->trait_values() << " " 
       << invasion_fitness(i->species())/i->species()->turnover_rate_r()
       << endl;
  }
}
void NewWeb::invasion_fitness_relief(const char * details){
  if(steady_state.begin()==steady_state.end()){
    FATAL_ERROR("Must compute steady_state first: do relax()");
  }
  
  //int focal_species=s.n_animals+random_integer(number_of_plants());
  int focal_species=random_integer(s.n_animals);
  NewSpecies::trait_sweep_t sweep=
    s[focal_species].trait_sweep("0");

  ofstream os("relief.dat");
  for(NewSpecies::trait_sweep_t::iterator i=sweep.begin();
      i!=sweep.end();
      ++i){
    cout << i->trait_values() << " ";
    NewSpecies focal2=*(i->species());
    NewSpecies::trait_sweep_t sweep2=
      focal2.trait_sweep("1");
    for(NewSpecies::trait_sweep_t::iterator j=sweep2.begin();
	j!=sweep2.end();
	++j){
      os << invasion_fitness(j->species())/j->species()->turnover_rate_r()
	 << " ";
    }
    os << endl;
  }
  cout << endl;
}

void NewWeb::invasion_fitness_samples(int n){
  ofstream of("fitness_samples.dat");
  average_meter p_invasion;
  if(n>0){
    for(int i=n;i-->0;){
      NewSpecies sample = NewSpecies(s,0,s.n_animals);
      double f = invasion_fitness(&sample)/sample.turnover_rate_r();
      p_invasion.sample(f>0);
      of << f << endl;
    }
  }else if(n < 0){
    for(int i=-n;i-->0;){
      NewSpecies sample = NewSpecies(s,s.n_animals,number_of_plants());
      double f = invasion_fitness(&sample)/sample.turnover_rate_r();
      p_invasion.sample(f>0);
      of << f << endl;
    }
  }
  REPORT(p_invasion);
}

void NewWeb::line_print(ODE_vector const & state, ostream &co){
  double log10=log(10);
  double log10_extinct_biomass=
    (minimum_biomass_considered?
     log(minimum_biomass_considered)/log10:
     -15 );
  for(int i=0;i<species_at.size();i++){
    if(species_at[i] == column_unused){
      co << log10_extinct_biomass;
    }else{
      co << state[species_at[i]]/log10;
    }
    co << " ";
  }
}

void NewWeb::write_link_coefficient_matrix(const char * filename){
  ofstream C(filename);
  C.precision(16);
  C.setf(std::ios::scientific);
  C.setf(std::ios::floatfield);
  bool no_carnivores=get_cfg_parameter("no_carnivores");

  for(int i=(no_carnivores ? s.n_animals : 0); i<s.size(); i++){
    for(int j=0;j<s.n_animals;j++){
      C << s[j].foraging_strength_c_on(s[i]) ;
      if(j+1 < s.n_animals)
	C << ",";
    }
    C << endl;
  }
}


// Helper for write_link_coefficient_matrix:
const string PGM_magic="P5"; 
class smaller_in {
  const NewVector & z;
public:
  smaller_in(NewVector & v):z(v){};
  bool operator()(int a,int b) const{
    return z[a]<z[b];
  }
};

void csv_write(NewMatrix M,const char * filename){
  ofstream os(filename);
  os.precision(16);
  os.setf(std::ios::scientific);
  //os.setf(std::ios::floatfield);

  for(int i=0;i<M.SIZE1();++i){
    for(int j=0;j<M.SIZE2();++j){
      os << M(i,j);
      os << ",";
    }
    os << endl;
  }
}

void critical_vector(const NewVector info, const NewMatrix & M,
		     const char * filename){
  ASSERT(M.SIZE1() == M.SIZE2());
  const int S = M.SIZE1();
  
  arma::cx_mat eigvec;
  arma::cx_vec eigval;
  arma::eig_gen(eigval, eigvec, M);

  int imin=0;
  for(int i=S;i-->0;){
    if(abs(eigvec(imin)) > abs(eigvec(i))){
      imin = i;
    }
  }

  arma::vec critical = arma::real(eigvec.col(imin));
  
  std::ofstream os(filename);
  for(int i=0;i<S;i++){
    os << info[i] << " " << critical[i] << std::endl;
  }
  return;
}

void alpha_image(const NewMatrix alpha, const char *filename)
{
  std::ofstream os(filename);
  
  const int S = alpha.SIZE1();
  
  os << PGM_magic << endl;
  os << alpha.SIZE2() << endl; //width
  os << alpha.SIZE1() << endl; //height
  os << 255 << endl; //maximum graylevel, i.e. black/white image
  
  for(int i=0; i<alpha.SIZE1(); ++i){
    for(int j=0; j<alpha.SIZE2(); ++j){
      int value=(-alpha(i,j)+1)*128;
      value = min(255,max(0,value));
      os.put(value);
    }
  }
}

template <typename T> int sgn(T & val) {
    return (0 < val) - (val < 0);
}

NewMatrix diagonal_one(const NewMatrix & m){
  return arma::diagmat(((0 < m.diag()) - (m.diag() < 0))/sqrt(abs(m.diag()))) *
    m *
    arma::diagmat(1/sqrt(abs(m.diag())));
}

// Compute competitive overlaps. Simple version, also allow for
// competition between plants and animals
void NewWeb::write_competition_matrix(const char * filename){

  int Sp=number_of_plants();
  int Sc=number_of_animals();
  int S=Sp+Sc;
  const double dx = 1e-5;

  // Prepares sorting
  permutation perm(S);
  for(int i=perm.size();i-->0;){
    perm[i]=i;
  }
#if 0
  // ... by trophic level
  // Compute trophic levels for later sorting
  NewVector level;
  NewMatrix F;
  {
    link_strength_matrix intake=intake_matrix();
    link_strength_matrix ifrac=in_fraction(intake);
    F=ifrac;
    NewVector ones(S);
    for(int i=S;i-->0;)
      ones[i]=1;
    level=solve(NewIdentityMatrix(S,S)-F,ones);
  }
  REPORT(level.size());
  sort(perm.begin(),perm.end(),smaller_in(level));
#else
  // ... by size
  NewVector body_size(S);
  for(int i=S;i-->0;){
    body_size[i]=s[i].bodymass();
  }
  REPORT(body_size.size());
  sort(perm.begin(),perm.end(),smaller_in(body_size));
#endif
  

  //// Linearize: generate A (=Aprime) epsAT and Cp matrix and sigma
  //vector
  NewMatrix iB(S,S);
  iB = arma::diagmat(iB);
  for(int i=S;i-->0;){
    iB(i,i)=1/s(i).biomass_abundance_B();
  }

  NewMatrix J=numerical_Jacobian(dx);

  // The following is the correct formula, because J is for
  // _logarithmic_ biomasses.
  NewMatrix G=prod(J,iB);

  
  // Split G according to signs only:
  NewMatrix C=-G;
  NewMatrix A=NewZeroMatrix(S,S);
  NewMatrix epsAT=NewZeroMatrix(S,S);
  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      if(i != j and C(i,j) >= 0 and C(j,i) <= 0){
	if(!(is_plant(i) and is_plant(j))){
	  A(i,j)=C(i,j);C(i,j)=0;
	  epsAT(j,i)=-C(j,i);C(j,i)=0;
	}
      }
    }
  }

#if 0
  csv_write(A,"Aprime.csv");
  csv_write(arma::trans(epsAT),"epsilonA.csv");
  csv_write(C,"C.csv");
#endif
  
  NewMatrix hatC, ihatC;

  // Compute a first guess for the effective competition matrix hatC:
  hatC=C;
  average_meter average_plant_competition;
  for(int i=S;i-->Sc;) average_plant_competition.sample(hatC(i,i));
  hatC.submat(0,0,Sc-1,Sc-1)=
    average_plant_competition.readout() * NewIdentityMatrix(Sc,Sc);
  // Compute the effective competition matrix hatC:

#if 1
  CompMatrixSchur(hatC, ihatC, epsAT, A, C);
#else 
  CompMatrix(hatC, ihatC, epsAT, A, C);
#endif
  
  NewMatrix alpha=NewZeroMatrix(0,0);
  NewMatrix norm(S,S);
  int negative_count=0;
  for(int j=S;j-->0;){
    if(hatC(j,j) < 0) negative_count++;
    double diag_entry=hatC(j,j);
    if(diag_entry==0){
      REPORT(Sc);
      REPORT(j);
      WARNING("zero self-competition, isolating species!!");
      hatC.row(j)*=0;
      hatC.col(j)*=0;
      A.row(j)*=0;
      A.col(j)*=0;
      epsAT.row(j)*=0;
      epsAT.col(j)*=0;
      C.row(j)*=0;
      C.col(j)*=0;
      C(j,j)=1;
      hatC(j,j)=1;
      norm(j,j)=1;
    }else{
      norm(j,j)=sgn(diag_entry)/sqrt(abs(diag_entry));
    }
  }
  norm = arma::diagmat(norm);
    
  REPORT(negative_count);
  WARNING("Getting alpha again");

  NewMatrix alpha_old = alpha;
  alpha = norm * hatC * abs(norm);

#if 0 // also use Schur method and compare
  alpha_image(alpha,"i.pgm");

  CompMatrixSchur(hatC, ihatC, epsAT, A, C);
  
  negative_count=0;
  for(int j=S;j-->0;){
    if(hatC(j,j) < 0) negative_count++;
    double diag_entry=hatC(j,j);
    if(diag_entry==0){
      REPORT(Sc);
      REPORT(j);
      WARNING("zero self-competition, isolating species!!");
      hatC.row(j)*=0;
      hatC.col(j)*=0;
      A.row(j)*=0;
      A.col(j)*=0;
      epsAT.row(j)*=0;
      epsAT.col(j)*=0;
      C.row(j)*=0;
      C.col(j)*=0;
      C(j,j)=1;
      hatC(j,j)=1;
      norm(j,j)=1;
    }else{
      norm(j,j)=sig(diag_entry)/sqrt(abs(diag_entry));
    }
  }
  norm = arma::diagmat(norm);
    
  REPORT(negative_count);
  WARNING("Getting alphaSchur");

  NewMatrix alphaIter = alpha;
  NewMatrix alphaSchur = norm * hatC * abs(norm);
  alpha = alphaSchur;
  
  alpha_image(alphaSchur,"s.pgm");
  alpha_image(arma::join_horiz(alphaIter,alphaSchur),filename);
#endif
  
  NewVector beta(S);
  for(int i=S;i-->0;){
    beta[i]=s[i].biomass_abundance_B()/abs(norm(i,i));
  }
  
  NewVector sigma_eff= prod(alpha, beta);

  {
    std::ofstream os("beta.dat");
    os << beta;
  }
  {
    std::ofstream os("B.dat");
    for(int i=0;i<S;i++){
      os << s[i].biomass_abundance_B() << endl;
    }
  }
  {
    std::ofstream os("sigma_eff.dat");
    os << sigma_eff;
  }
  {
    std::ofstream os("norm.dat");
    for(int i=0; i<S; i++){
      os << norm(i,i) << endl;
    }
  }
    

  // Deal with -1 on diagonal of alpha (needs improving)
  negative_count=0;
  for(int i=S;i-->0;){
    if(alpha(i,i)<0){
      if(isnan(alpha(i,i))){
	FATAL_ERROR("nan on diagonal of alpha: " << i);
      }
      negative_count++;
    }
  }
  if(negative_count){
    WARNING(negative_count << " times -1 on diagonal: dropping");
    //WARNING(negative_count << " times -1 on diagonal");
  }
  
  //negative_count=0;
  REPORT(negative_count);

  NewMatrix unordered_alpha(alpha);
  
  // Clean up an re-order alpha
  if(1){
    NewMatrix alpha_new(S-negative_count,S-negative_count);
    int ii=0;
    for(int i=0;i<S;i++){
      if(alpha(perm[i],perm[i])>0){
	int jj=0;
	for(int j=0;j<S;j++){
	  if(alpha(perm[j],perm[j])>0){
	    alpha_new(ii,jj)=alpha(perm[i],perm[j]);
	    jj++;
	  }
	}
	ii++;
      }
    }
    alpha=alpha_new;
    S-=negative_count;
  }else{
    negative_count=0;
  }
  
#if 0
  // Generate image of alpha
#if 0 // do not include plants
  alpha_image(alpha.submat(Sp,Sp,S-1,S-1),filename);
#else // do include plants
  alpha_image(alpha,filename);
#endif
#endif

    
//   {
//     std::ofstream os("foodweb.pgm");
    
// #if 0 // do not include plants
//     os << PGM_magic << endl;
//     os << Sc << endl; //width
//     os << Sc << endl; //height
//     os << 255 << endl; //maximum graylevel, i.e. black/white image
    
//     for(int i=Sc;i-->0;){
//       for(int j=Sc;j-->0;){
// 	int value=(-F(perm[j],perm[i])+1)*255;
// 	value = min(255,max(0,value));
// 	os.put(value);
//       }
//     }// for i
// #else // do include plants
//     os << PGM_magic << endl;
//     os << S << endl; //width
//     os << S << endl; //height
//     os << 255 << endl; //maximum graylevel, i.e. black/white image
    
//     for(int i=S;i-->0;){
//       for(int j=S;j-->0;){
// 	int value=(-F(perm[j],perm[i])+1)*255;
// 	value = min(255,max(0,value));
// 	os.put(value);
//       }
//     }// for i
// #endif
//   }

  // Compute and save spectrum
  eigen_report(alpha,"alpha_spectrum.dat");
  critical_vector(body_size,unordered_alpha,"alpha_critical.dat");
  // Compute and save animal spectrum
  {
    NewMatrix alphaC(unordered_alpha.submat(0,0,Sc-1,Sc-1));
    eigen_report(alphaC,"alphaC_spectrum.dat");
  }
  {
    NewMatrix alphaCS(diagonal_one(Schur_complement(unordered_alpha,0,Sc-1)));
    eigen_report(alphaCS,"alphaCS_spectrum.dat");
    critical_vector(body_size,alphaCS,"alphaCS_critical.dat");
  }
  {
    vector<int> fish(0);
    for(int i=S;i-->0;){
      if(!is_plant(i) && s(i).bodymass() > fishMlowerthreshold &&
	 s(i).bodymass() < fishMupperthreshold ){
	fish.push_back(i);
      }
    }
    arma::uvec sel(fish.size());
    if(fish.size()){
      copy(fish.begin(),fish.end(),sel.begin());
      NewMatrix alpha_fish=diagonal_one(Schur_complement(unordered_alpha,sel));
      eigen_report(alpha_fish,"alphaFS_spectrum.dat");
      critical_vector(body_size(sel),alpha_fish,"alphaFS_critical.dat");
      eigen_report(unordered_alpha(sel,sel),"alphaF_spectrum.dat");
    }
  }
  // Compute and save spectrum of Raising Operator
  eigen_report(epsAT*ihatC,"R_spectrum.dat");
  // Compute and save spectrum of Lowering Operator
  eigen_report(A*ihatC,"L_spectrum.dat");
}


// Explicitly define constructors in order to be able to specify
// csc_eating as symmetric:
NewWeb::precomputed_entry_t::
precomputed_entry_t():
  csc_eating(SortedMatrix::symmetric),
  csc_being_eaten(SortedMatrix::asymmetric),
  c(),
  cT(){
}


void NewWeb::initialize_precomputed(){
  precomputed.resize(0,s);
  precomputed.resize(s.size(),s);

  link_strength_matrix similarity;
  if(not adjust_prey_similarity_width){
    // Pre-pre compute similarities in order to precompute faster:
    similarity.resize(s.size());
    if(do_switching){
      for(int i=s.size();i-->0;){
	similarity[i][i]=
	  s(i).
	  switching_similarity_to(s(i),switching_similarity_width_w_s);
	for(int j=i;j-->0;){
	  similarity[i][j]=similarity[j][i]=
	    s(i).
	    switching_similarity_to(s(j),switching_similarity_width_w_s);
	}
      }
    }
  }
  // Now, do the sorted vectors first.  We cannot speed up much here:
  SortedVector dummy_vector;
  double vector_truncation_epsilon=dummy_vector.truncation_epsilon();
  WARN_IF(SortedVector().truncation_epsilon()==0,"Better set this to non-zero");
  for(int k=s.n_animals;k-->0;){
    for(int i=s.size();i-->0;){
      double c=s[k].foraging_strength_c_on(s[i]);
      precomputed[k].c[i]=c;
      if(!do_switching){
	precomputed[i].cT[k]=c;
      }
    }
  }
  for(int i=s.size();i-->s.n_animals;){
    for(int j=s.size();j-->s.n_animals;){
      precomputed[i].c[j]=s[i].plant_hampering_by(s[j]);
    }
  }

  /////// Test how many prey per predator are covered
  // for(int i=number_of_animals();i-->0;){
  //   REPORT(precomputed[i].c.size());
  // }
  // exit(1);

  if(!do_switching)
    return;

  // Now, do the sorted matrices, and do them smart:
  SortedMatrix dummy_matrix;
  double matrix_truncation_epsilon=dummy_matrix.truncation_epsilon();
  WARN_IF(SortedMatrix().truncation_epsilon()==0,"Better set this to non-zero");
  vector< int > prey(s.size()); // list of prey of a predator
  for(int k=s.n_animals;k-->0;){
    int Z=0;
    for(int i=s.size();i-->0;){
      if(precomputed[k].c[i] > matrix_truncation_epsilon){
	prey[Z++]=i;
      }
    }
    for(int ii=Z;ii-->0;){
      int i=prey[ii];
      double pre_k_i=precomputed[k].c[i];
      if(adjust_prey_similarity_width){
	for(int jj=ii+1;jj-->0;){
	  int j=prey[jj];
	  double entry=pre_k_i*
	    precomputed[k].c[j]*
	    s(i).switching_similarity_to(s(j),s(k).prey_similarity_width());
	  if(entry > matrix_truncation_epsilon){
	    precomputed[k].csc_eating[i][j]=entry;
	    precomputed[i].csc_being_eaten[j][k]=entry;
	    precomputed[j].csc_being_eaten[i][k]=entry;
	  }
	}
      }else{
	for(int jj=ii+1;jj-->0;){
	  int j=prey[jj];
	  double entry=pre_k_i*
	    precomputed[k].c[j]*
	    similarity[i][j];
	  if(entry > matrix_truncation_epsilon){
	    precomputed[k].csc_eating[i][j]=entry;
	    precomputed[i].csc_being_eaten[j][k]=entry;
	    precomputed[j].csc_being_eaten[i][k]=entry;
	  }
	}
      }
    }
  }
}


NewWeb::precomputed_t & 
NewWeb::precomputed_t::operator=(const NewWeb::precomputed_t & other){
  if (this != &other){
    // first, sloppy resize
    // instead of resize(size), do the following:
    while(size() < other.size()){
      push_back(new precomputed_entry_t());
    }
    while(size() > other.size()){
      delete pop_back();
    }

    // now, copy elementwise;
    copy(other.begin(),other.end(),begin());
  }
  return *this;
}

NewWeb::precomputed_t::precomputed_t(const precomputed_t & other){
  for(precomputed_container_t::const_iterator i=other.begin();
      i!=other.end();
      ++i){
    push_back(new precomputed_entry_t(*i));
  }
}

NewWeb::precomputed_t::~precomputed_t(){
  while(size() > 0){
    delete pop_back();
  }
}

void NewWeb::precomputed_t::resize(int size,species_list_t &s){
  // for convenience:
  precomputed_container_t & pre=*this;

  // instead of pre.resize(size), do the following:
  ASSERT(size >= 0);
  while(pre.size() < size){
    pre.push_back(new precomputed_entry_t());
  }
  while(pre.size() > size){
    delete pre.pop_back();
  }

  for(int i=size;i-->0;){
    precomputed_entry_t & pre_i=pre[i];
    pre_i.c.resize(size);
    if(do_switching){
      pre_i.csc_eating.resize(size);
      pre_i.csc_being_eaten.resize(size);
    }else{
      pre_i.cT.resize(size);
    }
  }
}
 
void NewWeb::precomputed_t::update(int at,species_list_t &s){
  ASSERT(at<s.size());

  // for convenience:
  precomputed_container_t & pre=*this;
  precomputed_entry_t & pre_at=pre[at];
  NewSpecies & s_at=s[at];

  pre_at.c.clear();
  pre_at.cT.clear();
  pre_at.csc_eating.clear();
  pre_at.csc_being_eaten.clear();

  double matrix_truncation_epsilon=SortedMatrix().truncation_epsilon();
  vector< int > prey(s.size()); // list of prey of at
  vector< int > predator(s.size()); // list of prey of at
  int GEN=0;
  int VUL=0;

  // Now, do the sorted vectors first.  We remember prey and predators:
  double vector_truncation_epsilon=SortedVector().truncation_epsilon();
  if(at<s.n_animals){
    for(int i=s.size();i-->0;){
      double c=
	pre_at.c[i]=s_at.foraging_strength_c_on(s[i]);

      // if(i>=s.n_animals){
      // 	for(int k=100000;k-->0;){
      // 	  NewSpecies p=NewSpecies(NewSpecies::plant);
      // 	  NewSpecies a=NewSpecies(NewSpecies::animal);
      // 	  p.set_unique_id();
      // 	  a.set_unique_id();
      // 	  //double c=a.foraging_strength_c_on(p);
      // 	  double lnc=a.log_niche_matching(p);
      // 	  REPORT(c);
      // 	  static average_meter lna;
      // 	  lna.sample(lnc);
      // 	  REPORT(lna);
      // 	  REPORT(lna.std());
      // 	}
      // }


      if(do_switching){
	if(c > matrix_truncation_epsilon ){
	  prey[GEN++]=i;
	}
      }else{
	pre[i].cT[at]=c;
      }
    }
  }else{
    for(int i=s.size();i-->s.n_animals;){
      pre[i].c[at]=s[i].plant_hampering_by(s_at);
      pre_at.c[i]=s_at.plant_hampering_by(s[i]);
    }
  }
  for(int k=s.n_animals;k-->0;){
    double c=
      pre[k].c[at]=s[k].foraging_strength_c_on(s_at);
    if(do_switching){
      if(c > matrix_truncation_epsilon ){
	predator[VUL++]=k;
      }
    }else{
      pre_at.cT[k]=c;
    }
  }

  if(!do_switching)
    return;

  // Now, do the sorted matrices, and do them smart:
  for(int ii=GEN;ii-->0;){
    int i=prey[ii];
    double pre_i=pre_at.c[i];
    for(int jj=ii+1;jj-->0;){
      int j=prey[jj];
      double entry=pre_i*
	pre_at.c[j]*
	s[i].switching_similarity_to(s[j],s_at.prey_similarity_width());
      if(entry > matrix_truncation_epsilon){
	pre_at.csc_eating[i][j]=entry;
	pre[i].csc_being_eaten[j][at]=entry;
	pre[j].csc_being_eaten[i][at]=entry;
      }
    }
  }
  for(int kk=VUL;kk-->0;){
    int k=predator[kk];
    double pre_k_at=pre[k].c[at];
    for(int i=s.size();i-->0;){
      double pre_k_i=pre[k].c[i];
      if(pre_k_i > matrix_truncation_epsilon){
	double entry=
	  pre_k_i*pre_k_at*
	  s[i].switching_similarity_to(s_at,s(k).prey_similarity_width());
	if(entry > matrix_truncation_epsilon){
	  pre[k].csc_eating[i][at]=entry;
	  pre[i].csc_being_eaten[at][k]=entry;
	  pre_at.csc_being_eaten[i][k]=entry;
	}
      }
    }
  }
}
 
void NewWeb::precomputed_t::move(int from,int to,species_list_t &s){
  // We can now be sure that move deletes all value assigned to the
  // "from" index, so we do not need to set value to zero in
  // update(...) anymore.

  // for convenience:
  precomputed_container_t & pre=*this;
  precomputed_entry_t & pre_from=pre[from];

  pre[to]=pre_from;
  pre_from.c.clear();
  pre_from.cT.clear();
  pre_from.csc_eating.clear();
  pre_from.csc_being_eaten.clear();

  for(int i=s.size();i-->0;){
    precomputed_entry_t & pre_i=pre[i];
    pre_i.c.move(from,to);
    pre_i.cT.move(from,to);
    pre_i.csc_eating.move(from,to);
    pre_i.csc_being_eaten.move(from,to);
  }
}

NewWeb::saved_state::saved_state():t(0){};
NewWeb::saved_state::
saved_state(double tt,const vector_with_max & B,const vector_with_max &cf):
  t(tt),biomass_B(B),common_factor(cf){};
NewWeb::saved_state::
saved_state(const saved_state & s):
  t(s.t),biomass_B(s.biomass_B),common_factor(s.common_factor){};

int number_of_CPUs(){
  int number_of_cpus;
#ifdef _SC_NPROCESSORS_ONLN
  number_of_cpus = sysconf(_SC_NPROCESSORS_ONLN); 
#elif defined(HW_NCPU)
  {
    int mib[2];
    size_t len;

    mib[0] = CTL_HW;
    mib[1] = HW_NCPU;
    len = sizeof(number_of_cpus);
    sysctl(mib, 2, &number_of_cpus, &len, NULL, 0);
  } 
#else
  WARNING("Could not detect number of CPUs.");
  number_of_cpus = 1;
#endif
  REPORT(number_of_cpus);
  return number_of_cpus;
}

#ifdef SET_CPU_AFFINITY
/* DEAD CODE (only needed with dead code below)
// kernel < 2.6.19 does not have getcpu, so we need to work around:
#ifndef KERNEL_VERSION
#define KERNEL_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))
#endif
#ifndef LINUX_VERSION_CODE
#define LINUX_VERSION_CODE KERNEL_VERSION(0,0,0)
#endif
#if LINUX_VERSION_CODE < KERNEL_VERSION(2,6,19)
#include <sys/syscall.h>
int getcpu(int * cpu_id, void*, void*){
  pid_t pid = syscall(SYS_getpid);
  std::stringstream stringbuffer;
  stringbuffer << "cat /proc/" << pid << "/stat | awk '{print $39}'";
  FILE* file = popen(stringbuffer.str().c_str(), "r");
  if (fscanf(file, "%d", cpu_id) == EOF){
    FATAL_ERROR("Could not find current cpu");
  }
  pclose(file);
  return 0;
}
#endif
 END OF DEAD CODE */
 
CPU_set_list_t find_my_CPUs(){
  int n_cpus=number_of_CPUs();
  vector<int> cpu_list(n_cpus);

  bool CPUs_assigned=false;

  /// First, use heuristics to find a good list of CPUs:
  if(getenv("PBS_VERSION") and
     strlen(getenv("PBS_VERSION"))>=6 and
     strncmp(getenv("PBS_VERSION"),"TORQUE",6) == 0 and
     getenv("PBS_JOBID")){
    // Looks like we can find the assigned CPU using a torque command:
    FILE * pipe=popen("qstat -f -1 ${PBS_JOBID}|grep exec_host","r");
    ALWAYS_ASSERT(pipe!=0);
    int i=0;
    while((!feof(pipe)) and !ferror(pipe) and i<n_cpus){
      char c=getc(pipe);
      if(c == '/'){
	fscanf(pipe,"%i",&cpu_list[i]);
	++i;
      }
      ALWAYS_ASSERT(ferror(pipe)==0);
    }
    cpu_list.resize(i);

    CPUs_assigned=true;
  }
  /* DEAD CODE: COULD NOT MAKE THIS WORK 
   else{
    // Find free CPUs manually:
    cpu_set_t all_cpus;
    for(int i = 0;i<n_cpus;++i){//try all cpus
      CPU_SET(i, &all_cpus);
    }

    cpu_set_t try_these_cpus;
    cpu_set_t all_cpus;

    for(int i = 0;i<n_cpus;++i){//try all cpus
      CPU_SET(i, &all_cpus);
      CPU_SET(i, &try_these_cpus);
    }
    // Now successively disallow those that we are using to prioritize those
    // that are available:
    int current_cpu;
    int i=0;
    while(true){
      getcpu(&current_cpu,0,0);
      cpu_list[i++]=current_cpu;
      if(!CPU_ISSET(current_cpu, &try_these_cpus)){
	FATAL_ERROR("Running on forbidden CPU");
      };
      if(i>=n_cpus) break;
      CPU_CLR(current_cpu, &try_these_cpus);
      if(sched_setaffinity(0, CPU_SETSIZE, &try_these_cpus)){
	REPORT(EFAULT);
	REPORT(EINVAL);
	REPORT(EPERM);
	REPORT(ESRCH);
	FATAL_ERROR("Problem setting CPU affinity (errno=" << errno <<")");
      };
    }
    for(int i=0;i<n_cpus-1;i++){
      getcpu(&current_cpu,0,0);
      cpu_list[i]=current_cpu;
      if(!CPU_ISSET(current_cpu, &try_these_cpus)){
	FATAL_ERROR("Running on forbidden CPU");
      }
      CPU_CLR(current_cpu,&try_these_cpus);
      if(sched_setaffinity(0, CPU_SETSIZE, &try_these_cpus)){
	FATAL_ERROR("Problem setting CPU affinity (errno=" << errno <<")");
      }
    }
    if(sched_setaffinity(0, CPU_SETSIZE, &all_cpus)){
      FATAL_ERROR("Problem setting CPU affinity (errno=" << errno <<")");
    }
  }
END DEAD CODE */

  max_num_threads=min((int)cpu_list.size(),max_num_threads);

  if(CPUs_assigned){
    cout << "Assigned CPUs: ";
    for(int i=0;i<cpu_list.size()-1;i++){
      cout << cpu_list[i] << ", ";
    }
    cout << *cpu_list.rbegin() << endl;
  }

  // We need this if we could not assign cpus.
  cpu_set_t all_cpus;
  for(int i = 0;i<n_cpus;++i){//try all cpus
    CPU_SET(i, &all_cpus);
  }

  /// cpu_list is now an array filled with different cpu ids.  Now we
  /// have to convert these into (trivial) cpu sets.
  CPU_set_list_t cpus(cpu_list.size());
    
  for(int i = 0;i<cpus.size();++i){
    CPU_ZERO(&cpus[i]); // clear cpu set
    if(CPUs_assigned){
      CPU_SET(cpu_list[i], &cpus[i]);
    }else{
      memcpy(&cpus[i],&all_cpus,sizeof(cpu_set_t));
    }
  }

  REPORT(CPUs_assigned);
  return cpus;
}
#endif

