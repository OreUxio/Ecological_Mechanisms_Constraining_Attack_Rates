// $Id: NewSpecies.cc 2462 2016-01-23 08:08:45Z axel $

#include "NewSpecies.h"
#include "random.h"
#include "evaluate.h"
#include "error.h"
#include "Statistics.h"
#include <float.h>
// #ifdef DEBUGGING
// #include <CLHEP/Matrix/Matrix.h> // for TRACING
// #endif
#include <fstream>

using namespace std;

static my_evaluator_t eval_first;


double NewSpecies::plant_fraction=0.5;
static double Mmax=0;//eval_first("1 * kilogram");//unused
static double Mopt=eval_first("1 * gram");
static double Mini=eval_first("10 * milligram");
static double Mmin=eval_first("1 * milligram");
static double body_mass_cutoff_exponent=2;
static int plant_size_evolves=true;

int niche_space_dimensions_D=5;
static double initial_number_of_individuals=1.0;
double minimum_biomass_considered=0;
static double typical_bodymass_ratio_d=1.1; ///< exp(size diffusion constant)
static double mu_V=0; ///< deprecate
static double plant_mu_V_boost=0; ///< deprecate
static double mu_V_plant=0.1;
static double mu_V_animal=0.1;
static double mu_F=0.1;
static double mu_G=0.1;
double predator_relative_size_Lambda=1; ///< tolerance to size inversion
double trophic_niche_width_w_t=0;
double trophic_niche_width_spread=0;
static double plant_niche_width_w_r=0;
static double plant_competition_drop_off_exponent=4;
static double plant_diffuse_competition=0;
static double plant_diffuse_competition_probability=1;
static double plant_diffuse_competition_symmetry=0;
static int unit_diagonal_plant_competition=0;
static int experimental=0;
int do_switching=1;
double switching_similarity_width_w_s=0;
int adjust_prey_similarity_width=0;
double relative_switching_similarity_width=0;
static double universal_switching_similarity=0;///< used when totally_random_link_strength_sigma is set
static double typical_relative_error_in_link_strength=0;
static double initial_aggressivity=1;
static double typical_aggressivity_ratio_da=1;
static double aggressivity_decay=1;
static int use_aggressivity=1;
static double carnivore_aggressivity_boost=1;
static double typical_niche_width_ratio=0;
static double niche_width_decay=0;
static double plant_vulnerability_separation=0;
static int vulnerability_separate_in_all_dimensions=0;
static double plant_animal_axis_offset=0;
static double plant_vulnerability_relative_weight=1;
double plant_shadowing=0;
double plant_hardening_time=0;
int big_plants_eaten=1;
int small_plants_win=false;
static double plant_speedup=1;
static int biomass_balancing=0; ///< deprecate, does nothing
static double foraging_variability_f0=1;
static int flat_foraging_distribution=0;
static int flat_vulnerability_distribution=0;
static int flat_plant_distribution=0;
static double plant_growth_rate=0; //zero means plant_growth_rate is variable
int emulate_carrying_capacity_bug=0;
int emulate_animal_body_mass_cutoff_bug=0;
int emulate_log_random_matching_bug=0;
static int no_carnivores=0;
static int no_plant_competition=0;
static int animal_physiology_version=1;
int plant_physiology_version=1;
static int trait_evolution_version=1;
static double plant_vulnerability_variability=1;
static double animal_vulnerability_variability=1;
static double plant_trait_variability=1;
static int trait_space_model=1;
static int use_Mopt=1;
static double log10_bodymass_ratio_window_center=4;
static double log10_bodymass_ratio_window_sigma=2;
static double cutoff_prey_predator_mass_ratio=0;
static double totally_random_link_strength_sigma=0;

/// Determines link-strength distributions through theoretical formulae
/** With value set to 2, also modify evolution of aggressivity accordingly.
 **/
static int controlled_random_link_strengths=0;
static double randomly_sampled_foraging_traits=0;
static int use_Mmax_fish=0;
static int further_change_trade_off_multiplier=1;
static double large_plant_trade_off_exponent=0;

// ALLOMETRY:

// Make use of allometric relations here: Robert Henry Peters: "The
// ecological implications of body size", McGill University, (1983)
static double allometric_base_unit=eval_first("1*kilogram");

// Peters (1993), p. 33:
static double energy_content_of_wet_biomass_q = eval_first("7e6*joule/kilogram");
// Peters (1983), p. 29:
static double metabolic_rate_allometric_exponent = 0.751;
// //homeotherms:
// double metabolic_rate_allometric_prefactor = eval_first("4.1*watt");
//poikilotherms:
static double metabolic_rate_allometric_prefactor = eval_first("0.14*watt");
// //unicellular:
// double metabolic_rate_allometric_prefactor = eval_first("0.018*watt");
// //average over all these:
// double metabolic_rate_allometric_exponent = eval_first("(3-(-12))/(3-(-13))");
// double metabolic_rate_allometric_prefactor =
// eval_first("10^(-12-(-13)*(3-(-12))/(3-(-13)))");

static double conversion_efficiency_epsilon=0.2;
// literature values, Peters (1983), p. 143: 0.02 (homeotherms) to 0.2
// (poikilotherms)

// from KOOTEN et al., Ecology, 85(7), 2004, pp. 1979-1991: !! I got
// this wrong, exponent should be 0.62-1, but left is here as is for
// backward compatibility of code.
static double
attack_rate_allometric_exponent = 0.62;
static double
attack_rate_allometric_prefactor = 1;

static double
small_prey_exponent_alpha = 0;

static double
big_prey_exponent_beta = 0.25;

// from Jonathan M. Jeschke & Ralph Tollrian, Ethology 111, 187-206 (2005)
static double
handling_time_allometric_prefactor = eval_first("21.49*minutes/gram*kilogram");
static double
handling_time_allometric_exponent = eval_first("-0.849+1");

// from Niklas & Enquist PNAS 98(5), 2922 (2001):
static double
growth_rate_allometric_prefactor=eval_first("0.208/year");
static double
growth_rate_allometric_exponent=eval_first("0.763-1");

static double growth_rate_mean_random_deviation=0.0;
static double growth_rate_std_random_deviation=0.0;


// from Tang & Peters Journal of Plankton Research 17(2), 303-315 (1995):
static double
plant_loss_rate_allometric_prefactor=eval_first("15.1/year");
static double
plant_loss_rate_allometric_exponent=eval_first("-0.07");

// guess:
double loss_rate_over_max_production_rate_r=0.2;

static double
body_size_allometric_prefactor=eval_first("2.58*meter");
static double
body_size_allometric_exponent=eval_first("0.264");


// Peters (1993) p. 170 mentiones a total net primary production of "1
// to 50 tons dry mass ha^-1 y^-1:
double area_per_compartment=eval_first("2*ha");
static double dominating_species_primary_production=
area_per_compartment*eval_first("50*tons/year/ha");

//Savage 2004:
static double max_growth_rate_allometric_prefactor=
  eval_first("2.6*10^(-8)/second");
static double max_growth_rate_allometric_exponent=
  eval_first("-1/4");
static double max_growth_rate_allometric_prefactor_zoo=
  eval_first("7.39*10^(-9)/second");
static double max_growth_rate_allometric_exponent_zoo=
  eval_first("-1/4");
// Schimtz and Lavigne (1984):
static double max_growth_rate_allometric_prefactor_marmam=
  eval_first("2.19*10^(-8)/second");
//Clarke and Johnston 1999:
static double respiration_rate_allometric_prefactor=
  eval_first("1.59*10^(-8)/second");

static double respiration_rate_mean_random_deviation=0.0;
static double respiration_rate_std_random_deviation=0.0;

static double respiration_rate_allometric_exponent=
  eval_first("-1/4");
//Lopez-Urrutia et al. 2006:
static double respiration_rate_allometric_prefactor_zoo=
  eval_first("1.23*10^(-8)/second");
static double respiration_rate_allometric_exponent_zoo=
  eval_first("-0.13");
//Lavigne et al. 1986
static double respiration_rate_allometric_prefactor_marmam=
  eval_first("4.76*10^(-7)/second");
//Peters:
static double production_over_respiration=0.8;
static double Mmax_fish=345*eval_first("1*kilogram");
static double marmam_pred_press_prefactor=1;
static double marmam_pred_press_exponent=1;
static double trade_off_multiplier_Mmin_prefactor=1;


// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] =
{
  {"plant_fraction", CFG_DOUBLE, &(NewSpecies::plant_fraction)},
  CFGDOUBLE(Mmin),
  CFGDOUBLE(Mopt),
  CFGDOUBLE(Mini),
  CFGDOUBLE(Mmax),
  CFGDOUBLE(body_mass_cutoff_exponent),
  CFGINT(plant_size_evolves),
  CFGINT(niche_space_dimensions_D),
  CFGDOUBLE(initial_number_of_individuals),
  CFGDOUBLE(minimum_biomass_considered),
  CFGDOUBLE(typical_bodymass_ratio_d),
  CFGDOUBLE(mu_V),
  CFGDOUBLE(mu_V_plant),
  CFGDOUBLE(mu_V_animal),
  CFGDOUBLE(mu_F),
  CFGDOUBLE(mu_G),
  CFGDOUBLE(predator_relative_size_Lambda),
  CFGDOUBLE(trophic_niche_width_w_t),
  CFGDOUBLE(trophic_niche_width_spread),
  CFGDOUBLE(plant_niche_width_w_r),
  CFGDOUBLE(plant_competition_drop_off_exponent),
  CFGDOUBLE(plant_diffuse_competition),
  CFGDOUBLE(plant_diffuse_competition_probability),
  CFGDOUBLE(plant_diffuse_competition_symmetry),
  CFGINT(unit_diagonal_plant_competition),
  CFGINT(experimental),
  CFGINT(do_switching),
  CFGDOUBLE(switching_similarity_width_w_s),
  CFGDOUBLE(relative_switching_similarity_width),
  CFGDOUBLE(universal_switching_similarity),
  CFGINT(adjust_prey_similarity_width),
  CFGDOUBLE(plant_vulnerability_separation),
  CFGINT(vulnerability_separate_in_all_dimensions),
  CFGDOUBLE(plant_animal_axis_offset),
  CFGDOUBLE(plant_vulnerability_relative_weight),
  CFGDOUBLE(plant_shadowing),
  CFGINT(big_plants_eaten),
  CFGINT(small_plants_win),
  CFGDOUBLE(plant_speedup),
  CFGDOUBLE(plant_growth_rate), //zero means plant_growth_rate is variable
  CFGDOUBLE(plant_hardening_time),
  CFGDOUBLE(foraging_variability_f0),
  CFGINT(flat_foraging_distribution),
  CFGINT(flat_vulnerability_distribution),
  CFGINT(flat_plant_distribution),
  CFGDOUBLE(allometric_base_unit),
  CFGDOUBLE(energy_content_of_wet_biomass_q),
  CFGDOUBLE(metabolic_rate_allometric_exponent),
  CFGDOUBLE(metabolic_rate_allometric_prefactor),
  CFGDOUBLE(conversion_efficiency_epsilon),
  CFGDOUBLE(attack_rate_allometric_exponent),
  CFGDOUBLE(attack_rate_allometric_prefactor),
  CFGDOUBLE(typical_aggressivity_ratio_da),
  CFGDOUBLE(typical_niche_width_ratio),
  CFGDOUBLE(niche_width_decay),
  CFGDOUBLE(totally_random_link_strength_sigma),
  CFGINT(controlled_random_link_strengths),
  CFGDOUBLE(randomly_sampled_foraging_traits),
  CFGDOUBLE(aggressivity_decay),
  CFGINT(use_aggressivity),
  CFGDOUBLE(carnivore_aggressivity_boost),
  CFGDOUBLE(handling_time_allometric_prefactor),
  CFGDOUBLE(handling_time_allometric_exponent),
  CFGDOUBLE(growth_rate_allometric_prefactor),
  CFGDOUBLE(growth_rate_allometric_exponent),
  CFGDOUBLE(growth_rate_mean_random_deviation),
  CFGDOUBLE(growth_rate_std_random_deviation),
  CFGDOUBLE(plant_loss_rate_allometric_prefactor),
  CFGDOUBLE(plant_loss_rate_allometric_exponent),
  CFGDOUBLE(loss_rate_over_max_production_rate_r),
  CFGDOUBLE(body_size_allometric_prefactor),
  CFGDOUBLE(body_size_allometric_exponent),
  CFGINT(biomass_balancing),
  CFGDOUBLE(area_per_compartment),
  CFGDOUBLE(dominating_species_primary_production),
  CFGDOUBLE(small_prey_exponent_alpha),
  CFGDOUBLE(big_prey_exponent_beta),
  CFGINT(emulate_carrying_capacity_bug),
  CFGINT(emulate_animal_body_mass_cutoff_bug),
  CFGINT(emulate_log_random_matching_bug),
  CFGINT(no_carnivores),
  CFGINT(no_plant_competition),
  CFGINT(animal_physiology_version),
  CFGINT(plant_physiology_version),
  CFGDOUBLE(max_growth_rate_allometric_prefactor_zoo),
  CFGDOUBLE(max_growth_rate_allometric_exponent_zoo),
  CFGDOUBLE(max_growth_rate_allometric_prefactor_marmam),
  CFGDOUBLE(max_growth_rate_allometric_prefactor),
  CFGDOUBLE(max_growth_rate_allometric_exponent),
  CFGDOUBLE(respiration_rate_allometric_prefactor_zoo),
  CFGDOUBLE(respiration_rate_allometric_exponent_zoo),
  CFGDOUBLE(respiration_rate_allometric_prefactor_marmam),
  CFGDOUBLE(respiration_rate_allometric_prefactor),
  CFGDOUBLE(respiration_rate_allometric_exponent),
  CFGDOUBLE(respiration_rate_mean_random_deviation),
  CFGDOUBLE(respiration_rate_std_random_deviation),
  CFGDOUBLE(production_over_respiration),
  CFGINT(trait_evolution_version),
  CFGDOUBLE(plant_vulnerability_variability),
  CFGDOUBLE(plant_mu_V_boost),
  CFGDOUBLE(animal_vulnerability_variability),
  CFGDOUBLE(plant_trait_variability),
  CFGINT(trait_space_model),
  CFGINT(use_Mopt),
  CFGDOUBLE(initial_aggressivity),
  CFGDOUBLE(log10_bodymass_ratio_window_center),
  CFGDOUBLE(log10_bodymass_ratio_window_sigma),
  CFGDOUBLE(cutoff_prey_predator_mass_ratio),
  CFGDOUBLE(Mmax_fish),
  CFGDOUBLE(marmam_pred_press_prefactor),
  CFGDOUBLE(marmam_pred_press_exponent),
  CFGDOUBLE(trade_off_multiplier_Mmin_prefactor),
  CFGINT(use_Mmax_fish),
  CFGINT(further_change_trade_off_multiplier),
  CFGDOUBLE(large_plant_trade_off_exponent),
  CFGDOUBLE(typical_relative_error_in_link_strength),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

double NewSpecies::log_DBL_MAX=5*log(DBL_MAX)/7.0; // we added some safety margin
double NewSpecies::log_DBL_MIN=5*log(DBL_MIN)/7.0; // we added some safety margin

int NewSpecies::next_unique_id=0;
int NewSpecies::next_sequential_id=0;

void NewSpecies::data_mapping(){
  if(the_remember->read_or_write()==rememberRead){
    int is_plant;
    if(PERMANENT(is_plant)==0){
      if(is_plant)
	the_taxon=plant;
    }
  }
  PERMANENT(the_taxon);
  PERMANENT(the_aggressivity_g);
  if(typical_niche_width_ratio){
    if(PERMANENT(the_niche_width_wt)>0){
      the_niche_width_wt=trophic_niche_width_w_t;
    }
  }
  PERMANENT(the_mean_bodymass_M);
  for(int i=0;i<niche_space_dimensions_D;i++){
#ifndef HAS_TYPEOF
    double & VComp=V[i]; //guessing type
#else
    typeof(V[0]) & VComp=V[i];
#endif
    if(PERMANENT(VComp)>0) VComp=0;
  }
  for(int i=0;i<niche_space_dimensions_D;i++){
#ifndef HAS_TYPEOF
    double & FComp=F[i];//guessing type
#else
    typeof(F[0]) & FComp=F[i];
#endif
    if(PERMANENT(FComp)>0) FComp=0;
  }
  if(the_remember->read_or_write()==rememberRead)
    set_parameters_given_body_mass();
  PERMANENT(the_biomass_abundance_B);
  PERMANENT(the_unique_id);
  PERMANENT(the_parent_unique_id);
  PERMANENT(the_clade_id);
  if(next_unique_id<=the_unique_id)
    next_unique_id=the_unique_id+1;
  if(PERMANENT(the_sequential_id)>0){
    the_sequential_id=the_unique_id;
  }
  next_sequential_id=max<int>(the_sequential_id+1,next_sequential_id);
  PERMANENT(the_random_tag);
  PERMANENT(the_fishing_mortality);
  PERMANENT(the_number_of_invading_offspring);
}

void NewSpecies::adjust_trait_vector_dimensions(){
  V.resize(niche_space_dimensions_D);
  F.resize(niche_space_dimensions_D);
}


void NewSpecies::set_unique_id(){
  the_unique_id=next_unique_id++;
  TRACE(the_unique_id,SPECIES);
  the_random_tag=random_int();
}

void NewSpecies::set_sequential_id(){
  the_sequential_id=next_sequential_id++;
}

double NewSpecies::get_trade_off_multiplier()const{

  // If trade_off_multiplier>=1, the loss rate is so large that
  // fitness is <=0, even at optimal conditions.
  double trade_off_multiplier=0;
  switch(trait_space_model){
  case 1:
    // don't use trade_off_multiplier
    break;
  case 2:
    {
      const double delta_log_M=
	log(the_mean_bodymass_M/Mopt)/
	log(Mmin/Mopt);
      trade_off_multiplier=pow(fabs(delta_log_M),body_mass_cutoff_exponent);
    }
    break;
  case 3:
    {
      if(use_Mopt){
	const double delta_log_M=
	  log(the_mean_bodymass_M/Mopt)/
	  log(Mmin/Mopt);
	trade_off_multiplier=
	  pow(fabs(delta_log_M),body_mass_cutoff_exponent);
 	if(emulate_animal_body_mass_cutoff_bug and not is_a_plant()){
	  trade_off_multiplier=
	    delta_log_M*delta_log_M;
	}
      }else{
	// if use_Mmax_fish = 1, then upper body mass limit imposed
	if(use_Mmax_fish){
	  if((the_mean_bodymass_M>Mmin)&&(the_mean_bodymass_M<Mmax_fish)){
	    trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
	  } else{
	  trade_off_multiplier=1e100;
	  }
	  // else use_Mmax_fish = 0
	} else{
	  if(the_mean_bodymass_M>Mmin){
	    trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
	  } else{
	    trade_off_multiplier=1e100;
	  }
	}
      }
      // only run below if further_change_trade_off_multiplier = 1
      if(further_change_trade_off_multiplier){
	if(is_a_plant() and not no_plant_competition){
	  trade_off_multiplier+=
	    dot(relative_vulnerability(),relative_vulnerability())/
	    (plant_vulnerability_variability*plant_vulnerability_variability)/
	    niche_space_dimensions_D;
	  trade_off_multiplier+=
	    dot(niche_G(),niche_G())/
	    (plant_trait_variability*plant_trait_variability)/
	    niche_space_dimensions_D;
	  trade_off_multiplier=
	    2*((1.0/3.0)*trade_off_multiplier-0.5);
	}else{//it's an animal or non-competing plant
	  const double variability=
	    ( is_a_plant() ?
	      plant_vulnerability_variability :
	      animal_vulnerability_variability );
	  trade_off_multiplier+=
	    dot(relative_vulnerability(),relative_vulnerability())/
	    (variability*variability)/
	    niche_space_dimensions_D;
	  trade_off_multiplier=
	    2*((1.0/2.0)*trade_off_multiplier-0.5);
	}
      }
    }
    break;
  case 4:
    {
      // Trade off only at small M to represent greater respiration costs due to loss of production apparatus
      // This multiplier approaches infinity at Mmin and is set to a very large number for M<=Mmin
      if(the_mean_bodymass_M>Mmin){
	trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
      } else{
	trade_off_multiplier=1e100;
      }
    }
    case 5:
    {
      // Same as case 4 except that for consumers, trade_off_multiplier is increased further by a term
      // representing predation by marine mammals (increasing function of body mass).
      if(is_a_plant()){
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
	} else{
	  trade_off_multiplier=1e100;
	}
      }else{//it's an animal
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=(1/((the_mean_bodymass_M/Mmin)-1))+((the_mean_bodymass_M/Mmax_fish)*marmam_pred_press_prefactor);
	} else{
	  trade_off_multiplier=1e100;
	}
      }
    }
    case 6:
    {
      // Same as case 5 but with a different term representing marine mammal predation.
      if(is_a_plant()){
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
	} else{
	  trade_off_multiplier=1e100;
	}
      }else{//it's an animal
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=(1/((the_mean_bodymass_M/Mmin)-1))
	    +(pow(log(the_mean_bodymass_M/Mmin)/log(Mmax_fish/Mmin),marmam_pred_press_exponent)*marmam_pred_press_prefactor);
	} else{
	  trade_off_multiplier=1e100;
	}
      }
    }
    case 7:
    {
      // Same as case 4 except that for consumers, trade_off_multiplier is increased to a very
      // high value if body mass is above maximum limit for fish.
      if(is_a_plant()){
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=1/((the_mean_bodymass_M/Mmin)-1);
	} else{
	  trade_off_multiplier=1e100;
	}
      }else{//it's an animal
	if((the_mean_bodymass_M>Mmin)&&(the_mean_bodymass_M<Mmax_fish || !use_Mmax_fish)){
	  trade_off_multiplier=(1/((the_mean_bodymass_M/Mmin)-1));
	} else{
	  trade_off_multiplier=1e100;
	}
      }
    }
    break;
  case 8:
    {
      // Same as case 4 except that there is an extra factor in trade_off_multiplier
      // controlling the increase as M approache Mmin
      if(the_mean_bodymass_M>Mmin){
	trade_off_multiplier=1/(trade_off_multiplier_Mmin_prefactor*((the_mean_bodymass_M/Mmin)-1));
      } else{
	trade_off_multiplier=1e100;
      }
    }
    break;
  case 9:
    {
      // Same as case 7 except that there is an extra factor in trade_off_multiplier
      // controlling the increase as M approache Mmin
      if(is_a_plant()){
	if(the_mean_bodymass_M>Mmin){
	  trade_off_multiplier=1/(trade_off_multiplier_Mmin_prefactor*((the_mean_bodymass_M/Mmin)-1));
	} else{
	  trade_off_multiplier=1e100;
	}
      }else{//it's an animal
	if((the_mean_bodymass_M>Mmin)&&
	   (the_mean_bodymass_M<Mmax_fish || !use_Mmax_fish) ){
	  trade_off_multiplier=1/(trade_off_multiplier_Mmin_prefactor*((the_mean_bodymass_M/Mmin)-1));
	} else{
	  trade_off_multiplier=1e100;
	}
      }
    }
    break;
  default:
    FATAL_ERROR("trait_space_model " << trait_space_model << " unknown");
  }

  if(is_a_plant()){
    trade_off_multiplier+=
      large_plant_trade_off_exponent*log(the_mean_bodymass_M/Mmin);
  }

  return trade_off_multiplier;
}

void NewSpecies::set_parameters_given_body_mass(){
  TRACE(the_mean_bodymass_M,SPECIES);

  the_log_mean_bodymass_M=
    log(the_mean_bodymass_M);
  double M=
    the_mean_bodymass_M/allometric_base_unit;


  double trade_off_multiplier=
    get_trade_off_multiplier();

  if(controlled_random_link_strengths){
    if(!use_aggressivity==5){
      FATAL_ERROR("controlled_random_link_strengths requires use_aggressivity==5");
    }
    if(!totally_random_link_strength_sigma){
      FATAL_ERROR("controlled_random_link_strengths requires totally_random_link_strength_sigma");
    }
  }

  switch(animal_physiology_version){
  case 1:
    {
      the_handling_time_T=handling_time_allometric_prefactor*
	pow(M,handling_time_allometric_exponent);
      double metabolic_rate_eta=
	metabolic_rate_allometric_prefactor*
	pow(M,metabolic_rate_allometric_exponent);
      double turnover_time=
	the_mean_bodymass_M*
	energy_content_of_wet_biomass_q/
	metabolic_rate_eta;
      //  cout << "turnover_time: " << turnover_time << std::endl;
      the_turnover_rate_r= 1/turnover_time;
      double raw_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon; // preliminary
      // e x/T_h - r x = raw_e x /T_h
      the_conversion_efficiency_epsilon =
	raw_conversion_efficiency_epsilon +
	the_handling_time_T*the_turnover_rate_r;
      ALWAYS_ASSERT(the_conversion_efficiency_epsilon < 1);
    }
    break;
  case 2:
    {
      double max_growth_rate=
	max_growth_rate_allometric_prefactor*
	pow(M,max_growth_rate_allometric_exponent);
      the_handling_time_T=
	conversion_efficiency_epsilon/
	max_growth_rate;
      the_turnover_rate_r=
	max_growth_rate/(production_over_respiration);
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon+
	the_turnover_rate_r*the_handling_time_T;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
//       REPORT(conversion_efficiency_epsilon);
//       REPORT(the_turnover_rate_r);
//       REPORT(max_growth_rate);
//       REPORT(1/the_handling_time_T);
//       REPORT(the_conversion_efficiency_epsilon);
//       exit(1);
    }
    break;
  case 3:
    {
      the_turnover_rate_r=
	respiration_rate_allometric_prefactor*
	pow(M,respiration_rate_allometric_exponent);
      double max_growth_rate=
	(production_over_respiration)*the_turnover_rate_r;
      if(max_growth_rate < 0){
	the_handling_time_T=0;//simple Type I response
      }else{
	the_handling_time_T=
	  conversion_efficiency_epsilon/
	  max_growth_rate;
      }
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon+
	the_turnover_rate_r*the_handling_time_T;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
//       REPORT(conversion_efficiency_epsilon);
//       REPORT(the_turnover_rate_r);
//       REPORT(max_growth_rate);
//       REPORT(1/the_handling_time_T);
//       REPORT(the_conversion_efficiency_epsilon);
//       exit(1);
    }
    break;
  case 4:
    // if M <= 10^(-6) kg, then both the max_growth_rate and the the_turnover_rate are
    // specified by allometric scaling laws for zooplankton
    {
      double max_growth_rate;
      if(M<=10^(-6)){
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor_zoo*
	  pow(M,respiration_rate_allometric_exponent_zoo);
	max_growth_rate=
	  max_growth_rate_allometric_prefactor_zoo*
	  pow(M,max_growth_rate_allometric_exponent_zoo);
      }else{
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor*
	  pow(M,respiration_rate_allometric_exponent);
	max_growth_rate=
	  (production_over_respiration)*the_turnover_rate_r;
      }
      the_handling_time_T=
	conversion_efficiency_epsilon/
	max_growth_rate;
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon+
	the_turnover_rate_r*the_handling_time_T;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
      //       REPORT(conversion_efficiency_epsilon);
      //       REPORT(the_turnover_rate_r);
      //       REPORT(max_growth_rate);
      //       REPORT(1/the_handling_time_T);
      //       REPORT(the_conversion_efficiency_epsilon);
      //       exit(1);
    }
    break;
  /*case 5:
    // same as case 3 except that conversion_efficiency_epsilon is now
    // the_conversion_efficiency_epsilon, since direct values for this
    // can be found from Hendriks (2007)
    {
      the_turnover_rate_r=
	respiration_rate_allometric_prefactor*
	pow(M,respiration_rate_allometric_exponent);
      double max_growth_rate=
	(production_over_respiration)*the_turnover_rate_r;
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon;
      double conversion_efficiency_epsilon_zero=
	the_conversion_efficiency_epsilon/
	(1+(the_turnover_rate_r/max_growth_rate));
      the_handling_time_T=
	conversion_efficiency_epsilon_zero/
	max_growth_rate;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon_zero);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
      //       REPORT(conversion_efficiency_epsilon);
      //       REPORT(the_turnover_rate_r);
      //       REPORT(max_growth_rate);
      //       REPORT(1/the_handling_time_T);
      //       REPORT(the_conversion_efficiency_epsilon);
      //       exit(1);
    }
    break;*/

  case 5:
    // same as case 3 except that conversion_efficiency_epsilon is now
    // the_conversion_efficiency_epsilon, since direct values for this
    // can be found from Hendriks (2007)
    // with optional random variation introduced in the_turnover_rate
    {
      the_turnover_rate_r=
	respiration_rate_allometric_prefactor*
        (1+respiration_rate_mean_random_deviation+respiration_rate_std_random_deviation*
        pseudo_standard_normal_random_with(*this,0))*
	pow(M,respiration_rate_allometric_exponent);
      double max_growth_rate=
	(production_over_respiration)*the_turnover_rate_r;
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon;
      double conversion_efficiency_epsilon_zero=
	the_conversion_efficiency_epsilon/
	(1+(the_turnover_rate_r/max_growth_rate));
      the_handling_time_T=
	conversion_efficiency_epsilon_zero/
	max_growth_rate;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon_zero);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
      //       REPORT(conversion_efficiency_epsilon);
      //       REPORT(the_turnover_rate_r);
      //       REPORT(max_growth_rate);
      //       REPORT(1/the_handling_time_T);
      //       REPORT(the_conversion_efficiency_epsilon);
      //       exit(1);
    }
    break;

  case 6:
    // same as case 4 except that conversion_efficiency_epsilon is now
    // the_conversion_efficiency_epsilon, since direct values for this
    // can be found from Hendriks (2007)
    {
      double max_growth_rate;
      if(M<=10^(-6)){
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor_zoo*
	  pow(M,respiration_rate_allometric_exponent_zoo);
	max_growth_rate=
	  max_growth_rate_allometric_prefactor_zoo*
	  pow(M,max_growth_rate_allometric_exponent_zoo);
      }else{
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor*
	  pow(M,respiration_rate_allometric_exponent);
	max_growth_rate=
	  (production_over_respiration)*the_turnover_rate_r;
      }
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon;
      double conversion_efficiency_epsilon_zero=
	the_conversion_efficiency_epsilon/
	(1+(the_turnover_rate_r/max_growth_rate));
      the_handling_time_T=
	conversion_efficiency_epsilon_zero/
	max_growth_rate;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon_zero);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
//       REPORT(conversion_efficiency_epsilon);
//       REPORT(the_turnover_rate_r);
//       REPORT(max_growth_rate);
//       REPORT(1/the_handling_time_T);
//       REPORT(the_conversion_efficiency_epsilon);
//       exit(1);
    }
    break;
  case 7:
    // if M > Mmax_fish, then both the max_growth_rate and the the_turnover_rate are
    // specified by allometric scaling laws for marine mammals
    {
      double max_growth_rate;
      if(M>(Mmax_fish/allometric_base_unit)){
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor_marmam*
	  pow(M,respiration_rate_allometric_exponent);
	max_growth_rate=
	  max_growth_rate_allometric_prefactor_marmam*
	  pow(M,max_growth_rate_allometric_exponent);
      }else{
	the_turnover_rate_r=
	  respiration_rate_allometric_prefactor*
	  pow(M,respiration_rate_allometric_exponent);
	max_growth_rate=
	  (production_over_respiration)*the_turnover_rate_r;
      }
      the_conversion_efficiency_epsilon=
	conversion_efficiency_epsilon;
      double conversion_efficiency_epsilon_zero=
	the_conversion_efficiency_epsilon/
	(1+(the_turnover_rate_r/max_growth_rate));
      the_handling_time_T=
	conversion_efficiency_epsilon_zero/
	max_growth_rate;
      if(the_conversion_efficiency_epsilon>1 && !no_carnivores){
	REPORT(conversion_efficiency_epsilon_zero);
	REPORT(the_turnover_rate_r);
	REPORT(max_growth_rate);
	REPORT(1/the_handling_time_T);
	REPORT(the_conversion_efficiency_epsilon);
	FATAL_ERROR("conversion_efficiency_epsilon is too large!");
      }
//       REPORT(conversion_efficiency_epsilon);
//       REPORT(the_turnover_rate_r);
//       REPORT(max_growth_rate);
//       REPORT(1/the_handling_time_T);
//       REPORT(the_conversion_efficiency_epsilon);
//       exit(1);
    }
    break;
  default:
    FATAL_ERROR("this animal_physiology_version does not exist yet");
  }
  // attack rate now defined
  switch(use_aggressivity){
  case 0:
   {
     the_attack_rate_a=attack_rate_allometric_prefactor*
       pow(M,attack_rate_allometric_exponent);
   }
   break;
  case 1: // attack rate scaled with respiration rate
   {
     the_attack_rate_a=the_aggressivity_g*the_turnover_rate_r/area_per_compartment;
   }
   break;
  case 2: // attack rate NOT scaled with respiration rate
   {
     the_attack_rate_a=the_aggressivity_g/area_per_compartment;
   }
   break;
  case 3: // attack rate scaled with respiration rate BUT with attack_rate_allometric_exponent as exponent
   {
     the_attack_rate_a=
       (the_aggressivity_g*(respiration_rate_allometric_prefactor*pow(M,attack_rate_allometric_exponent)))/
       area_per_compartment;
   }
   break;
  case 4: // attack rate just rescaled with aggressivity
   {
     the_attack_rate_a=
       attack_rate_allometric_prefactor*
       pow(M,attack_rate_allometric_exponent)*
       the_aggressivity_g;
   }
   break;
  case 5: // for controlled_random_link_strengths
   {
     double as_ratio=-log(the_aggressivity_g);
     double target_S=sqrt(2*M_PI)*as_ratio*exp(as_ratio*as_ratio/2);
     double sigma=sqrt(2*log(target_S))/the_niche_width_wt;
     const double assumed_plant_growth_rate=1;
     the_attack_rate_a=
       the_turnover_rate_r*assumed_plant_growth_rate/
       (dominating_species_primary_production*
	the_conversion_efficiency_epsilon)*
       exp(-sigma*as_ratio);
   }
   break;
  default:
    FATAL_ERROR("this use_aggressivity does not exist yet");
  }

  if(production_over_respiration > 0){
    // implement trade-offs:
    const double rmax=
      the_conversion_efficiency_epsilon/the_handling_time_T-
      the_turnover_rate_r;
    if(rmax > the_turnover_rate_r){
	static bool warned=false;
	if(not warned){
	  REPORT(the_turnover_rate_r);
	  REPORT(rmax);
	  WARNING("the_turnover_rate_r<rmax");
	  warned=true;
	}
	if(trade_off_multiplier<0){
	  the_turnover_rate_r+=the_turnover_rate_r*trade_off_multiplier;
	}else{
	  the_turnover_rate_r+=rmax*trade_off_multiplier;
	}
      }else{
	the_turnover_rate_r+=rmax*trade_off_multiplier;
    }
  }else{
    the_turnover_rate_r+=the_turnover_rate_r*trade_off_multiplier;
  }

  if(is_a_plant()){
    switch(plant_physiology_version){
    case 1: // (the old mess)
      the_plant_growth_rate_sigma=growth_rate_allometric_prefactor*
	pow(M,growth_rate_allometric_exponent);
      if(!emulate_carrying_capacity_bug){
	the_carrying_capacity=
	  dominating_species_primary_production/
	  the_plant_growth_rate_sigma;
      }else{
	// old buggy version
	double metabolic_rate_eta=
	  metabolic_rate_allometric_prefactor*
	  pow(M,metabolic_rate_allometric_exponent);
	double turnover_time=
	  the_mean_bodymass_M*
	  energy_content_of_wet_biomass_q/
	  metabolic_rate_eta;
	the_carrying_capacity=
	  dominating_species_primary_production*turnover_time;
      }
      //REPORT(turnover_time*the_plant_growth_rate_sigma);
      if(plant_growth_rate){
	the_turnover_rate_r=plant_growth_rate;
      }else{
	if(!emulate_carrying_capacity_bug){
	  the_turnover_rate_r=1/the_plant_growth_rate_sigma;
	}else{
	  // old buggy version
	  the_turnover_rate_r*=plant_speedup;
	}
      }
      break;
    case 2:
      // monoculture GPP, sigma_max allometry and production
      // efficiency are given
      the_loss_rate_over_max_production_rate_r=
	loss_rate_over_max_production_rate_r;
      {// implement trade-offs:
	if(trait_space_model == 3 &&
	   the_loss_rate_over_max_production_rate_r<0.5){
	  static bool warned=false;
	  if(not warned){
	    REPORT(the_loss_rate_over_max_production_rate_r);
	    WARNING("the_loss_rate_over_max_production_rate_r<0.5");
	    WARNING("This may lead to negative effective loss rates.");
	    warned=true;
	  }
	  the_loss_rate_over_max_production_rate_r+=
	    ( trade_off_multiplier>0 ?
	      (1-loss_rate_over_max_production_rate_r) :
	      loss_rate_over_max_production_rate_r
	      )
	    *
	    trade_off_multiplier;
	}else{
	  the_loss_rate_over_max_production_rate_r+=
	    (1-loss_rate_over_max_production_rate_r)*
	    trade_off_multiplier;
	}
      }
      the_plant_growth_rate_sigma=//this is the gross growth rate!!
	(growth_rate_allometric_prefactor*
	 pow(M,growth_rate_allometric_exponent))/
	(1-loss_rate_over_max_production_rate_r);
      the_carrying_capacity=
	1/(the_plant_growth_rate_sigma*
	   log(1/loss_rate_over_max_production_rate_r)/
	   dominating_species_primary_production );
      the_turnover_rate_r=the_plant_growth_rate_sigma*
	loss_rate_over_max_production_rate_r;
      break;
    /*case 3: // LV type plant interaction:
      // monoculture GPP, sigma_max allometry and production
      // efficiency are given
      the_loss_rate_over_max_production_rate_r=0; // this is not used!
      the_plant_growth_rate_sigma=//this is the max net growth rate
	(growth_rate_allometric_prefactor*
	 pow(M,growth_rate_allometric_exponent));
      the_carrying_capacity=dominating_species_primary_production/
	the_plant_growth_rate_sigma;
      the_turnover_rate_r=the_plant_growth_rate_sigma;
      {// implement trade-offs:
	the_plant_growth_rate_sigma*=(1-trade_off_multiplier);
      }
      break;*/
    case 3: // LV type plant interaction:
      // monoculture GPP, sigma_max allometry and production
      // efficiency are given
      // with random variation introduced in the_turnover_rate
      the_env_effect=
	(growth_rate_mean_random_deviation+growth_rate_std_random_deviation*
	 pseudo_standard_normal_random_with(*this,0));
      the_loss_rate_over_max_production_rate_r=0; // this is not used!
      the_plant_growth_rate_sigma=//this is the max net growth rate
	(growth_rate_allometric_prefactor*
	 pow(M,growth_rate_allometric_exponent));
      the_carrying_capacity=dominating_species_primary_production/
	the_plant_growth_rate_sigma;
      the_turnover_rate_r=the_plant_growth_rate_sigma;
      {// implement trade-offs:
	the_plant_growth_rate_sigma*=(1-trade_off_multiplier);
      }
      break;
    case 4:
      // same as case 2 except that the loss rate, i.e. the_turnover_rate_r, is specified as
      // an allometric scaling law
      the_turnover_rate_r=plant_loss_rate_allometric_prefactor*pow(M,plant_loss_rate_allometric_exponent);
      the_plant_growth_rate_sigma=//this is the gross growth rate!!
	(growth_rate_allometric_prefactor*
	 pow(M,growth_rate_allometric_exponent))+the_turnover_rate_r;
      loss_rate_over_max_production_rate_r=the_turnover_rate_r/the_plant_growth_rate_sigma;
      the_loss_rate_over_max_production_rate_r=
	loss_rate_over_max_production_rate_r;
      {// implement trade-offs:
	if(trait_space_model == 3 &&
	   the_loss_rate_over_max_production_rate_r<0.5){
	  static bool warned=false;
	  if(not warned){
	    REPORT(the_loss_rate_over_max_production_rate_r);
	    WARNING("the_loss_rate_over_max_production_rate_r<0.5");
	    WARNING("This may lead to negative effective loss rates.");
	    warned=true;
	  }
	  the_loss_rate_over_max_production_rate_r+=
	    ( trade_off_multiplier>0 ?
	      (1-loss_rate_over_max_production_rate_r) :
	      loss_rate_over_max_production_rate_r
	      )
	    *
	    trade_off_multiplier;
	}else{
	  the_loss_rate_over_max_production_rate_r+=
	    (1-loss_rate_over_max_production_rate_r)*
	    trade_off_multiplier;
	}
      }
      the_carrying_capacity=
	1/(the_plant_growth_rate_sigma*
	   log(1/loss_rate_over_max_production_rate_r)/
	   dominating_species_primary_production );
      break;
    default:
      FATAL_ERROR("this plant_physiology_version does not exist yet");
    }
  }
}

bool NewSpecies::increase_loss_rate_by(double rate){
  if(is_a_plant()){
    if(plant_physiology_version==3){
      FATAL_ERROR("cannot increase loss rate for plant_physiology_version"
		  << plant_physiology_version);
    }
    the_loss_rate_over_max_production_rate_r+=
      rate/the_plant_growth_rate_sigma;
    if(the_loss_rate_over_max_production_rate_r < 0){
      WARNING("negative loss rate");
      return false;
    }
  }else{ // it's an animal
    the_turnover_rate_r+=
      rate;
    if(the_turnover_rate_r < 0){
      WARNING("negative turnover rate");
      return false;
    }
  }
  return true;
}

// the default constructor is rather dumb for not confusing startup:
NewSpecies::NewSpecies():
#if 1 // this serves only to suppress some complaint from valgrind
  the_taxon(),
  the_mean_bodymass_M(),
  the_log_mean_bodymass_M(),
  the_biomass_abundance_B(),
  the_conversion_efficiency_epsilon(),
  the_attack_rate_a(),
  the_handling_time_T(),
  the_turnover_rate_r(),
  the_fishing_mortality(),
  the_plant_growth_rate_sigma(),
  the_unique_id(),
  the_parent_unique_id(),
  the_clade_id(),
  the_number_of_invading_offspring(),
#endif
  V(niche_space_dimensions_D),
  F(niche_space_dimensions_D)
{
  the_fishing_mortality=0;
}

NewSpecies::~NewSpecies(){}; //just don't try to inline

/*static*/
void NewSpecies::initialize() {
  // Weights of trophic dimensions:
  W=1.0/(2*trophic_niche_width_w_t*trophic_niche_width_w_t)*
    trait_t(niche_space_dimensions_D,1);
  W[0]*=plant_vulnerability_relative_weight;
  Wp=1.0/(2*plant_niche_width_w_r*plant_niche_width_w_r)*
    trait_t(niche_space_dimensions_D,1);
  TRACE(W,SPECIES);
  if(Mmax){
    WARNING("Mmax is set to " << Mmax << "but not used anymore!");
  }

  // Taxa:
#define ASSIGN_TAXON_NAME(name) the_taxon_name[name]=#name;
  ASSIGN_TAXON_NAME(plant);
  ASSIGN_TAXON_NAME(animal);
  // Check if all taxa have names assigned:
  for(int i=n_taxa;i-->0;){
    if(the_taxon_name[i]=="")
      FATAL_ERROR("Taxon with id " << i << " has no name assigned!");
  }
#undef ASSIGN_TAXON_NAME
}

std::string NewSpecies::the_taxon_name[n_taxa];
const std::string & NewSpecies::taxon_name() const{
  return the_taxon_name[the_taxon];
}


// Slow and simple implementation:
NewSpecies::taxon_t NewSpecies::taxon_with_name(const string & name){
  for(int i=n_taxa;i-->0;){
    if(the_taxon_name[i]==name){
      return taxon_t(i);
    }
  }
  FATAL_ERROR("No taxon with name " << name << " defined");
  return not_a_taxon;//this point is never reached.
}

double NewSpecies::V_center() const{
  return (is_a_plant() ? +1 : -1) * plant_vulnerability_separation +
    plant_animal_axis_offset;
}

void NewSpecies::separate(NewSpecies::trait_t & V) const{
  if(vulnerability_separate_in_all_dimensions){
    for(int i=niche_space_dimensions_D;i-->0;){
      V[i]+=V_center();
    }
  }else{
    V[0]+=V_center();
  }
}

void NewSpecies::unseparate(NewSpecies::trait_t & V) const {
  if(vulnerability_separate_in_all_dimensions){
    for(int i=niche_space_dimensions_D;i-->0;){
      V[i]-=V_center();
    }
  }else{
    V[0]-=V_center();
  }
}


NewSpecies::trait_t NewSpecies::add_reflecting(trait_t F,trait_t delta,
					       double radius){
  // add delta to F, but reflect on a hypersphere
  const double f02=radius*radius;
  if(dot(F,F) > f02){
    WARNING("Trait vector too long for hard limit.");
    WARNING("Re-scaling");
    F*=radius/sqrt(dot(F,F));
  }


  trait_t tentative_result=F+delta;
  while(dot(tentative_result,tentative_result)>=f02){
    //the fraction of delta until next reflection (solution of a
    //quadratic equation):
    double s0=
      (-dot(delta,F)+
       sqrt(dot(delta,F)*dot(delta,F)+(f02-dot(F,F))*dot(delta,delta) ))
      /
      dot(delta,delta);
    // the reflection point:
    F=F+s0*delta;
    // the normal to the sphere at reflection:
    trait_t normal=F/radius;
    // the way to go after reflection

    /////REPORT(abs(delta)*(1-s0)); //this...
    delta=(delta-normal*2*(dot(normal,delta)))*(1-s0);
    /////REPORT(abs(delta)); //...and this should be the same.
    tentative_result=F+delta;
  }
  return tentative_result;
}

NewSpecies::trait_t NewSpecies::random_in_sphere(double radius){
  // generates a random trait evenly distributed in the hypersphere
  trait_t normal;
  double a;
  do{
    for(int i=niche_space_dimensions_D;i-->0;){
      normal[i]=gaussian(0,1);
    }
    a=abs(normal);
  }while(a==0);//should never happen
  normal/=a;
  return normal*(pow(unirand(),1.0/niche_space_dimensions_D)*radius);
}

void NewSpecies::set_random_trait_maybe_flat(trait_t &t,double radius,
					     int do_reflect){
  if(do_reflect){
    t=random_in_sphere(radius);
  }else{
    for(int i=niche_space_dimensions_D;
	i-->0;){
      t[i]=gaussian(0,radius);
    }
  }
}

NewSpecies::NewSpecies(int taxon):
  // The next line needs to be fixed to include other taxa:
  the_taxon(taxon_t(taxon==-1 ? unirand()<plant_fraction : taxon)),
  the_mean_bodymass_M(the_taxon==plant?Mini:Mini*pow(10,log10_bodymass_ratio_window_center)),
  the_aggressivity_g(initial_aggressivity),
  the_niche_width_wt(trophic_niche_width_w_t+
		     (unirand()-0.5)*trophic_niche_width_spread),
  V(niche_space_dimensions_D),
  F(niche_space_dimensions_D)
{
  if(trait_space_model!=3){
    set_random_trait_maybe_flat(V,
				(is_a_plant()?
				 plant_vulnerability_variability:
				 animal_vulnerability_variability ),
				flat_vulnerability_distribution );

    // Foraging traits are defaulted with the same distribution as
    // vulnerability traits (!) so this looks a bit confusing:
    set_random_trait_maybe_flat(F,
				(is_a_plant()?
				 plant_trait_variability:
				 plant_vulnerability_variability ),
				(is_a_plant()?
				 flat_plant_distribution:
				 flat_vulnerability_distribution ));
  }else{
    set_random_trait_maybe_flat(V,
				(is_a_plant()?
				 plant_vulnerability_variability:
				 animal_vulnerability_variability ),
				false);
    if(randomly_sampled_foraging_traits &&
       sampling_species_list->size() > sampling_species_list->n_animals ){
      sample_F_randomly();
    }else{
      set_random_trait_maybe_flat(F,
				  plant_vulnerability_variability,
				  false);
    }
  }

  separate(V);

  TRACE(V,SPECIES);
  TRACE(F,SPECIES);

  the_fishing_mortality=0;
  set_parameters_given_body_mass();
  set_unique_id();
  initiate_biomass();
}


/// Gives a new NewSpecies with properties (except for switching
/// exponent) assembled randomly from the samples.  Animals from
/// animals, plants from plants.  Aborts if there are not
/// corresponding samples.
NewSpecies::NewSpecies(const species_list_t & samples,int offset,
		       int sample_size) :
  V(niche_space_dimensions_D),
  F(niche_space_dimensions_D),
  the_number_of_invading_offspring()
{
  ASSERT(sample_size>0);
  ASSERT(sample_size+offset <= samples.size());

  the_taxon=
    samples(offset+random_integer(sample_size)).the_taxon;
  //(but usually all is_plant should be the same)

  the_mean_bodymass_M=
    samples(offset+random_integer(sample_size)).the_mean_bodymass_M;

  the_aggressivity_g=
    samples(offset+random_integer(sample_size)).the_aggressivity_g;

  {
    const permutation p=random_permutation(niche_space_dimensions_D);
    const trait_t & V_sample=samples(offset+random_integer(sample_size)).V;
    for(int i=niche_space_dimensions_D;
	i-->0;){
      V[p[i]]=V_sample[i];
    }
  }

  {
    const permutation p=random_permutation(niche_space_dimensions_D);
    const trait_t & F_sample=samples(offset+random_integer(sample_size)).F;
    for(int i=niche_space_dimensions_D;
	i-->0;){
      F[p[i]]=F_sample[i];
    }
  }

  TRACE(the_taxon,SPECIES);
  TRACE(V,SPECIES);
  TRACE(F,SPECIES);

  set_unique_id();
  the_clade_id=the_unique_id; //clad id is unique id of clad founder
  the_parent_unique_id=the_unique_id;
  TRACE(the_parent_unique_id,SPECIES);
  set_parameters_given_body_mass();

  initiate_biomass();
}


NewSpecies::trait_t NewSpecies::W=NewSpecies::trait_t(0);
NewSpecies::trait_t NewSpecies::Wp=NewSpecies::trait_t(0);

void NewSpecies::add_random_to(trait_t &t,
			       double step ){
  for(int i=niche_space_dimensions_D;i-->0;){
    t[i]+=gaussian(0,step);
  }
}

void NewSpecies::add_random_maybe_reflecting_to(trait_t &t,
						double step,double radius,
						int do_reflect){
  trait_t delta;
  for(int i=niche_space_dimensions_D;i-->0;){
    delta[i]=gaussian(0,step);
  }
  if(do_reflect){
    t=add_reflecting(t,delta,radius);
  }else{
    switch(trait_evolution_version){
    case 1:
      t+=delta;
      t*=(1/sqrt(1+(step/radius)*(step/radius)));
      break;
    case 2:
      ALWAYS_ASSERT(step < radius);
      t*=sqrt(1-(step/radius)*(step/radius));
      t+=delta;
      break;
    default:
      FATAL_ERROR("trait_evolution_version " << trait_evolution_version
		  << " not implemented");
    }
  }
}

int NewSpecies::number_of_kinds_of_traits(taxon_t taxon){
  switch(taxon){
  case animal:
    return 4;
  case plant:
    return 3;
  }
  FATAL_ERROR("taxon " << taxon << " is unknown");
  return 0;//should never be reached.
}

const std::string vulnerability_trait_name="V";
const std::string foraging_trait_name="F";
const std::string plant_trait_name="G";
const std::string aggressivity_trait_name="A";
const std::string body_mass_trait_name="M";
const std::string niche_width_trait_name="W";

const std::string & NewSpecies::trait_name(int i){
  kind_of_trait_t t=kind_of_trait_t(i);

  switch(t){
  case foraging_traits:
    return (is_a_plant()?plant_trait_name:foraging_trait_name);
    break;
  case vulnerability_traits:
    return vulnerability_trait_name;
    break;
  case body_mass_trait:
    return body_mass_trait_name;
    break;
  case aggressivity_trait:
    return aggressivity_trait_name;
    break;
  case niche_width_trait:
    return niche_width_trait_name;
    break;
  default:
    FATAL_ERROR("no such kind of trait");
  }
}

int NewSpecies::trait_dimensions(int i){
  kind_of_trait_t t=kind_of_trait_t(i);

  switch(t){
  case foraging_traits:
    return niche_space_dimensions_D;
    break;
  case vulnerability_traits:
    return niche_space_dimensions_D;
    break;
  case body_mass_trait:
    return 1;
    break;
  case aggressivity_trait:
    return 1;
    break;
  default:
    FATAL_ERROR("no such kind of trait");
  }
}

NewSpecies::trait_t NewSpecies::get_trait_values(int i){
  kind_of_trait_t t=kind_of_trait_t(i);

  switch(t){
  case foraging_traits:
    return F;
    break;
  case vulnerability_traits:
    return V;
    break;
  case body_mass_trait:
    {
      trait_t M;
      M[0]=the_log_mean_bodymass_M;
      return M;
      break;
    }
  case aggressivity_trait:
    {
      trait_t a;
      a[0]=log(the_attack_rate_a);
      return a;
      break;
    }
  default:
    FATAL_ERROR("no such kind of trait");
  }
}

const species_list_t * NewSpecies::sampling_species_list=0;

void NewSpecies::mutate_selected_traits(trait_selection_t selection,
					double slowdown /*defaults to 1*/){

  if(mu_V){// old input data used, convert to new version:
    mu_V_animal=mu_V;
    mu_V_plant=mu_V*plant_mu_V_boost;
    mu_V=0;
    plant_mu_V_boost=0;
  }

  if(selection & (trait_selection_t(1)<<body_mass_trait)){
    // mutate body mass
    if(slowdown>=1){
      the_mean_bodymass_M*=
	pow(typical_bodymass_ratio_d,gaussian(0,1));
      // if(trait_space_model==4 && the_mean_bodymass_M < Mmin){
      // 	// reflect to value > Mmin on the log scale:
      // 	the_mean_bodymass_M=Mmin*Mmin/the_mean_bodymass_M;
      // }
    }
  }

  if(selection & (trait_selection_t(1)<<aggressivity_trait)){
    // mutate aggressivity
    if(slowdown>=1 && !is_a_plant()){
      if(use_aggressivity){
	if(controlled_random_link_strengths <= 1){
	  the_aggressivity_g*=
	    aggressivity_decay*
	    pow(typical_aggressivity_ratio_da,gaussian(0,1));
	}else{
	  double root=sqrt(-log(the_aggressivity_g));
	  root+=-log(aggressivity_decay)+
	    log(typical_aggressivity_ratio_da)*gaussian(0,1);
	  the_aggressivity_g=exp(-root*root);
	}
      }
    }
  }

  if(selection & (trait_selection_t(1)<<aggressivity_trait)){
    // mutate trophic niche width
    if(slowdown>=1 && !is_a_plant()){
      if(typical_niche_width_ratio){
	the_niche_width_wt*=
	  niche_width_decay*
	  pow(typical_niche_width_ratio,gaussian(0,1));
      }
    }
  }

  if(selection & (trait_selection_t(1)<<vulnerability_traits)){
    // mutate vulnerability traits
    if(trait_space_model!=3){
      unseparate(V);
      add_random_maybe_reflecting_to(V,
				     slowdown*
				     (is_a_plant()?
				      mu_V_plant:
				      mu_V_animal ),
				     (is_a_plant()?
				      plant_vulnerability_variability:
				      animal_vulnerability_variability ),
				     flat_vulnerability_distribution );
      separate(V);
    }else{//trait_space_model==3
      add_random_to(V,slowdown*( is_a_plant() ? mu_V_plant : mu_V_animal ));
    }
  }

  if(selection & (trait_selection_t(1)<<foraging_traits)){
    // mutate foraging or plant competition traits
    if(trait_space_model!=3){
      if(is_a_plant()){
	trait_t & G=F;
	add_random_maybe_reflecting_to(G,
				       slowdown*mu_G,
				       plant_trait_variability,
				       flat_plant_distribution );
      }else{
	add_random_maybe_reflecting_to(F,
				       slowdown*mu_F,
				       foraging_variability_f0,
				       flat_foraging_distribution );
      }
    }else{
      if(is_a_plant()){
	trait_t & G=F;
	add_random_to(G,slowdown*mu_G);
      }else{
	if(randomly_sampled_foraging_traits>0){
	  if(randomly_sampled_foraging_traits>=1){
	    sample_F_randomly();
	  }else{
	    double scale=randomly_sampled_foraging_traits*slowdown;
	    trait_t hold=F;
	    sample_F_randomly();
	    F*=scale;
	    switch(trait_evolution_version){
	    case 1:
	      F+=hold;
	      F*=(1/sqrt(1+scale*scale));
	      break;
	    case 2:
	      hold*=sqrt(1-scale*scale);
	      F+=hold;
	      break;
	    }
	    sample_F_randomly();
	  }
	}else{
	  add_random_to(F,slowdown*mu_F);
	}
      }
    }
  }

  set_parameters_given_body_mass();

  return;
}

void NewSpecies::sample_F_randomly(){
  const species_list_t & samples=*sampling_species_list;
  if(no_carnivores){
    int n_plants=samples.size()-samples.n_animals;
    for(int i=niche_space_dimensions_D;i-->0;){
      F[i]=samples[samples.n_animals+random_integer(n_plants)].
	V[random_integer(niche_space_dimensions_D)];
    }
  }else{
    int n_species=samples.size();
    for(int i=niche_space_dimensions_D;i-->0;){
      F[i]=samples[random_integer(n_species)].
	V[random_integer(niche_space_dimensions_D)];
    }
  }
}


double NewSpecies::get_threshold_biomass_B() const {
  return
    max<double>(2*the_mean_bodymass_M,
		minimum_biomass_considered);
}

void NewSpecies::initiate_biomass(){
  the_biomass_abundance_B=
    initial_number_of_individuals*
    get_threshold_biomass_B();
}



NewSpecies *
NewSpecies::speciate_as_new(double slowdown /*defaults to 1*/) const{
  NewSpecies * child_p=new NewSpecies(*this);
  NewSpecies & child=*child_p;  // for convenience

  trait_selection_t traits_to_mutate;
  if(is_a_plant() and !plant_size_evolves){
    traits_to_mutate=
      (trait_selection_t(-1) ^
       (trait_selection_t(1) << body_mass_trait) );
  }else{
    traits_to_mutate=trait_selection_t(-1);
  }

  child.mutate_selected_traits(traits_to_mutate,slowdown);

  TRACE(child.the_taxon,SPECIES);
  TRACE(child.V,SPECIES);
  TRACE(child.F,SPECIES);

  child.set_unique_id();
  child.the_parent_unique_id=child.the_unique_id;
  child.the_clade_id=child.the_unique_id;
  child.the_number_of_invading_offspring=0;

  child.initiate_biomass();

  return child_p;
}

NewSpecies * NewSpecies::speciate(double slowdown /*defaults to 1*/){
  NewSpecies *child_p=speciate_as_new(slowdown);
  child_p->the_clade_id=the_clade_id;
  child_p->the_parent_unique_id=the_unique_id;
  return child_p;
}



double NewSpecies::log_niche_matching(const NewSpecies & prey) const{
  if(totally_random_link_strength_sigma)
    return log_random_matching(prey);

  double dot_product=0;
  if(!typical_niche_width_ratio){
    for(int i=niche_space_dimensions_D;i-->0;){
      double diff=prey.V[i]-F[i];
      dot_product+=(diff*diff)*W[i];
    }
  }else{
    for(int i=niche_space_dimensions_D;i-->0;){
      double diff=prey.V[i]-F[i];
      dot_product+=(diff*diff);
    }
    dot_product/=2*the_niche_width_wt*the_niche_width_wt;
  }
  return -dot_product;
//   return -dot(prey.V-F,W*(prey.V-F));
}

double NewSpecies::log_random_matching(const NewSpecies & prey) const{
  double standard_normal=
    pseudo_standard_normal_random_with(prey);

  // From this, generate a "random" link strength:
  if(not controlled_random_link_strengths){
    // conventional case
    return (the_niche_width_wt/trophic_niche_width_w_t)*
      totally_random_link_strength_sigma * standard_normal;
  }else{
    // with controlled_random_link_strengths
    double as_ratio=-log(the_aggressivity_g);
    double target_S=sqrt(2*M_PI)*as_ratio*exp(as_ratio*as_ratio/2);
    double sigma=sqrt(2*log(target_S))/the_niche_width_wt;
    return sigma *
      ( do_switching ? 0.5*(1+universal_switching_similarity) : 1 ) *
      standard_normal;
  }
}

double
NewSpecies::pseudo_standard_normal_random_with(const NewSpecies & other,
					       const int choice) const{

  ALWAYS_ASSERT(0 <= choice && choice < 4);

  // First, compute the MD5 sum of the two random tags
  MD5 md5;
  random_tag_t in_data[]={the_random_tag,
			  other.the_random_tag,
			  the_unique_id,
			  other.the_unique_id};
  if(emulate_log_random_matching_bug){
    md5.update((unsigned char *)in_data, sizeof(in_data)/2);
  }else{
    md5.update((unsigned char *)in_data, sizeof(in_data));
  }
  md5.finalize();
  unsigned int * out_data=
    (unsigned int *) md5.digest_ptr(); // digest is worth 4 integers as least.

  // Now we have four "random" numbers depending on the random tags of
  // consumer and resource.  We convert these into a normally
  // distributed number using the Box–Muller method:

  const double int_norm=1/(((double) (1<<16))*((double) (1<<16)));

  double standard_normal;

  if(choice & 1){
    standard_normal=
      sqrt(-2*log((out_data[choice & 2]+0.5)*int_norm))*
      sin((2*M_PI*int_norm)*out_data[(choice & 2) + 1]);
  }else{
    standard_normal=
      sqrt(-2*log((out_data[choice & 2]+0.5)*int_norm))*
      cos((2*M_PI*int_norm)*out_data[(choice & 2) + 1]);
  }


  return standard_normal;
}

double NewSpecies::log_small_prey_matching(const NewSpecies & prey) const{
  return
    (small_prey_exponent_alpha *
     (prey.the_log_mean_bodymass_M-the_log_mean_bodymass_M) )
    ;
}


double NewSpecies::log_large_prey_matching(const NewSpecies & prey) const{
  return
    (big_plants_eaten && prey.is_a_plant() ?
     0 :
     -predator_relative_size_Lambda * prey.the_mean_bodymass_M / the_mean_bodymass_M )
    ;
}


double NewSpecies::log_mass_matching(const NewSpecies & prey) const{
  if(prey.the_mean_bodymass_M/the_mean_bodymass_M <
     cutoff_prey_predator_mass_ratio){
    return 2*log_DBL_MIN;
  }

  const double log_mass_ratio=
    the_log_mean_bodymass_M-prey.the_log_mean_bodymass_M;
  const double x=log_mass_ratio-
    log10_bodymass_ratio_window_center*log(10.0);
  if(x>0){
    return -small_prey_exponent_alpha*x;
  }else{// x<=0
    if(prey.is_a_plant() and big_plants_eaten){
      return -0.25*x;//emergency measure!!
    }else{
      return big_prey_exponent_beta*x;
    }
  }
}


double NewSpecies::log_other_matching(const NewSpecies & prey) const{
  double log_eatable_fraction=0;
  if(plant_hardening_time && prey.is_a_plant()){
    log_eatable_fraction=
      log(1-exp(-prey.the_plant_growth_rate_sigma*plant_hardening_time));
  }
  return log_eatable_fraction +
    ( prey.is_a_plant() ? 0 : log(carnivore_aggressivity_boost) );
}

double NewSpecies::log_foraging_strength_c_on(const NewSpecies & prey) const{
  if(no_carnivores && !prey.is_a_plant())
    return log_DBL_MIN-1;

  // REPORT(log_mass_matching(prey));
  // REPORT(log_niche_matching(prey));
  // REPORT(log_other_matching(prey));
  // static average_meter lnm;
  // lnm.sample(log_niche_matching(prey));
  // REPORT(lnm);
  // REPORT(lnm.std());

  double log_link_strength_error=0;
  if(typical_relative_error_in_link_strength){
    log_link_strength_error=
      typical_relative_error_in_link_strength*
      random_normal_with(prey);
  }

  return
    log_link_strength_error
    +
    log_mass_matching(prey)
    +
    log_niche_matching(prey)
    +
    log_other_matching(prey)
    ;
}

double NewSpecies::foraging_strength_c_on(const NewSpecies & prey) const{
  double log_c=log_foraging_strength_c_on(prey);
  if(log_c < log_DBL_MIN)
    return exp(log_DBL_MIN);
  else if(log_c > log_DBL_MAX)
    return exp(log_DBL_MAX);
  else
    return exp(log_c);
}

double NewSpecies::random_with(const NewSpecies & other,int which) const{
  // First, compute the MD5 sum of the two random tags
  MD5 md5;
  random_tag_t in_data[]={the_random_tag,
			  other.the_random_tag,
			  the_unique_id,
			  other.the_unique_id};
  if(emulate_log_random_matching_bug){
    md5.update((unsigned char *)in_data, sizeof(in_data)/2);
  }else{
    md5.update((unsigned char *)in_data, sizeof(in_data));
  }
  md5.finalize();
  unsigned int * out_data=
    (unsigned int *) md5.digest_ptr(); // digest is worth 4 integers as least

  // Now we have four "random" numbers depending on the random tags of
  // consumer and resource. We generate one evenly distributed random
  // number, depending on the value of the parameter 'which'.
  const double int_norm=1/(((double) (1<<16))*((double) (1<<16)));

  return out_data[which]*int_norm;
}

double NewSpecies::random_normal_with(const NewSpecies & other,int which) const{
  // First, compute the MD5 sum of the two random tags
  MD5 md5;
  // Use Box-Muller method to convert even random numbers into
  // standard normal ones:
  double U=random_with(other,1);
  double V=random_with(other,2);
  switch(which){
  case 0:
    return sqrt(-2*log(U)) * cos(2*M_PI*V);
  case 1:
    return sqrt(-2*log(U)) * sin(2*M_PI*V);
  default:
    FATAL_ERROR("random_normal_with can only return two different random numbers");
  }
}

double NewSpecies::plant_hampering_by(const NewSpecies & other) const{
  if(plant_shadowing && -the_mean_bodymass_M/other.the_mean_bodymass_M < log_DBL_MIN){
    return 0;
  }else{
    double dot_product=0;
    for(int i=niche_space_dimensions_D;i-->0;){
      double diff=F[i]-other.F[i];
      dot_product+=diff*diff*Wp[i];
    }
    double competition;
    if(no_plant_competition){
      competition=int(*this == other);
      //competition is now 1 if *this==other, 0 otherwise
      if(experimental && competition){
	competition *= exp(4*(random_with(other,1)-0.5));
      }
    }else{
      competition=
	exp(-pow(dot_product,plant_competition_drop_off_exponent/2)-
	    (plant_shadowing ?
	     the_mean_bodymass_M/other.the_mean_bodymass_M :
	     0  ));
      //competition is now 1 if *this==other, some smaller value otherwise
    }
    // The following if-clause adds plant_diffuse_competition to
    // competition with probability
    // plant_diffuse_competition_probability.  The result is then
    // divided by the carrying capacity of the other species and
    // returned.
    if(plant_diffuse_competition_symmetry){
      if(other.random_with(*this) <
	 plant_diffuse_competition_probability ){
	competition +=
	  plant_diffuse_competition *
	  plant_diffuse_competition_symmetry / 2;
      }
    }
    if(plant_diffuse_competition_probability==1 ||
       random_with(other)<plant_diffuse_competition_probability ){
      competition +=
	plant_diffuse_competition *
	(2-plant_diffuse_competition_symmetry) / 2;
    }
    if(unit_diagonal_plant_competition && *this == other)
      competition = 1;
    if(small_plants_win){
      return
	competition*the_plant_growth_rate_sigma/
	dominating_species_primary_production;
    }else{
      return
	competition/other.the_carrying_capacity;
      // competition*
      // other.the_plant_growth_rate_sigma/dominating_species_primary_production;
    }
  }
}


bool NewSpecies::outside_allowed_trait_space(){
  if(trait_space_model==3){
    return false;
  }

  if(is_a_plant() &&
     flat_plant_distribution &&
     dot(F,F) > plant_trait_variability*plant_trait_variability
     ){
    return true;
  }

  if(!is_a_plant() &&
     flat_foraging_distribution &&
     dot(F,F) > foraging_variability_f0*foraging_variability_f0 ){
    return true;
  }

  if(flat_vulnerability_distribution){
    unseparate(V);
    if(dot(V,V) >
       (is_a_plant() ?
	plant_vulnerability_variability*plant_vulnerability_variability :
	animal_vulnerability_variability*animal_vulnerability_variability ))
      return true;
    separate(V);
  }

  return false;
}

double NewSpecies::switching_similarity_to(const NewSpecies & other,
					   const double width) const{
  if(totally_random_link_strength_sigma){
    if(*this == other || !do_switching){
      return 1;
    }else{
      return universal_switching_similarity;
    }
  }else{
    /*
    // this turns out to be too slow, so we do it explicityly:
    return exp(-dot(other.V-V,Ws*(other.V-V)));
    */
    double sum=0;
    for(int i=niche_space_dimensions_D;i-->0;){
      double delta=V(i)-other.V(i);
      sum+=(delta*delta);
    }
    return exp(-sum*(1.0/(2*width*width)));
  }
}

double NewSpecies::foraging_distance_to(const NewSpecies & prey) const{
  if(is_a_plant())
    return 0;
  else
    return sqrt(dot(prey.V-F,W*(prey.V-F)));
}


void NewSpecies::fix_foraging_vector(){
  if(is_a_plant()){
    F/=foraging_variability_f0;
  }
}

void NewSpecies::perturb_physiological_parameters(double factor){
  if(is_a_plant()){
    the_plant_growth_rate_sigma*=factor;
    the_carrying_capacity*=factor; // this requires initialization of
				   // NewWeb afterwards!
  }else{
    // is animal
    the_turnover_rate_r/=factor;
  }
}


void NewSpecies::make_predictions(){
  cout << "************ STARTING PREDICTIONS **************" << endl;
  my_evaluator_t units;
  const double gross_conversion_efficiency=
    conversion_efficiency_epsilon*(1+1/production_over_respiration);

  REPORT(gross_conversion_efficiency);

  const double small_plant_growth_rate=growth_rate_allometric_prefactor*
    pow(Mmin/allometric_base_unit,growth_rate_allometric_exponent);
  REPORT(small_plant_growth_rate/units("day^(-1)"));

  cout << "log10 plant animal mass ratio:" << endl;
  REPORT(log10(pow(max_growth_rate_allometric_prefactor/
		   (growth_rate_allometric_prefactor*
		    conversion_efficiency_epsilon)
		   ,4)));
  const double respiration_allometric_prefactor=
    max_growth_rate_allometric_prefactor/production_over_respiration;

  const double max_plant_growth_rate_sigma=growth_rate_allometric_prefactor*
    pow(Mmin/allometric_base_unit,growth_rate_allometric_exponent);

  const double true_max_production=dominating_species_primary_production*
    log(1/get_cfg_parameter("loss_rate_over_max_production_rate_r"));

  const double min_carrying_capacity=true_max_production/
    max_plant_growth_rate_sigma;

  REPORT(max_plant_growth_rate_sigma/units("seconds^(-1)"));
  REPORT(min_carrying_capacity/area_per_compartment/units("kg/meter^2"));

  const double carrying_capacity_allometric_prefactor=
    true_max_production/growth_rate_allometric_prefactor;

  REPORT(log10(carrying_capacity_allometric_prefactor/area_per_compartment/units("kg/meter^2")));

  // output the PPMR window
  ofstream PPMR_window("PPMR_window.dat");
  NewSpecies predator(animal),prey(animal);
  predator.the_log_mean_bodymass_M=1;
  for(int i=-20;i<20;i++){
    prey.the_log_mean_bodymass_M=i;
    PPMR_window << i/log(10) << " " << predator.log_mass_matching(prey) << endl;
  }

  // attacke rate math:
  NewSpecies eater(animal),eaten(plant);
  double mass_matching_factor=exp(eater.log_mass_matching(eaten));
  double Ta=eater.attack_rate_a()*eater.handling_time_T();
  REPORT(eater.bodymass());
  REPORT(eaten.bodymass());
  REPORT(mass_matching_factor);
  REPORT(Ta*eval("kilogram"));
  double alpha=4.25,c=8.55;// power-law regression of FVdist.dat
  double expected_power_sum=
    alpha*pow(trophic_niche_width_w_t*sqrt(2),alpha)*tgamma(1+alpha/2)*
    c;
  REPORT(expected_power_sum);
  double predicted_mean_plant_biomass_density=
    1/
    (
     production_over_respiration*
     Ta*
     expected_power_sum*
     mass_matching_factor*
     area_per_compartment);
  REPORT(predicted_mean_plant_biomass_density*eval("meter^2/kilogram"));

#if 0
  cout << "*********** BEGIN Equivalent Miniweb parameters *******" << endl;
  average_meter log_link_strength;
  for(int i=10000;i-->0;){
    eaten.mutate_selected_traits(trait_selection_t(1)<<vulnerability_traits);
    eater.mutate_selected_traits(trait_selection_t(1)<<vulnerability_traits);
    eater.F=eater.V;
    log_link_strength.sample(eater.log_foraging_strength_c_on(eaten));
  }
  double unit_time=1/eaten.plant_growth_rate_sigma();
  double
    mu=log_link_strength.readout()+
    log(eater.attack_rate_a()*unit_time/eaten.plant_hampering_by(eaten)),
    sigma=log_link_strength.std(),
    epsilon=eater.conversion_efficiency_eps(),
    r=eater.turnover_rate_r()*unit_time,
    T=eater.handling_time_T()/unit_time,
    K=1,
    tau=1, // must be 1 to avoid some confusion in Miniweb
    animal_fraction=1/(get_cfg_parameter("plant_animal_ratio")+1);

  double true_carrying_capacity=
    1/eaten.plant_hampering_by(eaten)/area_per_compartment*
    eval("meter^2/kilogram");
  REPORT(true_carrying_capacity);
  REPORT(log10(true_carrying_capacity));

  double prediced_f=r*T/epsilon;
  REPORT(prediced_f);

  REPORT(mu);
  REPORT(sigma);
  REPORT(epsilon);
  REPORT(r);
  REPORT(T);
  REPORT(K);
  REPORT(tau);
  REPORT(animal_fraction);

  cout << "*********** END Equivalent Miniweb parameters *******" << endl;
#endif

  if(get_cfg_parameter("no_animals")){
    const double v=plant_trait_variability/plant_niche_width_w_r,
      D=niche_space_dimensions_D;
    const double predicted_plant_richness=1/
      (4*(exp(2*pow(v,4)/D)-1)*exp(-2*pow(v,2)+2*pow(v,4)/D));
    REPORT(predicted_plant_richness);
  }
}

std::ostream & operator<<(std::ostream & os,NewSpecies::taxon_t t){
  os << NewSpecies::the_taxon_name[t];
  return os;
}

std::istream & operator>>(std::istream & is,NewSpecies::taxon_t & t){
  string name;
  is >> name;
  t=NewSpecies::taxon_with_name(name);
  return is;
}

NewSpecies::trait_sweep_t
NewSpecies::trait_sweep(const char * details) const{
  int dim=atoi(details);
#define TRAIT F[dim]
  int n_points=100;
  NewSpecies model=*this;
  for(int i=2;i<niche_space_dimensions_D;++i){
    model.F[i]=0;
  }
  trait_sweep_t sweep;
  sweep.resize(n_points,trait_sweep_sample_t());
  double max_val=10*plant_trait_variability;
  double min_val=-max_val;

  sequence< double > delta, val, new_delta;
  delta[0]=model.TRAIT-min_val;
  for(int i=0;i<n_points;++i){
    val[0]=min_val+i*(max_val-min_val)/n_points;
    new_delta=model.TRAIT-val;
    if(abs(new_delta)<abs(delta)){
      delta=model.TRAIT-val;
    }
  }
  REPORT(model.TRAIT);


  for(int i=0;i<n_points;++i){
    val[0]=min_val+i*(max_val-min_val)/n_points;
    val+=delta;
    NewSpecies * s=new NewSpecies(*this);
    s->TRAIT=val[0];
    s->set_parameters_given_body_mass();
    sweep[i]=trait_sweep_sample_t(val,s);
  }
  return sweep;
}
