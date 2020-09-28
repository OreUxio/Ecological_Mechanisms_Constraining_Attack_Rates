// -*- mode: c++ -*-
// $Id: NewSpecies.h 2462 2016-01-23 08:08:45Z axel $
#ifndef _NEWSPECIES_H_
#define _NEWSPECIES_H_

#define NEW_POPULATION_DYNAMICS

#include <vector>
#include <algorithm>
#include <string>
#include "remember.h"
#include "sequence.h"
#include "md5.h"

// These externs are defined in NewSpecies.cc
extern double predator_relative_size_Lambda;
extern double plant_hardening_time;
extern int big_plants_eaten;
extern double plant_shadowing; // why double??? should be int
extern int niche_space_dimensions_D;
extern double area_per_compartment;
extern int emulate_carrying_capacity_bug;
extern int emulate_animal_body_mass_cutoff_bug;
extern int emulate_log_random_matching_bug;
extern int plant_physiology_version;
extern double loss_rate_over_max_production_rate_r;
extern double minimum_biomass_considered;
extern double trophic_niche_width_w_t;
extern double switching_similarity_width_w_s;
extern double relative_switching_similarity_width;
extern int adjust_prey_similarity_width;
extern int do_switching; ///< determines whether to include prey-switching

class NewSpecies; // Forward declaration.
/// A list of animal and plant species.
class species_list_t: public sequence<NewSpecies>{
public:
  /// Number of animals in web.
  /** This information is needed to organize the species list, because
      we keep first all animals, and then all plants in the list. */
  int n_animals;
  bool is_plant(int i) const {return i>= n_animals;}
  species_list_t():n_animals(0){};
  __attribute__ ((__always_inline__)) 
  NewSpecies & operator[](unsigned int i) {
    return sequence<NewSpecies>::operator[](i);
  }
  __attribute__ ((__always_inline__))
  const NewSpecies & operator[](unsigned int i) const{
    return sequence<NewSpecies>::operator[](i);
  }
};

class NewSpecies : public permanent
{
public:
  // This is perhaps the mimimal public interface to taxonomy I can think of:
  /// The most specific taxon of a species.
  const std::string & taxon_name() const;
  /// True if this species does not eat from other compartments.
  inline bool is_a_plant() const;
  enum taxon_t {animal=0,plant=1,n_taxa,not_a_taxon};
  friend std::ostream & operator<<(std::ostream & os,taxon_t t);
  friend std::istream & operator>>(std::istream & is,taxon_t & t);
private:
  static std::string the_taxon_name[n_taxa];
  /// Assigns indices to taxa.
  static taxon_t taxon_with_name(const std::string & name);
  /// The taxon of this species.
  taxon_t the_taxon;
public:
  enum kind_of_trait_t {
    body_mass_trait,
    vulnerability_traits,
    foraging_traits,
    aggressivity_trait,
    niche_width_trait,
    n_traits
  };
  static int number_of_kinds_of_traits(taxon_t taxon);
  const std::string & trait_name(int t);
  typedef unsigned int trait_selection_t;
public:
  /// A class representing trait vectors.
  /** Elementary arithmetic operations on vectors, and a few more
      specific members are defined. The user must make sure that the
      numbers of dimensions in arithmetic operations are compatible.*/
  class trait_t : public sequence<double>{
  public:
    trait_t():sequence<double>(niche_space_dimensions_D){};
    trait_t(const sequence<double> & s):sequence<double>(s){};
    trait_t(int i):sequence<double>(i){};
    trait_t(int i,double d):sequence<double>(i,d){};
    trait_t operator=(const trait_t & other){
      return trait_t(sequence<double>::operator=(other));
    }
    trait_t & operator+=(const trait_t & other){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)+=other(i);
      }
      return *this;
    }
    trait_t & operator-=(const trait_t & other){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)-=other(i);
      }
      return *this;
    }
    trait_t & operator*=(const trait_t & other){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)*=other(i);
      }
      return *this;
    }
    trait_t & operator/=(const trait_t & other){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)/=other(i);
      }
      return *this;
    }
    trait_t & operator*=(const double x){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)*=x;
      }
      return *this;
    }
    trait_t & operator/=(const double x){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)*=1/x;
      }
      return *this;
    }
    trait_t & operator+=(const double x){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)+=x;
      }
      return *this;
    }
    trait_t & operator-=(const double x){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)-=x;
      }
      return *this;
    }
    /// Square distance to another trait vector.
    double dist_square(const trait_t & other) const{
      double sum=0;
      double d;
      for(int i=niche_space_dimensions_D;i-->0;){
	d=(*this)(i)-other(i);
	sum+=d*d;
      }
      return sum;
    }
    /// Multipy-and-add operation:
    void add_multiple(double a,const trait_t other){
      for(int i=niche_space_dimensions_D;i-->0;){
	(*this)(i)+=a*other(i);
      }
    }
    //// Compiler warned this would never the be used! So, I took it out.
    // /// Convert trait_t vector to STL vector:
    // operator std::vector<double> () const {
    //   return 
    // 	std::vector<double>(begin(),end());
    // }

//     operator const HepVector () const {
//       HepVector v(niche_space_dimensions_D);
//       return 
// 	std::vector<double>(&(*this)[0],&(*this)[niche_space_dimensions_D]);
//     }
  };
  void adjust_trait_vector_dimensions();
private:
  /// Physiologically most easy value of first vulnerability dimension.
  /** In the Population-dynamical matching model, the vulnerabilities
      of plants and animals cover different, overlapping regions in
      trophic trait space.  This is achived by shifting the
      "preferred" value of the first component of the vulnerability
      trait vector away from the coordinate center by an amount
      depending on the type of species.  This amount is here
      computed. If vulnerability_separate_in_all_dimensions is
      non-zero, not only the first, but all components are shifted by
      this amount.*/
  double V_center() const;
  /// Separate the centers of vulerability depending on type of species
  void separate(trait_t & V) const;
  /// Undo separation of the centers of vulerability.
  void unseparate(trait_t & V) const;
private:
  // These members allow experimentation with different ways of
  // mutating traits.
  /// Add \a delta to \a F, reflecting the result back into a sphere.
  static trait_t add_reflecting(trait_t F,trait_t delta,double radius);
  /// Chose a trait randomly, evenly distributed in a sphere.
  static trait_t random_in_sphere(double radius);
  /// Mutate \a t by \a step
  static void add_random_to(trait_t &t, double step);
  /// Mutate \a t by \a step, with mechanism depending on \a do_reflect.
  static void add_random_maybe_reflecting_to(trait_t & t,
					     double step,double radius,
					     int do_reflect);
  /// Randomly choose \a t, by mechanism depending on \a do_reflect.
  static void set_random_trait_maybe_flat(trait_t &t,double radius,
					  int do_reflect);
public:
  bool outside_allowed_trait_space();
private:
  static double log_DBL_MIN; // the values for NewSpecies may be 
  static double log_DBL_MAX; // different from elsewhere
  static int next_unique_id;
  static int next_sequential_id; 

public:
  static double plant_fraction; ///< Default plant fraction.
private:
  // Several constants here that determine population dynamics of a
  // species.
  double the_mean_bodymass_M;
  double the_log_mean_bodymass_M;
  double the_biomass_abundance_B;
  double the_conversion_efficiency_epsilon;
  double the_attack_rate_a;
  double the_handling_time_T;
  double the_turnover_rate_r;
  double the_fishing_mortality;
  double the_plant_growth_rate_sigma;
  double the_carrying_capacity;
  double the_aggressivity_g;
  double the_loss_rate_over_max_production_rate_r;
  double the_niche_width_wt;
  double the_number_of_invading_offspring;
  double the_env_effect;
public:
  double the_GP;
  double the_saturation_strength;
  double the_light_strength;
  double the_top_down_strength;
  double the_competition_strength;
private:  
  int the_unique_id;  ///< Each species is assigned a unique number
  int the_sequential_id;  ///< Each species is assigned a unique number
  int the_parent_unique_id;  ///< Unique id of the parent species
  int the_clade_id; ///< Keeps track of clad
  typedef unsigned int random_tag_t;
  random_tag_t the_random_tag; ///< Used to link species randomly.
  trait_t V;  ///< Vulnerability traits
  trait_t F;  ///< Foraging traits or plant traits.
public:
  /// Returns vulnerbility trait vector.
  const trait_t & vulnerability_V() const {return V;}
  /// Returns vulnerbility trait vector relative to preferred value.
  const trait_t relative_vulnerability() const {
    trait_t Vr=V;
    unseparate(Vr);
    return Vr;
  }
  /// Returns foraging trait vector.
  const trait_t & foraging_F() const {return F;}
  /// Returns plant trait vector.
  const trait_t & niche_G() const {return F;}
  static trait_t W; ///< Weights of trophic trait dimensions.
  static trait_t Wp; ///< Weights of plant trait dimensions.
  static void initialize();///< Initialize static members, NewSpecies.

  ///A scaling factor of dimension TIME representative of
  ///population dynamics of the species.
  double slowness() const{
    return 
      (is_a_plant()?
       1/the_plant_growth_rate_sigma :
       the_handling_time_T/the_conversion_efficiency_epsilon);
  }
private:
  void set_unique_id();
public:
  void set_sequential_id();
  void set_parameters_given_body_mass();
  void fix_foraging_vector();///< Used when reading older webs.
private:
  virtual void data_mapping();///< For saving/loading species.
public:
  NewSpecies();
  ~NewSpecies();
  /// Sets species properties as if invading.  
  /** \a taxon==-1 means choose species type at random. */
  explicit NewSpecies(int taxon);
  /// Combine a random species from \a samples.
  /** Gives a new NewSpecies with properties assembled randomly from
      the samples.  Animals from animals, plants from plants.  Aborts
      if there are not corresponding samples. (Check code before
      usage.)*/
  NewSpecies(const species_list_t & samples,int offset, int sample_size); 
  /// Are the species identical?
  bool operator==(const NewSpecies & other) const{
    return the_unique_id==other.the_unique_id;
  }
  /// Generate a mutated version of this species.
  /** The parameter \a slowdown can be used to make all mutation steps
      smaller than the default, for example to search for optimal
      trait values in a given environment.*/
  void mutate_selected_traits(trait_selection_t selection,
			      double slowdown=1);
  void sample_F_randomly();
  /// Sample species to sample random traits:
  static const species_list_t * sampling_species_list;
  NewSpecies * speciate_as_new(double slowdown=1) const;
  /// Same as speciate_as_new(double), but assigns a clade_id.
  NewSpecies * speciate(double slowdown=1);
private:
  void initiate_biomass();
public:
  double get_threshold_biomass_B()const;
  double get_trade_off_multiplier()const;
  /// Component of niche overlap used by log_availability_structure().
  double log_niche_matching(const NewSpecies & prey) const;  
  /// Alternative to log_niche_matching(..).
  /** Called if totally_random_link_strength_sigma is != 0; */
  double pseudo_standard_normal_random_with(const NewSpecies & other,
					    const int choice=0) const;
  double log_random_matching(const NewSpecies & prey) const;  
  /// Component of niche overlap used by log_availability_structure().
  double log_small_prey_matching(const NewSpecies & prey) const;  
  /// Component of niche overlap used by log_availability_structure().
  double log_large_prey_matching(const NewSpecies & prey) const;  
  /// Component of niche overlap used by log_availability_structure().
  double log_mass_matching(const NewSpecies & prey) const;  
  /// Component of niche overlap used by log_availability_structure().
  double log_other_matching(const NewSpecies & prey) const;  
  /// Determines trophic affinities c to other species.
  double log_foraging_strength_c_on(const NewSpecies & prey) const;  
  /// Determines trophic affinities c to other species.
  double foraging_strength_c_on(const NewSpecies & prey) const;
  /// Generate a uniform "random number" between 0 and 1 determined by
  /// the species pair
  double random_with(const NewSpecies & other,int which=0) const;
  /// Generate a standard-normal "random number" between 0 and 1
  /// determined by the species pair
  double random_normal_with(const NewSpecies & other,int which=0) const;
  /// Determines interaction with other plants.
  double plant_hampering_by(const NewSpecies & prey) const;
  /// Distance to prey in trait space.
  double foraging_distance_to(const NewSpecies & prey) const;
  /// Determines switching similarity between two prey species.
  /** Switching similarities are assumed to be independent of the
      consumer. */
  double switching_similarity_to(const NewSpecies & other,
				 const double width) const;
  inline double prey_similarity_width() const;
  /// Returns total biomass of the species.
  double biomass_abundance_B() const{
    return the_biomass_abundance_B;
  }
  double get_fishing_mortality()const{
    return the_fishing_mortality;
  }
  void set_fishing_mortality(double m){
    the_fishing_mortality=m;
  }
  void set_biomass_abundance_B(double B){
    the_biomass_abundance_B=B;
  }
  void set_aggressivity_g(double g){
    the_attack_rate_a*=g/the_aggressivity_g;
    the_aggressivity_g=g;
  }
  void set_niche_width_wt(double wt){
    the_niche_width_wt=wt;
  }
  void set_conversion_efficiency_eps(double eps){
    the_conversion_efficiency_epsilon=eps;
  }
  void set_plant_growth_rate_sigma(double sig){
    the_plant_growth_rate_sigma=sig;
  }


  /// Returns a characteristic mean body mass of the species.
  /** This is sometimes assumed to be the mean body mass.*/
  double bodymass() const{
    return the_mean_bodymass_M;
  }
  /// Returns the species' unique id.
  int unique_id() const{
    return the_unique_id;
  }
  /// Returns the species' unique id.
  int sequential_id() const{
    return the_sequential_id;
  }
  /// Returns the species' parent's unique id.
  int parent_unique_id() const{
    return the_parent_unique_id;
  }
  /// Returns the species' clade id.
  int clade_id() const{
    return the_clade_id;
  }
  taxon_t taxon() const{
    return the_taxon;
  }
  double mean_bodymass_M() const{
    return the_mean_bodymass_M;
  }
  double log_mean_bodymass_M() const{
    return the_log_mean_bodymass_M;
  }

  double conversion_efficiency_eps() const{
    return the_conversion_efficiency_epsilon;
  }
  double attack_rate_a() const{
    return the_attack_rate_a;
  }
  double handling_time_T() const{
    return the_handling_time_T;
  }
  double turnover_rate_r() const{
    return the_turnover_rate_r+the_fishing_mortality;
  }
  double fishing_mortality() const{
    return the_fishing_mortality;
  }
  double fishing_yield() const{
    return the_fishing_mortality*the_biomass_abundance_B;
  }
  double plant_growth_rate_sigma() const{
    return the_plant_growth_rate_sigma;
  }
  double carrying_capacity() const{
    return the_carrying_capacity;
  }
  double aggressivity_g() const{
    return the_aggressivity_g;
  }
  double niche_width_wt() const{
    return the_niche_width_wt;
  }
  double get_loss_rate_over_max_production_rate_r() const{
    return the_loss_rate_over_max_production_rate_r;
  }
  bool increase_loss_rate_by(double rate);

  double env_effect() const{
    return the_env_effect;
  }

  /// Do perturbation of physiological parameters.
  /** The \a factor controlls the amount of perturbation */
  void perturb_physiological_parameters(double factor);

  /// Attempt some predictions based on parameters.
  static void make_predictions();

  // Trait sweeps:
  /// Entry in a trait sweep
  class trait_sweep_sample_t {
    sequence< double > the_trait_value;
    NewSpecies * the_species;
  public:
    const trait_sweep_sample_t & operator=(const trait_sweep_sample_t & other){
      if(the_species) delete the_species;
      the_trait_value=other.the_trait_value;
      the_species=(other.the_species?new NewSpecies(*(other.the_species)):0);
      return *this;
    }
    trait_sweep_sample_t(const trait_sweep_sample_t & other):
      the_trait_value(other.the_trait_value),
      the_species(other.the_species?new NewSpecies(*(other.the_species)):0){
    }
    trait_sweep_sample_t():the_trait_value(),the_species(0){}
    trait_sweep_sample_t(sequence< double > & v,NewSpecies * s):
      the_trait_value(v),the_species(s){}
    ~trait_sweep_sample_t(){if(the_species) delete the_species;}
    sequence< double > trait_values(){return the_trait_value;}
    NewSpecies * species(){return the_species;}
  };
  /// A trait sweep
  typedef std::vector< trait_sweep_sample_t > trait_sweep_t;
  /// Generates a trait sweep
  trait_sweep_t trait_sweep(const char * details) const;

  int trait_dimensions(int i);
  trait_t get_trait_values(int i);
  int notify_speciation(const NewSpecies & s){
    return ++the_number_of_invading_offspring;
  }
  int number_of_invading_offspring() const{
    return the_number_of_invading_offspring;
  }
  int age() const{
    return next_sequential_id-the_sequential_id;
  }
};

// some more inlines:

inline bool NewSpecies::is_a_plant() const{
  return the_taxon==plant;
}

/// A normal dot product between trait vectors.
inline double dot(const NewSpecies::trait_t & x,const NewSpecies::trait_t & y){
  double sum=0;
  for(int i=x.size();i-->0;){
    sum+=x(i)*y(i);
  }
  return sum;
}

inline double NewSpecies::prey_similarity_width() const{
  if(relative_switching_similarity_width){
    return the_niche_width_wt*switching_similarity_width_w_s;
  }else if(adjust_prey_similarity_width){
    return the_niche_width_wt*switching_similarity_width_w_s/
      trophic_niche_width_w_t;
  }else{
    return switching_similarity_width_w_s;
  }
}      

inline bool operator<(const NewSpecies & s1, const NewSpecies & s2){
  return s1.sequential_id() < s2.sequential_id();
}

#endif // _NEWSPECIES_H_
