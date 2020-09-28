// -*- mode: c++ -*-
// $Id: NewWeb.h 2505 2017-02-28 11:07:00Z axel $
#ifndef _NEWWEB_H_
#define _NEWWEB_H_


#include <list>
#include <stack>
#include <set>
#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/barrier.hpp>
#include <sched.h>
#include <hash_map>

#include "NetworkAnalysis.h"
#include "ODE.h"
#include "NewSpecies.h"
#include "link_strength_matrix.h"
#include "remember.h"
#include "snapshot.h"
#include "otherwebs.h"
#include "Statistics.h"
#include "SortedVector.h"
#include "SortedMatrix.h"
#include "relax.h"
#define NDEBUG // don't assert(..) things
#include "ptr_vector.h"
#include "CompMatrix.h"

/// Maximal number of allowed threads
extern int max_num_threads; 

/// Maximum size of stady state to use:
extern int max_steady_state_size; ///< 0 mean don't use

/// Strength of invasion pressure from other webs (obsolet).
extern double unified_invasion_pressure_kappa;

/// Strength of invasion pressure from other webs.
/**
   Probability by which a mutant of a species from another web is more
   likely to be added than an extant species.
*/
extern double any_invasion_pressure_q;

/// Implement switching?
extern int do_switching;

/// Use experimental optimized code?
extern int use_packed_simulation;

/// Write some more info about evolutionary process to stdout
extern int log_steady_state_stats;

#ifdef _GNU_SOURCE
#define SET_CPU_AFFINITY
#endif

#ifdef SET_CPU_AFFINITY
/// Linking threads to specific cores is supposed to reduce pages faults.
typedef std::vector<cpu_set_t> CPU_set_list_t;
extern CPU_set_list_t CPU_set_list;
#endif

/// Food web of the Population Dynamical Matching Model.
class NewWeb : public relaxing_dynamical_object, public permanent
{
public:
  /// Data used to compute population dynamics.
  /** Must be computed and updated upon insertion and deletion of
      species. */
  class precomputed_entry_t{
  public:
    /// Speed up prey switching calculations.
    /** For the consumer this matrix associated with, represents a
	product M_ij = c_i s_ij c_j, where i and j are two prey
	species, c_i and c_j are trophic link strengths to these
	species, and s_ij is the switching similarity between i and j
	in the eyes of the consumer.  Multiplying this matrix from the
	left and the right with vectors of abundances yields a number
	that can be used to calculate the total intake of the
	consumer.  See the draft manuscript vanLeeuwen2008.pdf and
	hotspots.cc for details.
    */
    SortedMatrix csc_eating; 
    /// Speed up prey switching calculations.
    /** For the resource species i this matrix is associated with,
	represents a matrix M_jk = c_ik s_ijk c_jk, where k is a
	consumer species, and j is another resource of this consumer,
	c_ik and c_jk are trophic link strengths to resources of k, and
	s_ijk represents the switching similarity between i and j in
	the eyes of k.  These matrices can be used to quickly compute
	the total biomass losses of i by consumption. See the draft
	manuscript vanLeeuwen2008.pdf and hotspots.cc for details.
    */
    SortedMatrix csc_being_eaten;
    /// Vector of trophic interaction strengths to resources.
    SortedVector c;
    /// Vector of trophic interaction strengths from consumers.
    /** Used only when not switching. */
    SortedVector cT;
    friend bool operator==(const precomputed_entry_t &,
			   const precomputed_entry_t &);
    precomputed_entry_t();
  };
  /// Compute a precomputed_entry_t for a test species.
  /** While \a pr is compute, the \a precomputed member of the \a
      NewWeb is not updated.  As a result, only effects on \a sp can
      be obtained from this, but not effects of \a sp on other species.
      Used to compute fitness.
  */
  void precompute(NewSpecies * sp,precomputed_entry_t & pr);

  typedef stdx::ptr_vector< precomputed_entry_t > precomputed_container_t;
  /// Stores and manages precomputed entries.
  class precomputed_t:public precomputed_container_t{
  public:
    precomputed_t():precomputed_container_t(){};
    precomputed_t(const precomputed_t & other);
    ~precomputed_t();
    precomputed_t & operator=(const precomputed_t & other);
    /// Re-compute all entries.
    /** The target species is \a k.*/ 
    void update(int k,species_list_t &s);
    /// Update entries when re-indexing species \a j as \a k.
    void move(int j,int k,species_list_t &s);
    /// Resize all members.
    void resize(int size,species_list_t &s);
  };
  precomputed_t precomputed;
  void initialize_precomputed();///< fast initialization of precomputed
    
  int speciation_counter;
  int invasion_counter;
  int max_ever_plant_richness;
  // the remaining data members are just helper variables
public:
  /// Is species \a i a plant species?
  int is_plant(int i) const {return i>=s.n_animals;}
  link_strength_matrix fx; ///< fx[i][j]=fx(j,i) (!!!!) is the flow from j to i
  /// Dynamically computed biomasses.
  /** Valid only after calling dynamics(...)!*/
  vector_with_max biomass_B; 
  /// Dynamically computed common factors in functional responses.
  /** In many formulae for functional responses, large parts are the
      same for all prey species.  The \a common_factor for a consumer
      species is this part times the consumer species' biomass. Valid
      only after calling dynamics(...)!*/
  vector_with_max common_factor; 
private:
  /// Data needed to easily replay dynamics for fitness calculations.
  class saved_state{
  public:
    double t;
    vector_with_max biomass_B;
    vector_with_max common_factor;
    saved_state();
    saved_state(const saved_state & s); ///< Copy constructor.
    /// Constructor recording a state.
    saved_state(double t,const vector_with_max & B,const vector_with_max &cf);
  };
  /// Represents a recorded steady state.
  class steady_state_t : public 
  stdx::ptr_vector< saved_state > {
  public:
    void clear(){
      while(begin()!=end()){
	delete detach(end()-1);
      }
    }
    ~steady_state_t(){clear();}
  };

  /// The last simulated steady state.
  steady_state_t steady_state;
  /// Records steady state during relaxation.
  virtual void record_for_steady_state();
  virtual void steady_state_is_fixed_point();

public:
  species_list_t s; ///< Vector of all extant species.
  void check_internal_consistency(int check_data_consistency);///< Debugging.
  // When species are added and deleted, indices of other species can
  // be shuffled around.  In order to be able to related to species
  // across addition and deletion events (mainly for plotting
  // population time series), each species is assigned a column.
  // Unused columns can later be reused.
  sequence<int> species_at; ///< Index of species assigned to column \a i.
  static const int column_unused; ///< Value representing an unused column.
  sequence<int> assigned_column; ///< Column a species \a i is assigned to. 
  static const int parent_not_extant; ///< Special species index.
  Otherwebs otherwebs;///< Other webs this one exchanges species with.

private:
  struct ancestral_record {
    int column;
    NewSpecies extinct_ancestor;
    int number_of_known_offspring;
    ancestral_record(){};
    ancestral_record(const ancestral_record & a):
      column(a.column), extinct_ancestor(a.extinct_ancestor),
      number_of_known_offspring(a.number_of_known_offspring) {};
    explicit ancestral_record(const NewSpecies & S) :
      column(parent_not_extant), extinct_ancestor(S),
      number_of_known_offspring(0) {};
    explicit ancestral_record(int c) :
      column(c),number_of_known_offspring(0) {};
  } ;
  class Ancestry : 
    public __gnu_cxx::hash_map< int, ancestral_record > {
  public:
    void note_speciation(const NewSpecies & offspring,const NewWeb & w);
    void note_impending_extinction(const NewSpecies & going_extinct);
  };
  
  Ancestry ancestors;
  // Manage column assignment:
private:
  int column_allocate(); ///< Get an unused column index.
  void column_free(int i); ///< Mark column \a i as unused.

  // Manage species indexing:
  /// Get a species index.
  /** Argument indicates if index is for a plant species. A newly
      added species can be assigned this index.  Then precomputed data
      needs to be updated.  Usually, you want to use the high-level
      function insert(...) instead.*/
  int get_unused_index_and_allocate_memory(int plant);
  /// Move all data pertaining to species \a old_i to position \a new_i.
  void move_species(int old_i,int new_i); 
  /// High-level species insertion.
  /** Inserts the species pointed to by \a sp into the list of extant
      species \a s, updates precomputed data and does other kinds of
      book-keeping.  The memory \a sp points to is freed.  Returns
      index assigned to species. */
public:
  int insert(const NewSpecies * sp); 
protected:
  /// High-level species deletion.
  /** This should be the direct opposite of insert(const NewSpecies * sp). 
      Return value is 0 if i was last species in s, 1 otherwise. */
  int delete_species(int i); 
  int delete_species_and_report(int i);
  /// Helper to find the parent of a species among extant species:
  int find_extant_parent(const NewSpecies & sp);

public:
  void delete_forbidden_species(); 

private:
  /// Manage permanent data.
  /** For saving to and loading from a file. */
  virtual void data_mapping();

public:
  /// Effective values of B_i*B_j, corrected for switching.
  link_strength_matrix biomass_action_products;
  /// Compute biomass_action_products (defunct).
  void compute_biomass_action_products();

  NewWeb(); ///< Construct an empty food web.
  virtual ~NewWeb();
  NewWeb(const NewWeb&);
  NewWeb& operator=(const NewWeb&);

  void clear(); ///< Delete all species.

  // Virtual functions from abstract base ODE_dynamical_object:
  int number_of_variables() const; 
  void write_state_to(ODE_vector &) const;
  void read_state_from(const ODE_vector &);

  /// Updates the matrix fx of biomass flows between species.
  /** Returns 0 if all went well, another value otherwise. */
  int compute_flows(ODE_vector const & state,bool initialize=true);
private:
  typedef std::vector< SortedVector > compute_flows_csc_matrix_t;
  std::vector< compute_flows_csc_matrix_t > compute_flows_csc;
public:
  void fix_fx(); ///< Make sure there are no flows to plants.
  double slowness(int i) const {return s(i).slowness();}
  ///< Returns a scaling factor of dimension TIME representative of
  ///population dynamics of species \a i.

  ///helpers for multi-threading;
private:
  class do_dynamics_manager_t{
  public:
    do_dynamics_manager_t();
    do_dynamics_manager_t(const do_dynamics_manager_t & other);
    do_dynamics_manager_t& operator=(const do_dynamics_manager_t& other);
    ~do_dynamics_manager_t();
    void initialize_threads_maybe(NewWeb * base_class);
    void start_threads(int nthreads,NewWeb * base_class);
    void stop_threads();
    int get_num_threads(){return _num_threads;}
    boost::thread_group  threads;
    struct task_t {
      NewWeb * dispatcher; 
      int index;
      void operator()();
    };
    std::vector< task_t > task;
    boost::mutex mutex;
    boost::condition condition;
    boost::barrier * barrier;
    std::vector< double > common_factor_max;
    static const double unset; // some value < 0;
    bool max_not_set(int i){return common_factor_max[i]<0;}
    const ODE_vector * state;
    ODE_vector * time_derivative;
    bool stop_now;
    int threads_ready;
  private:
    int _num_threads;
  };
  do_dynamics_manager_t do_dynamics_manager;
  friend class packed_simulation;
public:
  class packed_simulation * dynamics_accelerator;
  void prepare_for_integration();
  void cleanup_after_integration();
  void do_dynamics(int offset, 
		   ODE_vector const & state, 
		   ODE_vector & time_derivative);
  void do_eating(int begin, int end, ODE_vector & time_derivative);
  void do_plants(int begin, int end, ODE_vector & time_derivative);
  void do_being_eaten(int begin, int end, ODE_vector & time_derivative);
public:
  // More virtual functions from abstract base ODE_dynamical_object:
  /// Compute the time derivative given a \a state of log biomasses.
  int dynamics(ODE_vector const & state, 
	       ODE_vector & time_derivative);
  /// Experimental preconditioner for cvode ODE solver.
  void precondition(ODE_vector const & state,
		    ODE_vector const & in,
		    ODE_vector & out,
		    realtype gamma,
		    bool left_rather_than_right);
  bool has_preconditioner(){return true;};
  bool has_inherent_rates(){return true;};
  void get_inherent_rates(ODE_vector & rates);
  /// Print log_10 abundances of species in their assigned columns.
  /** Used to generate abundance time series */
  void line_print(ODE_vector const & state, std::ostream &co);

  /// Quick test if state contains large negative values.
  /** These would represent species that should be removed as extinct. */
  bool small_values_in(ODE_vector & state,
		       const species_set_t &conserved=species_set_t());
  /// Generate precomputed data from scratch.
  void initialize_for_integration();
  int number_of_species() const; ///< Number of species in this web.
  int number_of_plants() const; ///< Number of plant species in this web.
  int number_of_animals() const; ///< Number of animal species in this web.
  int number_of_fish() const; ///< Number of fish species in this web.
  double animal_biomass() const; ///< Total biomass of animals in web.
  double plant_biomass() const; ///< Total biomass of plants in web.
  double 
  biomass_larger(double threshold) const; ///< Total biomass of large things.

  // Miscellaneous:
  /// As dynamics(), but make sure old version is used to get all its
  /// side effects.
  int dynamics_for_side_effects(ODE_vector const & state, 
				ODE_vector & time_derivative);
  void delete_rare_fish_species(double percentage);
  void delete_rare_fish_species_5per_2(double slope);
  void delete_rare_fish_species_10per_2(double slope);
  void delete_random_animal_species();
  void delete_fish_species_at_random(double x);
  void delete_fish_species_by_size(double x);
  void delete_fish_species_by_size_reverse(double x);
  void delete_fish_species_by_biomass(double x);
  void delete_fish_species_by_biomass_reverse(double x);
  void delete_fish_species_by_TL(double x);
  void delete_fish_species_by_TL_reverse(double x);
  void delete_fish_species_by_Conn(double x);
  void delete_fish_species_by_Conn_reverse(double x);
  void delete_all_species_with_biomass_less_than(double m);
  void delete_fish_species_except_1(double x);
  /// This does what it says, too.  
  /** For \a conserved species the biomass has to be by a large value
      (\a conservation_tolerance) smaller than the body mass for
      deletion.  Returns set of deleted species.*/
  species_set_t 
  delete_all_species_with_less_than_one_individual(const species_set_t
						   &conserved=species_set_t() );
  /// Deletes species with biomass less than one individual.
  /** However, contrary to
      delete_all_species_with_less_than_one_individual(const
      species_set_t &conserved), \a ln \a biomass abundances are here
      passed as the \a si argument. For \a conserved species the
      biomass has to be by a large value (\a conservation_tolerance)
      smaller than the body mass for deletion.  Returns set of deleted
      species. */
  species_set_t delete_species_larger_than_exp(const sequence<double> & si,
					       const species_set_t
					       &conserved=species_set_t() );
  /// Delete species by trophic level.
  void delete_all_species_with_level_higher_than(double h);
  /// Delete all plants with negative growth rates
  void delete_all_plants_with_negative_growth_rates(species_set_t& inserted);
  /// Scale primary productivity
  void scale_primary_productivity(double x);
  /// See if there are plants larger than the largest animal
  bool large_plants();
  void forbid_fishing();
  void copy_fishing_mortalities_from(const NewWeb & web);
  /// Increase the mortality of a large abundant species by r.
  void fish(double r);
  void fishF(double r);
  void fishall(double r);
  void fishMrange(double r);
  void fishTLrange(double r);
  void fishMrangeconst(double r);
  void fishMrangeconst2(double r,double s);
  void fishTLrangeconst(double r);
  void fishMrangeline(double r);
  void fishTLrangeline(double r);
  void fishMrangeCS(double r);
  /// Delete proportion \a r of all plants.
  void delete_plant_fraction(double r);
  /// Delete species in upper \a r part of body mass range covered in web.
  void delete_large_species_fraction(double r);
  void delete_10_species_near_size(double size);
  void delete_plants_by_vulnerability_trait(double r);
  void delete_all_animals(); ///< Dito.
  void custom_perturbation(double x); ///< undefined
  void delete_rare_species_fraction(double r);
  void delete_lower_fraction_of_property(const sequence<double> & property,
					 double r);
  void delete_weeds();
  /// Reduce trophic interaction matrix to large-eats-small links. 
  Interaction_Matrix trivial_loop_removal(const Interaction_Matrix & im);

  // Fitness and species addition:

  /// Return mutant of an given species. 
  /** The argument \a parent is internal index of this species. */
  NewSpecies* speciate_selected(double parent);
  /// Return mutant of an extant species chosen at random. 
  /** The argument \a plant_fraction specifies the probability of
      choosing a plant.  It is prespected only if
      unified_invasion_pressure_kappa!=0. */
  NewSpecies* speciate(double plant_fraction=NewSpecies::plant_fraction);
  /// Return mutant of an extant species chosen with probability
  /// proportional to biomass.
  NewSpecies* speciate_by_biomass(int get_a_plant);

  /// Linear growth rate of test species \a sp in momentary state.
  /** Species \a sp is not yet inserted into web.  Assumes
      common_factor and biomass_B to be set! */
  double linear_fitness(NewSpecies * sp);
  /// Fast but numerically slighly different version of linear_fitness(sp).
  /** Species \a sp is not yet inserted into web.  Assumes
      common_factor and biomass_B to be set! */
  double fast_linear_fitness(NewSpecies * sp);
  /// Compute invasion fitness relative to saved steady_state.
  /** Species \a sp is not yet inserted into web. */
  double invasion_fitness(NewSpecies * sp);
private:
  /// Low-level function to compute linear growth rate of a test species.
  double linear_fitness(NewSpecies * sp,
			const vector_with_max &biomass_B,
			const vector_with_max &common_factor,
			const precomputed_entry_t & pr );
public:
  /// Momentary d ln B(t)/dt of species \a test temporarily inserted
  /// into web for this computation (obsolete).
  double fitness(NewSpecies & test);
  /// Momentary d ln B(t)/dt of an extant species (obsolete).
  double fitness(int test);

  /// Return a new species as a possible invader.
  /** If mean_field_threshold!=0, a "typical" species is possibly
      guessed from the overall distribution of the extant species (may
      not be working).  Otherwise, depending on the value of \a
      n_neutral and the status of \a otherwebs, either a mutant of a
      species from another web, or a species from a "neutral"
      distribution over niche space is returned.*/
  NewSpecies* invade(double plant_fraction);
  /// Return a new species as a possible invader.
  /** Returns a mutant of a species from \a otherwebs picked with
      probability proportional to biomass in otherwebs.  If no \a
      otherwebs are available, a random sepies is returned. */      
  NewSpecies* invade_by_biomass(double get_a_plant);
  /// Returns a random species obtained either through invade(double)
  /// or through speciate(double) (obsolete).
  /** Details depend on unified_invasion_pressure_kappa and other
      constants.*/
  NewSpecies* 
  invade_or_speciate(double plant_fraction=NewSpecies::plant_fraction);
  /// Returns a random species obtained either through invade(double)
  /// or through speciate(double).
  /** any_invasion_pressure_q controls what is done. */
  NewSpecies* invade_or_speciate_any(double dummy=0);
  /// Returns a random species obtained either through invade(double)
  /// or through speciate(double).
  /** any_invasion_pressure_q controls what is done. */
  NewSpecies* invade_or_speciate_any_by_biomass(double dummy=0);
  /// This macro iterates attempted additions until a species with
  /// positive invasion fitness is found.
  /** Fit candidates are inserted into the web, their index is
      returned. */
  typedef NewSpecies* (NewWeb::*finder_function_t)(double);
private:
  NewSpecies * find_only_fit(double plant_fraction, finder_function_t f);
public:
  NewSpecies * invade_only_fit(double plant_fraction);
  NewSpecies * speciate_selected_only_fit(double plant_fraction);
  NewSpecies * speciate_only_fit(double plant_fraction);
  NewSpecies * invade_or_speciate_only_fit(double plant_fraction);
  NewSpecies * invade_or_speciate_any_only_fit(double plant_fraction);
  NewSpecies * invade_or_speciate_any_by_biomass_only_fit(double plant_fraction);
//   /// This macro iterates attempted additions until a species with
//   /// positive invasion fitness is found.
// #define ONLY_FIT_AFTER_ANNEALING_DECLARE(X) int X##_only_fit_after_annealing(double plant_fraction)
//   ONLY_FIT_AFTER_ANNEALING_DECLARE(invade);
//   ONLY_FIT_AFTER_ANNEALING_DECLARE(speciate);
//   ONLY_FIT_AFTER_ANNEALING_DECLARE(invade_or_speciate);
// #undef ONLY_FIT_AFTER_ANNEALING_DECLARE
  // Outdated protocol for species addition.
  int invade_only_fit_or_speciate_only_fit
      (double plant_fraction=NewSpecies::plant_fraction);
  // Outdated protocol for species addition.
  int invade_only_fit_or_speciate_only_fit_after_annealing
      (double plant_fraction=NewSpecies::plant_fraction);
  // Outdated protocol for species addition.
  int invade_fit_or_speciate_fit
      (double plant_fraction=NewSpecies::plant_fraction);

  // Diagnostics:

  /// Returns precomputed matrix of trophic biomass flows.
  const link_strength_matrix & get_intake_matrix();
  /// Computes matrix of trophic biomass flows.
  link_strength_matrix intake_matrix(bool initialize=true);
  /// Compute diet of one extant predator species
  void get_diet(int predator,sequence< double > & diet);
  /// Computes matrix of trophic biomass flows, divided by resources
  /// turnover rate.
  link_strength_matrix scaled_intake_matrix(bool initialize=true);
  /// Compute matrix of trophic affinities. 
  link_strength_matrix matching_matrix();
  /// Enforce equilibrium by adjusing loss rates:
  void enforce_equilibrium();
  /// Jacobian at current state for log biomasses.
  link_strength_matrix numerical_Jacobian(double dx=0.001);
  /// Write a rank-abundances plot to \a filename.
  void rank_abundance_plot(const char * filename,bool with_plants);
  /// Write a size spectrum to \a filename.
  void size_spectrum(const char * filename,
		     bool with_plants, bool with_animals);
  /// Write a biomass spectrum to \a filename.
  void biomass_spectrum(const char * filename,
			bool with_plants, bool with_animals);
  /// Write various information on each species to \a filename
  /** Makes use of information in a Snapshop that was computed from this
      web before. */
  void species_table(const char * filename,bool with_plants,double threshold,
		     Snapshot & snap, bool jolly=false);
  void species_table2(const char * filename,bool with_plants,double threshold,
		     Snapshot & snap, link_strength_matrix & flows, int colno, bool jolly=false);
  void species_table3(const char * filename,bool with_plants,double threshold,
		     Snapshot & snap, link_strength_matrix & flows, bool jolly=false);
  /// Write various information on most important trophic links to \a
  /// filename
  void link_table(const char * filename,double threshold,
 		  link_strength_matrix & flows,
 		  const sequence<double> & biomass_B);
  void link_table2(const char * filename,double threshold,
 		  link_strength_matrix & flows,
 		  const sequence<double> & biomass_B);
  /// Generate matrix of distances in trophic niche space.
  link_strength_matrix trophic_distance_matrix();
  int get_invasion_counter();
  int get_speciation_counter();
  int get_max_ever_plant_richness();

  /// Compute momentary strengths top-down effects (defunct).
  /** Results are stored as parts of NewSpecies classes. */
  void get_top_down_strength();
  /// Compute Jacobian and write its eigenvalues to \a filename.
  void eigenvalue_list(const char * filename,double dx=0);
  /// Return vector of biomasses of all species.
  sequence<double> get_biomass_B() const;
  /// Return vector of body masses of all species.
  sequence<double> get_bodymass_M() const;
  /// Return vector of log body masses of all species.
  sequence<double> get_log_bodymass_M() const;
  /// Returns a vector flagging species that are plants.
  sequence<int> get_is_plant();
  /// Print a bunch of stats that are easily computed.
  void fast_stats();
  /// Whether or not to compute additional diagnostics while simulating.
  bool compute_diagnostics;
  /// Write histograms of distances in niche space to [ap]<\a filename>.
  /** Used to test hypothesis regarding the fractal distribution of
      traits in niche space */
  void distance_histogram(const char * filename);
  /// Write lists of distances in niche space to  *<\a filename>.
  void distance_lists(const char * filename);
  /// Analyze contributions to effective niche overlap (defunct).
  multinormal_distribution log_availability_structure();
  /// Meaning of the components in output of log_availability_structure().
  typedef enum 
    {log_av_large_prey, ///< cutoff for prey to large
     log_av_small_prey, ///< cutoff for prey to small
     log_av_niche,      ///< overlap in abstract niche space
     log_av_biomass,    ///< variation in prey biomass
     log_av_temporal,   ///< correlation of predator/prey abundances
     log_av_other,      ///< other contributions
     log_av_n_components} log_av_components;
  /// Perturb physiological parameters to model environmental change.
  void model_environmenatal_change(double strength);
  /// Compute the popular Living Plante Index of biodiversity decline.
  double Living_Planet_Index(NewWeb & baseline,
			     bool with_plants, bool with_animals);
  /// Compute the Biodiversity Intentness Index.
  double Biodiversity_Intactness_Index(NewWeb & baseline,
				       bool with_plants, bool with_animals);
  /// Returns the gross primary production of web.
  /** Requires GPP of each species to be computed first, by setting
      compute_diagnostics==true, simulating, and possibly averaging.*/
  double GPP();
  void report_evolutionary_pressures(int iterations=100);
  void invasion_fitness_curve(const char * details);
  void invasion_fitness_relief(const char * details);
  void invasion_fitness_samples(int n);
  void write_link_coefficient_matrix(const char * filename);
  void write_competition_matrix(const char * filename);
  double biomass_yield();
};


#ifdef SET_CPU_AFFINITY 
// minor support stuff
int number_of_CPUs();
CPU_set_list_t find_my_CPUs();
#endif

#endif // _NEWWEB_H_

