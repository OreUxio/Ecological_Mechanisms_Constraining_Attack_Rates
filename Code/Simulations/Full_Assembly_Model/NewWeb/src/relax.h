// -*- mode: c++ -*-
// $Id: relax.h 2452 2015-04-24 18:21:19Z axel $
#ifndef _RELAX_H_
#define _RELAX_H_

#include <set>
#include "ODE.h"
#include "sequence.h"
#include "error.h"

/// Abstract base class for systems to relax to a dynamic steady state.
class relaxing_dynamical_object : public ODE_dynamical_object{
public:
  typedef std::set<int> species_set_t;
  virtual ~relaxing_dynamical_object();
  // record_for_steady_state can optionally take two arguments:
  virtual void record_for_steady_state(const ODE_vector & state,
				       const ODE_vector & ddt){
    record_for_steady_state();
  };
  virtual void record_for_steady_state(){};
  virtual void steady_state_is_fixed_point(){};
  virtual bool small_values_in(ODE_vector & state,
			       const species_set_t&  conserved)=0;
  virtual species_set_t 
  delete_species_larger_than_exp(const sequence<double> & si,
				 const species_set_t&  conserved)=0;
  virtual species_set_t 
  delete_all_species_with_less_than_one_individual(const species_set_t
						   & conserved)=0;
  virtual bool is_active(int i);
  virtual double active_sum(ODE_vector &v);

  /// Simulates an ODE until it reaches a dynamic steady state.
  /** In particular, this member solves the problem to decide when the
      steady state has been reached.  Special functionalities handle
      extinct species. */
  double relax(double relaxation_time,
	       const species_set_t & newly_inserted=species_set_t(),
	       const bool do_auto_extinguish=true);
  double relax(double relaxation_time,
	       const species_set_t & newly_inserted,
	       species_set_t & deleted,
	       const bool do_auto_extinguish=true);
};

/// Makes one relaxing_dynamical_object out of two.
/** This is useful for comparing approximations of ODE. */
// Much of the clutter below should go into relax.cc!
class combined_relaxing_dynamical_object : public relaxing_dynamical_object{
private:
  relaxing_dynamical_object * const part1;
  relaxing_dynamical_object * const part2;
  void split_conserved(const species_set_t& conserved,
		       species_set_t &set1,
		       species_set_t &set2);
  void merge_into_deleted(const species_set_t &set1,
			  const species_set_t &set2,
			  species_set_t & deleted);
public:
  combined_relaxing_dynamical_object(relaxing_dynamical_object * p1,
				     relaxing_dynamical_object * p2):
    part1(p1),part2(p2){};

  // virtuals from ODE_dynamical_object
  virtual int dynamics(ODE_vector const & state, 
		       ODE_vector & time_derivative){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    const ODE_vector state1(&state[0],n1);
    const ODE_vector state2(&state[n1],n2);
    ODE_vector td1(&time_derivative[0],n1);
    ODE_vector td2(&time_derivative[n1],n2);

    part1->current_time=current_time;
    part2->current_time=current_time;
    
    int i=part1->dynamics(state1,td1);
    if(i){
      return i;
    }else{
      return part2->dynamics(state2,td2);
    }
  }

  virtual void write_state_to(ODE_vector & state) const{
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    ODE_vector state1(&state[0],n1);
    ODE_vector state2(&state[n1],n2);

    part1->write_state_to(state1);
    part2->write_state_to(state2);
  }
    
  virtual void read_state_from(const ODE_vector & state){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    ODE_vector state1(&state[0],n1);
    ODE_vector state2(&state[n1],n2);

    part1->read_state_from(state1);
    part2->read_state_from(state2);
  }
    
  virtual int number_of_variables() const {
    return part1->number_of_variables()+part2->number_of_variables();
  }
  virtual void line_print(ODE_vector const & state,std::ostream &co){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    ODE_vector state1(&state[0],n1);
    ODE_vector state2(&state[n1],n2);
    
    part1->line_print(state1,co);
    part2->line_print(state2,co);
  }

  // virtuals from relaxing_dynamical_object
  virtual void prepare_for_integration(){
    part1->prepare_for_integration();
    part2->prepare_for_integration();
  };
  virtual void record_for_steady_state(){
    part1->record_for_steady_state();
    part2->record_for_steady_state();
  };
  virtual bool small_values_in(ODE_vector & state,
			       const species_set_t&  conserved){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    ODE_vector state1(&state[0],n1);
    ODE_vector state2(&state[n1],n2);
    species_set_t set1,set2;
    split_conserved(conserved,set1,set2);
    bool result=part1->small_values_in(state1,set1)
      || part2->small_values_in(state2,set2);
    
    return result;
  }
  virtual species_set_t 
  delete_species_larger_than_exp(const sequence<double> & si,
				 const species_set_t&  conserved){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    species_set_t set1,set2;
    split_conserved(conserved,set1,set2);
    sequence<double> si1(n1),si2(n2);
    for(int i=n1;i-->0;){
      si1[i]=si[i];
    }
    for(int i=n2;i-->0;){
      si2[i]=si[i+n1];
    }
    
    set1= part1->delete_species_larger_than_exp(si1,set1);
    set2= part2->delete_species_larger_than_exp(si2,set2);

    species_set_t deleted;
    merge_into_deleted(set1,set2,deleted);
    
    return deleted;
  }

  virtual species_set_t 
  delete_all_species_with_less_than_one_individual(const species_set_t
						   & conserved){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    species_set_t set1,set2;
    split_conserved(conserved,set1,set2);
    
    set1= part1->delete_all_species_with_less_than_one_individual(set1);
    set2= part2->delete_all_species_with_less_than_one_individual(set2);

    species_set_t deleted;
    merge_into_deleted(set1,set2,deleted);
    
    return deleted;
  }
    
  virtual bool is_active(int i){
    int n1=part1->number_of_variables();
    return (i<n1 ? part1->is_active(i) : part2->is_active(i-n1));
  }
  virtual double active_sum(ODE_vector &v){
    int n1=part1->number_of_variables();
    int n2=part2->number_of_variables();
    ODE_vector v1(&v[0],n1);
    ODE_vector v2(&v[n1],n2);
    return part1->active_sum(v1)+part2->active_sum(v2);
  }
};

#endif // _RELAX_H_
