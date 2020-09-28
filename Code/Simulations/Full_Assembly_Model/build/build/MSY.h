// -*- mode: c++ -*-
// $Id: MSY.h 3 2005-12-01 07:13:32Z cvsrep $
#ifndef _MSY_H_
#define _MSY_H_

/// Experimental yield optimization interface

// Usage:
// MSY msy=MSY(web,fished);
// msy.yield();
// msy.F();

#include <nlopt.hpp>
#include <algorithm>
#include <string>
#include <iostream>

#include "NewMatrix.h"
#include "evaluate.h"
#include "NewWeb.h"

double equilibrium_yield(const std::vector<double> &x, 
			 std::vector<double> &grad, void* f_data);

double equilibrium_yield(const std::vector<double> &x, 
			 std::vector<double> &grad, void* f_data);

double isres_constraint(const std::vector<double> &x, 
			std::vector<double> &grad, void* f_data);

extern double MSY_maxE;
extern int MSY_penalize;
extern double MSY_decline_threshold;


class MSY_t 
{
private:
  void grid(std::ostream &os,int step,int depth,std::vector<double> logE);
  double saved_constraint;
  friend double 
  isres_constraint(const std::vector<double> &x, 
		   std::vector<double> &grad, void* f_data);
protected:
  const NewWeb & web;
  nlopt::opt opt;
  nlopt::result result;
  double best_Y_on_grid;
  std::vector<double> logE_on_grid;
  std::vector<double> logE;
  std::vector<double> E;
  double msy;
  double max_E;
  sequence<double> unperturbed_spectrum;
  virtual
  double worst_decline(const NewWeb & web);
  double simple_penalization(NewWeb & web,
			     double decline_threshold=
			     MSY_decline_threshold);
  double sustainability_penalization(NewWeb & web);
public:
  virtual 
  int n()=0;

  MSY_t(const NewWeb & w,double maxE=MSY_maxE);
  virtual 
  void find_MSY()=0;

  // The goal function:
  double operator()(const std::vector<double> &logE, 
		    std::vector<double> &grad, 
		    void* f_data,
		    bool penalize=MSY_penalize);

  // Modify web according to parameters
  virtual 
  double operator()(NewWeb & web,
		    const std::vector<double> &logE, 
		    std::vector<double> &grad);
  virtual
  void set_mortalities(NewWeb & web, 
		       const std::vector<double> &E)=0;

  double yield(){
    return opt.last_optimum_value();
  }
  
  void optimally_exploited_web(NewWeb & web);

  void grid(const char * filename,int steps=10);
};

class MSY_species : public MSY_t
{
  MSY_species();
protected:
  std::vector<int> fished;
  virtual
  double worst_decline(const NewWeb & web);
public:
  virtual 
  int n(){return fished.size();}

  MSY_species(const NewWeb & w,std::vector<int> & f);
  virtual 
  void find_MSY();    

  virtual
  void set_mortalities(NewWeb & web, const std::vector<double> &E);
  void get_mortalities(const NewWeb & web, std::vector<double> &E);
};

class MSY_fleets : public MSY_t 
{
  std::vector<std::string> catchability_function;
  std::vector<std::vector<double> > catchability;
  MSY_fleets();
  void compute_catchabilities();
public:
  virtual
  int n(){return catchability.size();}

  MSY_fleets(const NewWeb & w,std::vector<std::string> & cf);
  virtual 
  void find_MSY();

  virtual
  void set_mortalities(NewWeb & web, const std::vector<double> &E);
};


class harvest_controller_t : public MSY_species, public NewWeb
/// We need our own copy of NewWeb web here, because we want to modify it.
{
private:
  bool constructing;
protected:
  NewMatrix the_interaction_matrix;
  NewVector the_production_rate;
  NewVector the_old_B;
  NewVector the_old_F;
  void recompute_LV_approximation();
  NewMatrix
  Gtilde_from_G(const NewMatrix & G);
  NewMatrix
  Gbar_from_G(const NewMatrix & G);
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   )=0;
  void apply_fishing_mortalities(const ODE_vector & state, 
				 ODE_vector & time_derivative);
  void set_fishing_mortalities();
public:
  harvest_controller_t(const MSY_species & msy);
  virtual
  int dynamics(ODE_vector const & state, 
	       ODE_vector & time_derivative);
  virtual
  void recompute_strategy(){};
};

class transposed_interaction_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
public:
  transposed_interaction_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
  };
};

class productive_state_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_B;
public:
  productive_state_controller_t(const MSY_species & msy);
  virtual void recompute_strategy();
};

class target_pressure_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
public:
  target_pressure_controller_t(const MSY_species & msy);
  virtual void recompute_strategy();
};

class individual_productive_state_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  void recompute_target_B();
protected:
  NewVector the_target_B;
public:
  individual_productive_state_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_B();
  };
};

class soft_individual_productive_state_controller_t : 
  public individual_productive_state_controller_t{
  NewVector the_relaxation_rate;
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
public:
  soft_individual_productive_state_controller_t(const MSY_species & msy) :
    individual_productive_state_controller_t(msy) {};
  void recompute_strategy();
};


class individual_target_pressure_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
  void recompute_target_F();
public:
  individual_target_pressure_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_F();
  };
};

class individual_transposed_interaction_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
public:
  individual_transposed_interaction_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
  };
};


class individual_B_productive_state_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_B;
  void recompute_target_B();
public:
  individual_B_productive_state_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_B();
  };
};

class individual_B_target_pressure_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
  void recompute_target_F();
public:
  individual_B_target_pressure_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_F();
  };
};

class individual_B_transposed_interaction_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
public:
  individual_B_transposed_interaction_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
  };
};

class growth_rate_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
  void recompute_target_F();
public:
  growth_rate_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_F();
  };
};

class CFP_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
  void recompute_target_F();
public:
  CFP_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_F();
  };
};

// Data poor fishing:
class DPF_controller_t : public harvest_controller_t
{
  virtual
  void compute_fishing_mortalities(const ODE_vector & state, 
				   const ODE_vector & time_derivative,
				   NewVector & F
				   );
  NewVector the_target_F;
  void recompute_target_F();
public:
  DPF_controller_t(const MSY_species & msy);
  virtual
  void recompute_strategy(){
    recompute_LV_approximation();
    recompute_target_F();
  };
};



#endif // _MSY_H_
