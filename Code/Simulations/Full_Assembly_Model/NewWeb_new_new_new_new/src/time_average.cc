// -*- mode: c++ -*-
// $Id: time_average.cc 2175 2011-05-30 17:20:36Z axel $

#include "time_average.h"
#include "error.h"
#include "NewWeb.h"
#include "period_cutter.h"
#include "Integrator.h"

static int time_average_max_steps=1<<30;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGINT(time_average_max_steps),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);

time_average::time_average(NewWeb & web, double time_average_max_t){
  int S=web.number_of_species();
  int n=web.number_of_variables();
  ODE_vector start_state(n),derivative(n);

  // make sure all arrays are allocated:
  web.prepare_for_integration();
  web.write_state_to(start_state);
  web.dynamics(start_state,derivative);
  web.compute_biomass_action_products();

  // point to what to integrate:
  sequence< double* > locations;
  int k=0;
  for(int i=S;i-->0;){
    locations[k++]=&web.biomass_B(i).passivated_reference();
    locations[k++]=&web.s(i).the_top_down_strength;
    locations[k++]=&web.s(i).the_saturation_strength;
    locations[k++]=&web.s(i).the_light_strength;
    locations[k++]=&web.s(i).the_GP;
    for(int j=web.number_of_animals();j-->0;){
      locations[k++]=&web.fx(i,j);
      locations[k++]=&web.biomass_action_products(i,j);
    }
  }
  Integrator averager(locations);
  

  //integrating_dynamical_object iweb(web,locations);
  double &t=web.current_time;
  period_cutter_t period_cutter(t); 
  double t_start=t;
  double t_stop=t+time_average_max_t,dtsum=0;
  int nsteps=0;

  { //begin scope ode_state
    ODE_state ode_state(&web); 
    //initialize compute flows:
    web.compute_flows(ode_state,true);
		      
    while((!nsteps || t<t_stop) && nsteps < time_average_max_steps){

      // step integrator:
      int integrator_failure=ode_state.integrate_one_step(t_stop);
      if(integrator_failure){
	FATAL_ERROR("integrator failure");
      }
      nsteps++;

      web.compute_flows(ode_state,false); //called for side effects
      web.compute_biomass_action_products(); //called for side effects
      averager.sample(t);

      // compute distance:
      double s=0;
      for(int i=start_state.size();i-->0;)
	s+=(ode_state[i]-start_state[i])*(ode_state[i]-start_state[i]);
      double dist=sqrt(s);

      // ask period_cutter for interpretations:
      period_cutter.sample(dist,t);
      if(period_cutter.can_predict_end_time()){
	t_stop=period_cutter.predicted_end_time();
      }else if(period_cutter.noisy_steady_state()){
	WARNING("looks like a noisy steady state");
	t_stop=t;
	break;
      }else if(period_cutter.convergence_to_steady_state() &&
	       period_cutter.way_to_convergence() < 0.000003){
	WARNING("looks like convergence to some steady state");
	averager.disable_write_back();
	goto finish;
      }
    }// end of integrator loop
  } // end scope of ode_state
  if(nsteps >= time_average_max_steps || t - t_start > time_average_max_t){
    WARNING("failed to detect periodicity, convergence, or steady state");
  }
  if(period_cutter.can_predict_end_time()){
    WARNING("looks like an oscillation with period " << 
	    period_cutter.period_length() );
  }
  averager*=(1/(t-t_start));
 finish:
  REPORT(nsteps);
  REPORT(t-t_start);
  web.read_state_from(start_state);
  web.compute_flows(start_state); // sets fx as a side effect
  web.fix_fx();
  
}// end scope of averager, write back to web 

time_average::~time_average(){};
