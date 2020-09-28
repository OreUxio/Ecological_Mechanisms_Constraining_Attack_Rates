// -*- mode: c++ -*-
// $Id: relax.cc 2452 2015-04-24 18:21:19Z axel $

#include <fstream>
#include <algorithm>
#include <gsl/gsl_cdf.h>

#include "relax.h"
#include "polyfit.h"
#include "random.h"
#include "evaluate.h"


// #define RELAX_DEBUG

my_evaluator_t eval_now;
static int min_samples=5;
static double evenness_pickiness=1.5;
static double extrema_tolerance=1;
static double extrema_relaxer=1;
static double extinction_test_p=1e-3;
static double efficacious_amount_drifted_up=0.1;
static double fixed_point_wobble_tolerance=0.01;
static double waiting_r=0.5;
static int drift_accuracy_max_points=min_samples;
const double golden_ratio=0.5*(1+sqrt(5.0));
static double dephasing_time=eval_now("1*year")/golden_ratio;
static int record_relaxation_dynamics=true;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  {
    CFGINT(record_relaxation_dynamics),
    CFGDOUBLE(fixed_point_wobble_tolerance),
    CFGDOUBLE(extinction_test_p),
    CFGDOUBLE(waiting_r),
    {0, CFG_END, 0}
  };
static cfg_add cfg_dummy(cfg);


// Compute correlation between x and y
static double dummy;
static double correlation_coefficient(const double * x,const double * y,int N,
				      double & amount_drifted=dummy){
  double sum_sq_x = 0;
  double sum_sq_y = 0;
  double sum_coproduct = 0;
  double mean_x = x[0];
  double mean_y = y[0];
  for(int i=1;i<N;i++){
    double sweep = i/(i+1.0);
    double delta_x = x[i] - mean_x;
    double delta_y = y[i] - mean_y;
    sum_sq_x += delta_x * delta_x * sweep;
    sum_sq_y += delta_y * delta_y * sweep;
    sum_coproduct += delta_x * delta_y * sweep;
    mean_x += delta_x / i;
    mean_y += delta_y / i;
  }
  const double pop_sd_x = sqrt( sum_sq_x / N );
  const double pop_sd_y = sqrt( sum_sq_y / N );
  const double cov_x_y = sum_coproduct / N;
  const double correlation = cov_x_y / (pop_sd_x * pop_sd_y);
  amount_drifted=correlation*pop_sd_y;
  return correlation;
}

// Compute relative accuracy of drift coefficient in linear
// regression:
static double drift_accuracy(const double * x,const double * y,int N){
  int n_drift_points=std::min<int>(drift_accuracy_max_points,N);
  double step_size=N/double(n_drift_points);
  double square_sum=0;
  double drift=(y[N-1]-y[0])/(x[N-1]-x[0]);
  for(int i=1;i<n_drift_points;i++){
    const double delta_x=
      x[int(i*step_size)]-x[int((i-1)*step_size)];
    const double delta_y=
      y[int(i*step_size)]-y[int((i-1)*step_size)]-drift*delta_x;
    square_sum+=delta_y*delta_y/delta_x;
  }
  return (y[N-1]-y[0])/sqrt(square_sum/(N-2)*(x[N-1]-x[0]));
}

// Characterize shape of curve y(x):
typedef enum {linear_drift,exponential_relaxation,complex_shape} curve_shape_t;

static curve_shape_t shape(const double * x,const double * y,int preN,
		    double &fixed_point){
  const int N=std::min<int>(preN,20);
  sequence<double> xx(N);
  sequence<average_meter> y_av(N);
  const double x0=x[0];
  const int overloads=preN/N;
  const int rest=preN%N;
  
  // Resample on a new grid: distribute gridpoints x over N data
  // points, and get average of y around these points:
  int ii=0;
  for(int i=0;i<N;i++){
    average_meter xi;
    for(int j=0;j<overloads+(i<N-rest?0:1);j++){
      xi.sample(x[ii]-x0);
      y_av(i).sample(y[ii]+1);
      y_av(i).sample(y[ii]-1);
      ii++;
    }
    xx(i)=xi.readout();
  }
  // Fit y on grid:
  fitted_function f(xx,y_av);

  if(f.polynomial_order()<=1){
#ifdef RELAX_DEBUG
    WARNING("linear_drift (1)");
#endif
    return linear_drift;
  }

//   if(f(0)-f(x[preN-1]-x0) > log(10)){
// #ifdef RELAX_DEBUG
//     WARNING("linear_drift (2)");
// #endif
//     return linear_drift;
//   }
  fitted_function df=f.derivative();
#ifdef RELAX_DEBUG
  REPORT(f.polynomial_order());
  REPORT(f[0]);
  REPORT(f[1]);
  REPORT(f[2]);
  REPORT(df[0]);
  REPORT(df[1]);
  REPORT(df[2]);
#endif
  // fit a straight line to the derivative, using an evenly spaced
  // grid:
  const int n_grid_points=N;
  const double delta_x=(x[preN-1]-x0)/n_grid_points;
  sequence<average_meter> logdf_av(n_grid_points);
  xx.resize(n_grid_points);
  int df_initial_sign= ( df(0)>0 ? +1 : -1 );
#if 1
  // This usually not happen, as we only look at downward going
  // curves:
  if(df_initial_sign == +1){
#ifdef RELAX_DEBUG
    WARNING("complex_shape (1)");
#endif
    return complex_shape;
  }
#endif
  for(int i=n_grid_points;i-->0;){
    // generate evenly space grid:
    xx(i)=i*delta_x;
    const double dfx=df_initial_sign*df(xx(i));
    if(dfx<0) {
#ifdef RELAX_DEBUG
      REPORT(dfx);
#endif
#ifdef RELAX_DEBUG
    WARNING("complex_shape (2)");
#endif
      return complex_shape;
    }
    // Estimated error in estimation of derivative df/dx (no idea
    // where this formula came from). The fabs around the var_at
    // expression is needed because var_at is sometimes numerically
    // negative.
    double Delta_dfx=sqrt(fabs(df.var_at(xx(i))/df.var_at(xx(0))))*
      (1+fabs(df(xx(0)))+fabs(df(xx(n_grid_points-1))));
    //REPORT(Delta_dfx);
    if(Delta_dfx<=0) Delta_dfx=1;
    // Generate log(-df/dx) time series with corresponding sampling
    // errors:
    logdf_av(i).sample(log(dfx)+Delta_dfx/dfx);
    logdf_av(i).sample(log(dfx)-Delta_dfx/dfx);
  }
#ifdef RELAX_DEBUG
  REPORT(logdf_av);
#endif
  // Fit a line to log(-df/dx):
  fitted_function linear_fit(xx,logdf_av,2);
  
  // Test if data is compatible with an exponential relaxation:
  const double z=2;
  for(int i=n_grid_points;i-->0;){
    const double residual=(exp(linear_fit(xx(i)))-df_initial_sign*df(xx(i)));
#ifdef RELAX_DEBUG
    REPORT(residual/sqrt(df.var_at(xx(i))));
#endif
    // does this look like an exponential decay?
    if(fabs(residual/sqrt(df.var_at(xx(i)))) > z){
#ifdef RELAX_DEBUG
    WARNING("complex_shape (3)");
#endif
      return complex_shape;
    }
  }
  
  // looks like an exponential_relaxation, get the fixed point:
  const double decay_rate=-linear_fit[1];
#ifdef RELAX_DEBUG
  REPORT(decay_rate);
#endif
  if(decay_rate<-0.1*fabs(linear_fit[0]/(xx(n_grid_points-1)-xx(0)))){
#ifdef RELAX_DEBUG
    WARNING("complex_shape (4)");
#endif
    return complex_shape;
  }
  const double pre_factor=df_initial_sign*exp(linear_fit[0])/decay_rate;
  
  // make an estimate of the fixed point for each point of the time
  // series:
  average_meter FP;
  for(int i=n_grid_points;i-->0;){
    FP.sample(f(xx(i))+pre_factor*exp(-decay_rate*xx(i)));
  }
#ifdef RELAX_DEBUG
  REPORT(xx);
  REPORT(y_av);
  REPORT(FP);
#endif
  fixed_point=FP.readout();

#ifdef RELAX_DEBUG
    WARNING("exponential_relaxation");
#endif
  return exponential_relaxation;
}


// for debugging:
#ifdef RELAX_DEBUG
#include "NewWeb.h"
static sequence< int > candidates;
static NewWeb * web_for_debugging;
#endif

// A class that can be wrapped around an ODE_state of log abundances
// to watch out for extinctions:
class extinguisher_t {
  ODE_state & _state;
  relaxing_dynamical_object & _dynamical_object;
  typedef enum {max,min,n_return_cases} _return_case;
  static const _return_case _reference_return_case;
  std::simple_vector< sequence<double> > _t;
  std::simple_vector< sequence<ODE_vector> > _return_state;
  double _start_time;

  static const int _states_kept=3;
  long int _states_read;
  double _sum_history[_states_kept];
  double _t_history[_states_kept];
  ODE_vector _state_history[_states_kept];
  int current(){return (_states_read)%_states_kept;}
  int last(){return (_states_read+(_states_kept-1))%_states_kept;}
  int second_last(){return (_states_read+(_states_kept-2))%_states_kept;}
  double _max_sum;
  double _min_sum;
  bool _last_member_dummy;
  double time_since_start() const {
    return _state.time_since_start();
  }
    
  // This member is were the decision on continuing or not made:
  bool _internal_recommendation(const sequence<double>& t,
				const sequence<ODE_vector>& _return_state,
				sequence<double> & min_point) const {
    const int n=t.size();
    const double t_max=t[n-1];
    // first, split off the third quarter and fourth quarter of the
    // time period covered by _return_state(s):
    int upper_half_start=0,upper_quarter_start=0;
    for(int i=n;i-->0;){
      if(t[i]<t_max*(3.0/4.0)){
	if(!upper_quarter_start){
	  upper_quarter_start=i+1;
	}
	if(t[i]<t_max/2){
	  upper_half_start=i+1;
	  break;
	}
      }
    }
    
    // Test if samples are evenly distributed in second half, which
    // would indicate that we reached a steady state:
    const int samples_third_quarter=upper_quarter_start-upper_half_start;
    const int samples_fourth_quarter=n-upper_quarter_start;
    const int N=
      samples_third_quarter+
      samples_fourth_quarter;
    const int Delta=
      samples_third_quarter-
      samples_fourth_quarter;

    // According to Przyborowski and Wilenski (C-test: Homogeneity of
    // results in testing samples from Poisson series, Biometrica 31,
    // 313-323, 1940) samples_third_quarter is binomially distributed
    // within the N samples;
#ifdef RELAX_DEBUG
    REPORT(N);
    REPORT(Delta);
#endif

    if(N<=2 ||
       N < min_samples || Delta*Delta*evenness_pickiness > N){
      // Sampling did not reach a steady state.  Let's wait.
      return false;
    }
    
    // Sampling rate has reached steady state.  We are ready to make
    // recommendations:
#ifdef RELAX_DEBUG
    REPORT(t_max);
#endif

    // Prepare computing confidence intervals:
    double student_t=gsl_cdf_tdist_Pinv(extinction_test_p,N-2);
#ifdef RELAX_DEBUG
    REPORT(student_t);
#endif
    double student_t2=student_t*student_t;
    // A negative correlation that we would recognize with confidence:
    double target_r=-sqrt(student_t2/(N-2+student_t2));
#ifdef RELAX_DEBUG
    REPORT(target_r);
#endif
      
    // Prepare computing confidence intervals for waiting.  This
    // includes a Bonferroni correction, otherwise we would wait
    // forever with large numbers of variables.
    double corrected_student_t=gsl_cdf_tdist_Pinv(extinction_test_p/n,N-2);
#ifdef RELAX_DEBUG
    REPORT(corrected_student_t);
#endif
    double corrected_student_t2=corrected_student_t*corrected_student_t;
    // A negative correlation that we would recognize with confidence:
    double corrected_target_r=-sqrt(corrected_student_t2/(N-2+corrected_student_t2));
#ifdef RELAX_DEBUG
    REPORT(corrected_target_r);
#endif

    double z_error=gsl_cdf_ugaussian_Pinv(pow(0.5,1.0/n))/sqrt(N-3);
    double corrected_waiting_r=
      tanh(atanh(waiting_r)+z_error);
#ifdef RELAX_DEBUG
    REPORT(z_error);
    REPORT(tanh(waiting_r));
    REPORT(corrected_waiting_r);
#endif
      
    std::simple_vector<double> y(N);
#ifdef RELAX_DEBUG
    candidates.resize(0);
#endif
    // Now, look for each variable separately:
    for(int k=_return_state[0].size();k-->0;){
      if(_dynamical_object.is_active(k)){
	// Sample out second half of time series:
	for(int i=upper_half_start;i<n;i++){
	  y[i-upper_half_start]=_return_state[i][k];
	}
	double amount_drifted;
	double r=
	  correlation_coefficient(&t(upper_half_start),&y(0),N,amount_drifted);
	double z=
	  // relative accuracy of slope:
	  drift_accuracy(&t(upper_half_start),&y(0),N);
	if(r<target_r||z<student_t){
#ifdef RELAX_DEBUG
	  candidates[candidates.size()]=k;
	  REPORT(k);
	  REPORT(web_for_debugging->assigned_column(k));
	  REPORT(z);
	  REPORT(r);
#endif
	  // 	  for(int i=upper_half_start;i<n;i++)
	  // 	    std::cout << t[i]+ << " " << _return_state[i][k]/log(10) << std::endl;
	}
// 	if(amount_drifted >= efficacious_amount_drifted_up){
// 	  // This series is clearly going up.  Let's wait until it
// 	  // reaches steady state.
// #ifdef RELAX_DEBUG
// 	  WARNING(k << "'s going up");
// 	  REPORT(web_for_debugging->assigned_column(k));
// 	  REPORT(amount_drifted);
// #endif
// 	  return false;
// 	}
	if(amount_drifted < -fixed_point_wobble_tolerance and
	   r<corrected_target_r){
	  // This time series appears to be decaying.  Have a close
	  // look to recognize the shape of its curve:
	  double fixed_point;
	  curve_shape_t s;
// 	  double rr=r;
// 	  REPORT(rr);
	  try{
	    s=shape(&t(upper_half_start),&y(0),N,fixed_point);
	  }catch(polyfit_error e){
	    WARNING("polyfit_error " << e.message);
	    s=complex_shape;
	  }
#ifdef RELAX_DEBUG
	  REPORT(k);
	  REPORT(web_for_debugging->assigned_column(k));
	  REPORT(z);
	  REPORT(r);
	  REPORT(s);
	  REPORT(amount_drifted);
#endif
	  switch(s){
	  case linear_drift:
	    // This is a candidate for deletion.
	    min_point[k]=-INFINITY;
	    break;
	  case exponential_relaxation:
	    // This may be another candidate for deletion (depends on
	    // fixed point).
#ifdef RELAX_DEBUG
	    REPORT(log(web_for_debugging->s(k).bodymass()));
#endif
	    min_point[k]=fixed_point;
	    break;
	  case complex_shape:
	    // Something strange is happening.  Let's wait.
	    return false;
	  }
	}else{
	  // We could not decide to delete this one.  Does it require
	  // further observation?
	  if(fabs(amount_drifted) > fixed_point_wobble_tolerance and
	     fabs(r) > corrected_waiting_r){
	    // Result is not clear, let's wait a bit more.
#ifdef RELAX_DEBUG
	    WARNING(k << "'s not clear");
	    REPORT(web_for_debugging->assigned_column(k));
	    REPORT(z);
	    REPORT(r);
	    REPORT(amount_drifted);
#endif
	    return false;
	  }
	  // This variable will probably stay where it is right now.
#ifdef RELAX_DEBUG
	  REPORT(k);
	  REPORT(web_for_debugging->assigned_column(k));
	  REPORT(z);
	  REPORT(r);
#endif
	  min_point[k]=y[N-1];
	}
      }
    }
    // There was not reason not to make a recommendation, so we make
    // one.
    return true;
  }

  bool _at_fixed_point(double t_max){
    if(_return_state[_reference_return_case].size()<3){
      return false;
    }
    double t_min=t_max*2/3;
    for(int k=_state.size();k-->0;){
      if(_dynamical_object.is_active(k)){
	double mi=_return_state[_reference_return_case][_return_state[_reference_return_case].size()-1][k];
	double ma=mi;
	for(_return_case c=_return_case(0);
	    c<n_return_cases;
	    c=_return_case(c+1)){
	  const int n_max=_return_state[c].size();
	  int n=n_max-1;
	  for(; n>=0 and _t[c][n]>t_min; n--){
	    double val=_return_state[c][n][k];
	    mi=std::min(mi,val);
	    ma=std::max(ma,val);
	  }
	  if(n_max-n<=2){
	    return false;
	  }
	}
	if(fabs(mi-ma)> fixed_point_wobble_tolerance ){
	  return false;
	}
      }
    }
    WARNING("FIXED POINT");
    return true;
  }


public:
  extinguisher_t (ODE_state & state, 
		  relaxing_dynamical_object * obj)
 :_state(state),
  _dynamical_object(*obj),
  _t(n_return_cases),
  _return_state(n_return_cases,sequence<ODE_vector>(0,ODE_vector(state))),
  _start_time(state.start_time()),
  _max_sum(0),_min_sum(1),
  _states_read(0),
  _last_member_dummy(false){
    
  }

  
  bool recommendation(sequence<double> & min_point){
    // preprocess and save current state:
    _states_read++;
    const int current_state=current();
    const int last_state=last();
    // Compute the sum of all active variables:
    double sum=_dynamical_object.active_sum(_state);

    // Record current sum:
    _sum_history[current_state]=sum;
    _t_history[current_state]=time_since_start();
    _state_history[current_state]=_state;
    if(_states_read<3){
      // initialize max and min:
      if(_max_sum<_min_sum){
	_max_sum=_min_sum=sum;
      }
      return false;
    }
    
    // Look first for local maxima, then for local minima (which are
    // not used, after all)...
    for(_return_case c=_return_case(0);c<n_return_cases;c=_return_case(c+1)){
      bool catch_it=false;
      switch(c){
      case max:
	if((_sum_history[last_state] > _sum_history[current()]) &&
	   (_sum_history[last_state] > _sum_history[second_last()]) ){
	  // Look for a new "quite large" maximum, but be forgiving
	  // over time.
	  if(_sum_history[last_state] >
	     (_max_sum-=extrema_relaxer)-extrema_tolerance){
	    _max_sum=_sum_history[last_state];
	    catch_it=true;
	  }
	}
	break;
      case min:
	if(_sum_history[last_state] < _sum_history[current()] &&
	   _sum_history[last_state] < _sum_history[second_last()] ){
	  if(_sum_history[last_state] < 
	     (_min_sum+=extrema_relaxer)+extrema_tolerance){
	    _min_sum=_sum_history[last_state];
	    catch_it=true;
	  }
	}
      }
      if(catch_it){
#ifdef RELAX_DEBUG
	REPORT(c);
	REPORT(_min_sum);
	REPORT(_max_sum);
#endif
	// We found a new extremum, 
	int n=_t[c].size();
	_t[c][n]=_t_history[last_state];
	_return_state[c][n]=_state_history[last_state];
	// Decide whether to compute a recommendation.  The criterion
	// is to do it only when all but the first few digits of n are
	// zero.
	while(n>333){
	  if(n%10!=0)
	    return false;
	  n/=10;
	}
	if(c==_reference_return_case){
	  // Base on the time-series of extrema, compute whether we
	  // should stop here, and which species to delete then:
	  if(_internal_recommendation(_t[c],_return_state[c],min_point)){
	    return true;
	  }
	  //// This is supposed to stop relaxing if a fixed point is
	  //// reached.  However, I could not see any significant speed
	  //// up.
#if 1
	  if(_at_fixed_point(time_since_start())){
	    for(int i=_state.size();i-->0;){
	      min_point[i]=_state[i];
	    }
	    return true;
	  }
#endif
	}
      }
    }
    return false;
  }
};

const extinguisher_t::_return_case extinguisher_t::_reference_return_case=
  extinguisher_t::min;

  
relaxing_dynamical_object::~relaxing_dynamical_object(){};

bool relaxing_dynamical_object::is_active(int i){
  return true;
}

double relaxing_dynamical_object::active_sum(ODE_vector & s){
  double sum=0;
  for(int i=s.size();i-->0;){
    sum+=s[i];
  }
  return sum;
}

static bool integrator_trouble_last_step=false; //for integrator trouble

double relaxing_dynamical_object::relax(double relaxation_time,
					const species_set_t & newly_inserted,
					const bool do_auto_extinguish){
  species_set_t deleted;
  return
    relax(relaxation_time,newly_inserted,deleted,do_auto_extinguish);
}
  
// Iterate ODE while watching for recommendations by extinguisher_t to
// stop and extinguish rare species.  The rest is book keeping.
double relaxing_dynamical_object::relax(double relaxation_time,
					const species_set_t & newly_inserted,
					species_set_t & deleted,
					const bool do_auto_extinguish){
  
  REPORT(relaxation_time);

  bool recording_dynamics=record_relaxation_dynamics;
  std::ofstream record;
  if(recording_dynamics) record.open("dynamics.dat");

  double &t=current_time;
  double t_start=t;
  double true_t_stop=t+relaxation_time;
 restart_relaxation:
  double t_stop=true_t_stop;
  int nsteps=0;
  ODE_vector ddt(number_of_variables());

  bool auto_extinguish=do_auto_extinguish;
  sequence<double> min_point(0);

#ifdef RELAX_DEBUG
  web_for_debugging=(NewWeb *)this;
#endif

  while(t<t_stop){
    if(number_of_variables()==0){
      return t-t_start;
    }
    { //begin scope ode_state
      ODE_state ode_state(this); 
      extinguisher_t extinguisher(ode_state,this);

      do{
	if(exit_now){
	  throw terminal_condition("exit_now set during relax");
	}

	// step integrator:
	double last_t=t;
	int integrator_failure;
	integrator_failure=ode_state.integrate_one_step(t_stop);
	if(newton_failure){
	  ///!!!! actually, the previous ODE state should be restored here!!!!
	  newton_failure=0;
	  integrator_failure=1;
	}
	if(integrator_failure){
	  ode_state.diagnosis();
	  if(integrator_trouble_last_step){
	    WARNING("Integrator trouble again");
	    REPORT(ode_state);
	    throw terminal_condition("repeated integrator trouble");
	  }else{
	    WARNING("Integrator trouble");
	  }
	  integrator_trouble_last_step=true;
	  break;
	}else{
	  integrator_trouble_last_step=false;
	}
	nsteps++;
	
	// For debugging integrator:
	// 	static int round=0;
	// 	if(++round==1000){
	// 	  round=0;
	// 	  ode_state.diagnosis();
	// 	}
	
	// These two lines save time series of biomass_B and
	// common_factor for later computations of invasion fitness.
	// However, this is not the best way of doing this.  We'd
	// better save the ode state vector in each step, and later
	// compute biomass_B and common_factor from this if and when
	// we really need it, e.g. using calls to dynamics() as done
	// here.
	dynamics(ode_state,ddt);
	record_for_steady_state(ode_state,ddt);
	
	if(auto_extinguish){
	  bool made_recommendation=extinguisher.recommendation(min_point);
	  if(made_recommendation){
	    // move to a random phase before exiting:
	    t_stop=t_start+(t-t_start)*(1+unirand()/min_samples);
	    auto_extinguish=false; //stop extinguishing for now
#ifdef RELAX_DEBUG
	    REPORT(min_point);
#endif
	  }
	}

	// we keep this in to get some feedback of what's happening:
	// delete it if you want to:
	if(recording_dynamics and t != t_start){
	  int old_prec=std::cout.precision();
	  record.precision(19);
	  record << t-t_start;
	  record.precision(old_prec);
	  record << "  " 
		 << ode_state << std::endl;
	}
	// end of output


      }while(t<t_stop && !small_values_in(ode_state,newly_inserted));
      //ode_state.diagnosis();
    } //end scope of ode_state
    if(do_auto_extinguish && !auto_extinguish && !min_point.empty()){
#ifdef RELAX_DEBUG
//       std::cout << "candidates: ";
//       for(int i=candidates.size();i-->0;){
// 	std::cout << assigned_column(candidates(i)) << " ";
//       }
//       std::cout << std::endl;
#endif
      species_set_t just_deleted=
	delete_species_larger_than_exp(min_point,newly_inserted);
      if(!just_deleted.empty()){
	std::cout << "... at " << t-t_start << std::endl;
	deleted.insert(just_deleted.begin(),just_deleted.end());
	goto restart_relaxation;
      }
    }
    species_set_t just_deleted=
      delete_all_species_with_less_than_one_individual(newly_inserted);
    if(!just_deleted.empty()){
      std::cout << "... at " << t-t_start << std::endl;
      deleted.insert(just_deleted.begin(),just_deleted.end());
    }
    if(auto_extinguish && !newly_inserted.empty() && deleted==newly_inserted){
      // Exactly all inserted species went extinct.  We should be in a
      // steady state again.  Move to a random phase and stop:
      t_stop=
	std::min(t_stop,
		 std::max(t+dephasing_time,
			  t+(t-t_start)*unirand()) );
      auto_extinguish=false; // Here indicating that we don't need to look
		             // for further extinctions.
      min_point.resize(0);
    }
  }//while(t<t_stop)

  // write again to give some account of what we have lastly extinguished.
  
  if(number_of_variables()>0){//running this without species would lead to crash
    ODE_state ode_state(this);
    int old_prec=std::cout.precision();
    record.precision(19);
    record << t-t_start;
    record.precision(old_prec);
    record << "  " 
	   << ode_state << std::endl;
    // end of output
  }

  const double true_relaxation_time=t-t_start;
  bool is_fixedpoint=true;
  for(int i=ddt.size();i-->0;){
    if(fabs(ddt[i])*true_relaxation_time>fixed_point_wobble_tolerance){
      is_fixedpoint=false;
      break;
    }
  }
  if(is_fixedpoint){
    steady_state_is_fixed_point();
  }
  
  return true_relaxation_time;
}


void 
combined_relaxing_dynamical_object::
split_conserved(const species_set_t& conserved,
		species_set_t &set1,
		species_set_t &set2){
  for(species_set_t::const_iterator i=conserved.begin();
      i!=conserved.end();++i){
    if( (*i|1) == 0 ){
      set1.insert(*i/2);
    }else{
      set2.insert(*i/2);
    }
  }  
}

void 
combined_relaxing_dynamical_object::
merge_into_deleted(const species_set_t &set1,
		   const species_set_t &set2,
		   species_set_t & deleted){
  for(species_set_t::const_iterator i=set1.begin(); i!=set1.end(); ++i){
    deleted.insert(2*(*i));
  }
  for(species_set_t::const_iterator i=set2.begin(); i!=set2.end(); ++i){
    deleted.insert(2*(*i)+1);
  }
}
