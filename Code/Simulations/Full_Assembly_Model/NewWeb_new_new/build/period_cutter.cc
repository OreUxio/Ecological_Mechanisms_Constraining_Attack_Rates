// -*- mode: c++ -*-
// $Id: period_cutter.cc 1968 2010-11-12 17:30:43Z axel $

#include "period_cutter.h"
#include <algorithm>

static double fluctuation_tolerance=0.01*1.6;
static double convergence_detection_epsilon=0.02;
static double noisy_steady_state_detection_p=0.95;
static int replacements_of_max2_for_chaos=2;
static int maxima_kept=40;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(fluctuation_tolerance),
  CFGDOUBLE(convergence_detection_epsilon),
  CFGDOUBLE(noisy_steady_state_detection_p),
  CFGINT(replacements_of_max2_for_chaos),
  CFGINT(maxima_kept),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);


const int period_cutter_t::the_number_of_samples_kept=10;

period_cutter_t::period_cutter_t(double t){reset(t);};

period_cutter_t::~period_cutter_t(){};

void period_cutter_t::reset(double t){
  the_shoot_off_speed_v0=0;//required?

  the_largest_maximum=0;//required?
  the_last_d[the_number_of_samples_kept-1]=0;
  the_last_t[the_number_of_samples_kept-1]=0;
  the_phase=no_point_sampled_yet;
  the_start_time=t;
  the_number_of_maxima=0;
  the_number_of_samples=0;
  the_last_d*=0;
  the_last_t*=0;
  the_largest_maximum=0;
  the_number_of_replacements_of_second_largest_maximum=0;
  the_second_largest_maximum=0;
  the_largest_maxima.resize(0);
  the_times_of_largest_maxima.resize(0);
  looked_like_steady_state=false;
}

void period_cutter_t::sample(double dist, double t){
  //  REPORT(t-the_start_time);
  if(t<=the_last_t[0])
    return;
  switch(the_phase){
  case no_point_sampled_yet:
    if(t > the_start_time){
      the_shoot_off_speed_v0=dist/(t-the_start_time);
      the_phase=test_for_shooting_off;
    }
    break;
  case test_for_shooting_off:
    if(dist < the_last_d[0]){
      the_phase=not_shot_off;
    }
    the_phase=no_maximum_reached_yet;
    break;
  case no_maximum_reached_yet:
    if(dist < the_last_d[0]){
      the_largest_maximum=the_last_d[0];
      the_number_of_maxima=1;
      the_phase=some_maxima_reached;
    }
    break;
  case some_maxima_reached:
    if(dist > the_last_d[0]){
      the_phase=increasing_again;
      break;
    }
    { //(scope of transition_radius)
      double transition_radius;
      if(the_number_of_samples<4)
	break;
      transition_radius=the_largest_maximum/4;
      //REPORT(transition_radius);
      if(transition_radius>
	 fluctuation_tolerance*the_largest_maximum/the_number_of_maxima){
	transition_radius=
	  fluctuation_tolerance*the_largest_maximum/the_number_of_maxima;
	//REPORT(transition_radius);
      }
      if(dist < transition_radius){
	the_phase=near_zero;
      }
    }
    break;
  case increasing_again:
    if(dist < the_last_d[0]){
      // for fractality test
      if(the_last_d[0] > the_largest_maximum){
	// replace largest maximum
	the_second_largest_maximum=the_largest_maximum;
	the_largest_maximum=the_last_d[0];
	the_number_of_replacements_of_second_largest_maximum=0;
      }else if(the_number_of_maxima >=2 &&
	       the_last_d[0] > the_second_largest_maximum){
	// replace second largest maximum
	the_second_largest_maximum=the_last_d[0];
	the_number_of_replacements_of_second_largest_maximum++;
      }

      // for general steady-state test
      if(the_largest_maxima.size() >= maxima_kept && 
	 the_largest_maxima[maxima_kept-1]>0){
	looked_like_steady_state=
	  looks_like_steady_state_internal();
      }

      if(the_largest_maxima.size()){
	//REPORT(the_largest_maxima);
	//REPORT(the_last_d[0]);
	int i=the_largest_maxima.size()-1;
	while(i>=0 && the_largest_maxima[i]<the_last_d[0]) i--;
	i++;
	int j=(the_largest_maxima.size()<maxima_kept ?
	       the_largest_maxima.size() :
	       maxima_kept-1 );
	for(;j>i;j--){
	  the_largest_maxima[j]=
	    the_largest_maxima[j-1];
	  the_times_of_largest_maxima[j]=
	    the_times_of_largest_maxima[j-1];
	}
	the_largest_maxima[i]=the_last_d[0];
	the_times_of_largest_maxima[i]=t;
      }else{
	the_largest_maxima[0]=the_last_d[0];
	the_times_of_largest_maxima[0]=t;
      }
      
      //for periodicity test
      the_number_of_maxima++;
      the_phase=some_maxima_reached;
    }
    break;
  case near_zero:
    break;
  }
  the_number_of_samples++;
  the_last_d.prepend(dist);
  the_last_t.prepend(t);
  the_last_d.resize(the_number_of_samples_kept);
  the_last_t.resize(the_number_of_samples_kept);
  //     std::cout << "phase " << the_phase 
  // 	      << " d=" << dist 
  // 	      << " after " 
  // 	      << t-the_start_time << std::endl;
}

bool period_cutter_t::looks_like_steady_state_internal(){
  // general steady-state test: the times of the maxima_kept largest
  // maxima are evenly distributed and uncorrelated.
  
  // test for even distribution:
  typedef enum {point20,point15,point10,point05,point01,
		n_levels} significance_t;
  significance_t alpha=point20;
  const double KS_quantile[][n_levels]=
    {{.900,.925,.950,.975,.995},
     {.684,.726,.776,.842,.929},
     {.565,.597,.642,.708,.828},
     {.494,.525,.564,.624,.733},
     {.446,.474,.510,.565,.669},
     {.410,.436,.470,.521,.618},
     {.381,.405,.438,.486,.577},
     {.358,.381,.411,.457,.543},
     {.339,.360,.388,.432,.514},
     {.322,.342,.368,.410,.490},
     {.307,.326,.352,.391,.468},
     {.295,.313,.338,.375,.450},
     {.284,.302,.325,.361,.433},
     {.274,.292,.314,.349,.418},
     {.266,.283,.304,.338,.404},
     {.258,.274,.295,.328,.392},
     {.250,.266,.286,.318,.381},
     {.244,.259,.278,.309,.371},
     {.237,.252,.272,.301,.363},
     {.231,.246,.264,.294,.356} };
  const int KS_max_n=20;
  const double KS_asymptote[]={1.07,1.14,1.22,1.36,1.63};

  static sequence<double> sorted;

  // standardize:
  sorted=the_times_of_largest_maxima;
  int n=sorted.size();
  sorted-=the_start_time;
  sorted/=(the_last_t[0]-the_start_time);

  // compute Durbin-Watson statistic to test for correlations:
  double dwnum=0,dwden=
    (sorted[0]-0.5)*(sorted[0]-0.5);
  for(int i=1;i<n;i++){
    double diff=sorted[i]-sorted[i-1];
    dwnum+=diff*diff;
    dwden+=(sorted[i]-0.5)*(sorted[i]-0.5);
  }
  double dw=dwnum/dwden;
  //REPORT(n);
  //REPORT(sorted);
  //REPORT(dwnum);
  //REPORT(dwden);
  //REPORT(dw);
  if(dw<1.1) return false;
  
  // compute Kolmogorov-Smirnov statistic D
  double D=0;
  std::sort(sorted.begin(),sorted.end());
  for(int i=n;i-->0;){
    if(sorted[i]-i/double(n)>D){
      D=sorted[i]-i/double(n);
    }
    if(-(sorted[i]-(i+1)/double(n))>D){
      D=-(sorted[i]-(i+1)/double(n));
    }
  }
  
  double D_quantile=
    (n<=KS_max_n ? 
     KS_quantile[n-1][alpha] :
     KS_asymptote[alpha]/sqrt(n) );
  
  //REPORT(D);
  //REPORT(D_quantile);
  return D < D_quantile;
}
     


bool period_cutter_t::can_predict_end_time() const {
  return 
    the_phase==near_zero || 
    the_phase==not_shot_off;
}

double period_cutter_t::predicted_end_time() const {
  //ALWAYS_ASSERT(can_predict_end_time());
  if(the_phase==not_shot_off){
    return the_last_t[1];
  }else{
    // compute d d^2/dt:
    double dd2dt=
      (the_last_d[0]*the_last_d[0]-
       the_last_d[1]*the_last_d[1])/
      (the_last_t[0]-the_last_t[1]);
    return -dd2dt/(2*the_shoot_off_speed_v0*the_shoot_off_speed_v0)+
      (the_last_t[0]+the_last_t[1])/2;
  }
}

double period_cutter_t::period_length() const{
  return predicted_end_time()-the_start_time;
}

double period_cutter_t::looks_like_chaos() const{
  if(looked_like_steady_state)
    return 100;

  // fractality test:
  if(the_number_of_maxima>=10){
    return the_number_of_replacements_of_second_largest_maximum/
      double(replacements_of_max2_for_chaos);
  }else
    return 0;
}

double period_cutter_t::time_inspecting() const{
  return the_last_t(0)-the_start_time;
}

double period_cutter_t::dddt(int i){
  return (the_last_d[i+0]-the_last_d[i+1])/(the_last_t[i+0]-the_last_t[i+1]);
}

bool period_cutter_t::convergence_to_steady_state(){
  //This test works pretty independent of loop detection, but you
  //need to collect samples, too.

  //We just test for an exponential decay of differences, we might
  //actually be quite fare away from the steady state!;

  if(!(the_last_d[4] && the_last_t[4])){
    // insufficient samples
    return false;
  }
  return 
    fabs(dddt(0)/dddt(1)-dddt(1)/dddt(2))/(the_last_t[1]-the_last_t[3])
    <convergence_detection_epsilon &&
    fabs(dddt(1)/dddt(2)-dddt(2)/dddt(3))/(the_last_t[2]-the_last_t[4])
    <convergence_detection_epsilon &&
    dddt(0)/dddt(1) < 1;
}

double period_cutter_t::way_to_convergence(){
  return fabs((dddt(0)*dddt(0))/(the_last_d[0]*(dddt(0)-dddt(1))/(the_last_t[0]-the_last_t[1])));
}

bool period_cutter_t::noisy_steady_state(){
  //This test works pretty independent of loop detection, but you
  //need to collect samples, too.
  if(!(the_last_d[the_number_of_samples_kept-1] && 
       the_last_t[the_number_of_samples_kept-1])){
    // insufficient samples
    return false;
  }

  double mean_d=
    sum(the_last_d)/the_number_of_samples_kept;
  double var_d=
    sum(the_last_d*the_last_d)/the_number_of_samples_kept
    - mean_d*mean_d;

  //  REPORT(the_number_of_maxima);
  //  REPORT(the_number_of_samples);
  return 
    the_number_of_maxima > 4 &&
    the_number_of_maxima > 
    the_number_of_samples/10.0 && // the expectation value for
    // uncorrelated noise is
    // the_number_of_samples/4.0
    mean_d < 4*sqrt(var_d);
}

int period_cutter_t::number_of_samples(){
  return the_number_of_samples-1; // one step between two samples :)
}
