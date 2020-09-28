// -*- mode: c++ -*-
// $Id: period_cutter.h 1431 2009-05-04 12:22:40Z axel $
#ifndef _PERIOD_CUTTER_H_
#define _PERIOD_CUTTER_H_

#include "sequence.h"
#include <math.h>

//#include <iostream>
//#include "error.h"

/// Identifies periodicity and other kinds of dynamics.
/** Among others, this perdicts when one oscillation cycle is over.
    You have to provide samples of the distance d to the initial state
    in phase space every now and then.

    The time of reaching the minimum distance is to be predicted in
    last phase of an oscillation.
    
    d^2 ~ v0^2*(t-t0)^2 + d0
    
    This can be done by predicting the zero crossing of d d^2/dt, which
    requires only two points of d, since we can get v0 from the start phase.
    
    d d^2/dt = 2*v0^2*(t-t0)    
*/
class period_cutter_t 
{
 private:
  typedef enum 
    {no_point_sampled_yet,
     test_for_shooting_off,
     no_maximum_reached_yet,
     some_maxima_reached,
     increasing_again,
     near_zero,
     passing_by,
     not_shot_off,
     number_of_phases} phase_t;
  double the_shoot_off_speed_v0;
  double the_largest_maximum;
  static const int the_number_of_samples_kept;
  sequence<double> the_last_d;
  sequence<double> the_last_t;
  phase_t the_phase;
  double the_start_time;
  int the_number_of_maxima;
  int the_number_of_samples;
  int the_number_of_replacements_of_second_largest_maximum;
  double the_second_largest_maximum;
  sequence<double> the_largest_maxima;
  sequence<double> the_times_of_largest_maxima;
  bool looks_like_steady_state_internal();
  bool looked_like_steady_state;
 public:
  period_cutter_t(double t); ///< \a t is start time.
  ~period_cutter_t();
  void reset(double t); ///< Dito; \a t is new start time.
  void sample(double dist, double t); ///< \a dist: distance to initial state.
  bool can_predict_end_time() const; ///< Is end of cycle in sight?
  double predicted_end_time() const; ///< Predicted end of cycle.
  double period_length() const; ///< Predicted oscillation period.
  double looks_like_chaos() const; ///< >=1 if state looks like chaos
  double time_inspecting() const; ///< Time since initial sample.
private:
  double dddt(int i=0);
public:
  bool convergence_to_steady_state(); ///< Approaching a fixed point?
  double way_to_convergence(); ///< How far is the fixed point?
  bool noisy_steady_state(); ///< Is this a fixed point + noise?
  int number_of_samples(); ///< Total number samples taken so far.
};

#endif // _PERIOD_CUTTER_H_
