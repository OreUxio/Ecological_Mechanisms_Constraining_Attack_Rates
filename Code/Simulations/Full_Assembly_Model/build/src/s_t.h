// -*- c++ -*-
//$Id: s_t.h 2014 2010-12-07 10:36:30Z axel $

#ifndef __S_T_H__
#define __S_T_H__

#include "props_t.h"

extern double lambda;
extern double plant_fraction; // -1 mean no distiction between animals and plants

class s_t {
private:
  s_t & operator=(const s_t & other); //this yields an error (why should it?)
  static const int unset;
  static int link_threshold;
  static double true_C0;
  static void set_link_threshold(double C0);
  // helper type to pass parameters below:
  struct param_t {double p_break;};
  static double 
    delta_p_break(double switching_probability,void * vparam);
  static double 
    compute_swapping_probability(double beta);

public:
  double s;
  bool is_plant;
  props_t weapons;
  props_t armor;
  static props_t::swapping_probability swap_weapons;
  static props_t::swapping_probability swap_armor;
  static void set_C0(double C0){
    set_link_threshold(C0);
  }
  static void set_beta_r(double beta_r){
    swap_weapons=
      compute_swapping_probability(beta_r);
  }
  static void set_beta_c(double beta_c){
    swap_armor=
      compute_swapping_probability(beta_c);
  }
  static void set_swap_armor(double p){
    swap_armor=p;
  }
  static void set_swap_weapons(double p){
    swap_weapons=p;
  }
public:
  s_t mutate(){
    weapons.swap_a_fraction(swap_weapons);
    armor.swap_a_fraction(swap_armor);
    return (*this);
  }
  bool eats(s_t & other){
    if(other.s+lambda > (*this).s && !is_plant){
      return 
	(this->weapons.matches_with(other.armor) >= link_threshold);
    }else{
      return false;
    }
  }
  void set_random_props();
  s_t(double ss):s(ss){
    set_random_props();
  }
  bool operator<(s_t & other){
    return s > other.s;
  }
  static double get_true_C0(){return true_C0;};
};

#endif // __S_T_H__