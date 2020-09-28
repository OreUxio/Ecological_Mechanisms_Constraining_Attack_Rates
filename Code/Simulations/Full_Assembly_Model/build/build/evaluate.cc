// $Id: evaluate.cc 2305 2013-06-01 07:57:00Z axel $
#include <iostream>
#include <stdio.h>
#include "evaluate.h"
#include "error.h"
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>

static double eval_expintEn(double a, double b) { 
  return gsl_sf_expint_En(a,b); 
} 
static double eval_gamma(double a) { 
  return tgamma(a); 
} 
static double eval_gamma_inc(double a,double x) { 
  return gsl_sf_gamma_inc(a,x); 
} 
static double eval_floor(double a) { 
  return floor(a); 
} 


my_evaluator_t::my_evaluator_t() : HepTool::Evaluator() {
  HepTool::Evaluator::setStdMath();
  //HepTool::Evaluator::setSystemOfUnits(); // SI Units
  //  HepTool::Evaluator::
  //  setSystemOfUnits(100.,1000,1.0,1.0,1.0,1.0,1.0); // cgs Units
  HepTool::Evaluator::
    setSystemOfUnits(100.,1000,1.0/(86400*360),1.0,1.0,1.0,1.0); //cga Units ([T]=Year)
  setVariable("kilo_", 1.e+03); // chilioi (Greek) "thousand"
  setVariable("mega_", 1.e+06); // megas (Greek) "large"
  setVariable("giga_", 1.e+09); // gigas (Greek) "giant"
  setVariable("tera_", 1.e+12); // teras (Greek) "monster"
  setVariable("peta_", 1.e+15); // pente (Greek) "five"
  setVariable("exa_", 1.e+18);
  setVariable("zetta_", 1.e+21);
  setVariable("yotta_", 1.e+24);
  setVariable("deci_", 1.e-01); // decimus (Latin) "tenth"
  setVariable("centi_", 1.e-02); // centum  (Latin) "hundred"
  setVariable("milli_", 1.e-03); // mille   (Latin) "thousand"
  setVariable("micro_", 1.e-06); // micro (Latin) or mikros (Greek) "small"
  setVariable("nano_", 1.e-09); // nanus (Latin) or nanos  (Greek) "dwarf"
  setVariable("pico_", 1.e-12); // pico (Spanish) "bit"
  setVariable("femto_", 1.e-15);
  setVariable("atto_", 1.e-18); 
  setVariable("zepto_", 1.e-21); 
  setVariable("yocto_", 1.e-24); 

  setVariable("sec",evaluate("1 * s"));
  setVariable("second",evaluate("1 * s"));
  setVariable("seconds",evaluate("1 * s"));
  setVariable("min",evaluate("60 * s"));
  setVariable("minute",evaluate("60 * s"));
  setVariable("minutes",evaluate("60 * s"));
  setVariable("h",evaluate("60 * min"));
  setVariable("hour",evaluate("60 * min"));
  setVariable("hours",evaluate("60 * min"));
  setVariable("day",evaluate("24 * h"));
  setVariable("days",evaluate("24 * h"));
  setVariable("month",evaluate("30 * days"));
  setVariable("months",evaluate("30 * days"));
  setVariable("week",evaluate("7 * days"));
  setVariable("weeks",evaluate("7 * days"));
  setVariable("year",evaluate("12 * month"));
  setVariable("years",evaluate("12 * month"));
  setVariable("ha",evaluate("100*100*meter*meter"));
  setVariable("microgram",evaluate("1e-6 * gram"));
  setVariable("ton",evaluate("1e6 * gram"));
  setVariable("tons",evaluate("1e6 * gram"));

  setVariable("true", double(1)); 
  setVariable("false", double(0)); 

  setFunction("expintEn", eval_expintEn);
  setFunction("gamma", eval_gamma);
  setFunction("gamma_inc", eval_gamma_inc);
  setFunction("floor", eval_floor);

  struct timeval t;
  if(gettimeofday(&t,0)) FATAL_ERROR("problem with system clock");
  srand(t.tv_usec); // a new seed every microsecond
}

bool my_evaluator_t::
if_error_hint_at_position(const char * const &  value){
  if(status() !=  HepTool::Evaluator::OK){
    if(status() != HepTool::Evaluator::WARNING_BLANK_STRING){
      // indicate error position:
      std::cout << value << std::endl;
      for(int i=0;i<error_position();i++) std::cout << "-";
      std::cout << "^" << std::endl;
      std::cout.flush();
      print_error();
    }
    return 1; //error
  }
  return 0; // OK
}

double my_evaluator_t::operator()(const char * const &  value){
  setVariable("random",rand());
  setVariable("unirand",rand()/double(RAND_MAX));
  double v=evaluate(value);
  if(if_error_hint_at_position(value))
    FATAL_ERROR("SYNTAX_ERROR");
  return v;
}
    

my_evaluator_t eval;    

const char * format_string1 = "%.0f";

// Use some pre-processor trickery to get correct accuracy into format string:
#ifndef DBL_MANT_DIG
#error DBL_MANT_DIG not defined in float.h!
#endif
#define FORMAT_STRING2a(DIGa) "%." #DIGa "g";
#define FORMAT_STRING2(DIG) FORMAT_STRING2a(DIG)
const char * format_string2 = FORMAT_STRING2(DBL_MANT_DIG);
#undef FORMAT_STRING2
#undef FORMAT_STRING2a

// only this one can be a friend of eval:
int evaluate_expression(char ** value){
  eval.setVariable("random",rand());
  double v=eval.evaluate(*value);

  // error handling:
  if(eval.if_error_hint_at_position(*value))
    return 1;//error
      
  const char * format_string =
    //try to force printing integer values as integers without
    //scientific 'E' notation:
    (fabs(v)+2*DBL_EPSILON >= 1 and fabs(int(fabs(v)+2*DBL_EPSILON)-fabs(v))<2*DBL_EPSILON ? 
     format_string1 : format_string2);

  int str_length;
  //snprintf behaves a bit differently (see man page)
#if defined(ON_SX5FSV) || defined(SX)
  str_length=128;
  char * new_value = (char *) malloc((str_length+1)*sizeof(char));
  while(snprintf(new_value,str_length,format_string,v) >= str_length-1){
    free(new_value);
    str_length*=2;
    new_value = (char *) malloc((str_length+1)*sizeof(char));
  }
#else
  str_length=snprintf(0,0,format_string,v); // count length of new value
  if(str_length<0)
    return 1; //error
  char * new_value = (char *) malloc((str_length+1)*sizeof(char));

  if(!new_value)
    return 1; //error
  sprintf(new_value,format_string,v); // write new value
#endif

  free(*value); 
  *value=new_value;
  return 0; // all right
}
    
void my_evaluator_t::set_variable(const char * name,double value){
  setVariable(name,value);
}

void set_evaluator_variable(const char * name,double value){
  eval.setVariable(name,value);
  //  #ifdef DEBUGGING
  std::cout << name << " = " << value << std::endl;
  //  #endif
}
