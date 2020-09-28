// -*- c++ -*-
//$Id: polyfit.h 2386 2014-10-08 10:13:51Z axel $

#ifndef __POLYFIT_H__
#define __POLYFIT_H__

#include <string>
#include "gsl/gsl_multifit.h"
#include "Statistics.h"
#include "sequence.h"

class polyfit_error{
public:
  const char * message;
  polyfit_error(const char * m): 
    message(m){};
};

/// A polynomial fit to the data provided at construction.
class fitted_function {
  int _n;
  gsl_vector *_a;
  gsl_matrix *_COV;
  explicit fitted_function(int n);
 public:
  ~fitted_function();
  fitted_function();
  fitted_function(const fitted_function & other);
  fitted_function const & operator= (fitted_function const& other);
  // polynomial fit, order automatically determined:
  fitted_function(sequence<double>& x,sequence<average_meter>& y,
		  bool re_estimate_co_variances=true);
  // multivariable fit:
  fitted_function(sequence< sequence<double> >& x,sequence<average_meter>& y,
		  bool re_estimate_co_variances=true);
  // (n-1)-th order polynomial fit:
  fitted_function(sequence<double>& x,sequence<average_meter>& y,int n,
		  bool re_estimate_co_variances=true); 
  // get value as predicted:
  double operator()(double x) const;
  double operator()(sequence<double> x) const;
  double cov_at(double x1, double x2) const;
  double cov_at(const sequence<double> &x1, const sequence<double> &x2) const;
  double var_at(double x) const;
  double var_at(const sequence<double> &x) const;
  double operator[](int n) const{
    ASSERT(n>=0);
    if(n<_n) return gsl_vector_get(_a,n);
    else return 0;
  }
  double var_of(int n) const{
    ASSERT(n>=0);
    if(n<_n) return gsl_matrix_get(_COV,n,n);
    else return 0;
  }
  std::string operator()(const std::string s) const;
  fitted_function derivative() const;
  int polynomial_order() const {return _n-1;}
};
    

#endif // __POLYFIT_H__
