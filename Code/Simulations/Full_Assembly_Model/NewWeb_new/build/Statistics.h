// -*- c++ -*-
// $Id: Statistics.h 2466 2016-05-01 23:27:44Z axel $

/// \file Classes useful for statistical analyses.

#ifndef __STATISTICS__
#define __STATISTICS__

// #include <CLHEP/Matrix/Matrix.h>
// #include <CLHEP/Matrix/Vector.h>

#include "NewMatrix.h"
#include "simple_vector.h"
#include "error.h"

/// Abstract base class to estimate distributions from samples.
/** This is technical overkill!?! Perhaps.  But very often you start
    doing a historgram, and later decide you might try something
    better.*/
class Distribution_Estimator{
 protected:
  int N; // number of samples
 public:
  Distribution_Estimator():N(0){};
  virtual void sample(double x)=0;
  // !! this gives a density in units of 1/[unit of x] !!
  virtual double estimate_density(double x)=0;
  virtual ~Distribution_Estimator(){};
};

/// Compute a histogram from samples.
class Histogram_Estimator : public Distribution_Estimator {
  std::simple_vector<int> the_histogram;
  double the_lowest_x_included;
  double the_bin_width;
 public:
  Histogram_Estimator(double min, double max, int N);
  Histogram_Estimator(double min, double width);
  int the_number_of_bins(){
    return the_histogram.size();
  }
  void sample(double x);
  void save(const char * filename);
  double estimate_density(double x);
  virtual ~Histogram_Estimator();
};

/* class average_meter{ */
/*   int samples; */
/*   double sum; */
/*   double square_sum; */
/* public: */
/*   average_meter() { */
/*     samples=0; */
/*     sum=0; */
/*     square_sum=0; */
/*   }; */
/*   void sample(double x){ */
/*     sum+=x; */
/*     square_sum+=x*x; */
/*     samples++; */
/*   } */
/*   double readout(){ */
/*     return sum/samples; */
/*   } */
/*   double var(){ */
/*     return (square_sum-sum*sum/samples)/(samples-1); */
/*   } */
/*   int n(){ */
/*     return samples; */
/*   } */
/* }; */

/// Collects samples (with weights) to compute their average.
class weighted_average_meter{
 protected:
  long double samples;
  long double sum;
  long double square_sum;
 public:
  weighted_average_meter():samples(0),sum(0),square_sum(0){};
  weighted_average_meter(const double x):samples(1),sum(x),square_sum(x*x){};
  weighted_average_meter(const double x,const double varx):
    samples(1),sum(x),square_sum(varx+x*x){};
  void sample(double x,double weight);
  double readout() const;
  double sample_var() const;
  double sample_std() const;
  weighted_average_meter  & 
  operator+=(const weighted_average_meter & other){
    samples+=other.samples;
    sum+=other.sum;
    square_sum+=other.square_sum;
    return *this;
  }
  weighted_average_meter  & 
  operator+=(const double y){
    square_sum+=2*y*sum+y*y*samples;
    sum+=y*samples;
    return *this;
  }
  weighted_average_meter  & 
  operator-=(const double y){
    return operator+=(-y);
  }
  operator double(){
    return readout();
  }
};

/// Collects samples to compute their average.
class average_meter : public weighted_average_meter {
 public:
  average_meter():weighted_average_meter(){};
  average_meter(const double x):weighted_average_meter(x){};
  average_meter(const weighted_average_meter m):weighted_average_meter(m){};
  void sample(double x);
  void sample(double x,double weight);// produces error message
  double var() const;
  double std() const;
  double error() const;
  double error_var() const;
  int n();
  average_meter & operator+=(const double x){
    sample(x);
    return *this;
  }
  average_meter  & 
  operator+=(const average_meter & other){
    this->weighted_average_meter::operator+=(other);
    return *this;
  }
};

std::ostream & operator<<(std::ostream &stream, 
			  const average_meter & av);

inline average_meter operator+(average_meter & a, double x){
  average_meter b;
  b+=a;
  b+=x;
  return b;
}

/// Samples vectors to estimate their multivariate normal distribution.
class multinormal_distribution {
public:
  NewVector mean;
  NewMatrix cov;
  multinormal_distribution():mean(),cov(){};
  multinormal_distribution(NewVector &m, NewMatrix &c);
  std::simple_vector<double> main_axis()const;
  std::simple_vector<double> raw_main_axis()const;
  void save(const char * name)const;
  void load(const char * name);
  double var_of_sum();
};
  

/// Used to test if empirical data is consistent with Monte-Carlo simulations.
class chi_square_meter {
  const std::simple_vector<bool> _selection;
  int _size;
  int _n_samples;
  NewVector _sum;
  NewMatrix _square_sum;
  void _initialize_selection_given();
  void postfix(std::simple_vector<double> const & data,std::simple_vector<bool> const & fix,
	       NewVector & data_star, NewVector & m_star, 
	       NewMatrix & cov_star, NewMatrix & icov_star) const;
public:
  explicit chi_square_meter(std::simple_vector<bool> sel) :
    _selection(sel), _n_samples(0){
    _initialize_selection_given();
  }
  explicit chi_square_meter(int size) :
    _selection(size,true), _n_samples(0){
    _initialize_selection_given();
  }
  void sample(std::simple_vector<double> data);
  std::simple_vector<double> mean() const;
  double chi_square_old(std::simple_vector<double> data) const;
  double chi_square(std::simple_vector<double> const & data,
		    std::simple_vector<bool> const & fix=std::simple_vector<bool>(0)) const;
  double chi_square(double * log_det_cov,
		    std::simple_vector<double> const & data,
		    std::simple_vector<bool> const & fix=std::simple_vector<bool>(0)) const;
  std::simple_vector<weighted_average_meter>
  mean_and_var(std::simple_vector<double> const & data,
	       std::simple_vector<bool> const & fix=std::simple_vector<bool>(0)) const;
  std::simple_vector<double> 
  deviation_old(std::simple_vector<double> const & data) const;
  std::simple_vector<double> 
  deviation(std::simple_vector<double> const & data,
			      std::simple_vector<bool> const & fix
			      =std::simple_vector<bool>(0)) const;
  multinormal_distribution
  estimate(std::simple_vector<double> const & data=std::simple_vector<double>(0),
	   std::simple_vector<bool> const & fix=std::simple_vector<bool>(0));
};

// save/load empirical data:
void save(std::simple_vector<double> const & data,
	  std::simple_vector<bool> const & select,
	  std::simple_vector<bool> const & fix,
	  const char * name);

void load(NewVector & data,
	  const char * name);

#endif // __STATISTICS__
