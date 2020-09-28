// -*- mode: c++ -*-
// $Id: vector_with_max.h 2150 2011-05-23 08:43:17Z axel $
#ifndef _VECTOR_WITH_MAX_H_
#define _VECTOR_WITH_MAX_H_

#include "vector"
#include "math.h"

/// Vector that keeps track of the max absolute value of its entries.
/** For use with SortedMatrix and SortedVector. */ 
class vector_with_max : private std::vector< double >
{
  double max;
  class active_reference_t;

public:
  vector_with_max():max(0){};
  vector_with_max(const vector_with_max &other):
    std::vector<double>(other),max(other.max){};
  vector_with_max(double const* b,double const* e):
    std::vector<double>(b,e){fix_max();};
  vector_with_max(float const* b,float const* e):
    std::vector<double>(b,e){fix_max();};
  explicit vector_with_max(size_t s):
    std::vector<double>(s,0),max(0){};
  const double * get_data() const {return &std::vector<double>::operator[](0);}
  double get_max() const {return max;}
  void reset_max(){max=0;};
  inline active_reference_t operator()(unsigned int i);
  inline double operator()(unsigned int i) const;
  size_type size() const {return std::vector<double>::size();};
  void resize(size_t n){std::vector<double>::resize(n);};
  const_iterator begin() const {return std::vector<double>::begin();}
  const_iterator end() const {return std::vector<double>::end();}
  vector_with_max& operator=(const vector_with_max& other);

  // thread-safe interface:
  reference operator[](size_type i){
    return std::vector<double>::operator[](i);
  }
  const_reference operator[](size_type i) const {
    return std::vector<double>::operator[](i);
  }
  void fix_max();
  void set_max(double value) {max=value;}
};

class vector_with_max::active_reference_t{
  double & v;
  vector_with_max & origin;

public:
  active_reference_t(double & x, vector_with_max & m):v(x),origin(m){};
  operator double () const{
    return v;
  }
  double operator=(double x) const{
    if(fabs(x) > origin.max)
      origin.max=x;
    v=x;
    return x;
  }
  double & passivated_reference(){return v;}
};

inline vector_with_max::active_reference_t 
vector_with_max::operator()(unsigned int i){
  return 
    active_reference_t(std::vector<double>::operator[](i),*this);
}

inline double
vector_with_max::operator()(unsigned int i) const{
  return std::vector< double >::operator[](i);
}

#endif // _VECTOR_WITH_MAX_H_
