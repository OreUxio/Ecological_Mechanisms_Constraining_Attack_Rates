// -*- mode: c++ -*-
// $Id: SortedVector.cc 2011 2010-11-30 23:45:14Z axel $

#include "SortedVector.h"
#include <iostream>
#include <algorithm>
#include <math.h>
#include <boost/preprocessor/repetition.hpp>
#include "error.h"
#include "Statistics.h"

using namespace std;

static double vector_accuracy=DBL_EPSILON;
static double vector_truncation_epsilon=0;

//#define SOVE_MEASURE_EFFICIENCY
#define SOVE_UNROLL_LOOP

#ifndef SOVE_ITERATIONS_UNROLLED // number of iterations to unroll
#define SOVE_ITERATIONS_UNROLLED 8
#endif

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  CFGDOUBLE(vector_accuracy),
  CFGDOUBLE(vector_truncation_epsilon),
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);


double & SortedVector::default_accuracy=vector_accuracy;  
double & SortedVector::default_truncation_epsilon=
  vector_truncation_epsilon;

const SortedVector::entry_t SortedVector::unused_entry=
  entry_t(location_t(-1),0);


SortedVector::SortedVector():
  _size(0),
  _accuracy(default_accuracy),
  _truncation_epsilon(default_truncation_epsilon),
  _clean(true)
{
}

double SortedVector::get_value(const location_t & loc) const{
  locator_t::iterator li=locator.find(loc);
  if(li==locator.end()){
    return 0;
  }else{
    return _entries[li->second].value;
  }
}

void SortedVector::set_value(double v,const location_t & loc){
  locator_t::iterator li=locator.find(loc);
  if(li!=locator.end()){
    entry_t & e=_entries[li->second];
    if(e.value==v)
      return;
    // delete entry if it was there already
    locator.erase(li);
    e=unused_entry;
    _clean=false;
  }
  if(fabs(v)>_truncation_epsilon){
    _entries.push_back(entry_t(loc,v));
    locator[loc]=_entries.size()-1;
    _clean=false;
  }
}

double SortedVector::get(int i){
  return operator[](i);
}

// inline
// bool 
// SortedVector::larger_value_than::operator()(const entry_t e1,const entry_t e2){
//   if(fabs(e1.value) != fabs(e2.value))
//     return e1.value>e2.value;
//   return e1.row < e2.row;
// }
  
void SortedVector::cleanup_helper() const {
  sort(_entries.begin(),_entries.end(),larger_value_than());
  
  //erase unused entries
  if(_entries.begin()!=_entries.end()){
    if(_entries.begin()->unused()){
      _entries.clear();
      locator.clear();
    }else{
      _container::iterator l=_entries.end();
      while((--l)->unused() ) 
	locator.erase(locator.find(*l));
      
      _entries.erase(++l,_entries.end());
    }
  }
  
  //     locator.clear(); // This line is expensive.  If this helps to
  // 		     // remove the sudden death phenomenon, try to do
  // 		     // this more efficient.
  for(int i=0;i<_entries.size();++i){
    locator[_entries[i]]=i;
  }
  _clean=true;
}


#define SOVE_ONE_ITERATION(z, n, unused)	\
sum+=i->value*v1[i->row];			\
++i;


double 
SortedVector::dot(const double* v1, double max_v1) const {
#ifdef SOVE_MEASURE_EFFICIENCY
  static average_meter av_turns;
  static average_meter av_speedup;
  static average_meter av_error;
#endif

  cleanup();
  double sum=0;
  const double max=max_v1;
  const double factor=max*(1.0/_accuracy);

  _container::const_iterator i=_entries.begin();
  const _container::const_iterator end=_entries.end();
#ifdef SOVE_MEASURE_EFFICIENCY
  int t=0;
#endif
  while(i!=end){
#ifdef SOVE_UNROLL_LOOP
    //Unroll loop heavily:
    if(i+SOVE_ITERATIONS_UNROLLED<end){
#ifdef SOVE_MEASURE_EFFICIENCY
      t+=SOVE_ITERATIONS_UNROLLED;
#endif
      // Repeated insertion of one iteration at preprocessor stage:
      BOOST_PP_REPEAT(SOVE_ITERATIONS_UNROLLED, SOVE_ONE_ITERATION, ~);
    }
#endif
#ifdef SOVE_MEASURE_EFFICIENCY
    t++;
#endif
    const entry_t & e=*i;
    const double val=e.value;
    sum+=val*v1[e.row];
    if(val*factor<sum)
      break;
    ++i;
  }
#ifdef SOVE_MEASURE_EFFICIENCY
  double error=0;
  while(i!=end){
    const entry_t & e=*i;
    const double val=e.value;
    error+=val*v1[e.row];
    ++i;
  }
  if(sum)
    av_error.sample(error/sum);
  av_turns.sample(t);
  if(t)
    av_speedup.sample(_entries.size()/double(t));
  REPORT(_entries.size());
  REPORT(t);
  REPORT(av_turns);
  REPORT(av_speedup);
  REPORT(av_error);
#endif
  return sum;
}

//// inlined
// double 
// SortedVector::dot(const vector_with_max & v1)const{
//   return dot(v1.get_data(),v1.get_max());
// }

bool SortedVector::operator==(const SortedVector & other){
  cleanup();
  other.cleanup();
  return this->_entries==other._entries;
}

void SortedVector::resize(int new_size){
  if(new_size<_size){
    for(_container::iterator i=_entries.begin();
      i!=_entries.end();++i){
      if(i->row >= new_size){
	locator.erase(locator.find(*i));
	*i=unused_entry;
	_clean=false;
      }
    }
  }
  _size=new_size;
  return;
}

void SortedVector::move(int from,int to){ 
  if(from==to) return;
  for(_container::iterator i=_entries.begin();
      i!=_entries.end();++i){
    if(i->row == to){
      locator.erase(locator.find(*i));
      *i=unused_entry;
      _clean=false;
    }else{ 
      if(i->row == from){
	locator.erase(locator.find(*i));
	i->row=to;
	locator[*i]=i-_entries.begin();
      }
    }
  }
}

std::ostream & operator<<(std::ostream &stream, const SortedVector &m){
  m.cleanup();
  stream << m._entries.size() << " entries" << std::endl;
  for(SortedVector::_container::iterator i=m._entries.begin();
      i!=m._entries.end();){
    stream << i->value << " " 
	   << i->row << " ";
    ++i;
    stream << std::endl;
  }
  return stream;
}

