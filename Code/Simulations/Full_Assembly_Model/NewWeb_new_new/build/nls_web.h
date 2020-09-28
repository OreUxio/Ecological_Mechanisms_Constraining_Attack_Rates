// -*- c++ -*-
//$Id: nls_web.h 247 2006-02-09 05:28:03Z cvsrep $

#ifndef __NLS_WEB_H__
#define __NLS_WEB_H__

#include <fstream>
#include <errno.h>
#include "sequence.h"
#include "error.h"

// vectors

class nls_vector : public sequence<double> {
public:
  nls_vector & operator+=(nls_vector & other){
    for(int i=(size()>other.size()?size():other.size());i-->0;){
      (*this)[i]+=other[i];
    }
    return *this;
  }
};

std::istream & operator>>(std::istream &stream, 
			  nls_vector &v);
    
// matrices

class nls_matrix : public sequence< sequence<double> > {
public:
  nls_matrix & operator+=(nls_matrix & other){
    for(int i=(size()>other.size()?size():other.size());i-->0;){
      (*this)[i]+=other[i];
    }
    return *this;
  }
};

std::istream & operator>>(std::istream &stream, 
			  nls_matrix &v);
    

class nls_web {
public:
  std::string header;
  int size;
  int number_of_living_compartments;
  sequence<std::string> name;
  nls_vector biomass;
  nls_vector input;
  nls_vector output;
  nls_vector respiration;
  nls_matrix flow;
  friend std::istream & operator>>(std::istream &stream, 
				   nls_web &w);
public:
  nls_web(){};
  nls_web(const char * filename){
    std::ifstream s(filename);
    if(!s){
      REPORT(strerror(errno));
      FATAL_ERROR("error while reading from file");
    }
    s >> *this;
  }
  nls_web & operator+=( nls_web & other);
  sequence<double> preference(int i);
  sequence<double> disposition(int i);
  sequence<double> inflows(int i);
  sequence<double> outflows(int i);
  bool is_lumped(int i);
  void strength_distribution(double th, const char * filename);
  nls_matrix trophic_intake_matrix();
};

// std::ostream & operator<<(std::ostream &stream, 
// 			  const nls_web &w);

std::istream & operator>>(std::istream &stream, 
			  nls_web &w);


#endif // __NLS_WEB_H__
