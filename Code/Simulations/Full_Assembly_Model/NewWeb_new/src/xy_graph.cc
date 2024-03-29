// -*- mode: c++ -*-
// $Id: xy_graph.cc 564 2006-07-30 00:00:49Z cvsrep $

#include <fstream>
#include <math.h>
#include "xy_graph.h"
#include "error.h"

void xy_graph::save(const char * filename){
  std::ofstream os(filename);
  os << *this;
}
void xy_graph::load(const char * filename){
  std::ifstream is(filename);
  is >> *this;
}

std::ostream & operator<<(std::ostream &stream, const xy_graph &g){
  stream << 
    g.the_x.format("%5g ")+
    g.the_y.format("%5g ");
  return stream;
}

std::istream & operator>>(std::istream &stream, xy_graph &g){
  int i=0;
  do{
    stream >> g.the_x[i] >> g.the_y[i++];
  }while(!stream.fail());
  g.the_x.resize(i-1);
  g.the_y.resize(i-1);
  return stream;
}

int xy_graph::size() const{
  return ( the_x.size() > the_y.size() ? the_x.size() : the_y.size() );
}

log_spectrum::log_spectrum(const xy_graph & data, double bin_factor){
  if(data.size()==0)
    return;

  //get absolute bin locations and their minimum:
  sequence<int> abs_bin;
  ALWAYS_ASSERT(data.get_x(0)>0);
  int min_abs_bin=int(floor(log(data.get_x(0))/log(bin_factor)));
  abs_bin[0]=min_abs_bin;
  for(int i=data.size();i-->1;){
    ALWAYS_ASSERT(data.get_x(i)>0);
    abs_bin[i]=int(floor(log(data.get_x(i))/log(bin_factor)));
    if(abs_bin[i]<min_abs_bin){
      min_abs_bin=abs_bin[i];
    }
  }

  // fill the bins:
  for(int i=data.size();i-->0;){
    the_y[abs_bin[i]-min_abs_bin]+=data.get_y(i);
  }
  
  int n=the_y.size();

  // normalize and compute bin boundaries:
  the_lower_end=pow(bin_factor,min_abs_bin);
  double lower_bound=the_lower_end;
  double fac=bin_factor-1;
  for(int i=0;i<n;i++){
    the_x[i]=lower_bound;
    the_y[i]/=lower_bound*fac;
    lower_bound*=bin_factor;
  }
  
  // double the last point to get nice histograms:
  the_x[n]=lower_bound; // now it's the upper bound
  the_y[n]=the_y[n-1];
  
  return;
}

#ifndef ON_SX5FSV
xy_graph xy_graph::log_xy(){
  return xy_graph(Map(::log10,the_x),Map(::log10,the_y));
}
#endif
