// $Id: random.cc 2345 2013-12-13 10:52:12Z axel $
#include <stdlib.h>
#include <math.h>

#include <algorithm>

#include "random.h"
#include "error.h"

#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include "boost/random/variate_generator.hpp"

using namespace boost;
using namespace boost::random;

// a random interger uniformly distributed between 0 and upper_bound-1
int random_integer(int upper_bound){
  int i = int(upper_bound*unirand());
  TRACE(i,RANDOM);
  return i;
}

double unirand(){
//  return double(random())/RAND_MAX;
  static uniform_01<> u;
  double uni=u(myRandomEngine);
  TRACE(uni,RANDOM);
  return uni;
}


// a random double uniformly distributed between 0 and upper_bound
double random_double(double upper_bound){
  double r=upper_bound*unirand();
  TRACE(r,RANDOM);
  return r;
}

double gaussian(double mean, double std){
  static normal_distribution<> dist;
  static 
    variate_generator<myRandomEngine_t&,
		      boost::normal_distribution<> > 
    n(myRandomEngine, dist);
  double gaussian = mean + std * n();
  TRACE(gaussian,RANDOM);
  return gaussian;
}

int poisson(double mu){
  poisson_distribution<> dist(mu);
  variate_generator<myRandomEngine_t&,
		    boost::poisson_distribution<> > 
    p(myRandomEngine, dist);
  int poisson = p();
  TRACE(poisson,RANDOM);
  return poisson;
}

double random_exponential(double mean){
  exponential_distribution<> dist(1/mean);
  variate_generator<myRandomEngine_t&,
		    boost::exponential_distribution<> > 
    e(myRandomEngine, dist);
  double exponential = e();
  TRACE(exponential,RANDOM);
  return exponential;
}

myRandomEngine_t myRandomEngine;


void set_random_seed(long int seed){
  myRandomEngine.seed(seed);
  TRACE(seed,RANDOM);
};

class random_ranking {
private:
  std::vector<double> ranking;
public:
  random_ranking(int size):ranking(std::vector<double>(size)){
    for(int i=0;i<size;i++){
      ranking[i]=random_double(1);
    }
  }
  bool operator()(int s1, int s2){
    return 
      (ranking[s1] >= ranking[s2] );
  }
};

permutation random_permutation(int s){
  // sort by mass
  permutation perm(s);
  for(int i=0;i<s;i++){
    perm[i]=i;
  }
  stable_sort(perm.begin(),perm.end(),random_ranking(s));
  return perm;
}

permutation permutation::inverse(){
  int s=this->size();
  permutation p2(s);
  for(int i=0;i<s;i++)
    p2[(*this)[i]]=i;
  return p2;
}

