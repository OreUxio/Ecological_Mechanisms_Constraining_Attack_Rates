// $Id: random.h 2345 2013-12-13 10:52:12Z axel $

/// \file The project's interfact to random-number generators.

#ifndef __RANDOM_H__
#define __RANDOM_H__

// a fast random 32 bit interger
inline unsigned int random_int();

#include <vector>
#include <boost/random/mersenne_twister.hpp>

typedef boost::mt19937 myRandomEngine_t;
extern myRandomEngine_t myRandomEngine;

// a random 32 bit random integer 
inline unsigned int random_int(){
  return myRandomEngine();
}

// a uniformly distributed random number:
double unirand();

// a random interger uniformly distributed between 0 and upper_bound-1
int random_integer(int upper_bound);

// a random double uniformly distributed between 0 and upper_bound
double random_double(double upper_bound);

//**// starting form here we swich to using libCLHEP:

double gaussian(double mean, double std);

int poisson(double mu);

double random_exponential(double mean);

// sets the seed for all random number generators:
void set_random_seed(long int i);

/// member randint(i) generates random numbers < i
class randint{
public:
  unsigned operator()(unsigned i){
    return unsigned(i*unirand());
  }
};

/// Permutation, defined here for random permutations
class permutation : public std::vector<int> {
public:
  permutation(int size=0):std::vector<int>(size){};
  permutation inverse();
};

permutation random_permutation(int s);

#endif // __RANDOM_H__
