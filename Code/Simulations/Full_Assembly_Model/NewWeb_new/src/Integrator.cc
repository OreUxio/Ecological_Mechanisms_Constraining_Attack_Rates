// -*- mode: c++ -*-
// $Id: Integrator.cc 2067 2011-01-26 19:22:44Z tak $

#include "Integrator.h"


// This FORTRAN code should give some inspiration on how to do the
// integration better.

// subroutine simpne ( x, y, num, result )
// !
// !***********************************************************************
// !
// !! SIMPNE approximates the integral of unevenly spaced data.
// !
// !
// !  Discussion:
// !
// !    The routine repeatedly interpolates a 3-point Lagrangian polynomial 
// !    to the data and integrates that exactly.
// !
// !  Reference:
// !
// !    Philip Davis and Philip Rabinowitz,
// !    Methods of Numerical Integration,
// !    Blaisdell Publishing, 1967.
// !
// !  Modified:
// !
// !    30 October 2000
// !
// !  Parameters:
// !
// !    Input, real X(NUM), contains the X values of the data, in order.
// !
// !    Input, real Y(NUM), contains the Y values of the data.
// !
// !    Input, integer NUM, number of data points.  NUM must be at least 3.
// !
// !    Output, real RESULT.
// !    RESULT is the approximate value of the integral.
// !
//   implicit none
// !
//   integer num
// !
//   real del(3)
//   real e
//   real f
//   real feints
//   real g(3)
//   integer i
//   integer n
//   real pi(3)
//   real result
//   real sum1
//   real x(num)
//   real x1
//   real x2
//   real x3
//   real y(num)
// !
//   result = 0.0E+00
 
//   if ( num <= 2 ) then
//     write ( *, '(a)' ) ' '
//     write ( *, '(a)' ) 'SIMPNE - Fatal error!'
//     write ( *, '(a)' ) '  NUM <= 2.'
//     stop
//   end if
 
//   n = 1
 
//   do
 
//     x1 = x(n)
//     x2 = x(n+1)
//     x3 = x(n+2)
//     e = x3*x3-x1*x1
//     f = x3*x3*x3-x1*x1*x1
//     feints = x3-x1
//     del(1) = x3-x2
//     del(2) = x1-x3
//     del(3) = x2-x1
//     g(1) = x2+x3
//     g(2) = x1+x3
//     g(3) = x1+x2
//     pi(1) = x2*x3
//     pi(2) = x1*x3
//     pi(3) = x1*x2
 
//     sum1 = 0.0E+00
//     do i = 1, 3
//       sum1 = sum1 + y(n-1+i)*del(i)*(f/3.0E+00-g(i)*0.5E+00*e+pi(i)*feints)
//     end do
//     result = result - sum1 / ( del(1) * del(2) * del(3) )
 
//     n = n+2

//     if ( n + 1 >= num ) then
//       exit
//     end if

//   end do
 
//   if ( mod(num,2) /= 0 ) then
//     return
//   end if

//   n = num-2
//   x3 = x(num)
//   x2 = x(num-1)
//   x1 = x(num-2)
//   e = x3*x3-x2*x2
//   f = x3*x3*x3-x2*x2*x2
//   feints = x3-x2
//   del(1) = x3-x2
//   del(2) = x1-x3
//   del(3) = x2-x1
//   g(1) = x2+x3
//   g(2) = x1+x3
//   g(3) = x1+x2
//   pi(1) = x2*x3
//   pi(2) = x1*x3
//   pi(3) = x1*x2
 
//   sum1 = 0.0E+00
//   do i = 1, 3
//     sum1 = sum1 + y(n-1+i) * del(i) * &
//       ( f / 3.0E+00 - g(i) * 0.5E+00 * e + pi(i) * feints )
//   end do
 
//   result = result - sum1 / ( del(1) * del(2) * del(3) )
 
//   return
// end

#ifdef USE_SIMPLE_INTEGRATOR

Integrator::Integrator(const sequence< double * > & l):
  the_locations(l),
  write_back(true),
  started(false)
{
  the_integrals.resize(the_locations.size());
  y1.resize(the_locations.size());
}

void Integrator::sample(double t){
  if(started){
    const double dt=t-last_t;
    for(int i=the_locations.size();i-->0;){
      double y2=*the_locations(i);
      the_integrals(i)+=(y1(i)+y2)*dt/2;
      y1(i)=y2;
    }
    last_t=t;
  }else{
    started=true;
    start_t=t;
    last_t=t;
    for(int i=the_locations.size();i-->0;){
      y1(i)=*the_locations(i);
    }
  }
}

#else // use complex integrator

Integrator::Integrator(const sequence< double * > & l):
  the_locations(l),
  write_back(true),
  num(0)
{
  the_integrals.resize(the_locations.size());
  y1.resize(the_locations.size());
  y2.resize(the_locations.size());
}

void Integrator::sample(double t){
  ++num;
  if(num==1){
    start_t=x1=t;
    for(int j=the_locations.size();j-->0;){
      y1(j)=*the_locations(j);
    }
    return;
  }

  if(num%2==0){//num is even
    x2=t;
    for(int j=the_locations.size();j-->0;){
      y2(j)=*the_locations(j);
      the_integrals(j) += (x2-x1)*(y1(j)+y2(j))/2;//preliminary guess
    }
  }else{//num is odd
    const double x3 = t;
    const double e = x3*x3-x1*x1;
    const double f = x3*x3*x3-x1*x1*x1;
    const double feints = x3-x1;
    sequence<double> del(4),g(4),pi(4);
    del(1) = x3-x2;
    del(2) = x1-x3;
    del(3) = x2-x1;
    g(1) = x2+x3;
    g(2) = x1+x3;
    g(3) = x1+x2;
    pi(1) = x2*x3;
    pi(2) = x1*x3;
    pi(3) = x1*x2;
    
    for(int j=the_locations.size();j-->0;){
      the_integrals(j) -= (x2-x1)*(y1(j)+y2(j))/2;//clear preliminary guess

      double y3=*the_locations(j);
      double sum1 = 0;
      sum1 += y1(j)*del(1)*(f/3-g(1)/2*e+pi(1)*feints);
      sum1 += y2(j)*del(2)*(f/3-g(2)/2*e+pi(2)*feints);
      sum1 += y3   *del(3)*(f/3-g(3)/2*e+pi(3)*feints);
      the_integrals(j) -= sum1 / ( del(1) * del(2) * del(3) );

      y1(j)=y3;
    }
    x1=t;
  }
  last_t=t;
}

#endif

void Integrator::operator*=(double x){
  the_integrals*=x;
}

void Integrator::disable_write_back(){
  write_back=false;
}

Integrator::~Integrator(){
  if(write_back){
      for(int i=the_locations.size();i-->0;){
	*the_locations[i]=the_integrals[i];
      }
  }
}

sequence< double > Integrator::get_integrals(){
  return the_integrals;
}

void Integrator::reset(){
  for(int i=the_integrals.size();i-->0;){
    the_integrals[i]=0;
    started=false;
  }
}



