// -*- c++ -*-
//$Id: sequence.h 2478 2016-10-23 11:40:59Z axel $

#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>
#include "simple_vector.h"

template<typename T>
std::string format(const char * fmt, T x);

/// A slow but fool-proof extension of the STL vector.
/** Observe: if \a v is a simple_vector, than v[n] (with n>=0) is
    always OK, but v(n) only if v.size()>n . */
template <typename T>
class sequence : public std::simple_vector<T> {
private:
  static T default_element;
  typedef typename std::simple_vector<T>::iterator iter;
public:
  sequence<T>() : std::simple_vector<T>(){};
  sequence<T>(const sequence<T> & s) : std::simple_vector<T>(s){};
  sequence<T>(const std::simple_vector<T> & v) : std::simple_vector<T>(v){};
  sequence<T>(int i,T x=T()) : std::simple_vector<T>(i,x){/*shorten() suppressed*/};
  T & operator[](unsigned int i) {
    if(i >= std::simple_vector<T>::size())
      std::simple_vector<T>::resize(i+1);
    return std::simple_vector<T>::operator[](i);
  }
  const T & operator[](unsigned int i) const{
    if(i >= std::simple_vector<T>::size()){
      return default_element;
    }else{
      return std::simple_vector<T>::operator[](i);
    }
  }
  T & operator()(unsigned int i) {
    ASSERT(i < std::simple_vector<T>::size());
    return std::simple_vector<T>::operator[](i);
  }
  const T & operator()(unsigned int i) const{
    ASSERT(i < std::simple_vector<T>::size());
    return std::simple_vector<T>::operator[](i);
  }
  void shorten(){
    while( std::simple_vector<T>::rbegin()!=std::simple_vector<T>::rend() && 
	   *std::simple_vector<T>::rbegin()==default_element ) 
      std::simple_vector<T>::pop_back();
  }
  sequence<T> & operator=(const sequence<T> & v);
  template<typename D>
  sequence<T> & operator=(const sequence<D> & v);
  template<typename D>
  sequence<T> & operator+=(const sequence<D> & v);
  template<typename D>
  sequence<T> & operator*=(const sequence<D> & v);
  template<typename D>
  sequence<T> & operator-=(const sequence<D> & v);
  template<typename D>
  sequence<T> & operator/=(const sequence<D> & v);
  template<typename D>
  sequence<T> & operator+=(D v);
  template<typename D>
  sequence<T> & operator-=(D v);
  template<typename D>
  sequence<T> & operator*=(D v);
  template<typename D>
  sequence<T> & operator/=(D v);
  sequence<T> operator+(T y) const {
    return sequence<T>(*this)+=y;
  }
  sequence<T> operator-(T y) const {
    return sequence<T>(*this)-=y;
  }
  sequence<T> operator*(T y) const {
    return sequence<T>(*this)*=y;
  }
  sequence<T> operator/(T y) const{
    return sequence<T>(*this)/=y;
  }
  sequence<T> operator+(sequence<T> y) const {
    return sequence<T>(*this)+=y;
  }
  sequence<T> operator-(sequence<T> y) const {
    return sequence<T>(*this)-=y;
  }
  sequence<T> operator*(sequence<T> y) const {
    return sequence<T>(*this)*=y;
  }
  sequence<T> operator/(sequence<T> y) const{
    return sequence<T>(*this)/=y;
  }
  sequence<T> reverse_cumulative_sum() const ;
  template<typename D>
  operator sequence<D> () const;
  sequence<std::string> format(const char * fmt) const;
  void push_back(const T & y){
    (*this)[this->std::simple_vector<T>::size()]=y;
  }
  T last(){return std::simple_vector<T>::operator[](std::simple_vector<T>::size()-1);}
};

template<typename T>
T sequence<T>::default_element=T();


template<typename T>
template<typename D>
sequence<T>::operator sequence<D> () const{
  sequence<D> d;
  for(int i=std::simple_vector<T>::size()-1;i>=0;i--) // doing backwards saves mallocs
    d[i]=D(std::simple_vector<T>::operator[](i));
  return d;
}

template<typename T>
sequence<std::string> sequence<T>::format(const char * fmt) const{
  sequence<std::string> s;
  for(int i=std::simple_vector<T>::size()-1;i>=0;i--) // doing backwards saves mallocs
    s[i]=::format(fmt,std::simple_vector<T>::operator[](i));
  return s;
}

template<typename T>
std::string format(const char * fmt, T x){
  const int default_limit=32;
  static char buffer[default_limit]; 
  int retval=snprintf(buffer,default_limit,fmt,x);
  if(retval<0) return std::string("ERROR: sprintf trouble");
  if(retval<default_limit)
    return buffer;
  // buffer was too small:
  char *buffer2;
  buffer2 = new char[retval+1]; // one more for the trailing zero
  retval=sprintf(buffer2,fmt,x);
  std::string buffer2_string(buffer2);
  delete buffer2;
  if(retval<0) return std::string("ERROR: sprintf trouble");
  return buffer2_string;
}

template<typename T>
sequence<T> & sequence<T>::operator=(const sequence<T> & v){
  if(std::simple_vector<T>::size()!=v.std::simple_vector<T>::size())
    this->resize(v.std::simple_vector<T>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<T>::const_iterator j=v.std::simple_vector<T>::begin(),e=v.std::simple_vector<T>::end();
  while(j!=e) {
    (*i++)=(*j++);
  }
  //shorten();
  return *this;
}

template<typename T> template<typename D>
sequence<T> & sequence<T>::operator=(const sequence<D> & v){
  if(std::simple_vector<T>::size()!=v.std::simple_vector<D>::size())
    resize(v.std::simple_vector<D>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<D>::const_iterator j=v.std::simple_vector<D>::begin(),e=v.std::simple_vector<D>::end();
  while(j!=e) {
    (*i++)=(*j++);
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator+=(const sequence<D> & v){
  if(std::simple_vector<T>::size()<v.std::simple_vector<D>::size())
    this->resize(v.std::simple_vector<D>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<D>::const_iterator j=v.std::simple_vector<D>::begin(),e=v.std::simple_vector<D>::end();
  while(j!=e) {
    (*i)=T(*i + (*j++));
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator*=(const sequence<D> & v){
  if(std::simple_vector<T>::size()<v.std::simple_vector<T>::size())
    this->resize(v.std::simple_vector<T>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<D>::const_iterator j=v.std::simple_vector<T>::begin(),e=v.std::simple_vector<T>::end();
  while(j!=e) {
    (*i)=T(*i * (*j++));
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator-=(const sequence<D> & v){
  if(std::simple_vector<T>::size()<v.std::simple_vector<T>::size())
    this->resize(v.std::simple_vector<T>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<D>::const_iterator j=v.std::simple_vector<T>::begin(),e=v.std::simple_vector<T>::end();
  while(j!=e) {
    (*i)=T(*i - (*j++));
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator/=(const sequence<D> & v){
  if(std::simple_vector<T>::size()<v.std::simple_vector<T>::size())
    this->resize(v.std::simple_vector<T>::size());
  iter i=std::simple_vector<T>::begin();
  typename std::simple_vector<D>::const_iterator j=v.std::simple_vector<T>::begin(),e=v.std::simple_vector<T>::end();
  while(j!=e) {
    (*i)=T(*i / (*j++));
    i++;
  }
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator+=(D v){
  iter i=std::simple_vector<T>::begin(),e=std::simple_vector<T>::end();
  while(i!=e) {
    (*i)=T(*i + v);
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator-=(D v){
  iter i=std::simple_vector<T>::begin(),e=std::simple_vector<T>::end();
  while(i!=e) {
    (*i)=T(*i - v);
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator*=(D v){
  iter i=std::simple_vector<T>::begin(),e=std::simple_vector<T>::end();
  while(i!=e) {
    (*i)=T(*i * v);
    i++;
  }
  //shorten();
  return *this;
}
template<typename T> template<typename D>
sequence<T> & sequence<T>::operator/=(D v){
  iter i=std::simple_vector<T>::begin(),e=std::simple_vector<T>::end();
  while(i!=e) {
    (*i)=T(*i / v);
    i++;
  }
  return *this;
}

template<typename T>
sequence<T> sequence<T>::reverse_cumulative_sum() const{
  sequence<T> s;
  T sum=0;
  for(int i=std::simple_vector<T>::size()-1;i>=0;i--) 
    s[i]=(sum+=std::simple_vector<T>::operator[](i));
  return s;
}
  

template<typename T> 
inline sequence<T> operator+(T x, const sequence<T> & y){
  return sequence<T>(y.std::simple_vector<T>::size(),x)+=y;
}

template<typename T> 
inline sequence<T> operator-(T x, const sequence<T> & y){
  return sequence<T>(y.std::simple_vector<T>::size(),x)-=y;
}

template<typename T> 
inline sequence<T> operator*(T x, const sequence<T> & y){
  return sequence<T>(y.std::simple_vector<T>::size(),x)*=y;
}

template<typename T> 
inline sequence<T> operator/(T x, const sequence<T> & y){
  return sequence<T>(y.std::simple_vector<T>::size(),x)/=y;
}

template<typename T>
std::ostream & operator<<(std::ostream &stream, const sequence<T> &s){
  typename std::simple_vector<T>::const_iterator i=s.std::simple_vector<T>::begin(),e=s.std::simple_vector<T>::end();
  while(i!=e) {
    stream << *i << " " ;
    i++;
  }
  return stream;
}

template<typename T>
std::istream & operator>>(std::istream &stream, sequence<T> &s){
  int i=0;
  do{
    stream >> s[i++];
  }while(!stream.fail());
  s.resize(i-1);
//   std::string line;
//   stream >> line; //read one line
//   std::istringstream is(line);
//   int i=0;
//   std::cout << line << " " << is.eof() << std::endl;
//   while(!is.eof()){
//     s.std::simple_vector<T>::resize(i+1);
//     is >> s.std::simple_vector<T>::operator[](i);
//     std::cout << i << " " << s[i] << std::endl;
//     i++;
//   }
  return stream;
}

// strings get vertical output:
std::ostream & operator<<(std::ostream &stream, 
			  const sequence<std::string> &s);

template<typename T, typename D>
sequence<D> Map(D f(T),  const sequence<T> t){
  sequence<D> d;
  for(int i=t.std::simple_vector<T>::size()-1;i>=0;i--) 
    d[i]=f(t[i]);
  return d;
}

#include<numeric>

template<typename T>
T sum(sequence<T> x){
  return std::accumulate(x.begin(), x.end(), T());
}

template<typename T>
T mean(sequence<T> x){
  return accumulate(x.begin(), x.end(), T())/x.std::simple_vector<T>::size();
}

inline sequence<std::string> operator+(const char c[], const sequence<std::string> & y){
  return std::string(c)+y;
}

// Some simple abbreviations.  Should not really be here, but then
// again, it's ok:

template< typename T1, typename T2 > double dot(const T1 &x,const T2 &y){
  return std::inner_product(x.begin(),x.end(),y.begin(),0.0);
}
template< typename T1 > double abs2(const T1 &x){
  return dot(x,x);
}

#include <math.h>
template< typename T1 > double abs(const sequence< T1 > &x){
  return sqrt(abs2(x));
}



#endif // __SEQUENCE_H__
