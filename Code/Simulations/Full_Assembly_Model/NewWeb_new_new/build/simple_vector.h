// -*- mode: c++ -*-
// $Id: simple_vector.h 2319 2013-08-04 16:58:58Z axel $
#ifndef _SIMPLE_VECTOR_H_
#define _SIMPLE_VECTOR_H_

// this vector should be simple enought to admit automatic parallelization,
// but complete enough to allow implementation of sequence<> templates

// currently the memory management is VERY inefficient!!

#include "error.h"

#undef simple_

namespace std {

#if 0
template<typename T>
/// Variant of STL vector more suitable for parallel computing.
class simple_vector{
public:
  typedef int index_t;
private:
  T * data;
  index_t the_size;
  index_t the_allocated_size;
public:
  simple_vector(const simple_vector<T> &other){
    index_t s=other.the_size;
    the_size=the_allocated_size=s;
    data = new T[s];
    ASSERT(data);
    for(index_t i=s;i-->0;){
      data[i]=other.data[i];
    }
  }
  simple_vector(const index_t s=0){
    the_size=the_allocated_size=s;
    data = new T[s];
    ASSERT(data);
    //if(the_size>0) ASSERT(data[0]==T()); //this fails for T=int!
    //for classes T this amounts to a double initialization.  Is there
    //a work-around?
    for(index_t i=the_size;i-->0;)
      data[i]=T();
  }
  simple_vector(const index_t s,const T & d){
    the_size=the_allocated_size=s;
    data = new T[s];
    ASSERT(data);
    //for classes T this amounts to a doulbe initialization.  Is there
    //a work-around?
    for(index_t i=the_size;i-->0;)
      data[i]=d;
  }
  ~simple_vector(){
    delete[] data;
  }
  simple_vector<T> & resize(index_t s);
  simple_vector<T> & operator=(const simple_vector<T> & other){
    if(&other != this){
      this->resize(other.the_size);
      for(index_t i=other.the_size;i-->0;){
	(*this)[i]=other[i];
      }
    }
    return * this;
  };
  bool operator==(const simple_vector<T> & other) const{
    if(the_size!=other.the_size)
      return false;
    for(index_t i=the_size;i-->0;){
      if(!((*this)[i]==other[i]))
	return false;
    }
    return true;
  };
  simple_vector<T> & prepend(const T & x){
    index_t s=the_size+1;
    if(s>the_allocated_size){
      // realloc should do it here, too.
      T * new_data = new T[s];
      for(index_t i=the_size;i-->0;)
	new_data[i+1]=data[i];
      delete[] data;
      data=new_data;
      the_allocated_size=s;
    }else{
      for(index_t i=the_size;i-->0;)
	data[i+1]=data[i];
    }
    the_size=s;
    data[0]=x;
    return *this;
  }
  T & operator[](index_t i){
    ASSERT(i<the_size);
    return data[i];
  }
  const T & operator[](index_t i) const{
    ASSERT(i<the_size);
    return data[i];
  }
  T & operator()(index_t i){
    ASSERT(i<the_size);
    return data[i];
  }
  const T & operator()(index_t i) const{
    ASSERT(i<the_size);
    return data[i];
  }
  index_t size() const {return the_size;};
  typedef T * iterator;
  typedef const T * const_iterator;
  T * begin() const {return data;};
  T * end() const {return data+the_size;};
  const T * begin() const {return data;};
  const T * end() const {return data+the_size;};
  bool empty() const {return the_size==0;};
};
  
template <typename T>
simple_vector<T> & simple_vector<T>::resize(index_t s){
  if(s>the_allocated_size){
    the_allocated_size=2*s;
    // realloc should do it here, too.
    T * new_data = new T[the_allocated_size];
    //for classes T this amounts to a doulbe initialization.  Is there
    //a work-around?
    for(index_t i=the_size;i-->0;)
      new_data[i]=data[i];
    delete[] data;
    data=new_data;
  }
  if(s>the_size){
    // fill with default element:
    //for classes T this amounts to a doulbe initialization.  Is there
    //a work-around?
    for(index_t i=s;i-->the_size;){
      data[i]=T();
    }
  }
  the_size=s;
  return *this;
}

#else
}
#include<vector>
namespace std {
  template<typename T>
  /// Variant of STL vector more suitable for parallel computing.
  class simple_vector : public vector<T>{
  public:
    simple_vector():vector<T>(){};
    simple_vector(size_t i):vector<T>(i){};
    simple_vector(size_t i,const T & x):vector<T>(i,x){};
    simple_vector<T> & prepend(const T & x){
      this->insert(this->begin(),x);
      return *this;
    };
    T & operator()(size_t i){
      return vector<T>::operator[](i);
    };
    const T & operator()(size_t i) const{
      return vector<T>::operator[](i);
    };
  };
#endif
  
}; // namespace std

#endif // _SIMPLE_VECTOR_H_


