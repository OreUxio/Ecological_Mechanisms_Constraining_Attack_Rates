// -*- mode: c++ -*-
// $Id: SortedMatrix.h 2108 2011-05-17 00:02:55Z axel $
#ifndef _SORTEDMATRIX_H_
#define _SORTEDMATRIX_H_

//#include <set>
#include <vector>
#include <algorithm>
#include <hash_map>
#include <float.h>

#include <boost/serialization/split_member.hpp>

#include "vector_with_max.h"

/// Matrix class optimized for situation where most entries are small.
/** This allows a fast computation of sandwich products of the form
    a^T*M*b, where a and b are vectors, and M is a SortedMatrix. This
    works by representing the matrix as a list of entries sorted by
    falling magnitude. When computing a^T*M*b, only the largest
    contributions are taken into account. */ 
class SortedMatrix 
{
public:  
  typedef short int index_t;
private:

  /// Row and column of a matrix entry.
  class __attribute__ ((packed)) location_t{
  public:
      index_t row,column;
    location_t(index_t r,index_t c):row(r),column(c){};
    bool operator==(const location_t & other) const{
      return row==other.row && column==other.column;
    }
    bool unused(){return row<0;};
  };
  struct location_t_hash{
    size_t operator()( const location_t& l ) const
    {
      const size_t someprime=3571;
      return size_t(l.row)^(size_t(l.column)*someprime);
      //return size_t(l.row)+(size_t(l.column)<<(8*sizeof(index_t)));
    }
  };

  /// A matrix entry (value and location)
  class entry_t:public location_t{
  public:
    double value; 
    entry_t(const location_t & l,double v):location_t(l),value(v){};
    entry_t(index_t r, index_t c,double v):location_t(r,c),value(v){};
  };
  
  class larger_value_than{
  public:
    inline bool operator()(const entry_t e1,const entry_t e2);
  };


  class has_row_or_colum_equal_to_either_of{
    typedef std::pair<SortedMatrix::location_t, int>  arg_t;
  private:
    index_t a,b;
  public:
    has_row_or_colum_equal_to_either_of(index_t A,index_t B):
      a(A),b(B) {};
    bool operator()(arg_t lp) const{
      return 
	(lp.first.row==a) || (lp.first.column==a) || 
	(lp.first.row==b) || (lp.first.column==b);
    }
  };

  class has_location{
    const location_t loc;
  public:
    has_location(location_t l):loc(l){};
    bool operator()(entry_t other){
      return location_t(other)==loc;
    }
  };
    
  //typedef std::set< entry_t, larger_value_than > _container;
  typedef std::vector< entry_t > _container;
  mutable _container _entries;
  
  typedef __gnu_cxx::hash_map< location_t, int, location_t_hash > locator_t;
  mutable locator_t locator;

  int _size;
  double _accuracy; // controls when to stop summations
  double _truncation_epsilon; // controls which matrix elements to keep
public:
  typedef enum {asymmetric,symmetric} symmetry_t;
private:
  symmetry_t _symmetry;
  mutable bool _clean;
  
private:
  friend class boost::serialization::access;
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    cleanup();
    ar & _accuracy & _truncation_epsilon & _size;
    ar & _entries;
  }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar & _accuracy & _truncation_epsilon & _size;
    ar & _entries;
    _clean=false;
    cleanup();
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()


private:
  void cleanup_helper() const;
public:
  inline void cleanup() const;
  static double &default_accuracy,&default_truncation_epsilon;
  static const entry_t unused_entry; 
  class active_reference_t;
  class const_active_reference_t;
  class incomplete_location_t;
public:
  SortedMatrix(symmetry_t s=asymmetric);
  // access to individual matrix entries is expensive for
  // this class...
  inline incomplete_location_t operator[](index_t r);
  int size() const {return _size;};
  void resize(int size);
  void move(index_t from, index_t to);
  void clear(){_entries.clear(); locator.clear();_clean=true;}
  // ...but multiplication of vectors from both sides can be very fast:
  double sandwich_product(const double* v1, double max_v1, 
			  const double* v2, double max_v2)const;
  inline 
  double sandwich_product(const vector_with_max & v1,
			  const vector_with_max & v2)const;
  bool operator==(const SortedMatrix & other);
  double truncation_epsilon(){return _truncation_epsilon;}
  double accuracy() const {return _accuracy;}
  friend std::ostream & operator<<(std::ostream &stream, const SortedMatrix &m);
  friend class packed_simulation;
};

// Inline implementation details:

/// Updates a SortedMatrix upon assignment of entries.
class SortedMatrix::active_reference_t{
    const location_t loc;
    SortedMatrix & origin;
    double get_value() const;
    void set_value(double v) const;
  public:
    active_reference_t(const location_t l, SortedMatrix & m):loc(l),origin(m){};
    operator double () const{
      return get_value();
    }
    double operator=(double v) const{
      set_value(v);
      return v;
    }
  };
class SortedMatrix::incomplete_location_t{//only the row is known
  const index_t row;
  SortedMatrix & origin;
public:
  explicit incomplete_location_t(index_t r,SortedMatrix & o):row(r),origin(o){};
  active_reference_t operator[](index_t c) const{
    return active_reference_t(location_t(row,c),origin);
  }
};


inline SortedMatrix::incomplete_location_t SortedMatrix::operator[](index_t r){
  return incomplete_location_t(r,*this);
}
  
inline bool 
SortedMatrix::larger_value_than::operator()(const entry_t e1,const entry_t e2){
  if(fabs(e1.value) != fabs(e2.value))
    return e1.value>e2.value;
  if(e1.row != e2.row)
    return e1.row < e2.row;
  return e1.column < e2.column;
}
inline void SortedMatrix::cleanup() const{
  if(!_clean) cleanup_helper();
}

inline double 
SortedMatrix::sandwich_product(const vector_with_max & v1,
			       const vector_with_max & v2)const{
  return sandwich_product(v1.get_data(),v1.get_max(),
			  v2.get_data(),v2.get_max());
}


void SortedMatrix_tester();

#endif // _SORTEDMATRIX_H_
