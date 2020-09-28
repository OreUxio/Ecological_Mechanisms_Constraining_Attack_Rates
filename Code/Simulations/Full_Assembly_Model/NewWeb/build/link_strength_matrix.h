// -*- mode: c++ -*-
// $Id: link_strength_matrix.h 2466 2016-05-01 23:27:44Z axel $
#ifndef _LINK_STRENGTH_MATRIX_H_
#define _LINK_STRENGTH_MATRIX_H_

// /// a matrix of double which automatically adjusts its size:
// class link_strength_matrix : public sequence< sequence<double> > {
// public:
//   link_strength_matrix(){};
//   link_strength_matrix(const sequence< sequence<double> > &m):
//     sequence< sequence<double> >(m){};
//   link_strength_matrix(const Interaction_Matrix &m);
//   operator const HepMatrix ();
//   link_strength_matrix operator=(const link_strength_matrix & m){
//     return sequence< sequence<double> >::operator=(m);
//   }
// };

#include "NetworkAnalysis.h"
#include "sequence.h"

/// a matrix of double optimized for speed
typedef CMatrix<double,Container2DRow<double> > link_strength_matrix_base ;
class link_strength_matrix : public link_strength_matrix_base {
  int the_size;
  int make_odd(int i){
    return i+1-(i&1);
  }
public:
  link_strength_matrix();
  link_strength_matrix(const sequence< sequence<double> > &m);
  link_strength_matrix(const link_strength_matrix_base &m);
  link_strength_matrix(const Interaction_Matrix &m);
  int size() const{return the_size;}
  void resize(int requested_size);
  operator const NewMatrix () const;
  link_strength_matrix operator+(const link_strength_matrix & y)const;
  link_strength_matrix operator*(const double y)const;
  double & operator()(int i,int j){
    // !! indices are exchanged to get correspondence with paper.
    // There the indices are like this because energy flows go from i
    // to j.
    return (*this)[j][i];
  }
  const double & operator()(int i,int j) const {
    // !! indices are exchanged to get correspondence with paper.
    // There the indices are like this because energy flows go from i
    // to j.
    return (*this)[j][i];
  }
};

#endif // _LINK_STRENGTH_MATRIX_H_
