// -*- mode: c++ -*-
// $Id: matrix_transformers.h 1790 2010-05-11 17:22:36Z tak $
#ifndef _MATRIX_TRANSFORMERS_H_
#define _MATRIX_TRANSFORMERS_H_

#include "link_strength_matrix.h"
#include "NetworkAnalysis.h"

/// Threshold \a l at \a th and make it an Interaction_Matrix.
Interaction_Matrix 
threshold(const link_strength_matrix & l,double th);

/// Computes diet fractions from a link_strength_matrix.
link_strength_matrix 
in_fraction(link_strength_matrix fx);

/// Compute DPF from a diet fraction matrix, and write it to \a filename.
void
strength_distribution(const link_strength_matrix & l,double th,
		      const char * filename);
/// Compute DPF from a diet fraction matrix, and write it to \a filename.
/** Computes also Zc(1%) and an estimate of the diet-partitioning
    exponent nu.*/
void
strength_distribution(const link_strength_matrix & l,double th,
		      const char * filename,double &Zc01,double &nu);

// Following two functions are the same as the previous two except that
// they only consider species above a certain body mass, which is why
// there is an extra argument for the list of body masses.
void
strength_distribution_new(const link_strength_matrix & l,double th,
		      const char * filename,const sequence<double> &M);

void
strength_distribution_new(const link_strength_matrix & l,double th,
	const char * filename,const sequence<double> &M,
	double &Zc01,double &nu);

#endif // _MATRIX_TRANSFORMERS_H_
