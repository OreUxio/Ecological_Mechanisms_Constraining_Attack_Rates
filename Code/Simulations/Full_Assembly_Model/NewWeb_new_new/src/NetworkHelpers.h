// -*- c++ -*-
// $Id: NetworkHelpers.h 1431 2009-05-04 12:22:40Z axel $

#ifndef _NETWORKHELPERS_H_
#define _NETWORKHELPERS_H_

#include "NetworkAnalysis.h"

Interaction_Matrix forgiving_tsort(Interaction_Matrix &	m);
Interaction_Matrix the_largest_connected_subweb(const Interaction_Matrix & m);
Interaction_Matrix::histogram loop_free_chain_hist(Interaction_Matrix & m);
double loop_species_fraction(Interaction_Matrix & m);
double omnivore_fraction(Interaction_Matrix & m);
bool is_connected(Interaction_Matrix & m);
class graphs_stats_data {
public:
  graphs_stats_data(){};
  double oChnLg;
  double oChnSD;
  double oChnNo;
  double oLoop;
  double oOmniv;
  Interaction_Matrix sorted;
};
graphs_stats_data graph_analyze(Interaction_Matrix & m, 
				const bool go_through_multiplicities=true);
void find_lowest_level(Interaction_Matrix & m,std::simple_vector<bool> & lowest);
Interaction_Matrix::distribution trophic_height_vector(const Interaction_Matrix & m);

int Cy4(const Interaction_Matrix & m);

int niche_overlap_graph_analysis(const Interaction_Matrix & m);

#endif // _NETWORKHELPERS_H_
