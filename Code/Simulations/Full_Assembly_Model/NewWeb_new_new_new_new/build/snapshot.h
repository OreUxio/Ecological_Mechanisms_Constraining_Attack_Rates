// -*- mode: c++ -*-
// $Id: snapshot.h 2318 2013-07-15 00:53:10Z tak $
#ifndef _SNAPSHOT_H_
#define _SNAPSHOT_H_

#include "sequence.h"
#include "link_strength_matrix.h"
#include "NetworkAnalysis.h"

// generic description of a food-web state:

class Snapshot 
{
  int S;
  sequence<int> the_bottom;
  sequence<double> the_biomass_B;
  sequence<double> the_bodymass_M;
  sequence<double> the_level;
  sequence<double> the_level0;
  sequence<double> the_level1;
  sequence<double> the_cStar;
  link_strength_matrix the_flow;
  link_strength_matrix the_ifrac;
  Interaction_Matrix the_link;
  double the_area_per_compartment;
  void compute_level_if_required();
  void compute_ifrac_if_required();
  void compute_level0_if_required();
  void compute_level1_if_required();
  void compute_cStar_if_required();
  bool is_plant(int i);
  permutation the_rectification;
 public:
  Snapshot();
  void set_number_of_compartments(int s);
  void set_biomasses(const sequence<double> & B);
  void set_bodymasses(const sequence<double> & M);
  sequence< double > get_biomasses();
  sequence< double > get_bodymasses();
  void set_flows(const link_strength_matrix & flow);
  void set_links(const Interaction_Matrix & link);
  void set_bottom(const sequence<int> & bottom);
  void set_area(double area);
  void adjust_links_given_flows(double threshold);
  void adjust_bottom_given_links();
  const link_strength_matrix & ifrac();
  link_strength_matrix get_ifrac();
  Interaction_Matrix get_im();
  int number_of_species();
  int number_of_plants();
  int number_of_animals();
  int number_of_species_Mthreshold(double Mlowerthreshold, double Mupperthreshold);
  int number_of_species_TLthreshold(double TLlowerthreshold, double TLupperthreshold);
  const link_strength_matrix & get_intake_matrix();  //trash? 
  link_strength_matrix intake_matrix();              //trash?
  void rank_abundance_plot(const char * filename,bool with_plants);
  void size_spectrum(const char * filename,
		     bool with_plants, bool with_animals);
  void biomass_spectrum(const char * filename,
			bool with_plants, bool with_animals);
  void trophic_level_structure(const char * filename);
  void pajek_graph(const std::string & filename, double threshold);
  void species_table(const char * filename,bool with_plants,double threshold,
		     const sequence<std::string> & additional_headings=sequence<std::string>(),
		     const sequence<std::string> & additional_data=sequence<std::string>());
  void species_table2(const char * filename,bool with_plants,double threshold,
		     const sequence<std::string> & additional_headings=sequence<std::string>(),
		     const sequence<std::string> & additional_data=sequence<std::string>());
  void species_table3(const char * filename,bool with_plants,double threshold,
		     const sequence<std::string> & additional_headings=sequence<std::string>(),
		     const sequence<std::string> & additional_data=sequence<std::string>());
  void link_table(const char * filename, double threshold);
  void link_table2(const char * filename, double threshold);
  void fast_stats();
  permutation rectification();
  double Omega(permutation p,double level_threshold=0);
private:
  double Delta_Omega(const permutation & p,int i,int j); //change permuting i, j
public:
  double Living_Planet_Index(Snapshot & baseline,
			     bool with_plants, bool with_animals);
  double Shanon_diversity_H(bool with_plants, bool with_animals);
  double Shanon_diversity_H_Mthreshold(double Mlowerthreshold, double Mupperthreshold);
  double Shanon_diversity_H_TLthreshold(double TLlowerthreshold, double TLupperthreshold);
  double Simpson_s_diversity_D(bool with_plants, bool with_animals);
  double Simpson_s_diversity_D_Mthreshold(double Mlowerthreshold, double Mupperthreshold);
  double Simpson_s_diversity_D_TLthreshold(double TLlowerthreshold, double TLupperthreshold);
  double MTI_Mthreshold(double Mlowerthreshold, double Mupperthreshold);
  double MTI_TLthreshold(double TLlowerthreshold, double TLupperthreshold);
  double max_level();
  double biomass_weighted_mean_level();
  double mean_level();
};

#endif // _SNAPSHOT_H_
