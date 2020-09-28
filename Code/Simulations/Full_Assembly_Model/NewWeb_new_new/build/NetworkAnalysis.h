//$Id: NetworkAnalysis.h 2466 2016-05-01 23:27:44Z axel $

#ifndef __NETWORK_ANALYSIS__
#define __NETWORK_ANALYSIS__

#define NEW_POPULATION_DYNAMICS
#define NO_POPULATION_DYNAMICS

#include "error.h"
#include "random.h"

#include "CMatrix.h"
#include <iostream>
#include <iomanip>
#include <iosfwd>
#include <vector>
#if !defined(NO_POPULATION_DYNAMICS) && !defined(NEW_POPULATION_DYNAMICS)
#include "lists.h"
#elif !defined(NO_POPULATION_DYNAMICS) && defined(NEW_POPULATION_DYNAMICS)
#include "NewSpecies.h"
#include <list>
typedef std::list<NewSpecies *> Species_List ;
#endif
#include "sequence.h"
#include "NewMatrix.h"

// This should perhaps be a class member rather than a namespace?
namespace NetworkAnalysis
{
  /// Kinds of interactions. Currently only eats and none are used.
  typedef enum {none=0,eats,shadows,hampers,other} Interaction;
}

/// Base type for Interaction_Matrix
typedef CMatrix<NetworkAnalysis::Interaction> IMatrix;

/// Binary interaction matrix with many fancy operations defined.
class Interaction_Matrix : private IMatrix  {
  int the_size;
 public:
  const IMatrix & CMatrix(){return *this;}
#if !defined(NO_POPULATION_DYNAMICS)
  std::simple_vector< Species_List > the_species;
#else
  std::simple_vector<bool> the_species; //dummy
#endif
  void label_species_by_index();
  typedef NetworkAnalysis::Interaction entry;
  //Constructor & Copy Constructor
  Interaction_Matrix(int size=0):
    IMatrix(size,size),the_size(size),the_species(size){};
  Interaction_Matrix(const IMatrix& matrix):
    IMatrix(matrix),the_size(matrix.GetYSize()),the_species(the_size)
    {
      ASSERT(matrix.GetXSize()==matrix.GetYSize());
      the_size=matrix.GetYSize();
    };
  Interaction_Matrix(const Interaction_Matrix& matrix);
  Interaction_Matrix const & operator= (Interaction_Matrix const& matrix);
#if 0 //defined(DEBUGGING) & !defined(PARALLEL)
  const Container2DRow<NetworkAnalysis::Interaction> & operator [] (int i) 
    const
  {
    return IMatrix::operator[](i);
  };
  accountingContainer2DRow & operator [] (int i)  {
    return IMatrix::operator[](i);
  };
#else
  const NetworkAnalysis::Interaction * operator [] (int i) 
    const
    {
      return IMatrix::operator[](i);
    };
  NetworkAnalysis::Interaction * operator [] (int i)  
    {
      return IMatrix::operator[](i);
    };
#endif //DEBUGGING

public:
  // helper functions:
  Interaction_Matrix select(std::simple_vector<bool> &sel) const;
  Interaction_Matrix permute(permutation new_pos) const;
  int size() const {return the_size;};

public:
  // here comes the actual analysis stuff:
  bool eats(int i, int j) const{
    return (int)(*this)[i][j]==(int) NetworkAnalysis::eats;
  }
  bool connected(int i, int j){
    return eats(i,j) || eats(j,i);
  }
  int Number_of_Species_S() const;
  int Number_of_Links_L();
  void Print(std::ostream &os = std::cout) const;
  void PPrint(std::ostream &os = std::cout) const; 
  void tPrint(std::ostream &os = std::cout) const;
  void pgm_write(const char * filename);
  double connectance_C();
  double links_per_species_Z();
  Interaction_Matrix msort(); //sort by mass
  Interaction_Matrix tsort(); // topologial sort
  Interaction_Matrix tsort2();
  Interaction_Matrix random_shuffle();
  Interaction_Matrix trophic();
  Interaction_Matrix randomize();
  Interaction_Matrix remove_lowest_level();
  Interaction_Matrix lump_lowest_level();
  Interaction_Matrix largest_connected_subweb() const;
  void cumulative_prey_hist(std::ostream & co,
			    std::string prefix="");
  void cumulative_predator_hist(std::ostream & co,
				std::string prefix="");
  void prey_hist(std::ostream & co,
			    std::string prefix="");
  void predator_hist(std::ostream & co,
				std::string prefix="");
  //typedef sequence<long long int> histogram;
  typedef sequence<double> histogram;
  typedef sequence<double> distribution;
  histogram cumulative_prey_hist();
  histogram cumulative_predator_hist();
  histogram prey_hist();
  histogram predator_hist();
  distribution theoretical_cumulative_prey_dist(int n, double Z=-1);
  distribution theoretical_cumulative_predator_dist(int n, double Z=-1);
  distribution theoretical_prey_dist(int n, double Z=-1);
  distribution theoretical_predator_dist(int n, double Z=-1);
  bool connected(); // web is connected?
  double prop_T();//fraction of top predators
  double prop_I();//fraction of intermediate species
  double prop_B();//fraction of bottom species
  double prop_GenSD();//std of prey count normalized to (L/S)
  double prop_VulSD();//std of predator count normalized to (L/S)
 private:
  double prop_T(histogram & h);//fraction of top predators
  double prop_I(histogram & h1, histogram & h2);
  //fraction of intermediate species
  double prop_B(histogram & h);//fraction of bottom species
  double prop_GenSD(histogram & h);//std of prey count normalized to (L/S)
  double prop_VulSD(histogram & h);//std of predator count normalized to (L/S)
 public:
  /* class similarity_matrix_t :  */
  /* public ublas::symmetric_matrix<double>{ */
  /* public: */
  /* similarity_matrix_t(int l) :  */
  /*   symmetric_matrix(NewZeroMatrix(l,l)){}; */
  /* similarity_matrix_t() : symmetric_matrix(0){} */
  /* }; */
  typedef NewMatrix similarity_matrix_t; // must enforce symmetry
  similarity_matrix_t
    similarity_s();//fraction of prey and predators shared by a pair
  double prop_MxSim();//web average over max s
  typedef sequence<int> food_chain_t;
  histogram chain_hist(); 
  //all loop-free directed chains starting from bottom
  double prop_ChnLg(histogram & h); //average chain length
  double prop_ChnSD(histogram & h); //std chain length
  double prop_ChnNo(histogram & h); //log # of chains
  double prop_Cannib(); //cannibal fraction
  double prop_Loop(); //fraction of species involved in loops
  double prop_Omniv(); // fraction of species that have food chains of
		       // different length
  double prop_Clust(); // Clustering Coefficient
  double prop_Ddiet(bool FrenchVariant=false);
  int prop_Cy4() const;
  double prop_Nest(); // degree of nestedness
  enum {pS,pC,pZ,pT,pI,pB,pGenSD,pVulSD,pMxSim,
	pChnLg, pChnSD,pChnNo,pLoop,pCannib,pOmniv,
	poChnLg, poChnSD,poChnNo,poLoop,poOmniv,
	pDdiet,
	pfDdiet,
	pCy4,
	pNest,
	pClust,
	pSStab,
	pend} property_t;
  typedef sequence<double> prop_vec_t;
  prop_vec_t props(); // all properties.
  static sequence<const char *> prop_names(); // the names of these properties.
  void two_column_write(const std::string & filename,std::string deliminter=",");
  distribution trophic_height_vector();
  histogram shortest_path_level_vector();
  void dot_graph(const std::string & filename);
  double structural_stability();
  bool has_consecutive_ones();
  bool is_interval();
  bool is_chordal();
};

bool has_consecutive_ones(const CMatrix<NetworkAnalysis::Interaction> & m);

std::istream & operator>>(std::istream &is, Interaction_Matrix & im);

#include "NetworkHelpers.h" //only to force linking

#endif //__NETWORK_ANALYSIS__


