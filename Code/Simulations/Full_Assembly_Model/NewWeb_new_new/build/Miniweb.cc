// $Id: MoranI.cc 1155 2008-10-09 16:02:35Z axel $
#include <string>
const std::string version("$Revision: 1155 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/MoranI.cc,v $ $Date: 2008-10-09 17:02:35 +0100 (Thu, 09 Oct 2008) $");

#define LOGNORMALS
#define SORTED_VECTORS
#define USE_INVASION_FITNESS

const bool multitrophic_web=false;

#include <cmath>
#include <cstring>
#include <fstream>
#include <boost/range.hpp>

#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>

#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/DiagMatrix.h>

#include "sequence.h"
#include "Statistics.h"
#include "relax.h"
#include "random.h"
#include "snapshot.h"
#include "matrix_transformers.h"
#include "linpack_eigen.h"
#ifdef SORTED_VECTORS
#include "SortedVector.h"
#endif

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
#include "linpack_eigen.h"
#include "period_cutter.h"
#include "xy_graph.h"
#include "polyfit.h"
#include "Integrator.h"
#endif

using namespace std;

int steps_transient=0;
int requested_number_of_inactive_animals=10;
int requested_number_of_inactive_plants=10;
int min_number_of_plants=10;
double invasion_biomass=1e-5;
double log10_extinction_biomass=-100;

// adjustable parameters:
static double mu,sigma,epsilon,r,T,tau,K,C,a0,animal_fraction;
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGDOUBLE(mu),
    CFGDOUBLE(sigma),
    CFGDOUBLE(epsilon),
    CFGDOUBLE(r),
    CFGDOUBLE(T),
    CFGDOUBLE(tau),
    CFGDOUBLE(K),
    CFGDOUBLE(C),
    CFGDOUBLE(a0),
    CFGDOUBLE(animal_fraction),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);

bool has_suffix(const char * suff,const char * filename){
  return 0==strcmp(suff,filename+strlen(filename)-strlen(suff));
}

namespace io=boost::iostreams ;

struct opipestream : io::stream< io::file_descriptor_sink >
{
  typedef io::stream< io::file_descriptor_sink > base ;
  explicit opipestream( const char* command )
  : base( fileno( pipe = popen( command, "w" ) ) ) {}
  ~opipestream() { close() ; pclose( pipe ) ; };
private : 
  FILE* pipe ;
};

struct ipipestream : io::stream< io::file_descriptor_source >
{
  typedef io::stream< io::file_descriptor_source > base ;
  explicit ipipestream( const char* command )
  : base( fileno( pipe = popen( command, "r" ) ) ) {}
  ~ipipestream() { close() ; pclose( pipe ) ; };
private : 
  FILE* pipe ;
};

ostream & operator<<(ostream & os,CLHEP::HepGenMatrix & m){
  os << endl;
  for(int i=0;i<m.num_col();i++){
    for(int j=0;j<m.num_row();j++){
      os << m[i][j] << "\t";
    }
    os << endl;
  }
}

// ostream & operator<<(ostream & os,sequence<sor){
//   os << endl;
//   for(int i=0;i<m.num_col();i++){
//     for(int j=0;j<m.num_row();j++){
//       os << m[i][j] << "\t";
//     }
//     os << endl;
//   }
// }

void matrix_pgm(CLHEP::HepMatrix & m, const char * filename, 
		double mi, double ma){
  ofstream pic(filename);
  const int max_gray=256;
  pic << "P5" << endl;
  pic << m.num_col() << endl;
  pic << m.num_row() << endl;
  pic << max_gray-1 << endl;
  for(int i=0;i<m.num_row();i++){
    for(int j=0;j<m.num_col();j++){
      double val=m[i][j];
      val=1-(val-mi)/(ma-mi);
      if(val >= 1) val=0.9999999999;
      if(val < 0) val=0;
      pic << (unsigned char)(val*max_gray);
    }
  }
}
  

class miniweb_t : public relaxing_dynamical_object{
  template<typename ARRAY>
  inline void delete_replacing(int i,ARRAY & a){
    int last=a.size()-1;
    a[i]=a[last];
    a.resize(last);
  }

  template<typename ARRAY>
  inline void roll_up(int i,ARRAY & a){
    int last=a.size();
    a.resize(last+1);
    a[last]=a[i];
  }

  template<typename ARRAY>
  inline void double_roll_down(int i,ARRAY & a){
    int new_size=a.size()-1;
    a[i]=a[n_animals-1];
    a[n_animals-1]=a[new_size];
    a.resize(new_size);
  }

private:
  double _mu,_sigma,_epsilon,_animal_respiration,_T,_tau,_K,_C,_a0,
    _animal_fraction;
  int _max_print_length;
  sequence<double> available;
  sequence<double> fac;
#ifdef SORTED_VECTORS
  vector< SortedVector > a;
  vector< SortedVector > aT;
#else
  vector< vector<double> > a;
#endif
  double plant_sum;
  typedef vector<double> species_container_t;
  typedef boost::iterator_range<species_container_t::iterator> species_range;
  species_container_t species;
  vector<bool> active;
  int n_animals,n_inactive_animals,n_inactive_plants;
public:
  int step_count;

  // SERIALIZATION:
private:
  friend class boost::serialization::access;
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    ar 
      & _mu & _sigma & _epsilon & _animal_respiration
      & _T & _tau & _K & _C & _a0 & _animal_fraction 
      & _max_print_length
      & a
      & species
      & active
      & n_animals & n_inactive_animals & n_inactive_plants
      & step_count;
  }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    ar 
      & _mu & _sigma & _epsilon & _animal_respiration
      & _T & _tau & _K & _C & _a0 & _animal_fraction 
      & _max_print_length
      & a
      & species
      & active
      & n_animals & n_inactive_animals & n_inactive_plants
      & step_count;

    update_ranges();
#ifdef SORTED_VECTORS
    aT.resize(n_animals);
    for(int i=species.size();i-->0;){
      for(int j=n_animals;j-->0;){
	aT[j].resize(species.size());
	aT[j][i]=a[i][j];
      }
    }
#endif
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
  species_range plant,animal;
private:
  void update_ranges(){
    species_container_t::iterator boundary=
      species.begin()+n_animals;
    animal=species_range(species.begin(),boundary);
    plant=species_range(boundary,species.end());
  }

public:
  miniweb_t(){};
  miniweb_t(double mu, double sigma,double epsilon, double r,double T,
	    double tau, double K, double C, double a0, 
	    double animal_fraction):
    _mu(mu),
    _sigma(sigma),
    _epsilon(epsilon),
    _animal_respiration(r),
    _T(T),
    _tau(tau),
    _K(K),
    _C(C),
    _a0(a0),
    _animal_fraction(animal_fraction),
    _max_print_length(0),
    n_animals(0),
    n_inactive_animals(0),n_inactive_plants(0),
    species(0),
    active(0),
    animal(species),
    plant(species),
    step_count(0)
  {};
  ~miniweb_t(){};
  virtual void read_state_from(const ODE_vector & state){
    for(int i=species.size();i-->0;)
      species[i]=( active[i] ? exp(state[i]) : state[i] );
  }
  virtual void write_state_to(ODE_vector & state) const{
    for(int i=species.size();i-->0;)
      state[i]=( active[i] ? log(species[i]) : species[i] );
  }
  virtual int number_of_variables() const{
    return species.size();
  }
  int number_of_species(){
    return species.size()-n_inactive_animals-n_inactive_plants;
  }
  int number_of_animals(){
    return n_animals-n_inactive_animals;
  }
  int number_of_plants(){
    return species.size()-n_animals-n_inactive_plants;
  }
  void pre_dynamics(){
    available.resize(animal.size());
    fac.resize(animal.size());
    plant_sum=0;
    if(species.empty()){
      // there are no species...
      for(int i=animal.size();i-->0;){
	available(i)=0;
	fac(i)=1;
      }
    }else{
      // there are some species...
#ifdef SORTED_VECTORS
      double species_max=0;
      for(int j=species.size();j-->0;){
	if(species[j]>species_max)
	  species_max=species[j];
      }
#endif
      for(int i=animal.size();i-->0;){
#if defined(SORTED_VECTORS)
	available(i)=aT[i].dot(&(species[0]),species_max);
#else
	available(i)=0;
	for(int j=species.size();j-->0;){
	  available(i)+=species[j]*a[j][i];
	}
#endif
	fac(i)=1/(1+_T*available(i));
      }
      for(int j=plant.size();j-->0;){
	plant_sum+=plant[j];
      }
    }    
  }

  virtual int dynamics(ODE_vector const & state,
		       ODE_vector & time_derivative){

    read_state_from(state);
    for(int i=species.size();i-->0;){
      if(!active[i])
	species[i]=0;
    }

    pre_dynamics();

    for(int i=animal.size();i-->0;){
      
      time_derivative[i]=
	fac(i)*
	_epsilon*
	available(i)-_animal_respiration;
    }
    for(int j=plant.size();j-->0;){
      time_derivative[j+animal.size()]=
	(1-plant[j]/_K)/_tau;
    }
    if(!animal.empty()){
      // there are animals...
#ifdef SORTED_VECTORS
      double animal_max=0;
      for(int i=animal.size();i-->0;){
	if(animal[i]>animal_max)
	  animal_max=animal[i];
      }
      vector_with_max hungry_animal(animal.size());
      for(int i=animal.size();i-->0;){
	hungry_animal(i)=
	  fac(i)*animal[i];
      }
#endif
      for(int j=species.size();j-->(multitrophic_web?0:n_animals);){
#if defined(SORTED_VECTORS)
	double being_eaten;
	being_eaten=a[j].dot(hungry_animal);
#else
	double being_eaten=0;
	for(int i=animal.size();i-->0;){
	  being_eaten+=a[j][i]*animal[i]
	    *fac(i);
	}
#endif
	time_derivative[j]-=being_eaten;
      }
    }
    return 0;
  }
  double active_sum(ODE_vector & s){
    double sum=0;
    for(int i=s.size();i-->0;){
      if(active[i]){
	sum+=s[i];
      }
    }
    return sum;
  }
  bool is_active(int i){
    return active[i];
  }
private:
  double start_recording_time;
  double last_recording_time;
public:
  vector<double> abundance_integral;
  vector<double> effective_abundance_integral;

public:
  void add_inactive_species_if_required(){
#ifdef USE_INVASION_FITNESS
    while(n_inactive_animals<requested_number_of_inactive_animals){
      add_animal();
    }
    while(n_inactive_plants<requested_number_of_inactive_plants){
      add_plant();
    }
#endif
  }

private:
  void prepare_for_integration(){
    abundance_integral.resize(species.size());
    effective_abundance_integral.resize(n_animals);
    last_recording_time=DBL_MAX;
    
    for(int i=species.size();i-->0;){
      if(!active[i]){
	species[i]=0;
      }
    }
  }
  void record_for_steady_state(){
    if(current_time < last_recording_time){
      last_recording_time=current_time;
      start_recording_time=current_time;

      abundance_integral.resize(species.size());
      effective_abundance_integral.resize(n_animals);

      for(int i=abundance_integral.size();i-->0;){
	abundance_integral[i]=0;
      }
      for(int i=effective_abundance_integral.size();i-->0;){
	effective_abundance_integral[i]=0;
      }
    }else{
      double delta_t=current_time-last_recording_time;
      for(int i=species.size();i-->n_animals;){
	abundance_integral[i]+=species[i]*delta_t;
      }
      for(int i=n_animals;i-->0;){
	abundance_integral[i]+=species[i]*delta_t;
	effective_abundance_integral[i]+=species[i]*fac(i)*delta_t;
      }
      last_recording_time=current_time;
    }
  }

  double random_link_strength(){
#ifdef LOGNORMALS
    double normal=gaussian(_mu,_sigma);
    if(normal > -log(DBL_MAX)/2){
      return exp(normal);
    }else{
      return 0;
    }
#else
    return ( unirand()<_C ? _a0 : 0 );
#endif
  }

public:
  int add_animal(){
    int i=n_animals;
    roll_up(i,species);
    species[i]=0;
    roll_up(i,active);
    active[i]=false;
    n_inactive_animals++;
    roll_up(i,a);
    for(int ii=i;ii-->0;){// Determine predators of i
      // No roll up in second index of a needed here, because this animal
      // is added in the last column.
      a[i][ii]=(multitrophic_web?random_link_strength():0);
    }
#ifdef SORTED_VECTORS
    aT.resize(i+1);
    aT[i].resize(species.size());
    for(int ii=i;ii-->0;){
      roll_up(i,aT[ii]);
      aT[ii][i]=a[i][ii];
    }
#endif
    n_animals++;
    update_ranges();
    for(int j=a.size();j-->0;){// Determine prey of i
      a[j].resize(i+1);
#ifdef SORTED_VECTORS
      aT[i][j]=
#endif
	a[j][i]=(j<n_animals && !multitrophic_web ?
		 0 : random_link_strength() );
    }
#if DEBUGGING
    WARNING("added animal " << i);
#endif
    return i;
  }
  int add_plant(){
    int j=species.size();
    species.resize(j+1);
    species[j]=0;
    update_ranges();
    active.resize(j+1);
    active[j]=false;
    n_inactive_plants++;
    a.resize(j+1);
    a[j].resize(n_animals);
    
    for(int i=n_animals;i-->0;){// Determine grazers of j
#ifdef SORTED_VECTORS
      aT[i].resize(j+1);
      aT[i][j]=
#endif
	a[j][i]=random_link_strength();
    }
#if DEBUGGING
    WARNING("added plant " << j);
#endif
    return j;
  }
  void save_invasion_fitnesses(){
    static ofstream ra("ra.dat",ios_base::out | ios_base::app);
    static ofstream rp("rp.dat",ios_base::out | ios_base::app);
    bool do_animal=true,do_plant=true;
    for(int i=species.size();i-->0;){
      if(!active[i]){
	double fitness= species[i]/(current_time-start_recording_time);
	if(i<n_animals){
	  if(do_animal){
	    ra << fitness << endl;
	    delete_species(i);
	    do_animal=false;
	  }
	}else{
	  if(do_plant){
	    rp << fitness << endl;
	    delete_species(i);
	    do_plant=false;
	  }
	}
      }
    }
  }
  void activate_species(int i){
    species[i]=invasion_biomass;
    active[i]=true;
    if(i<n_animals){
      n_inactive_animals--;
      WARNING("activated animal " << i );
    }else{
      n_inactive_plants--;
      WARNING("activated plant " << i );
    }
  }

    
  bool next_species_is_animal(){
    if(number_of_plants() < min_number_of_plants){
      return false;
    }
    if(_animal_fraction){
      return unirand() < _animal_fraction;
    }else{
      return unirand() < (/*true ||*/ step_count < steps_transient ? 
			  0.5/*_animal_fraction*/ :
			  double(n_animals)/double(species.size()) );
    }
  }
  int add_species(){
    if(next_species_is_animal()){
      return add_animal();
    }else{
      return add_plant();
    }
  }
  void extinguish_species(int i){
    WARNING("species " << i << " went extinct");
    delete_species(i);
  }
  void delete_species(int i){
    if(i<animal.size()){
      // this was an animal:
      double_roll_down(i,species);
      if(!active[i]) n_inactive_animals--;
      double_roll_down(i,active);
      double_roll_down(i,a);
#ifdef SORTED_VECTORS
      delete_replacing(i,aT);
#endif
      for(int j=a.size();j-->0;){
	delete_replacing(i,a[j]);
      }
#ifdef SORTED_VECTORS
      for(int ii=aT.size();ii-->0;){
	double_roll_down(i,aT[ii]);
      }
#endif
      n_animals--;
      update_ranges();
#if DEBUGGING
      WARNING("deleted animal " << i);
#endif
    }else{
      // this was a plant
      delete_replacing(i,species);
      if(!active[i]) n_inactive_plants--;
      delete_replacing(i,active);
      delete_replacing(i,a);
#ifdef SORTED_VECTORS
      for(int ii=aT.size();ii-->0;){
	delete_replacing(i,aT[ii]);
      }
#endif
      update_ranges();
#if DEBUGGING
      WARNING("deleted plant " << i);
#endif
    }
  }
  void check_internal_consistency(){
    ALWAYS_ASSERT(species.size()==a.size());
    for(int j=species.size();j-->0;){
      ALWAYS_ASSERT(a[j].size()==n_animals);
      if(active[j])
	ALWAYS_ASSERT(species[j]>0);
    }
    int aacount=0;
    for(int i=n_animals;i-->0;){
      if(active[i]) aacount++;
    }
    ALWAYS_ASSERT(aacount==number_of_animals());
    int apcount=0;
    for(int i=species.size();i-->n_animals;){
      if(active[i]) apcount++;
    }
    ALWAYS_ASSERT(apcount==number_of_plants());

    if(!multitrophic_web){
      for(int i=n_animals;i-->0;){
	for(int j=n_animals;j-->0;){
	  ALWAYS_ASSERT(a[i][j]==0);
	}
      }
    }
    ALWAYS_ASSERT(animal.size()==n_animals);
    ALWAYS_ASSERT(plant.size()==(species.size()-n_animals));
#ifdef SORTED_VECTORS
    ALWAYS_ASSERT(aT.size()==n_animals);
    for(int i=n_animals;i-->0;){
      ALWAYS_ASSERT(aT[i].size()==species.size());
      for(int j=species.size();j-->0;){
	if(a[j][i]!=aT[i][j]){
	  REPORT(a[j][i]);
	  REPORT(aT[i][j]);
	  REPORT(aT[i]);
	  ALWAYS_ASSERT(a[j][i]==aT[i][j]);
	}
      }
    } 
#endif
  }

  int add_fit_species(){
    bool add_an_animal_now=
      next_species_is_animal();
    
//     REPORT(n_animals);
//     if(number_of_species()>0){
//       ODE_vector s(number_of_variables());
//       write_state_to(s);
//       REPORT(s);
//     }
//     REPORT(a);
//     REPORT(aT);

    while(true){
#ifdef USE_INVASION_FITNESS
      while(n_inactive_animals && n_inactive_plants){
	if(add_an_animal_now){
	  int i=0; // pointer to next animal to test
	  while(active[i]) i++;
	  //REPORT(species[i]);
	  if(species[i]<0){
	    delete_species(i);
	  }else{
	    activate_species(i);
#ifdef DEBUGGING
	    check_internal_consistency();
#endif
	    add_inactive_species_if_required();
	    return i;
	  }
	}else{
	  int j=n_animals;
	  while(active[j]) j++;
	  //REPORT(species[j]);
	  if(species[j]<0){
	    delete_species(j);
	  }else{
	    activate_species(j);
#ifdef DEBUGGING
	    check_internal_consistency();
#endif
	    add_inactive_species_if_required();
	    return j;
	  }
	}
      }     
      // we ran out of inactive species
      WARNING("consider using more inactive species");
      add_inactive_species_if_required();
      relax(1000);
    }
    FATAL_ERROR("This point should never be reached");
    return 0;
#else  // not using invasion fitness
    int S=number_of_variables()+1;
    ODE_vector l(S),dl(S);
    while(true){
      int i=(add_an_animal_now ? 
	     add_animal() :
	     add_plant() );
      ALWAYS_ASSERT(!active[i]);
      active[i]=true;
      if(i<n_animals){
	n_inactive_animals--;
      }else{
	n_inactive_plants--;
      }
      species[i]=1e-5;

#ifdef DEBUGGING
      check_internal_consistency();
#endif

      write_state_to(l);
      dynamics(l,dl);
    
      if(dl[i]<0){
	delete_species(i);
#ifdef DEBUGGING
	check_internal_consistency();
#endif
      }else{
	return i;
      }
    }
#endif
  }
  bool small_values_in(ODE_vector & state,
		       const species_set_t&  conserved){
    double log_threshold=log10_extinction_biomass*M_LN10-1;
    for(int i=state.size();i-->0;){
      if(active[i] && state[i]<log_threshold) return true;
    }
    return false;
  }

  virtual species_set_t 
  delete_species_larger_than_exp(const sequence<double> & si,
				 const species_set_t&  conserved){
    species_set_t deleted;
    double log_threshold=log10_extinction_biomass*M_LN10;
    for(int i=number_of_variables();i-->0;){
      if(active[i] && si(i)<log_threshold){
	extinguish_species(i);
	cout << "P1" << endl;
#ifdef DEBUGGING
	check_internal_consistency();
#endif
	deleted.insert(i);
      }
    }
    return deleted;
  }

  virtual species_set_t 
  delete_all_species_with_less_than_one_individual(const species_set_t
						   & conserved){
    species_set_t deleted;
    for(int i=number_of_variables();i-->0;){
      if(active[i] && species[i] <= pow(10,log10_extinction_biomass)){
	extinguish_species(i);
	cout << "P2" << endl;
#ifdef DEBUGGING
	check_internal_consistency();
#endif
	deleted.insert(i);
      }
    }
    return deleted;
  }    

  void delete_all_inactive_species(){
    for(int i=number_of_variables();i-->0;){
      if(!active[i]){
	delete_species(i);
      }
    }
  }

  void set_abundances(vector<double> & n){
    REPORT(n.size());
    REPORT(number_of_species());
    REPORT(number_of_variables());
    ALWAYS_ASSERT(n.size()==number_of_variables());
    copy(n.begin(),n.end(),species.begin());
  }

  void scale_abundance(int n, double f){
    species[n]*=f;
  }

  void evaluate(){
    //return;

    int n_inactive_species=n_inactive_animals+n_inactive_plants;

    Snapshot web_data;
    link_strength_matrix ls;

    for(int i=abundance_integral.size();i-->0;){
      abundance_integral[i]/=current_time-start_recording_time;
    }
    for(int i=n_animals;i-->0;){
      effective_abundance_integral[i]/=current_time-start_recording_time;
    }
    
    ls.resize(number_of_variables());
    ls=ls*0;
    for(int i=n_animals;i-->0;){
      if(active[i]){
	for(int j=species.size();j-->0;){
	  if(active[j]){
	    // Observe: This is still an estimate of average flows, but
	    // perhaps a rather good one.
	    ls[i][j]=a[j][i]*abundance_integral[j]*
	      effective_abundance_integral[i];
	  }
	}
      }
    }

//     web_data.set_number_of_compartments(x.size());
//     web_data.set_biomasses(x);
//     web_data.set_bodymasses(x*0+1);
//     web_data.set_flows(ls);
//     web_data.set_area(1);
//     web_data.adjust_links_given_flows(0.01);

    double nu_pred=
      sqrt(2*log(plant.size()))/_sigma;
    double Zc01_pred=
      sin(M_PI*nu_pred)/(M_PI*nu_pred)*pow(0.01/0.99,-nu_pred);
    double magic_number=3*_animal_respiration/(_epsilon*_a0*_K);
    double Zc01,nu;
    strength_distribution(in_fraction(ls),//web_data.get_ifrac(),
			  1e-8,
			  "strength.dat",
			  Zc01,nu);
    REPORT(Zc01);
    REPORT(Zc01_pred);
    REPORT(nu);
    REPORT(nu_pred);
    average_meter animal_dist;
    for(int i=animal.size();i-->0;){
      if(active[i]){
	animal_dist.sample(log(abundance_integral[i]));
      }
    }
    double sigma_animal=sqrt(animal_dist.var());
    REPORT(sigma_animal);
    double sigma_animal_pred=
      _sigma/sqrt(1+2*log(plant.size()));
    REPORT(sigma_animal_pred);
    REPORT(animal_dist);
    static average_meter plant_sat,Z,avul,pvul,satZ,f,avaT,favaT,S;
    if(step_count>steps_transient){
      for(int i=animal.size();i-->0;){
	if(active[i]){
	  avaT.sample(_T*available(i));
	  f.sample(fac(i));
	  favaT.sample(fac(i)*_T*available(i));
	}
      }
      Z.sample(Zc01);
      S.sample(number_of_variables());
      for(int j=plant.size();j-->0;){
	if(active[j]+n_animals){
	  plant_sat.sample(plant[j]/_K);
	}
      }
      for(int j=species.size();j-->0;){
	if(active[j]){
	  int v=0;
	  for(int i=animal.size();i-->0;){
	    if(active[i] && a[j][i])
	      v++;
	  }
	  if(j<n_animals){
	    avul.sample(v);
	  }else{
	    pvul.sample(v);
	    satZ.sample(v*species[j]/_K);
	  }
	}
      }
    }
    REPORT(plant_sat);
    REPORT(Z);
    REPORT(avul);
    REPORT(pvul);
    REPORT(satZ);
    REPORT(f);
    REPORT(avaT);
    REPORT(favaT);
    REPORT(S);

    {
#ifdef DEBUGGING
      check_internal_consistency();
#endif
      ofstream a("animals.dat"),p("plants.dat");
#if 1
      vector<double> sa(number_of_animals());
      vector<double> sp(number_of_plants());
      vector<double>::iterator j;
      j=sa.begin();
      for(int i=n_animals;i-->0;){
	if(active[i]){
	  *j = abundance_integral[i];
	  j++;
	}
      }
      ALWAYS_ASSERT(j==sa.end());
      
      j=sp.begin();
      for(int i=species.size();i-->n_animals;){
	if(active[i]){
	  *j = abundance_integral[i];
	  j++;
	}
      }
      ALWAYS_ASSERT(j==sp.end());
#else
      vector<double> sa(abundance_integral.begin(),
			abundance_integral.begin()+n_animals);
      vector<double> sp(abundance_integral.begin()+n_animals,
			abundance_integral.end());
#endif
      sort(sa.begin(),sa.end());
      sort(sp.begin(),sp.end());
      for(int i=0;i<sa.size();i++){
	if(sa[i]!=0){
	  a << sa[i] << endl;
	}
      }
      for(int j=0;j<sp.size();j++)
	if(sp[j]!=0){
	  p << sp[j] << endl;
	}
    }

  }
  void line_print(ODE_vector const & state, ostream &co){
    double log10=log(10);
    if(_max_print_length<number_of_species()){
      _max_print_length=number_of_species();
    }
    for(int i=0;i<species.size();i++){
      if(active[i]){
	co << state[i]/log(10.0) << " ";
      }
    }
    for(int i=_max_print_length-number_of_species();i-->0;){
      co << log10_extinction_biomass << " ";
    }
  }
  void competition_report(){
    int n=number_of_variables();
    CLHEP::HepMatrix A(n,n_animals);
    CLHEP::HepMatrix feed(n_animals,n);
    CLHEP::HepDiagMatrix norm(n_animals);
    pre_dynamics();

    for(int c=n_animals;c-->0;){
      double ss=0;
      for(int r=n;r-->0;){
	double f=species[r]*a[r][c]*fac(c);
	ss+=f*f;
	feed[c][r]=f;
	A[r][c]=a[r][c];
      }
      norm[c][c]=1/sqrt(ss);
    }
    if(!multitrophic_web){
      A=A.sub(n_animals+1,n,1,n_animals);
    }
    CLHEP::HepMatrix comp=norm*(feed*feed.T())*norm;
    matrix_pgm(comp,"comp.pgm",0,1);
    eigen_report(comp,"compev.dat");

    average_meter comp_av;
    for(int i=n_animals;i-->0;){
      for(int j=i;j-->0;){
	comp_av.sample(comp[i][j]);
      }
    }
    REPORT(animal.size());
    REPORT(plant.size());
    REPORT(comp_av);
    REPORT(comp_av.var());
    REPORT(number_of_animals()*comp_av);
    REPORT(number_of_animals()*comp_av.var());
    REPORT(comp_av.var()/(comp_av));
//     REPORT(comp);
//     REPORT(norm);

    int error;
    CLHEP::HepMatrix D=A.T()*A;
    eigen_report(D,"Wev.dat");
    CLHEP::HepDiagMatrix H(n_animals);
    for(int i=n_animals;i-->0;){
      ASSERT(D[i][i]>0);
      H[i][i]=1/sqrt(D[i][i]);
    }
    CLHEP::HepMatrix B=H*D*H;
    eigen_report(B,"Bev.dat");
    matrix_pgm(B,"B.pgm",0,1);
    
    average_meter B_av;
    for(int i=n_animals;i-->0;){
      for(int j=i-1;j-->0;){
	B_av.sample(B[i][j]);
      }
    }
    REPORT(B_av);
    REPORT(B_av.var()*plant.size());
    REPORT(B_av.var()*B.num_col());

    D=D.inverse(error);
    ASSERT(!error);
    D=A*D*A.T();
    matrix_pgm(D,"D.pgm",-0.03,1);    
    //eigen_report(D,"Dev.dat"); // Boring: 1 for each animal, 0 elsexv

    if(false){// some neutral theory
      int m=1000;
      n=int(2*m);
      double sigma=2*sqrt(2*log(double(n)));
      CLHEP::HepMatrix C(n,m);
      for(int i=n;i-->0;){
	for(int j=m;j-->0;){
	  C[i][j]=exp(gaussian(0,sigma));
	}
      }
      CLHEP::HepMatrix CC=C.T()*C;
      CLHEP::HepDiagMatrix DD(m);
      for(int j=m;j-->0;){
	DD[j][j]=1/sqrt(CC[j][j]);
      }
      CC=DD*CC*DD;
      average_meter CCav;
      for(int i=m;i-->0;){
	for(int j=i;j-->0;){
	  CCav.sample(CC[i][j]);
	}
      }
      REPORT(CCav.var()*n);
      matrix_pgm(CC,"CC.pgm",0,1);    
      eigen_report(CC,"CCev.dat");
    }
  }
  void eigenvalue_report(){
    int n=number_of_variables();
    CLHEP::HepMatrix J(n,n);
    numerical_Jacobian(J);
//     for(int i=n;i-->0;){
//       for(int j=n;j-->0;){
// 	J[i][j]*=species[i]/species[j];
//       }
//     }
    eigen_report(J,"eigenvalues.dat");
    const CLHEP::HepMatrix J_pa=J.sub(1+n_animals,n,1,n_animals);
    const CLHEP::HepMatrix J_ap=J.sub(1,n_animals,1+n_animals,n);
    const CLHEP::HepMatrix J_pp=J.sub(1+n_animals,n,1+n_animals,n);
    const CLHEP::HepMatrix J_aa=J.sub(1,n_animals,1,n_animals);
    REPORT(J_pa);
    for(int i=n_animals;i<n;i++){
      cout << species[i] << "\t";
    }
    cout << endl;
    REPORT(J_ap);
    REPORT(J_pp);
    REPORT(J_aa);
    int ierr;
    const CLHEP::HepMatrix C=
      -J_aa/*==0 ? */+J_ap*J_pp.inverse(ierr)*J_pa;
    if(ierr) WARNING("matrix inversion failed");
    eigen_report(C,"C1.dat");
    CLHEP::HepDiagMatrix D(n_animals);
    for(int i=n_animals;i-->0;){
      D[i][i]=1.0/animal[i];
    }
    CLHEP::HepMatrix X;
    X=C*D; // observe: J is Jacobian for logarithmic abundances!
//     for(int i=0;i<n_animals;i++){
//       cout << animal[i] << ":\t";
//       for(int j=0;j<n_animals;j++){
// 	cout << X[i][j] << "\t";
//       }
//       cout << endl;
//     }
    eigen_report(X,"C2.dat");
    bool failure=false;
    for(int i=n_animals;i-->0;){
      for(int j=n_animals;j-->0;){
	X[i][j]/= /*(X[i][i]>0?+1:-1) * */ sqrt(fabs(X[i][i]*X[j][j]));
	if(!finite(X[i][j])) failure=true;
      }
    }
    //X=D*C; // observe: J is Jacobian for logarithmic abundances!
    if(failure) return;
    cout << "****************************" << endl;
//     for(int i=0;i<n_animals;i++){
//       cout << animal[i] << ":\t";
//       for(int j=0;j<n_animals;j++){
// 	cout << X[i][j] << "\t";
//       }
//       cout << endl;
//     }
    eigen_report(X,"C3.dat");
  }
  double time_scale(){
    return _animal_respiration;
  }
private:
  void eigen_report(const CLHEP::HepMatrix & J,const char * filename){
    std::ofstream os(filename);
    
    ALWAYS_ASSERT(J.num_row()==J.num_col());
    CLHEP::HepVector ddr(J.num_row()),ddi(J.num_row());
    NewEigen(J,ddr,ddi);
    
    int npos=0;
    int old_prec=std::cout.precision();
    os.precision(4);  
    for(int i=0;i<J.num_row();i++){
      os << ddr[i] << " " << ddi[i] << std::endl;
      if(ddr[i]>0)
	npos++;
    }
    std::cout.precision(old_prec);

    if(npos){
//       WARNING(npos <<  (npos==1 ? 
// 			" eigenvalue has positive real part" : 
// 			" eigenvalues have positive real part" ));
    }
    return;
  }
};

void read(miniweb_t & web,int argc,char ** argv){
  std::string input_file_name;
  if(argc>2){
    input_file_name=argv[2];
  }else{
    input_file_name="Miniweb.dat";
  }
  
  if(has_suffix(".bz2", input_file_name.c_str())){
#if 1 // Alternative but equivalent formulations
    std::stringstream input_command(stringstream::out);
    input_command << "bunzip2 -c " << input_file_name;
    ipipestream in(input_command.str().c_str());
#else
    using namespace boost::iostreams;
    ifstream file(input_file_name.c_str(), ios_base::in | ios_base::binary);
    filtering_istream in;
    in.push(bzip2_decompressor());
    in.push(file);
#endif
    boost::archive::text_iarchive ia(in);
    ia >> web;
  }else{
    std::ifstream ifs(input_file_name.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> web;
  }
}


///////////////////////////////////////////////////////////////  
int main(int argc,char** argv){
  signal_handling();
  set_cfg_parameter("vector_accuracy","1e-5");
  set_cfg_parameter("vector_truncation_epsilon","1e-10");
  set_cfg_parameter("DEFAULT_RELATIVE_TOLERANCE","0.001");

  set_random_seed(getpid());

  average_meter mean_size;
  bool skip_averaging=true;
  int previous_size=0;
  
  double target_S=100;
  K=1e0;a0=1e0;epsilon=0.6;
  double feeding=0.6;
  const double plant_abunance_drop=0.5;

#ifdef LOGNORMALS
  double target_nu=0.5;
  REPORT(target_nu);
  sigma=sqrt(2*log(target_S))/target_nu;
  REPORT(sigma);
  mu=log(a0)-3.0/2.0*sigma*sigma;  //The ratio of E[X^2] to E[X]
				      //for this lognormal is exactly
				      //a0.  Should make sure that SS
				      //abundance of animals is about
				      //1/a0.
  REPORT(mu);
  double inflection_point=mu+sigma*sqrt(2*log(target_S));
  //double mean_sum=target_S*exp(mu+sigma*sigma/2);
  REPORT(exp(inflection_point)/a0);
  double effective_r=epsilon*K*exp(inflection_point)*(1-feeding)*
    plant_abunance_drop;
  C=1;//dummy;
  double target_Z=sin(M_PI*target_nu)/(M_PI*target_nu)*pow(0.01/0.99,-target_nu);
  effective_r/=exp(inflection_point);
  mu-=inflection_point;
  REPORT(mu);
  REPORT(effective_r);
#else
  double target_Z=10;
  double effective_r=epsilon*a0*K*target_Z*(1-feeding)*plant_abunance_drop;
  C=2*target_Z/target_S;
#endif

  T=feeding*epsilon/effective_r;
  r=effective_r-(multitrophic_web?K*plant_abunance_drop:0);
  //  double T=0.3/(a0*target_Z*K),r=0.2*epsilon/T;
  tau=1;
  T*=100*r;
  effective_r/=100*r;
  r/=100*r;
  
#ifdef LOGNORMALS
  animal_fraction=0.5;
#else
  animal_fraction=1/(1+target_Z/2);
#endif
  steps_transient=10;


//   mu = -392.2664;
//   sigma = 16.1713;
//   epsilon = 0.1;
//   r = 0.2/10;
//   T = 0.2*10;
//   K = 1;
//   tau = 1;
//   animal_fraction = 0.5;

//   effective_r=1;

  miniweb_t web;

  if( argc>1 && has_suffix(".cfg",argv[1]) ){
    read_parameters_from_file(argv[1]);
    ++argv;
    --argc;
  }

  if(argc>2) read(web,argc,argv);

  std::string activity;
  if(argc<=1){
    activity="run";

    REPORT(a0);
    REPORT(K);
    REPORT(epsilon);
    REPORT(T);
    REPORT(r);
    REPORT(tau);
    REPORT(C);
    REPORT(animal_fraction);
    REPORT(steps_transient);
    ALWAYS_ASSERT(r>=0);
    
    web=miniweb_t(mu,sigma,epsilon,r,T,tau,K,C,a0,animal_fraction);

    for(int i=2;i-->0;){
      web.activate_species(web.add_plant());
    }
  }else{
    activity=argv[1];
  }

  if(activity=="template"){
    write_cfg_template("tempate.cfg");
    exit(0);
  }
  
  if(activity=="evaluate"){
    web.relax(10000/web.time_scale());
    web.evaluate();
    exit(0);
  }

  if(activity=="eigen"){
    web.delete_all_inactive_species();
    web.eigenvalue_report();
    exit(0);
  }

  if(activity=="competition"){
    web.delete_all_inactive_species();
    web.competition_report();
    exit(0);
  }


  if(activity!="run"){
    FATAL_ERROR("unknown activity: " << activity);
  }

  web.add_inactive_species_if_required();
  while(++web.step_count<=100000 && !exit_now){
    REPORT(web.step_count);
    int new_species=
      web.add_fit_species();
    try{
      web.current_time=0;
      web.relax(1000/web.time_scale());
      web.save_invasion_fitnesses();
    }catch(int){
      web.delete_species(new_species);
    }
    REPORT(web.number_of_animals());
    REPORT(web.number_of_plants());
    if(skip_averaging){
      static int max_S=0;
      if(web.number_of_species()>max_S) max_S=web.number_of_species();
      REPORT(max_S);
      REPORT(web.step_count);
      if(web.step_count > max_S)
	skip_averaging=false;
      previous_size=web.number_of_species();
    }
    else
      mean_size.sample(web.number_of_species());
    REPORT(mean_size);
    REPORT(web.number_of_species());
    if(web.number_of_animals() && web.number_of_plants()){
      web.evaluate();
      if(false && web.abundance_integral.size()==
	 web.number_of_variables())
	{
	  miniweb_t web_copy(web);
	  web_copy.
	    set_abundances(web.abundance_integral);//must evaluate() web
	                                           //first!
	  
	  web_copy.delete_all_inactive_species();
	  
	  fixed_point_analyzer fpa(&web_copy);
	  fpa.set_time_scale_to(1/web.time_scale());
	  int fp_error=fpa.snap_to_fixed_point();
 	  web_copy.competition_report();
	  if(!fp_error){
	    //	    web_copy.eigenvalue_report();
	    web_copy.scale_abundance(0,1.001);
	    web_copy.relax(10000/web_copy.time_scale());
	    system("mv dynamics.dat recover.dat");
	    web_copy=web;
	    web_copy.delete_all_inactive_species();
	    web_copy.relax(10000/web_copy.time_scale());
	    system("mv dynamics.dat original.dat");
	  }
// 	  char c;
// 	  cin >> c;
      }
      if(web.step_count<=100 || 
	 0==(web.step_count%int(0.5+pow(10,-1+int(log10(double(web.step_count)))))) ){
	// save data to archive
	std::stringstream output_command(stringstream::out);
	output_command << "bzip2 -c > " << "Miniweb" << web.step_count 
		       << ".dat.bz2";
	opipestream ops(output_command.str().c_str());
	boost::archive::text_oarchive oa(ops);
	const miniweb_t cweb(web);
	oa << cweb;
      }
    }
  }
}
