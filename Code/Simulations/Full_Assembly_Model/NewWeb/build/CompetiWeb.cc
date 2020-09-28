// $Id: MoranI.cc 1155 2008-10-09 16:02:35Z axel $
#include <string>
const std::string version("$Revision: 1155 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/MoranI.cc,v $ $Date: 2008-10-09 17:02:35 +0100 (Thu, 09 Oct 2008) $");

#define LOGNORMALS
#define USE_INVASION_FITNESS

#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <complex>

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/DiagMatrix.h>
namespace CLHEP{
  double dot(const HepVector &v1, const HepVector &v2);
}
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <sys/time.h>

#include "sequence.h"
#include "Statistics.h"
#include "relax.h"
#include "random.h"
#include "snapshot.h"
#include "matrix_transformers.h"
#include "SortedVector.h"
#include "linpack_eigen.h"

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
#include "linpack_eigen.h"
#include "period_cutter.h"
#include "xy_graph.h"
#include "polyfit.h"
#include "Integrator.h"
#endif

#define SQR(X) ((X)*(X))

using namespace std;

int step_count=0;  
int steps_transient=0;
int requested_number_of_inactive_animals=50;
static int sampling_random_seed;
static int perturbing_random_seed;
static double typical_abundance=1e0;
static double log_typical_abundance=log(typical_abundance);

vector<double> vec(ODE_vector& v){
  return vector<double>(&v[0],&v[0]+v.size());
}

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);


class CompetiWeb_t : public relaxing_dynamical_object{
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


private:
  double _birth,
    _juvenile_interaction,_adult_interaction,
    _juvenile_i,_adult_i, 
    _juvenile_C, _adult_C, 
    _maturation_rate,
    _net_juvenile_population_growth,_adult_death,
    _these_are_all_population_dynamical_parameters;
  vector< SortedVector > _juvenile_a;
  vector< SortedVector > _adult_a;
  int _max_print_length;
  typedef vector<double> stage_container_t;
  stage_container_t juvenile,adult;
  vector<bool> active;
  int n_inactive_animals;
  bool _delete_by_deactivating;
  typedef enum {no=0,mean_field,local} _qna_t;
  _qna_t _is_qna;
  int _fishing_target;
  double _fishing_pressure;
  double juvenile_weight,adult_weight;// set by MF QNA

  // SERIALIZATION:
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & _birth &
      _juvenile_interaction & _adult_interaction &
      _juvenile_i & _adult_i & 
      _juvenile_C & _adult_C & 
      _maturation_rate &
      _net_juvenile_population_growth & _adult_death &
      _these_are_all_population_dynamical_parameters;
    ar & _juvenile_a & _adult_a;
    ar & juvenile & adult;
    ar & active;
    ar & n_inactive_animals;
    ar & _is_qna;
  }
  // not serialized members must be set in default constructor!
  CompetiWeb_t():
    _max_print_length(0),
    n_inactive_animals(0),
    _delete_by_deactivating(false),
    _fishing_target(0),
    _fishing_pressure(0)
  {};
  int _S() const {return juvenile.size();}
public:
  CompetiWeb_t(double birth, double maturation, 
	       double juvenile_interaction, double adult_interaction, 
	       double juvenile_i, double adult_i, 
	       double juvenile_C, double  adult_C, 
	       double juvenile_death, double adult_death, 
	       double these_are_all_population_dynamical_parameters):
    _birth(birth),
    _juvenile_interaction(juvenile_interaction),
    _adult_interaction(adult_interaction),
    _juvenile_i(juvenile_i),
    _adult_i(adult_i), 
    _juvenile_C(juvenile_C), 
    _adult_C(adult_C),
    _maturation_rate(maturation),
    _net_juvenile_population_growth(-_maturation_rate-juvenile_death),
    _adult_death(adult_death),
    _these_are_all_population_dynamical_parameters(1),
    _max_print_length(0),
    _juvenile_a(0),_adult_a(0),
    n_inactive_animals(0),
    juvenile(0),adult(0),
    active(0),
    _delete_by_deactivating(false),
    _is_qna(no),
    _fishing_pressure(0),
    _fishing_target(0)
  {};
  ~CompetiWeb_t(){};
  
  double eigenvalue_equation(double CC){
    // must be an increasing function of CC:
    return -_birth*_maturation_rate*_maturation_rate
      +(_adult_death+CC)
      *(-_net_juvenile_population_growth*_maturation_rate+
	(_adult_death+CC)*CC*
	_juvenile_interaction/_adult_interaction);
  }

  double guess_adult_nonlinearity(){
    // do a simple binary search to solve eigenvalue_equation:
    double upper=1;
    double lower=0;
    double middle;

    ALWAYS_ASSERT(eigenvalue_equation(lower)<0);

    while(eigenvalue_equation(upper)<0) upper*=2;

    middle=(upper+lower)/2;
    double ee=eigenvalue_equation(middle);
    double last_ee=2*ee;

    while(ee!=last_ee){
      if(ee<0)
	lower=middle;
      else 
	upper=middle;
      middle=(upper+lower)/2;
      last_ee=ee;
      ee=eigenvalue_equation(middle);
    }
    return middle;
  }

  CompetiWeb_t quasi_neutral_approximation(bool good=true){
    ALWAYS_ASSERT(_juvenile_C==_adult_C);
    ALWAYS_ASSERT(_juvenile_i==_adult_i);
    int S=_S();


    CompetiWeb_t QNA;
    QNA._is_qna=mean_field;
    QNA._birth=0;
    QNA._maturation_rate=0;
    QNA._adult_interaction=0;
    QNA._adult_death=0;
    QNA._adult_a.resize(0);
    QNA._adult_i=0;
    QNA._adult_C=0;

    double CC=guess_adult_nonlinearity();
    REPORT(CC);
    REPORT(CC/_adult_death);

    const double qna_growth_rate=
      (std::pow(_adult_death + CC,2)*_net_juvenile_population_growth + 
       _birth*(_adult_death + 2*CC)*_maturation_rate)/
      (std::pow(_adult_death + CC,2) + _birth*_maturation_rate);
    REPORT(qna_growth_rate);

    const double linear_growth_rate=
      (-_adult_death + _net_juvenile_population_growth +
       sqrt(4*_birth*_maturation_rate + 
	    std::pow(_adult_death + _net_juvenile_population_growth,2) ))
      /2;
    REPORT(linear_growth_rate);

    const double structure_relaxation_rate=
      +_net_juvenile_population_growth-_adult_death+CC;
    REPORT(structure_relaxation_rate);

    if(good){
      QNA._net_juvenile_population_growth=
	qna_growth_rate;
    }else{//bad: linear growth rate
      QNA._net_juvenile_population_growth=
	linear_growth_rate;
    }
    QNA._juvenile_a.resize(S);
    QNA._adult_a.resize(S);

    double denominator=(_adult_death+CC+_maturation_rate)*
      (std::pow(_adult_death+CC,2)+_birth*_maturation_rate);
    ALWAYS_ASSERT(denominator!=0);
    juvenile_weight=std::pow(_adult_death+CC,3)/denominator;
    adult_weight=_birth*std::pow(_maturation_rate,2)/denominator;

    if(!good){
      juvenile_weight=linear_growth_rate+_adult_death;
      adult_weight=_maturation_rate;
      denominator=juvenile_weight+adult_weight;
      juvenile_weight/=denominator;
      adult_weight/=denominator;
      juvenile_weight*=juvenile_weight;
      adult_weight*=adult_weight;
    }

    REPORT(juvenile_weight);
    REPORT(adult_weight);
    
    for(int i=S;i-->0;){
      QNA._juvenile_a[i].resize(S);
      QNA._adult_a[i].resize(S);
      for(int j=S;j-->0;){
	QNA._juvenile_a[i][j]=
	  juvenile_weight*_juvenile_a[i][j]+
	  adult_weight*_adult_a[i][j];
      }
    }
    QNA._juvenile_interaction=QNA._juvenile_a[0][0];

    // Generate underlying eigenvectors:
    QNA.juvenile.resize(S);
    QNA.adult.resize(S);
    QNA.ev_adult.resize(S);
    QNA.ev_juvenile.resize(S);
    QNA.adev_adult.resize(S);
    QNA.adev_juvenile.resize(S);

    const double adev_denom=std::pow(_adult_death+CC,2)+_birth*_maturation_rate;
    ALWAYS_ASSERT(adev_denom!=0);
    const double adev_j=(good?
			 (_adult_death+CC)*
			 (_adult_death+CC+_maturation_rate)/
			 adev_denom
			 :
			 1
			 );
    const double adev_a=(good?
			 _birth*(_adult_death+CC+_maturation_rate)/
			 adev_denom
			 :
			 1
			 );
    REPORT(adev_j);
    REPORT(adev_a);
    const double ad_ratio=adev_j/adev_a;
    const double ratio=(good?
			(_adult_death+CC)/_maturation_rate
			:
			(_adult_death+linear_growth_rate)/_maturation_rate
			);
    const double ev_j=ratio/(1+ratio);
    const double ev_a=1/(1+ratio);
    REPORT(ev_j);
    REPORT(ev_a);
    for(int i=S;i-->0;){
      QNA.adev_juvenile[i]=adev_j;
      QNA.adev_adult[i]=adev_a;
      QNA.ev_adult[i]=ev_a;
      QNA.ev_juvenile[i]=ev_j;
      
      QNA.juvenile[i]=adev_j*juvenile[i]+adev_a*adult[i];
      QNA.adult[i]=1;
    }
    double cc=CC*_juvenile_interaction/_adult_interaction*ratio;

    cout << "L = " << endl;
    cout << _net_juvenile_population_growth- cc<< "\t";
    cout << _birth << endl;
    cout << _maturation_rate << "\t";
    cout << -_adult_death-CC << endl;
    
    double lambda2=_net_juvenile_population_growth- cc-_adult_death-CC;
    REPORT(lambda2);
    REPORT(qna_growth_rate/lambda2);

    QNA._fishing_target=
      ( _fishing_target<S ? _fishing_target : _fishing_target-S );
    QNA._fishing_pressure=_fishing_pressure*
      ( _fishing_target<S ? adev_j : adev_a );
    
    // this is not part of the quasi neutral approximation anymore:
    QNA._juvenile_i=_juvenile_i;
    QNA._juvenile_C=_juvenile_C;

    QNA.active=active;
    QNA.n_inactive_animals=n_inactive_animals;

    return QNA;
  }

private:
  vector<double> ev_adult,ev_juvenile,adev_adult,adev_juvenile;
  vector<double> qna_growth;

  void compute_local_null_spaces(CompetiWeb_t & QNA,bool good){
    int S=_S();

    double max_juvenile=0,max_adult=0;
    for(int i=S;i-->0;){
      //ALWAYS_ASSERT(active[i]);
      if(active[i]){
	if(juvenile[i]>max_juvenile)
	  max_juvenile=juvenile[i];
	if(adult[i]>max_adult)
	  max_adult=adult[i];
      }
    }

    QNA.ev_adult.resize(S);
    QNA.ev_juvenile.resize(S);
    QNA.adev_adult.resize(S);
    QNA.adev_juvenile.resize(S);
    if(good){
      for(int i=S;i-->0;){
	double ratio=
	  (_adult_death+_adult_a[i].dot(&adult[0],max_adult))/
	  _maturation_rate;
	QNA.ev_adult[i]=1/(1+ratio);
	QNA.ev_juvenile[i]=ratio/(1+ratio);
	double ad_ratio=
	  _birth/(ratio*_maturation_rate);
	QNA.adev_adult[i]=ad_ratio/
	  (ad_ratio*QNA.ev_adult[i]+QNA.ev_juvenile[i]);
	QNA.adev_juvenile[i]=1/
	  (ad_ratio*QNA.ev_adult[i]+QNA.ev_juvenile[i]);
      }      
    }else{// not good
      for(int i=S;i-->0;){
 	double ratio=juvenile[i]/adult[i];
 	QNA.ev_adult[i]=1/(1+ratio);
 	QNA.ev_juvenile[i]=ratio/(1+ratio);
	QNA.adev_adult[i]=1;
	QNA.adev_juvenile[i]=1;
      }      
    }
  }


public:
  CompetiWeb_t local_quasi_neutral_approximation(bool good=true){
    ALWAYS_ASSERT(_juvenile_C==_adult_C);
    ALWAYS_ASSERT(_juvenile_i==_adult_i);
    int S=_S();

    CompetiWeb_t QNA;
    QNA._is_qna=local;
    QNA._birth=0;
    QNA._maturation_rate=0;
    QNA._adult_interaction=0;
    QNA._adult_death=0;
    QNA._adult_a.resize(0);
    QNA._adult_i=0;
    QNA._adult_C=0;

    QNA.qna_growth.resize(S);
    QNA.juvenile.resize(S);
    QNA.adult.resize(S);
    QNA._juvenile_a.resize(S);
    QNA._adult_a.resize(S);

    compute_local_null_spaces(QNA,good);


    for(int i=S;i-->0;){
      QNA.qna_growth[i]=
	(_birth*QNA.ev_adult[i]+
	 _net_juvenile_population_growth*QNA.ev_juvenile[i] )*
	QNA.adev_juvenile[i]+
	(_maturation_rate*QNA.ev_juvenile[i]-
	 _adult_death*QNA.ev_adult[i] )*
	QNA.adev_adult[i];
    }
    for(int i=S;i-->0;){
      for(int j=S;j-->0;){
	QNA._juvenile_a[i][j]=
	  (_juvenile_a[i][j]*QNA.ev_juvenile[j]*QNA.ev_juvenile[i] )*
	    QNA.adev_juvenile[i]+
	  (_adult_a[i][j]*QNA.ev_adult[j]*QNA.ev_adult[i] )*
	  QNA.adev_adult[i];
      }
    }
    const double linear_growth_rate=
      (-_adult_death + _net_juvenile_population_growth +
       sqrt(4*_birth*_maturation_rate + 
	    std::pow(_adult_death + _net_juvenile_population_growth,2) ))
      /2;
    REPORT(linear_growth_rate);
    double qna_growth_sum=0;
    for(int i=S;i-->0;){
      qna_growth_sum+=QNA.qna_growth[i];
    }
//     double factor=linear_growth_rate/(qna_growth_sum/S);
//     if(!good){
//       for(int i=S;i-->0;){
// 	for(int j=S;j-->0;){
// 	  QNA._juvenile_a[i][j]=QNA._juvenile_a[i][j]*factor;
// 	}
// 	QNA.qna_growth[i]*=factor;
//       }
//     }
    for(int i=S;i-->0;){
      QNA.juvenile[i]=
	juvenile[i]*QNA.adev_juvenile[i]+
	adult[i]*QNA.adev_adult[i];
      QNA.adult[i]=1;
    }

    QNA._fishing_target=
      ( _fishing_target<S ? _fishing_target : _fishing_target-S );
    QNA._fishing_pressure=_fishing_pressure*
      ( _fishing_target<S ? QNA.adev_juvenile[QNA._fishing_target]
	: QNA.adev_adult[QNA._fishing_target] );

    // this is not part of the quasi neutral approximation anymore:
    QNA._juvenile_i=_juvenile_i;
    QNA._juvenile_C=_juvenile_C;

    QNA.active=active;
    QNA.n_inactive_animals=n_inactive_animals;

    return QNA;
  }

  void delete_all_inactive_species(){
    for(int i=_S();i-->0;){
      if(!active[i]){
	delete_species(i);
      }
    }
  }

    
  void delete_by_deactivating(){
    _delete_by_deactivating=true;
  }
  double abundance_variation(const CompetiWeb_t &previous) const{
    ALWAYS_ASSERT(_S() == previous._S());
    double square_sum=0;
    int n=0;
    for(int i=_S();i-->0;){
      if(previous.active[i]){
	if(!active[i]){
	  //species i went extinct
	  square_sum+=std::pow(previous.juvenile[i],2);
	}else{
	  square_sum+=std::pow(previous.juvenile[i]-juvenile[i],2);
	}
	n++;
      }
    }
    return sqrt(square_sum/n);
  }

  virtual void read_state_from(const ODE_vector & state){
    int S=_S();
    for(int i=S;i-->0;)
      juvenile[i]=( active[i] ? exp(state[i]) : state[i] );
    for(int i=S;i-->0;)
      adult[i]=( active[i] ? exp(state[i+S]) : state[i+S] );
  }
  virtual void write_state_to(ODE_vector & state) const{
    int S=_S();
    for(int i=S;i-->0;)
      state[i]=( active[i] ? log(juvenile[i]) : juvenile[i] );
    for(int i=S;i-->0;)
      state[i+S]=( active[i] ? log(adult[i]) : adult[i] );
  }
  virtual int number_of_variables() const{
    return 2*_S();
  }
  int number_of_species(){
    return _S()-n_inactive_animals;
  }
  virtual int dynamics(ODE_vector const & state,
		       ODE_vector & time_derivative){
    int S=_S();
    read_state_from(state);

    if(_is_qna==local){
      return local_qna_dynamics(time_derivative);
    }

    double max_juvenile=0,max_adult=0;

    for(int i=S;i-->0;){
      if(!active[i]){
	juvenile[i]=0;
	adult[i]=0;
      }else{
	if(juvenile[i]>max_juvenile)
	  max_juvenile=juvenile[i];
	if(adult[i]>max_adult)
	  max_adult=adult[i];
      }
    }

    for(int i=S;i-->0;){
      double adult_juvenile_ratio=
	( active[i] ? adult[i]/juvenile[i] : exp(state[i+S]-state[i]) );

      if(!adult_juvenile_ratio){
	adult_juvenile_ratio=1;//catch some crap when _delete_by_deactivating
      }
      
      time_derivative[i]=
	_birth*adult_juvenile_ratio;
      time_derivative[i]+=
	_net_juvenile_population_growth;
      time_derivative[i]-=
	_juvenile_a[i].dot(&juvenile[0],max_juvenile);
      time_derivative[i+S]=
	_maturation_rate/adult_juvenile_ratio;//=maturation
      time_derivative[i+S]-=
	_adult_death;
      time_derivative[i+S]-=
	_adult_a[i].dot(&adult[0],max_adult);
    }
    time_derivative[_fishing_target]-=_fishing_pressure/
      (_fishing_target<S?juvenile[_fishing_target]:adult[_fishing_target-S]);
    return 0;
  }

  int local_qna_dynamics(ODE_vector & time_derivative){
    int S=_S();

    double max_juvenile=0;

    for(int i=S;i-->0;){
      if(!active[i]){
	juvenile[i]=0;
	adult[i]=0;
      }else{
	if(juvenile[i]>max_juvenile)
	  max_juvenile=juvenile[i];
      }
    }

    for(int i=S;i-->0;){
      time_derivative[i]=qna_growth[i];
      time_derivative[i]-=
	_juvenile_a[i].dot(&juvenile[0],max_juvenile);
      time_derivative[i+S]=0;
    }
    time_derivative[_fishing_target]-=_fishing_pressure/
      (juvenile[_fishing_target]);
    return 0;
  }

  template<typename MATRIX>
  int Jacobian(ODE_vector const & state,
	       ODE_vector const & dynamics,
	       MATRIX & jac){
    // ATTENTION: this is the Jacobian in linear, not in logarithmic
    // coordinates.  So this is not the Jacobian needed by
    // ODE_dynamical_object !
    int S=_S();
    read_state_from(state);
    
    double max_juvenile=0,max_adult=0;
    for(int i=S;i-->0;){
      if(!active[i]){
	juvenile[i]=0;
	adult[i]=0;
      }else{
	if(juvenile[i]>max_juvenile)
	  max_juvenile=juvenile[i];
	if(adult[i]>max_adult)
	  max_adult=adult[i];
      }
    }
    
    if(_is_qna!=local){
      for(int i=S;i-->0;){
	double adult_juvenile_ratio=
	  ( active[i] ? adult[i]/juvenile[i] : exp(state[i+S]-state[i]) );
	
	if(!adult_juvenile_ratio){
	  adult_juvenile_ratio=1;//catch some crap when _delete_by_deactivating
	}
	
	for(int j=S;j-->0;){
	  jac[i][j]=-_juvenile_a[i][j]*juvenile[i];
	  jac[i+S][j+S]=-_adult_a[i][j]*adult[i];
	  jac[i][j+S]=0;
	  jac[i+S][j]=0;
	}
	
	jac[i][i]+=_net_juvenile_population_growth -
	_juvenile_a[i].dot(&juvenile[0],max_juvenile);
	jac[i][i+S]=_birth;
	jac[i+S][i]=_maturation_rate;
	jac[i+S][i+S]-=_adult_death +
	  _adult_a[i].dot(&adult[0],max_adult);
	
      }
    }else{ //local qna:
      for(int i=S;i-->0;){
	for(int j=S;j-->0;){
	  jac[i][j]=-_juvenile_a[i][j]*juvenile[i];
	  jac[i+S][j+S]=0;
	  jac[i][j+S]=0;
	  jac[i+S][j]=0;
	}
	
	jac[i][i]+=qna_growth[i] -
	_juvenile_a[i].dot(&juvenile[0],max_juvenile);
	
      }      
    }
    return 0;
  }
    

  double active_sum(const ODE_vector & s) const{
    int S=_S();
    double sum=0;
    for(int i=S;i-->0;){
      if(active[i]){
	sum+=s[i];
	sum+=s[i+S];
      }
    }
    return sum;
  }
  bool is_active(int i) const {
    int S=_S();
    return active[i%S];
  }
private:
  double start_recording_time;
  double last_recording_time;
  vector<double> abundance_integral;
  vector<double> effective_abundance_integral;

public:
  void add_inactive_species_if_required(){
#ifdef USE_INVASION_FITNESS
    while(n_inactive_animals<requested_number_of_inactive_animals){
      add_animal();
    }
#endif
  }

private:
  void prepare_for_integration(){
    abundance_integral.resize(number_of_variables());
    last_recording_time=DBL_MAX;
    
    for(int i=_S();i-->0;){
      if(!active[i]){
	juvenile[i]=0;
	adult[i]=0;
      }
    }
  }
  void record_for_steady_state(){
    int S=_S();
    
    if(current_time < last_recording_time){
      last_recording_time=current_time;
      start_recording_time=current_time;
      for(int i=abundance_integral.size();i-->0;){
	abundance_integral[i]=0;
      }
      for(int i=effective_abundance_integral.size();i-->0;){
	effective_abundance_integral[i]=0;
      }
    }else{
      double delta_t=current_time-last_recording_time;
      for(int i=S;i-->0;){
	abundance_integral[i]+=juvenile[i]*delta_t;
	abundance_integral[i+S]+=adult[i]*delta_t;
      }
    }
  }

public:
  int add_animal(){
    const double asym=1;

    int i=_S();
    roll_up(i,juvenile);
    juvenile[i]=0;
    roll_up(i,adult);
    adult[i]=0;
    roll_up(i,active);
    active[i]=false;
    n_inactive_animals++;

    roll_up(i,_juvenile_a);
    _juvenile_a[i].resize(i+1);
    for(int ii=i;ii-->0;){
      roll_up(i,_juvenile_a[ii]);
      double base_interaction_strength=
	(unirand() < _juvenile_C ? _juvenile_i*_juvenile_interaction : 0);
      _juvenile_a[i][ii]=base_interaction_strength*
	(1+asym*(1-2*unirand()));
      _juvenile_a[ii][i]=base_interaction_strength*
	(1+asym*(1-2*unirand()));
    }
    _juvenile_a[i][i]=_juvenile_interaction;


    roll_up(i,_adult_a);
    _adult_a[i].resize(i+1);
    for(int ii=i;ii-->0;){
      roll_up(i,_adult_a[ii]);
      double base_interaction_strength=
	(unirand() < _adult_C ? _adult_i*_adult_interaction : 0);
      _adult_a[i][ii]=base_interaction_strength*
	(1+asym*(1-2*unirand()));
      _adult_a[ii][i]=base_interaction_strength*
	(1+asym*(1-2*unirand()));
    }
    _adult_a[i][i]=_adult_interaction;
#if DEBUGGING
    WARNING("added animal " << i);
#endif
    return i;
  }
  int add_species(){
    return add_animal();
  }
  int add_active_species(){
    int i=add_species();
    juvenile[i]=1e-5*typical_abundance;
    adult[i]=1e-5*typical_abundance;
    active[i]=true;
    n_inactive_animals--;
    return i;
  }

  void delete_species(int i){
#ifdef DEBUGGING
    check_internal_consistency();
#endif
    if(!_delete_by_deactivating){
      delete_replacing(i,juvenile);
      delete_replacing(i,adult);
      if(!active[i]) n_inactive_animals--;
      delete_replacing(i,active);
      delete_replacing(i,_juvenile_a);
      delete_replacing(i,_adult_a);
      for(int j=_S();j-->0;){
	delete_replacing(i,_juvenile_a[j]);
	delete_replacing(i,_adult_a[j]);
      }
    }else{// delete by deactivating:
      if(active[i]){
	n_inactive_animals++;
	REPORT(juvenile[i]);
	WARNING("deactivated animal " << i);
      }
      if(!active[i])
	WARNING("species " << i << " was not active");
      active[i]=false;
      juvenile[i]=0;
      adult[i]=0;
      _juvenile_a[i].resize(0);
      _adult_a[i].resize(0);
    }
#if DEBUGGING
    WARNING("deleted animal " << i);
#endif
  }
  void check_internal_consistency(){
    int S=_S();
    ALWAYS_ASSERT(_juvenile_a.size()==S);
    ALWAYS_ASSERT(_adult_a.size()==S);
    for(int i=S;i-->0;){
      if(active[i]){
	ALWAYS_ASSERT(_juvenile_a[i].size()==S);
	ALWAYS_ASSERT(_adult_a[i].size()==S);
      }
    }
  }
  int add_fit_species(){
#ifdef USE_INVASION_FITNESS
    while(n_inactive_animals){
      int i=0; // pointer to next animal to test
      while(active[i]) i++;
//       REPORT(juvenile[i]);
//       REPORT(adult[i]);
      if(juvenile[i]+adult[i]<0){
	delete_species(i);
      }else{
	juvenile[i]=1e-5*typical_abundance;
	adult[i]=1e-5*typical_abundance;
	active[i]=true;
	n_inactive_animals--;
#ifdef DEBUGGING
	WARNING("activated animal " << i );
	check_internal_consistency();
#endif
	add_inactive_species_if_required();
	return i;
      }
    }     
    // we ran out of inactive species
    WARNING("consider using more inactive species");
    add_inactive_species_if_required();
    relax(1000);
    return add_fit_species();
#else  // not using invasion fitness
    int n=number_of_variables()+1;
    int S=_S();
    ODE_vector l(n),dl(n);
    while(true){
      int i=add_species();
      ALWAYS_ASSERT(!active[i]);
      active[i]=true;
      n_inactive_animals--;
      juvenile[i]=1e-5*typical_abundance;
      adult[i]=1e-5*typical_abundance;
#ifdef DEBUGGING
      check_internal_consistency();
#endif
      write_state_to(l);
      dynamics(l,dl);
      
      if(dl[i]+dl[i+S]<0){
	delete_species(i);
      }else{
	return i;
      }
    }
#endif
  }
  void invasion_fitness_distribution(const char * filename){
    ofstream f(filename);
    const int S=_S();
    int i=0; // pointer to next animal to test
    while(i<S){
      while(active[i]) ++i;
      f << step_count << " " 
	<< juvenile[i]/(current_time-start_recording_time) << endl;
      ++i;
    }
  }
  void invasion_fitness_distribution(){
    std::stringstream output_file_name(stringstream::out);
    output_file_name << "fitness" << step_count << ".dat";
    invasion_fitness_distribution(output_file_name.str().c_str());
  }

  bool small_values_in(ODE_vector & state,
		       const species_set_t&  conserved){
    int S=_S();
    for(int i=state.size();i-->0;){
      if(active[i%S] && state[i]<-20.0*M_LN10+log_typical_abundance) 
	return true;
    }
    return false;
  }

  virtual species_set_t 
  delete_species_larger_than_exp(const sequence<double> & si,
				 const species_set_t&  conserved){
    species_set_t deleted;
    int S=_S();

    for(int i=S;i-->0;){
      if(active[i] && 
	 (si(i)<-15.0*M_LN10+log_typical_abundance || 
	  si(i+S)<-15.0*M_LN10+log_typical_abundance) ){
	delete_species(i);
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
    int S=_S();
    for(int i=S;i-->0;){
      if(active[i] && 
	 (juvenile[i] < exp(-15.0*M_LN10)*typical_abundance ||
	  adult[i] < exp(-15.0*M_LN10)*typical_abundance )){
	delete_species(i);
#ifdef DEBUGGING
	check_internal_consistency();
#endif
	deleted.insert(i);
      }
    }
    return deleted;
  }    
  void evaluate(){
    REPORT(number_of_species());

    int n_inactive_species=n_inactive_animals;

    for(int i=abundance_integral.size();i-->0;){
      abundance_integral[i]/=current_time-start_recording_time;
    }
    average_meter animal_dist;
    for(int i=_S();i-->0;){
      if(active[i]){
	animal_dist.sample(log(abundance_integral[i]));
      }
    }
    double sigma_animal=sqrt(animal_dist.var());
    REPORT(sigma_animal);
    REPORT(animal_dist);
    {
#ifdef DEBUGGING
      check_internal_consistency();
#endif
      ofstream a("animals.dat");
      vector<double> sa(number_of_species());
      vector<double>::iterator j;
      j=sa.begin();
      for(int i=_S();i-->0;){
    	if(active[i]){
    	  *j++ = juvenile[i]+adult[i];
    	}
      }
      ALWAYS_ASSERT(j==sa.end());
      
      sort(sa.begin(),sa.end());
      for(int i=0;i<sa.size();i++)
	a << sa[i] << " " << double(sa.size()-i)/sa.size() << endl;
    }
    
  }

  int most_abundant_species(){
    int S=_S();
    double max_abunance=0;
    int maxi=0;
    for(int i=S;i-->0;){
      if(active[i] && 
	 (max_abunance < juvenile[i]+(_is_qna?0:adult[i]))){
	max_abunance=juvenile[i]+(_is_qna?0:adult[i]);
	maxi=i;
      }
    }
    return maxi;
  }

  class less_abundant{
    const CompetiWeb_t & obj ;
  public:
    less_abundant(const CompetiWeb_t & o): obj(o){};
    bool operator()(int i,int j) const {
      return obj.population(i) < obj.population(j);
    }
  };

  int nth_abundant_species(int i){
    int S=_S();
    ALWAYS_ASSERT(i>0);
    ALWAYS_ASSERT(i<=S);
    vector<int> order(S);
    for(int j=S;j-->0;) order[j]=j;
    sort(order.begin(),order.end(),less_abundant(*this));
    return order[S-i];
  }

  void line_print(ODE_vector const & state, ostream &co){
    int S=_S();
    double log10=log(10);
    if(_max_print_length<S){
      _max_print_length=S;
    }
    for(int i=0;i<S;i++){
      if(active[i]){
	if(_is_qna){
	  // print only juveniles:
	  co << exp(state[i]) << " ";
	}else{
	  co << exp(state[i])+exp(state[i+S]) << " ";
	}
      }else if(_delete_by_deactivating){
	co << 00 << " ";
      }
    }
    for(int i=_max_print_length-S;i-->0;){
      co << 0.0 << " ";
    }
  }
  double adult_population_sum(){
    double s=0;
    for(int i=_S();i-->0;){
      if(active[i]) 
	s+=adult[i];
    }
    return s;
  }
  double juvenile_population_sum(){
    double s=0;
    for(int i=_S();i-->0;){
      if(active[i]) 
	s+=juvenile[i];
    }
    return s;
  }

  double quick_norm(const double &v,int n){
    double sum=0;
    for(int i=n;i-->0;){
      sum+=(&v)[i]*(&v)[i];
    }
    return sqrt(sum);
  }

  class ddr_larger {
    const CLHEP::HepVector &_ddr;
  public:
    ddr_larger(CLHEP::HepVector &dr):_ddr(dr){};
    bool operator()(int i,int j){return _ddr[i]>_ddr[j];}
  };

  void spectral_radius(){
    // This actually does something quite different form computing a
    // spectral radius.
    ALWAYS_ASSERT(!_is_qna);
    int S=_S();
    
    CompetiWeb_t qna=local_quasi_neutral_approximation(true);
    CompetiWeb_t mf_qna=quasi_neutral_approximation(true);
    int n=qna.number_of_variables();
    ALWAYS_ASSERT(n==2*S);
    
    CLHEP::HepMatrix J0D(n,n,0),A(n,n,0);
    CLHEP::HepDiagMatrix N(n);
    
    for(int i=S;i-->0;){
      N[i][i]=N[i+S][i+S]=juvenile[i]+adult[i];
      for(int j=S;j-->0;){
	A[i][j]=-_juvenile_a[i][j]*juvenile[i]/N[i][i];
	A[i+S][j+S]=-_adult_a[i][j]*adult[i]/N[i][i];
      }
    }	
    double max_juvenile=0,max_adult=0;
    for(int i=S;i-->0;){
      if(juvenile[i]>max_juvenile)
	max_juvenile=juvenile[i];
      if(adult[i]>max_adult)
	max_adult=adult[i];
    }
    CLHEP::HepVector Lambda0_inverse_Diag(S);
    average_meter tau_intr;
    for(int k=S;k-->0;){
      Lambda0_inverse_Diag[k]=
	1/(_net_juvenile_population_growth
	   -_juvenile_a[k].dot(&juvenile[0],max_juvenile)
	   -_adult_death
	   -_adult_a[k].dot(&adult[0],max_adult)
	   );
      tau_intr.sample(-Lambda0_inverse_Diag[k]);
    }
    for(int i=S;i-->0;){
      J0D[i][i]=qna.adev_adult[i]*qna.ev_adult[i]*
	Lambda0_inverse_Diag[i];
      J0D[i+S][i+S]=qna.adev_juvenile[i]*qna.ev_juvenile[i]*
	Lambda0_inverse_Diag[i];
      J0D[i][i+S]=-qna.adev_adult[i]*qna.ev_juvenile[i]*Lambda0_inverse_Diag[i];
      J0D[i+1][i]=-qna.adev_juvenile[i]*qna.ev_adult[i]*Lambda0_inverse_Diag[i];
    }
    
    CLHEP::HepMatrix vv(S,n,0),ww(n,S,0);
    for(int i=S;i-->0;){
      vv[i][i]=qna.adev_juvenile[i];
      vv[i][i+S]=qna.adev_adult[i];
      ww[i][i]=qna.ev_juvenile[i];
      ww[i+S][i]=qna.ev_adult[i];
    }
    
    CLHEP::HepMatrix Pi(S,S,-1.0/S);//Projector onto non-constant subspace.
    for(int i=S;i-->0;){
      Pi[i][i]+=1;
    }
    
    CLHEP::HepMatrix x=(vv*N*A)*J0D*(N*(A*ww))*Pi;
    CLHEP::HepMatrix a=vv*N*A*ww*Pi;
    x=x.T()*x;
    a=a.T()*a;
    double Delta_theo=sqrt(x.trace()/a.trace());
    REPORT(Delta_theo);
    double z=juvenile_weight*_juvenile_interaction/
      adult_weight/_adult_interaction;
    double Delta_simp=
      mf_qna._net_juvenile_population_growth*_adult_i*tau_intr*2/
      (z+1/z);
    REPORT(Delta_simp);
//     eigen_report(x,"x.dat");
//     eigen_report(a,"a.dat");
  }
  
  typedef complex<double> cmplx;
  double abs2(cmplx x){return x.real()*x.real()+x.imag()*x.imag();}
  double abs(cmplx x){return sqrt(abs2(x));}

  void eigenvalue_perturbation_analysis_2(int n_values){
    ALWAYS_ASSERT(!_is_qna);
    int S=_S();
    
    CompetiWeb_t qna=local_quasi_neutral_approximation(true);
    CompetiWeb_t mf_qna=quasi_neutral_approximation(true);
    int n=qna.number_of_variables();
    ALWAYS_ASSERT(n==2*S);
    
    CLHEP::HepMatrix J0D(n,n,0),AA(n,n,0);
    CLHEP::HepDiagMatrix NN(n);
    CLHEP::HepDiagMatrix N(S);
    
    for(int i=S;i-->0;){
      NN[i][i]=NN[i+S][i+S]=juvenile[i]+adult[i];
      N[i][i]=NN[i][i];
      for(int j=S;j-->0;){
	AA[i][j]=-_juvenile_a[i][j]*juvenile[i]/NN[i][i];
	AA[i+S][j+S]=-_adult_a[i][j]*adult[i]/NN[i][i];
      }
    }	
    double max_juvenile=0,max_adult=0;
    for(int i=S;i-->0;){
      if(juvenile[i]>max_juvenile)
	max_juvenile=juvenile[i];
      if(adult[i]>max_adult)
	max_adult=adult[i];
    }
    CLHEP::HepVector Lambda0_inverse_Diag(S);
    average_meter tau_intr;
    for(int k=S;k-->0;){
      Lambda0_inverse_Diag[k]=
	1/(_net_juvenile_population_growth
	   -_juvenile_a[k].dot(&juvenile[0],max_juvenile)
	   -_adult_death
	   -_adult_a[k].dot(&adult[0],max_adult)
	   );
      tau_intr.sample(-Lambda0_inverse_Diag[k]);
    }
    for(int i=S;i-->0;){
      J0D[i][i]=qna.adev_adult[i]*qna.ev_adult[i]*
	Lambda0_inverse_Diag[i];
      J0D[i+S][i+S]=qna.adev_juvenile[i]*qna.ev_juvenile[i]*
	Lambda0_inverse_Diag[i];
      J0D[i][i+S]=-qna.adev_adult[i]*qna.ev_juvenile[i]*Lambda0_inverse_Diag[i];
      J0D[i+1][i]=-qna.adev_juvenile[i]*qna.ev_adult[i]*Lambda0_inverse_Diag[i];
    }
    
    CLHEP::HepMatrix vv(S,n,0),ww(n,S,0);
    for(int i=S;i-->0;){
      vv[i][i]=qna.adev_juvenile[i];
      vv[i][i+S]=qna.adev_adult[i];
      ww[i][i]=qna.ev_juvenile[i];
      ww[i+S][i]=qna.ev_adult[i];
    }
    
    //CLHEP::HepMatrix A=vv*AA*ww;
    CLHEP::HepMatrix M=vv*NN*AA*ww;
    REPORT(M.num_row());
    REPORT(M.num_col());

    CLHEP::HepVector lambda_r(S),lambda_i(S);
    vector< vector<cmplx> > eigenvectors(S,vector<cmplx>(S));
    vector< vector<cmplx> > adeigenvectors(S,vector<cmplx>(S));
    vector<int> eval_order(S);
    ODE_vector fixed_point(S);
    CLHEP::HepMatrix eigenv(S,S);
    CLHEP::HepMatrix adeigenv(S,S);
    eigen(M,lambda_r,lambda_i,adeigenv,eigenv); 

    for(int i=0;i<S;i++){
      for(int u=0;u<S;u++){
	if(lambda_i[u]==0){
	  eigenvectors[u][i]=eigenv[u][i];
	  adeigenvectors[u][i]=adeigenv[u][i];
	}else{//undo lapack convention:
	  ALWAYS_ASSERT(u<S-1);
	  eigenvectors[u][i]=cmplx(eigenv[u][i],eigenv[u+1][i]);
	  eigenvectors[u+1][i]=cmplx(eigenv[u][i],-eigenv[u+1][i]);
	  adeigenvectors[u][i]=
	    cmplx(adeigenv[u][i],-adeigenv[u+1][i]);
	  adeigenvectors[u+1][i]=
	    cmplx(adeigenv[u][i],+adeigenv[u+1][i]);
	  u++;
	}
      }
    }

    for(int u=S;u-->0;){
      cmplx sum=0;
      for(int i=S;i-->0;){
	sum+=adeigenvectors[u][i]*eigenvectors[u][i];
      }
      double sum2=0;
      for(int i=S;i-->0;){
	adeigenvectors[u][i]/=sum;
	sum2+=(eigenvectors[u][i]*conj(eigenvectors[u][i])).real();
      }
    }

    for(int i=S;i-->0;){
      eval_order[i]=i;
    }

    stable_sort(eval_order.begin(),eval_order.end(),ddr_larger(lambda_r));
    for(int uu=0;uu < min<int>(4,S);uu++){
      int u=eval_order[uu];
      REPORT(lambda_r[u]);
    }

    CLHEP::HepMatrix
      eigenvec_r(S,S),eigenvec_i(S,S),
      adeigenvec_r(S,S),adeigenvec_i(S,S);

    for(int u=S;u-->0;){
      for(int i=S;i-->0;){
	eigenvec_r[u][i]=eigenvectors[u][i].real();
	eigenvec_i[u][i]=eigenvectors[u][i].imag();
	adeigenvec_r[u][i]=adeigenvectors[u][i].real();
	adeigenvec_i[u][i]=adeigenvectors[u][i].imag();
      }
    }

    // TEST:
//     CLHEP::HepMatrix ONE=
//       -eigenvec_r*adeigenvec_r.T()
//       +eigenvec_i*adeigenvec_i.T();
//     eigen_report(ONE,"one.dat");


    CLHEP::HepMatrix J1=NN*AA;
    CLHEP::HepDiagMatrix II(n,1);
    CLHEP::HepMatrix R=vv*J1*J0D*J1*ww;
    //CLHEP::HepMatrix R=vv*J1*(II-ww*vv)*J1*ww*tau_intr;
    //CLHEP::HepMatrix R=vv*J1*J1*ww*tau_intr;

    weighted_average_meter Delta_lambda;
    for(int uu=0;uu<S;uu++){
      int u=eval_order[uu];
      int up=u+1;
      CLHEP::HepVector a_r=eigenvec_r.sub(up,up,1,S).T();
      CLHEP::HepVector a_i=eigenvec_i.sub(up,up,1,S).T();
      CLHEP::HepVector b_r=adeigenvec_r.sub(up,up,1,S).T();
      CLHEP::HepVector b_i=adeigenvec_i.sub(up,up,1,S).T();
      CLHEP::HepVector Rar=R*a_r;
      CLHEP::HepVector Rai=R*a_i;
      cmplx lambda=cmplx(lambda_r[u],lambda_i[u]);
      cmplx correction=
	cmplx(-CLHEP::dot(b_r,Rar)+CLHEP::dot(b_i,Rai),
	      +CLHEP::dot(b_r,Rai)+CLHEP::dot(b_i,Rar));
      //REPORT(lambda);
      REPORT(abs(correction/lambda/lambda));
      double excitability2=0;
      for(int i=S;i-->0;){
	excitability2+=abs2(adeigenvectors[u][i]*N[i][i]/lambda);
      }
      Delta_lambda.sample(abs2(correction/lambda),1/abs2(lambda));
    }
    REPORT(sqrt(double(Delta_lambda))/2);
  }
  
  void eigenvalue_perturbation_analysis(int n_values){
    ALWAYS_ASSERT(!_is_qna);
    int S=_S();
    n_values=min<int>(S,n_values);
    
    CompetiWeb_t qna=local_quasi_neutral_approximation(false);
    int n=qna.number_of_variables();
    ALWAYS_ASSERT(n==2*S);

    ODE_vector state(n),dynamics(n);
    qna.write_state_to(state);
    CLHEP::HepMatrix preJ(n,n);
    qna.Jacobian(state,dynamics,preJ);
    CLHEP::HepMatrix J=preJ.sub(1,S,1,S);
    CLHEP::HepMatrix a(S,S),b(S,S);
    CLHEP::HepVector ddr(J.num_row()),ddi(J.num_row()),NN(J.num_row());
    eigen(J,ddr,ddi,b,a); // b and a are orthogonal
    for(int u=S;u-->0;){ // ... but need to be made orthonormal
      double sum=0;
      for(int i=S;i-->0;){
	sum+=b[u][i]*a[u][i];
      }
      for(int i=S;i-->0;){
	b[u][i]/=sum;
      }
    }
    
    vector<int> eval_order(S),rank(S);
    for(int i=S;i-->0;){
      eval_order[i]=i;
      rank[i]=i;
      NN[i]=qna.juvenile[i];
    }
    stable_sort(eval_order.begin(),eval_order.end(),ddr_larger(ddr));
    stable_sort(rank.begin(),rank.end(),ddr_larger(NN));

    double max_juvenile=0,max_adult=0;
    for(int i=S;i-->0;){
      if(!active[i]){
	juvenile[i]=0;
	adult[i]=0;
      }else{
	if(juvenile[i]>max_juvenile)
	  max_juvenile=juvenile[i];
	if(adult[i]>max_adult)
	  max_adult=adult[i];
      }
    }

    CLHEP::HepVector Lambda0_inverse_Diag(S);
    for(int k=S;k-->0;){
      Lambda0_inverse_Diag[k]=
	1/(_net_juvenile_population_growth
	   -_juvenile_a[k].dot(&juvenile[0],max_juvenile)
	   -_adult_death
	   -_adult_a[k].dot(&adult[0],max_adult)
	   );
    }

    CLHEP::HepMatrix Law(S,2*S); // Lambda*a*w vectors for each eigenvalue.
    CLHEP::HepMatrix LiLaw(S,2*S); 
    CLHEP::HepMatrix LLiLaw(S,2*S);
    for(int uu=0;uu<n_values;uu++){
      int u=eval_order[uu];
      int n=rank[S-1-uu];
      for(int k=S;k-->0;){
	double sum=0;
	for(int l=S;l-->0;){
	  sum+=_juvenile_a[k][l]*a[u][l]*qna.ev_juvenile[l];
	}
	Law[u][k]=sum*juvenile[k];
      }
      for(int k=S;k-->0;){
	double sum=0;
	for(int l=S;l-->0;){
	  sum+=_adult_a[k][l]*a[u][l]*qna.ev_adult[l];
	}
	Law[u][S+k]=sum*adult[k];
      }
      cout << "ev " << uu << ": " << ddr[u] << endl;
      cout << "est" << uu << ": " << -NN[n]*qna._juvenile_a[n][n] << endl;
      cout << "Law   " << ": " << quick_norm(Law[u][0],2*S) << endl;

      for(int j=S;j-->0;){
	double projection=
	  qna.ev_adult[j]*Law[u][j]-qna.ev_juvenile[j]*Law[u][S+j];
	projection*=Lambda0_inverse_Diag[j];
	LiLaw[u][j]=qna.adev_adult[j]*projection;
	LiLaw[u][j+S]=-qna.adev_juvenile[j]*projection;
      }
      cout << "LiLaw " << ": " << quick_norm(LiLaw[u][0],2*S) << endl;

      for(int k=S;k-->0;){
	double sum=0;
	for(int l=S;l-->0;){
	  sum+=_juvenile_a[k][l]*LiLaw[u][l];
	}
	LLiLaw[u][k]=sum*juvenile[k];
      }
      for(int k=S;k-->0;){
	double sum=0;
	for(int l=S;l-->0;){
	  sum+=_adult_a[k][l]*LiLaw[u][S+l];
	}
	LLiLaw[u][S+k]=sum*adult[k];
      }
      cout << "LLiLaw" << ": " << quick_norm(LLiLaw[u][0],2*S) << endl;

      double lambda2=0,dumb_lambda2=0;
      for(int i=S;i-->0;){
	lambda2-=b[u][i]*
	  (qna.adev_juvenile[i]*LLiLaw[u][i]+
	   +qna.adev_adult[i]*LLiLaw[u][S+i] );
	dumb_lambda2-=b[u][i]*
	  (qna.adev_juvenile[i]*juvenile[i]+
	   +qna.adev_adult[i]*adult[i] );
      }
      REPORT(lambda2);
      REPORT(lambda2/ddr[u]);
      REPORT(dumb_lambda2);
      ostringstream bname;
      bname << "b" << uu << ".dat";
      ofstream bb(bname.str().c_str());
      for(int i=0;i<S;i++){
	bb << b[u][rank[i]] << endl;
      }
      ostringstream aname;
      aname << "a" << uu << ".dat";
      ofstream aa(aname.str().c_str());
      for(int i=0;i<S;i++){
	aa << a[u][rank[i]] << endl;
      }
      double test_norm=0;
      for(int i=0;i<S;i++){
	test_norm+=a[u][rank[i]]*b[u][rank[i]];
      }
      REPORT(test_norm);
    }
    ofstream softness("softness.dat");
    ofstream estimation("estimation.dat");
    ofstream focus("focus.dat");
    for(int uu=0;uu<S;uu++){
      int u=eval_order[uu];
      int n=rank[S-1-uu];
      double sensitivity2_sum=0;
      for(int i=S;i-->0;){
	double sensitivity=b[u][i]*NN[i];
	sensitivity2_sum+=sensitivity*sensitivity;
      }
      double abs_lambda=sqrt(ddr[u]*ddr[u]+ddi[u]*ddi[u]);
      if(true or ddi[u]==0){
	softness << -ddr[u] << " "
		 << sqrt(sensitivity2_sum/S)/abs_lambda << endl;
	estimation  << -ddr[u] << " " 
		    << NN[n]*qna._juvenile_a[n][n] << endl;
	focus << -ddr[u] << " " << b[u][n]*a[u][n] << endl;
      }
    }
  }

  void eigenvalue_singular_analysis(int n_values){
    ALWAYS_ASSERT(!_is_qna);
    int S=_S();
    n_values=min<int>(S,n_values);
    
    CompetiWeb_t qna=local_quasi_neutral_approximation(true);
    int n=qna.number_of_variables();
    ALWAYS_ASSERT(n==2*S);

    ODE_vector state(n),dynamics(n);
    qna.write_state_to(state);
    CLHEP::HepMatrix preJ(n,n);
    qna.Jacobian(state,dynamics,preJ);
    CLHEP::HepMatrix J=preJ.sub(1,S,1,S);
    CLHEP::HepMatrix modified_J(S,S);
    CLHEP::HepMatrix a(S,S),b(S,S);
    CLHEP::HepVector ddr(J.num_row()),ddi(J.num_row()),NN(J.num_row());
    eigen(J,ddr,ddi,b,a); // b and a are orthogonal
    for(int u=S;u-->0;){ // ... but need to be made orthonormal
      double sum=0;
      for(int i=S;i-->0;){
	sum+=b[u][i]*a[u][i];
      }
      for(int i=S;i-->0;){
	b[u][i]/=sum;
      }
    }
    
    vector<int> eval_order(S),rank(S);
    for(int i=S;i-->0;){
      eval_order[i]=i;
      rank[i]=i;
      NN[i]=qna.juvenile[i];
    }
    stable_sort(eval_order.begin(),eval_order.end(),ddr_larger(ddr));
    stable_sort(rank.begin(),rank.end(),ddr_larger(NN));
    
    ofstream singular("singular.dat");
    for(double p=0;p>-4.1;p-=0.5){
      cout << p << " ";
      cout.flush();
      for(int i=S;i-->0;){
	for(int j=S;j-->0;){
	  modified_J[i][j]=J[i][j]*pow(NN[j],p);
	}
      }
      NewEigen(modified_J,ddr,ddi);
      stable_sort(eval_order.begin(),eval_order.end(),ddr_larger(ddr));
      
      singular << -p;
      for(int uu=0;uu<n_values;uu++){
	int u=eval_order[uu];
	int n=rank[S-1-uu];
	singular << " " << -ddr[u];
      }
      singular << endl;
    }
    cout << endl;
  }

  void eigenvalue_list(const char* filename){
    const int n=number_of_variables();
    CLHEP::HepMatrix J(n,n);
    ODE_vector state(n),ddt(n);
    
    write_state_to(state);
    
    Jacobian(state,ddt,J);
    
    if(_is_qna){
      CLHEP::HepMatrix K=J.sub(1,n/2,1,n/2);
      eigen_report(K,filename);
    }else{
      eigen_report(J,filename);
    }
  }
  void eigenvalue_list_of_a(const char* filename){
    ALWAYS_ASSERT(_is_qna);
    const int S=_S();
    CLHEP::HepMatrix A(S,S,0);
    for(int i=S;i-->0;){
      for(int j=S;j-->0;){
	A[i][j]=-(_juvenile_a[i][j]);
      }
    }
    
    eigen_report(A,filename);
  }
private:
  void eigen_report(CLHEP::HepMatrix & J,const char * filename){
    std::ofstream os(filename);

    cout << "Reporting for " << filename << endl;
    
    ALWAYS_ASSERT(J.num_row()==J.num_col());
    CLHEP::HepVector ddr(J.num_row()),ddi(J.num_row());
    NewEigen(J,ddr,ddi);
    
    int npos=0;
    int old_prec=std::cout.precision();
    os.precision(4);  
    for(int i=0;i<J.num_row();i++){
      os << -ddr[i] << " " << ddi[i] << std::endl;
      if(ddr[i]>0)
	npos++;
    }
    std::cout.precision(old_prec);
    if(npos){
      WARNING(npos <<  (npos==1 ? 
			" eigenvalue has positive real part" : 
			" eigenvalues have positive real part" ));
    }
    
    int S=J.num_row();
    double rho2=0;
    double median_r2=ddr[S/2]*ddr[S/2]+ddi[S/2]*ddi[S/2];
    double median_r=sqrt(median_r2);
    double r2sum=0;
    double Re_sum=0;
    average_meter mean_r,mean_r2;
    for(int i=S;i-->0;){
      double r2=ddr[i]*ddr[i]+ddi[i]*ddi[i];
      if(r2>rho2){
	rho2=r2;
      }
      mean_r2.sample(r2);
      r2sum+=r2;
      mean_r.sample(sqrt(r2));
      Re_sum+=ddr[i];
    }

    double rho=sqrt(rho2);
    REPORT(median_r);
    REPORT(mean_r);
    REPORT(sqrt(mean_r2.readout()));
    REPORT(rho);
    REPORT(sqrt(r2sum));
    REPORT(sqrt(r2sum-rho2));
    REPORT(Re_sum);
    return;
  }
public:
  void set_fishing_policy(int target,double pressure){
    int S=_S();
    if(_is_qna){
      pressure*=(target<S?
		 adev_juvenile[target]:
		 adev_adult[target-S]);
      if(target>=S) target-=S;
    }
    _fishing_target=target;
    _fishing_pressure=pressure;
  }
  double get_fishing_pressure(){
    return _fishing_pressure;
  }
  void scale_population(int i,double x){
    if(!_is_qna) adult[i]*=x;
    juvenile[i]*=x;
  }
  double population(int i) const {
    if(is_active(i)){
      if(_is_qna){
	return juvenile[i];
      }else{
	return juvenile[i]+adult[i];
      }
    }else{
      return 0;
    }
  }
  void report_parameters(){
    REPORT(_birth);
    REPORT(_juvenile_interaction);
    REPORT(_adult_interaction);
    REPORT(_juvenile_i);
    REPORT(_adult_i);
    REPORT(_juvenile_C);
    REPORT(_adult_C);
    REPORT(_maturation_rate);
    REPORT(_net_juvenile_population_growth);
    REPORT(_adult_death);
    REPORT(_these_are_all_population_dynamical_parameters);
    REPORT(_is_qna);
    REPORT(_fishing_pressure);
    REPORT(_fishing_target);
  }
  double msy(const int target_species,double standard_relaxation_time){
    CompetiWeb_t dummy;
    return msy(target_species,standard_relaxation_time,dummy);
  }
  double msy(const int target_species,double standard_relaxation_time,
	     CompetiWeb_t & saved_web,
	     const char * report_file_name=0,
	     const double initial_pressure=1){
    int S=_S();

    std::ofstream report;
    if(report_file_name){
      report.open(report_file_name);
    }

    delete_by_deactivating();
    relax(standard_relaxation_time);

    const double pressure_increment=1.05;
    double pressure=initial_pressure;
    double msy,first_pop,last_pop,pop_ratio;
      
    first_pop=population(target_species);
    if(report_file_name){
      report << 0 << " " << first_pop << endl;
    }

    bool was_stable=false;
    do{
      try{
	was_stable=true;
	CompetiWeb_t tester=*this;
	tester.set_fishing_policy(target_species+S,pressure);
	tester.relax(standard_relaxation_time/10);
      }catch(int){
	was_stable=false;
	pressure/=3;
	cout << "Restarting msy search with lower pressure." << endl;
      }
    }while(!was_stable);

    pressure/=pressure_increment;
    set_fishing_policy(target_species+S,pressure);
    relax(standard_relaxation_time);

    try{
      while(true){
	last_pop=population(target_species);
	set_fishing_policy(target_species+S,pressure);
	REPORT(pressure);
	saved_web=*this;
	relax(standard_relaxation_time/10);
	if(report_file_name){
	  report << pressure << " " 
		 << population(target_species) << endl;
	}
	pressure*=pressure_increment;
      }
    }catch(int){
      msy=pressure/pressure_increment;
    }
    pop_ratio=last_pop/first_pop;

    const double effective_growth_rate=2*msy/last_pop;
    const double effective_K=first_pop;

    REPORT(msy);
    REPORT(pop_ratio);
    REPORT(effective_growth_rate);
    REPORT(effective_K);
    REPORT(effective_growth_rate/effective_K);

    return msy;
  }
  double total_abundance_in(const CompetiWeb_t & other,int i){
    return adev_adult[i]*other.adult[i]+adev_juvenile[i]*other.juvenile[i];
  }
  double total_abundance_in(const CompetiWeb_t & other){
    int S=_S();
    double abundace_sum=0;
    ALWAYS_ASSERT(_is_qna && !other._is_qna);
    for(int i=S;i-->0;){
      abundace_sum+=total_abundance_in(other,i);
    }
    return abundace_sum;
  }
  double predict_species_richness_S(){
    delete_all_inactive_species();
    int S=_S();
    CompetiWeb_t qna=quasi_neutral_approximation(true);
    cout << "********************************************" << endl;
    ALWAYS_ASSERT(_adult_C==_juvenile_C);
    average_meter G_ij;
    average_meter CJ,CA,IJ,IA;
    for(int i=S;i-->0;){
      //cout << "";// suppress a compiler optimization bug
      for(int j=S;j-->0;){
	if(i!=j){ 
	  G_ij.sample(qna._juvenile_a[i][j]);
	  CJ.sample(_juvenile_a[i][j]>0);
	  CA.sample(_adult_a[i][j]>0);
	  if(_juvenile_a[i][j]>0)
	    IJ.sample(_juvenile_a[i][j]/_juvenile_a[0][0]);
	  if(_adult_a[i][j]>0)
	    IA.sample(_adult_a[i][j]/_adult_a[0][0]);
	}
      }
    }
    REPORT(_adult_C);
    REPORT(CJ);
    REPORT(CA);
    REPORT(_adult_i);
    REPORT(IJ);
    REPORT(IA);
    double C=_adult_C;
    //    C=(CJ.readout()+CA.readout())/2;
    REPORT(C*S);
    REPORT(CJ*S);
    REPORT(CA*S);
    ALWAYS_ASSERT(_adult_i==_juvenile_i);
    double I=_adult_i;
    //    I=(IA.readout()+IJ.readout())/2;
    REPORT(I);
    double G_MF=
      juvenile_weight*_juvenile_interaction+
      adult_weight*_adult_interaction;
    REPORT(juvenile_weight);
    REPORT(adult_weight);
    REPORT(_juvenile_interaction);
    REPORT(_adult_interaction);
    REPORT(G_MF);
    REPORT(qna._juvenile_a[0][0]);
    double var_G_ij=C*(4-3*C)*SQR(I)*
      (SQR(juvenile_weight*_juvenile_interaction)+
       SQR(adult_weight*_adult_interaction) )
      /3;
    REPORT(var_G_ij);
    REPORT(G_ij.var());
    double mean_G_ij=C*I*G_MF;
    REPORT(mean_G_ij);
    REPORT(G_ij.readout());
    double tau=1-1/(4-3*C);
    REPORT(tau);
    REPORT(1-(tau*var_G_ij+SQR(I*G_MF))/SQR(G_MF));
    double aspect_ratio=(1-tau)/(1+tau);
    REPORT(G_MF);
    REPORT(G_MF*aspect_ratio);
    double S_pred=SQR(G_MF)/(SQR(1+tau)*var_G_ij);
    double G_MF_CORRECTED=G_MF*(1-(1+(S-1.0)*C*I)/S);
    REPORT(G_MF_CORRECTED);
    double S_corrected=SQR(G_MF_CORRECTED)/(SQR(1+tau)*var_G_ij);
//     double xmin=G_MF-(1+tau)*sqrt(var_G_ij*S);
//     double xmax=G_MF+(1+tau)*sqrt(var_G_ij*S);
//     double ymax=(1-tau)*sqrt(var_G_ij*S);
//     REPORT(xmin);
//     REPORT(xmax);
//     cout << "Y=" << ymax << 
//       "*SQRT(1-(X-" << G_MF << ")^2/" <<  (1+tau)*(1+tau)*var_G_ij*S << ")" 
// 	 << endl;
    REPORT((1+tau)*sqrt(var_G_ij*S));
    REPORT((1-tau)*sqrt(var_G_ij*S));
    REPORT(S_pred);
    REPORT(S_pred*SQR(1-I*C));
    REPORT(S_corrected);
    REPORT(S);
    return S_pred;
  }

  void experimental(){
    int S=1000;
    double C=0.005;
    double I=1/sqrt(S*C);
    ALWAYS_ASSERT(I<1);
    CLHEP::HepMatrix M(S,S);
    for(int i=S;i-->0;){
      for(int j=S;j-->0;){
	M[i][j]=(unirand()<C ? I : 0);
      }
    }
    M=M+CLHEP::HepDiagMatrix(S,1);
    int ierr;
    M.invert(ierr);
    CLHEP::HepVector positive(S);
    for(int i=S;i-->0;){
      positive[i]=unirand();
    }
    positive=M*positive;
    ofstream pos("positive.dat");
    average_meter p;
    for(int i=S;i-->0;){
      pos << positive[i] << endl;
      p.sample(positive[i]);
    }
    REPORT(p);
    REPORT(p.var());
  }

  void logistic(){
    // get a logistic equation for each species:
    ALWAYS_ASSERT(_is_qna);
    delete_all_inactive_species();
    const int S=_S();
    REPORT(S);

    // Populate LV representation using CLHEP data types.
    CLHEP::HepMatrix a(S,S);
    CLHEP::HepVector r(S);
    for(int i=_S();i-->0;){
      for(int j=_S();j-->0;){
	a[i][j]=_juvenile_a[i][j];
      }
      switch(_is_qna){
      case local:
	r[i]=qna_growth[i];
	break;
      case mean_field:
	r[i]=_net_juvenile_population_growth;
	break;
      default:
	FATAL_ERROR("unknown qna type");
      }
    }
    // compute a logistic representation for each species:
    
    int ierr;
    CLHEP::HepMatrix ia=a.inverse(ierr);
    if(ierr){
      FATAL_ERROR("Matrix inversion failed");
    }
    CLHEP::HepVector N_star=ia*r;
    ofstream lo("logistic.dat");
    REPORT(_juvenile_a[1][1]);
    for(int i=0;i<S;i++){
      lo << N_star[i] << " "
	 << N_star[i]/ia[i][i] << endl;
    }
  };

  void print_a(){
    const int S=_S();
    cout << "AAAAAAAAAAAA" << endl;
    for(int i=0;i<S;i++){
      if(is_active(i)){
	for(int j=0;j<S;j++){
	  if(is_active(j)){
	    cout << _juvenile_a[i][j] << "\t";
	  }
	}
	cout << endl;
      }
    }
    cout << "AAAAAAAAAAAA" << endl;
  }

  friend class accuracy1_object;
  friend class decomposition_object;
};

class decomposition_object: public CompetiWeb_t{
  int n;
  typedef complex<double> cmplx;
  CLHEP::HepVector ddr,ddi;
  vector< vector<cmplx> > eigenvectors;
  vector< vector<cmplx> > adeigenvectors;
  vector<int> eval_order;
  ODE_vector fixed_point;
  std::ofstream modes,damped_modes,principal_dynamics;
  class ddr_larger {
    const CLHEP::HepVector &_ddr;
  public:
    ddr_larger(CLHEP::HepVector &dr):_ddr(dr){};
    bool operator()(int i,int j){return _ddr[i]>_ddr[j];}
  };
public:
  decomposition_object(CompetiWeb_t & o1):
    CompetiWeb_t(o1),
    n(number_of_variables()),
    ddr(n), ddi(n),
    eigenvectors(n,vector<cmplx>(n)),
    adeigenvectors(n,vector<cmplx>(n)),
    eval_order(n),
    fixed_point(n),
    modes("modes.dat"),
    damped_modes("damped_modes.dat"),
    principal_dynamics("principal_dynamics.dat")
  {
    CLHEP::HepMatrix J(n,n);
    ODE_vector ddt(n);
    write_state_to(fixed_point);
    Jacobian(fixed_point,ddt,J);
    for(int i=n;i-->0;){
      for(int j=n;j-->0;){
	J[i][j]*=exp(fixed_point[j]-fixed_point[i]);
      }
    }
    CLHEP::HepMatrix eigenv(n,n);
    CLHEP::HepMatrix adeigenv(n,n);
    eigen(J,ddr,ddi,adeigenv,eigenv); 
//     for(int u=0;u<n;u++){
//       cout << ddr[u] << " + I* " << ddi[u] << endl;
//     }
      
    for(int i=0;i<n;i++){
      for(int u=0;u<n;u++){
	if(ddi[u]==0){
	  eigenvectors[u][i]=eigenv[u][i];
	  adeigenvectors[u][i]=adeigenv[u][i];
	}else{//undo lapack convention:
	  eigenvectors[u][i]=cmplx(eigenv[u][i],eigenv[u+1][i]);
	  eigenvectors[u+1][i]=cmplx(eigenv[u][i],-eigenv[u+1][i]);
	  adeigenvectors[u][i]=
	    cmplx(adeigenv[u][i],-adeigenv[u+1][i]);
	  adeigenvectors[u+1][i]=
	    cmplx(adeigenv[u][i],+adeigenv[u+1][i]);
	  u++;
	}
      }
    }
    
    for(int u=n;u-->0;){
      cmplx sum=0;
      for(int i=n;i-->0;){
	sum+=adeigenvectors[u][i]*eigenvectors[u][i];
      }
      double sum2=0;
      for(int i=n;i-->0;){
	adeigenvectors[u][i]/=sum;
	sum2+=(eigenvectors[u][i]*conj(eigenvectors[u][i])).real();
      }
      cout << sum2 << endl;
//       sum=0;// test orthogonality
//       if(u>=2){
// 	for(int i=n;i-->0;){
// 	  sum+=adeigenvectors[u-1][i]*eigenvectors[u][i];
// 	}
//       }
//       if(ddi[u]) cout << "* ";
//       REPORT(sum);
    }
    for(int i=n;i-->0;){
      eval_order[i]=i;
    }
    stable_sort(eval_order.begin(),eval_order.end(),ddr_larger(ddr));
    for(int uu=0;uu < min<int>(4,n);uu++){
      int u=eval_order[uu];
      REPORT(ddr[u]);
    }
    ofstream mode_list("mode_list.dat");
    for(int uu=0;uu < n;uu++){
      int u=eval_order[uu];
      mode_list << uu << " " << ddr[u] << endl;
    }
    current_time=0;
//     ODE_vector state(n);  // test if these are really linear modes:
//     write_state_to(state);
//     for(int i=n;i-->0;){
//       state[i]+=0.01*eigenvectors[eval_order[1]][i];
//     }
//     read_state_from(state);
  }
  virtual void record_for_steady_state(const ODE_vector & state,
				       const ODE_vector & ddt){
    modes << current_time;
    vector<cmplx> amplitude(n);
    for(int uu=0;uu<n;uu++){
      int u=eval_order[uu];
      cmplx sum=0;
      for(int i=n;i-->0;){
	sum+=adeigenvectors[u][i]*(state[i]-fixed_point[i]);
      }
      amplitude[u]=sum;
      modes << " " << sqrt(sum.real()*sum.real()+sum.imag()*sum.imag());
    }
    modes << endl;

    const double cutoff_time=0.001;
    damped_modes << current_time;
    for(int uu=0;uu<n;uu++){
      int u=eval_order[uu];
      double x=cutoff_time*ddr[u];
      //amplitude[u]*=exp(-x*x/2);
      amplitude[u]*=exp(-ddr[u]*current_time);
      damped_modes << " " <<
	1/sqrt((amplitude[u]*conj(amplitude[u])).real());
    }
    damped_modes << endl;
    
    int S=number_of_species();
    principal_dynamics << current_time;
    for(int i=0;i<S;i++){
      double sum_J=fixed_point[i];
      double sum_A=fixed_point[i+S];
      for(int uu=0;uu<n;uu++){
	int u=eval_order[uu];
	sum_J+=(amplitude[u]*eigenvectors[u][i]).real();
	sum_A+=(amplitude[u]*eigenvectors[u][i+S]).real();
      }
      principal_dynamics << " " << exp(sum_J)+exp(sum_A) ;
    }
    principal_dynamics << endl;
  }
};

class accuracy1_object: public combined_relaxing_dynamical_object{
  CompetiWeb_t *_o1,*_o2;
  vector<double> max_diff,max_diff_corrected;
  std::ofstream _os,_ddtos,_rvsum,_sample_dynamics,_excitations;
  vector<double> _initial_pop_diff;
  int _suggested_target;
  vector<bool> _is_sample;
public:
  accuracy1_object(CompetiWeb_t & o1,
		   CompetiWeb_t & o2,
		   int suggested_target=-1):
    combined_relaxing_dynamical_object(&o1,&o2),
    _o1(&o1),_o2(&o2),max_diff(_o1->_S()),max_diff_corrected(_o1->_S()),
    _os("deviations.dat"),_ddtos("ddt_errors.dat"),
    _rvsum("rvsum.dat"),_sample_dynamics("sample_dynamics.dat"),
    _excitations("excitations.dat"),
    _initial_pop_diff(_o1->_S()),
    _suggested_target(suggested_target),
    _is_sample(_o1->_S(),false)
  {
    int S=_o1->_S();
    ALWAYS_ASSERT(o1._S()==o2._S());
    ALWAYS_ASSERT(!_o1->_is_qna);
    ALWAYS_ASSERT(_o2->_is_qna);
    for(int i=S;i-->0;){
      double o1_pop=
	//_o1->population(i);
	_o1->juvenile[i]*_o2->adev_juvenile[i]+_o1->adult[i]*_o2->adev_adult[i];
      _initial_pop_diff[i]=o1_pop-_o2->population(i);
    }
    set_random_seed(sampling_random_seed);
    if(suggested_target==-1){
      int n_samples=std::min(S,10);
      while(n_samples){
	int i=random_integer(S);
	if(not _is_sample[i]){
	  _is_sample[i]=true;
	  n_samples--;
	}
      }
    }else{
      _is_sample[suggested_target]=true;
    }
  }
  virtual void record_for_steady_state(const ODE_vector & state,
				       const ODE_vector & ddt){
    int S=_o1->_S();
    double d2sum=0,d2sum_corrected=0,max_d=0;
    double popsum=0;
    _sample_dynamics << current_time;
    _excitations << current_time;
    for(int i=S;i-->0;){
      double o1_pop=
	//_o1->population(i);
	_o1->juvenile[i]*_o2->adev_juvenile[i]+_o1->adult[i]*_o2->adev_adult[i];
      popsum+=o1_pop;
      double diff=fabs(o1_pop-_o2->population(i));
      d2sum+=diff*diff;
      max_diff[i]=
	max<double>(max_diff[i],diff);
      double diff_corrected=
	fabs(o1_pop-_o2->population(i)-_initial_pop_diff[i]);
      d2sum_corrected+=diff_corrected*diff_corrected;
      max_diff_corrected[i]=
	max<double>(max_diff_corrected[i],diff_corrected);
      max_d=max<double>(max_d,diff);
      if(_is_sample[i]){
	_sample_dynamics << " " << o1_pop << " " << _o2->population(i);
	_excitations << " " <<
	  _o1->juvenile[i]*_o2->ev_adult[i]-_o1->adult[i]*_o2->ev_juvenile[i];
      }
    }
    _sample_dynamics << endl;
    _excitations << endl;

    double vector_diff2_sum=0;
    double vector2_sum=0;
    for(int i=S;i-->0;){
      double ddt_web=ddt[i]*_o1->juvenile[i]*_o2->adev_juvenile[i]+ddt[i+S]*_o1->adult[i]*_o2->adev_adult[i];
      //double ddt_web=ddt[i]*_o1->juvenile[i]+ddt[i+S]*_o1->adult[i];
      double ddt_qna=ddt[i+_o1->number_of_variables()]*_o2->juvenile[i];
      double diff=(ddt_web-ddt_qna);///ddt_web;
      vector_diff2_sum+=diff*diff;
      vector2_sum+=ddt_web*ddt_web;
    }

    _os << current_time 
	<< " " << sqrt(d2sum/S)/(popsum/S) 
	<< " " << max_d 
	<< endl;
    _ddtos << current_time 
	   << " " << sqrt(vector_diff2_sum/vector2_sum) 
	   << " " << sqrt(vector2_sum/S) 
	   << endl;
    _rvsum << current_time 
	   << " " << popsum
	   << endl;
  }
  double max_deviation(int i){
    return max_diff[i];
  }
  double max_deviation_corrected(int i){
    return max_diff_corrected[i];
  }
  double rms_max_diff(){
    int S=_o1->_S();
    double sum2=0;
    for(int i=S;i-->0;){
      sum2+=max_diff[i]*max_diff[i];
    }
    return sqrt(sum2/S);
  }
  double rms_relative_max_diff(){
    int S=_o1->_S();
    double sum2=0;
    for(int i=S;i-->0;){
      double x=max_diff[i]/_o1->population(i);
      sum2+=x*x;
    }
    return sqrt(sum2/S);
  }
  double rms_max_diff_corrected(){
    int S=_o1->_S();
    double sum2=0;
    for(int i=S;i-->0;){
      sum2+=max_diff_corrected[i]*max_diff_corrected[i];
    }
    return sqrt(sum2/S);
  }
};

void read(CompetiWeb_t & web,int argc,char ** argv){
  std::string input_file_name;
  if(argc>2){
    input_file_name=argv[2];
  }else{
    input_file_name="CompetiWeb.dat";
  }
  std::ifstream ifs(input_file_name.c_str());
  boost::archive::text_iarchive ia(ifs);
  ia >> web;
}

///////////////////////////////////////////////////////////////  
static int use_true_qna=true;
static int use_local_qna=false;
static double standard_relaxation_time=1e5;
static int random_seed=-1;

void read_arguments(int &argc,char **&argv){
  char c;
  const char * options="hvblmr:";

  while(1){
    c=getopt(argc, argv, options);
    if (c == -1)
      break;

    switch(c){
    case 'v':
      printf("%s\n",version.c_str());
      exit(0);
      break;
    case 'b':
      use_true_qna=false;
      break;
    case 'l':
      use_local_qna=true;
      break;
    case 'm':
      use_local_qna=false;
      break;
    case 'r':
      random_seed=atoi(optarg);
      break;
    case 'h':
    default:
      fprintf(stdout,"\nusage: %s -%s\n",argv[0],options);
        exit(c=='h'?0:1);
      }
  }
  argc-=optind-1;
  argv+=optind-1;

}


int main(int argc,char** argv){
  set_cfg_parameter("vector_accuracy","1e-5");
  set_cfg_parameter("vector_truncation_epsilon","1e-25");
  set_cfg_parameter("DEFAULT_ABSOLUTE_TOLERANCE","1e-6");

  read_arguments(argc,argv);

  if(random_seed<0){
    struct timeval t;
    if(gettimeofday(&t,0)) FATAL_ERROR("problem with system clock");
    set_random_seed(t.tv_sec+t.tv_usec);
  }else{
    set_random_seed(random_seed);
  }      
  sampling_random_seed=random_integer(1000);
  perturbing_random_seed=random_integer(1000);

  average_meter mean_size;
  bool skip_averaging=true;
  int previous_size=0;
  int max_size=0;
  int focal_species=0; 

  const double I=0.1,C=0.5;
  const double Cnew=0.6;

  const double CONST=SQR(1+1-1.0/(4-3*C))*C*(4.0/3.0-C)*SQR(I);
  
  const double Inew=
    sqrt(CONST)/sqrt(SQR(1+1-1.0/(4-3*Cnew))*Cnew*(4.0/3.0-Cnew));

  CompetiWeb_t web(6 /* birth */, 
		   1 /* maturation */,
		   1.0 /*juvenile_interaction*/, 1.0 /*adult_interaction*/, 
		   I /*juvenile_i*/, I /*adult_i*/, 
		   C /*juvenile_C*/, C /* adult_C*/, 
		   1 /*juvenile_death*/, 1 /*adult_death*/, 
		   1 /*these_are_all_population_dynamical_parameters*/);
  
  std::string activity;
  if(argc<=1){
    activity="run";
  }else{
    activity=argv[1];
  }
  
  if(activity=="run"){

    if(argc>2) read(web,argc,argv);

    ofstream richness("rich.dat"),focal_species_abundance("focal.dat");
    CompetiWeb_t weight_web=web;
    weight_web.add_species();
    weight_web=weight_web.quasi_neutral_approximation();

    web.add_inactive_species_if_required();

    while(step_count<=500000){
      REPORT(step_count);
      int new_species=
	web.add_fit_species();
      try{
	web.delete_all_inactive_species();
	web.add_inactive_species_if_required();
	web.relax(standard_relaxation_time);
	step_count++;
	richness << step_count << " " 
		 << web.number_of_species() << endl;
	max_size=max<int>(max_size,web.number_of_species());

	if(web.number_of_species() < focal_species || 
	   !web.is_active(focal_species)){
	  focal_species=new_species;
	}
	focal_species_abundance 
	  << step_count << " " 
	  << weight_web.total_abundance_in(web,focal_species) << endl;
      
	if(skip_averaging){
	  if(web.number_of_species()<previous_size)
	    skip_averaging=false;
	  previous_size=web.number_of_species();
	}
	else
	  mean_size.sample(web.number_of_species());
	REPORT(mean_size);
	REPORT(web.number_of_species());
	if(step_count%20==0){//save eigenvectors:
	}
	web.evaluate();
	if(step_count<=100 || 
	   0==(step_count%int(0.5+pow(10,-1+int(log10(double(step_count)))))) ){
	  // save data to archive
	  std::stringstream output_file_name(stringstream::out);
	  output_file_name << "CompetiWeb" << step_count << ".dat";
	  std::ofstream ofs(output_file_name.str().c_str());
	  boost::archive::text_oarchive oa(ofs);
	  const CompetiWeb_t cweb(web);
	  oa << cweb;

	  web.invasion_fitness_distribution();
	}
	if(step_count>50*max_size){
	  exit(0);
	}
      }catch(int){
	web.delete_species(new_species);
      }

    }
  } else if(activity=="qna"){

    read(web,argc,argv);
    
    web.delete_all_inactive_species();
    CompetiWeb_t qna=web.quasi_neutral_approximation(use_true_qna);
    const CompetiWeb_t start_qna=qna;
      
    REPORT(web.juvenile_population_sum()+web.adult_population_sum());
    REPORT(qna.juvenile_population_sum());
    REPORT(qna.number_of_species());
      
    qna.delete_by_deactivating();
    qna.relax(standard_relaxation_time);
      
    REPORT(qna.juvenile_population_sum());
    REPORT(qna.number_of_species());
      
    REPORT(qna.abundance_variation(start_qna));
      
    if(0){
      cout << "***** perturbation experiment ****" << endl;
	
      int most_abundant=web.most_abundant_species();
      //while(!web.is_active(most_abundant)) most_abundant++;
      web.delete_by_deactivating();
      web.delete_species(most_abundant);
      cout << "**** relaxing perturbed original:" << endl;
      web.relax(standard_relaxation_time);
	
      cout << "**** make qna3 thereof:" << endl;
      const CompetiWeb_t qna3=web.quasi_neutral_approximation(use_true_qna);
	
      cout << "**** post-relaxing this to qna2:" << endl;
      CompetiWeb_t qna2=qna3;
      qna2.delete_by_deactivating();
      qna2.relax(standard_relaxation_time);
	
      cout << "**** relaxing perturbed qna:" << endl;
      qna.delete_species(most_abundant);
      qna.relax(standard_relaxation_time);
	
      REPORT(qna.abundance_variation(qna2));
      REPORT(qna.abundance_variation(qna3));
      REPORT(qna2.abundance_variation(qna3));
      REPORT(qna.abundance_variation(start_qna));
    }
    exit(0);
  }
  else if(activity=="accuracy1" || activity=="accuracy2"){
    
    //set_cfg_parameter("default_absolute_tolerance","1e-8");
    

    read(web,argc,argv);
    
    web.delete_all_inactive_species();

    CompetiWeb_t qna=
      (activity=="accuracy1" ?
       web.quasi_neutral_approximation(use_true_qna) :
       web.local_quasi_neutral_approximation(use_true_qna) );
    web.delete_by_deactivating();
    qna.delete_by_deactivating();
    const CompetiWeb_t prerelaxed_qna=qna;
//     if(activity=="accuracy1"){
//       qna.relax(standard_relaxation_time);
//       cout << "finished relaxation of qna" << endl;
//     }
    const CompetiWeb_t start_qna=qna;

    accuracy1_object a1o(web,qna);
    
    int number_of_species_to_reduce=int(web.number_of_species()*1.0/3.0);
    vector<bool> picked(web.number_of_species(),false);
    double reduction=1.0/3.0;
    set_random_seed(perturbing_random_seed);
    while(number_of_species_to_reduce>0){
      int s=unirand()*web.number_of_species();
      if(web.is_active(s) and not picked[s]){
	web.scale_population(s,reduction);
	qna.scale_population(s,reduction);
	picked[s]=true;
	number_of_species_to_reduce--;
      }
    }

    std::cout << "starting relaxation" << endl;
    a1o.relax(standard_relaxation_time);
    //web.relax(standard_relaxation_time);

    const CompetiWeb_t post_qna=
      (activity=="accuracy1" ?
       web.quasi_neutral_approximation(use_true_qna) :
       web.local_quasi_neutral_approximation(use_true_qna) );
    REPORT(qna.abundance_variation(start_qna));
    REPORT(qna.abundance_variation(post_qna));
    const double mean_abundance=
      (web.juvenile_population_sum()+web.adult_population_sum())/
      web.number_of_species();
    REPORT(start_qna.abundance_variation(prerelaxed_qna)/mean_abundance);
    REPORT(a1o.rms_max_diff());
    REPORT(a1o.rms_max_diff()/mean_abundance);
    if(false && activity=="accuracy1"){
      REPORT(a1o.rms_max_diff_corrected()/mean_abundance);
    }
    REPORT(a1o.rms_relative_max_diff());
    {
      std::stringstream output_file_name(stringstream::out);
      output_file_name << "relaxed_web.dat";
      std::ofstream ofs(output_file_name.str().c_str());
      boost::archive::text_oarchive oa(ofs);
      const CompetiWeb_t cweb(web);
      oa << cweb;
    }
  }
  else if(activity == "neutral"){

    read(web,argc,argv);
    web.delete_all_inactive_species();
    int S=web.number_of_species();

    CompetiWeb_t qna=
      (use_local_qna ?
       web.local_quasi_neutral_approximation(use_true_qna) :
       web.quasi_neutral_approximation(use_true_qna) );
       
    web.delete_by_deactivating();
    qna.delete_by_deactivating();

    const int target_species=web.nth_abundant_species(S/4+1);
    const double start_abundance=web.population(target_species);
    const double qna_msy=
      CompetiWeb_t(qna).msy(target_species,standard_relaxation_time);
    REPORT(qna_msy);
    CompetiWeb_t saved_web(web);
    const double msy=
      CompetiWeb_t(web).
      msy(target_species,standard_relaxation_time,saved_web,0,
	  qna_msy/pow(1.05,4));

    set_cfg_parameter("default_absolute_tolerance","1e-8");

    const double pressure=0.9*msy;
    
    if(!use_local_qna){
      qna.relax(standard_relaxation_time);
    }

    CompetiWeb_t saved_qna=qna; saved_web=web;
    const double web_start_abundance=qna.total_abundance_in(web);
    const double qna_start_abundance=qna.juvenile_population_sum();
    const double web_start_abundance1=
      qna.total_abundance_in(web,target_species);
    const double qna_start_abundance1=qna.population(target_species);

    //accuracy1_object a1o(web,qna,target_species);
    accuracy1_object a1o(web,qna);

    web.set_fishing_policy(S+target_species,pressure);
    qna.set_fishing_policy(S+target_species,pressure);

    web.relax(standard_relaxation_time);
    qna.relax(standard_relaxation_time);

    const double web_final_abundance=qna.total_abundance_in(web);
    const double qna_final_abundance=
      qna.juvenile_population_sum();
    const double web_final_abundance1=
      qna.total_abundance_in(web,target_species);
    const double qna_final_abundance1=
      qna.population(target_species);

    const double web_neutrality=1-
      (web_start_abundance-web_final_abundance)/
      (web_start_abundance1-web_final_abundance1);
    const double qna_neutrality=1-
      (qna_start_abundance-qna_final_abundance)/
      (qna_start_abundance1-qna_final_abundance1);

    web.set_fishing_policy(0,0);
    qna.set_fishing_policy(0,0);
    
    std::cout << "starting relaxation" << endl;
    a1o.relax(standard_relaxation_time);
    REPORT(target_species);
    REPORT(a1o.rms_max_diff());
    const double mean_abundance=
      (web.juvenile_population_sum()+web.adult_population_sum())/
      web.number_of_species();
    REPORT(a1o.rms_max_diff()/mean_abundance);
    REPORT(fabs(a1o.max_deviation(target_species)));
    REPORT(start_abundance);
    REPORT(fabs(a1o.max_deviation(target_species))/start_abundance);
    if(!use_local_qna){
      REPORT(a1o.max_deviation_corrected(target_species)/mean_abundance);
      REPORT(fabs(a1o.max_deviation_corrected(target_species))/start_abundance);
    }
    
    REPORT(web_start_abundance);
    REPORT(web_final_abundance);
    REPORT(web_start_abundance1);
    REPORT(web_final_abundance1);
    REPORT(web_neutrality);
    REPORT(qna_start_abundance);
    REPORT(qna_final_abundance);
    REPORT(qna_start_abundance1);
    REPORT(qna_final_abundance1);
    REPORT(qna_neutrality);
  }
  else if(activity=="msy"){
    
    read(web,argc,argv);
    
    web.delete_all_inactive_species();

    int S=web.number_of_species();

    CompetiWeb_t qna=web.local_quasi_neutral_approximation(use_true_qna);
    qna.set_fishing_policy(0,0);
    web.delete_by_deactivating();
    qna.delete_by_deactivating();
    qna.relax(standard_relaxation_time);
    cout << "finished relaxation of qna" << endl;

    const int target_species=S*unirand();//web.most_abundant_species();

    const double web_start_abundance=qna.total_abundance_in(web);
    const double qna_start_abundance=qna.juvenile_population_sum();
    const double web_start_abundance1=
      qna.total_abundance_in(web,target_species);
    const double qna_start_abundance1=qna.population(target_species);

    CompetiWeb_t saved_web(web),saved_qna(qna);

    const double qna_msy=
      qna.msy(target_species,standard_relaxation_time,saved_qna,"response.dat");
    const double web_msy=
      web.msy(target_species,standard_relaxation_time,saved_web,
	      "qna_respose.dat",
	      0.8*qna_msy);

    const double web_final_abundance=saved_qna.total_abundance_in(saved_web);
    const double qna_final_abundance=
      saved_qna.juvenile_population_sum();
    const double web_final_abundance1=
      saved_qna.total_abundance_in(saved_web,target_species);
    const double qna_final_abundance1=
      saved_qna.population(target_species);

    const double web_neutrality=1-
      (web_start_abundance-web_final_abundance)/
      (web_start_abundance1-web_final_abundance1);
    const double qna_neutrality=1-
      (qna_start_abundance-qna_final_abundance)/
      (qna_start_abundance1-qna_final_abundance1);
    
    REPORT(web_start_abundance);
    REPORT(web_final_abundance);
    REPORT(web_start_abundance1);
    REPORT(web_final_abundance1);
    REPORT(web_neutrality);
    REPORT(qna_start_abundance);
    REPORT(qna_final_abundance);
    REPORT(qna_start_abundance1);
    REPORT(qna_final_abundance1);
    REPORT(qna_neutrality);
    REPORT(qna_msy/web_msy);
    exit(0);
  }
  else if(activity == "evs"){

    read(web,argc,argv);
    
    web.relax(10);
    web.evaluate();
    CompetiWeb_t web2=web;
    web2.delete_all_inactive_species();
    web2.eigenvalue_list("evs.dat");
    CompetiWeb_t qna=web2.local_quasi_neutral_approximation(use_true_qna);
    qna.eigenvalue_list("local_qna_evs.dat");
    qna.eigenvalue_list_of_a("local_a_evs.dat");
    qna=web2.quasi_neutral_approximation(use_true_qna);
    qna.eigenvalue_list("qna_prerelax_evs.dat");
    CompetiWeb_t qna2(qna);
    fixed_point_analyzer fpa(&qna2);
    fpa.snap_to_fixed_point();
    qna2.eigenvalue_list("qna_fixedpoint_evs.dat");
    int pre_relaxation_S=qna.number_of_species();
    //    qna.relax(standard_relaxation_time/10);
    int post_relaxation_S=qna.number_of_species();
    int species_lost_relaxing=pre_relaxation_S-post_relaxation_S;
    REPORT(species_lost_relaxing);
    qna.eigenvalue_list("qna_evs.dat");
    qna.eigenvalue_list_of_a("a_evs.dat");
  }
  else if(activity == "eval"){

    read(web,argc,argv);
    
    web.delete_by_deactivating();
    web.relax(1);
    web.evaluate();
  }
  else if(activity == "quicky"){
    
    int i=web.add_active_species();
    web.relax(1e3);
    CompetiWeb_t qna=
      (use_local_qna?
       web.local_quasi_neutral_approximation(use_true_qna) :
       web.quasi_neutral_approximation(use_true_qna) );
    web.scale_population(i,0.1);
    qna.scale_population(i,0.1);
    accuracy1_object a1o(web,qna);
    std::cout << "starting relaxation" << endl;
    a1o.relax(standard_relaxation_time);
    REPORT(a1o.rms_max_diff());
    REPORT(a1o.rms_max_diff());
  }
  else if(activity == "parameters"){
    read(web,argc,argv);
    web.report_parameters();
  }
  else if(activity == "perturb"){
    read(web,argc,argv);
    web.relax(10);
    web.delete_all_inactive_species();
    web.eigenvalue_perturbation_analysis(10);
  }
  else if(activity == "singular"){
    read(web,argc,argv);
    web.relax(10);
    web.delete_all_inactive_species();
    web.eigenvalue_singular_analysis(1000);
  }
  else if(activity == "modes"){
    set_cfg_parameter("default_absolute_tolerance","1e-8");
    read(web,argc,argv);
    web.relax(standard_relaxation_time/10);
    cout << "Finished relaxing" << endl;
    web.delete_all_inactive_species();
    decomposition_object dob(web);
    int target_number_of_species=int(dob.number_of_species()*2.0/3.0);
    vector<bool> picked(dob.number_of_species(),false);
    while(dob.number_of_species()>target_number_of_species){
      int s=unirand()*dob.number_of_species();
      if(dob.is_active(s) and not picked[s]){
	target_number_of_species++;
	dob.scale_population(s,1.0/3.0);
	picked[s]=true;
      }
    }
    dob.relax(2000);
  }
  else if(activity == "rho"){
    read(web,argc,argv);
    web.relax(standard_relaxation_time/100);
    cout << "Finished relaxing" << endl;
    web.delete_all_inactive_species();
    web.spectral_radius();
  }
  else if(activity == "S"){
    if(argc>2) 
      read(web,argc,argv);
    else{
      web.add_inactive_species_if_required();
      while(++step_count<=2){
	REPORT(step_count);
	int new_species=
	  web.add_fit_species();
	try{
	  web.relax(standard_relaxation_time);
	}catch(int){
	  web.delete_species(new_species);
	}
      }
    }
    web.predict_species_richness_S();
  }
  else if(activity == "S-line"){
    ofstream os("S-line.dat");
    double delta=0.1;
    for(double x=0.0;x<10;x+=delta){
      CompetiWeb_t 
	webx(6 /* birth */, 
	     1 /* maturation */,
	     1.0 /*juvenile_interaction*/, 1.0 /*adult_interaction*/, 
	     0.1 /*juvenile_i*/, 0.1 /*adult_i*/, 
	     0.5 /*juvenile_C*/, 0.5 /* adult_C*/, 
	     x /*juvenile_death*/, 1 /*adult_death*/, 
	     1 /*these_are_all_population_dynamical_parameters*/);
      step_count=0;
      webx.add_inactive_species_if_required();
      while(++step_count<=2){
	REPORT(step_count);
	int new_species=
	  webx.add_fit_species();
	try{
	  webx.relax(standard_relaxation_time);
	}catch(int){
	  webx.delete_species(new_species);
	}
      }
      os << x << " " << webx.predict_species_richness_S() << endl;
    }
  }
  else if(activity == "logistic"){
    read(web,argc,argv);
    web.delete_all_inactive_species();
    CompetiWeb_t qna=web.quasi_neutral_approximation(use_true_qna);
    qna.logistic();
  }
  else if(activity == "experimental"){
    web.experimental();
  }
  else {
    FATAL_ERROR("unreconized activity: " << activity);
  }
}
  
