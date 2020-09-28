// $Id: MoranF.cc 1431 2009-05-04 12:22:40Z axel $
#include <string>
const std::string version("$Revision: 1431 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/MoranF.cc,v $ $Date: 2009-05-04 13:22:40 +0100 (Mon, 04 May 2009) $");

#include <signal.h>
#include "sequence.h"
#include "random_pick_field.h"
#include <iostream>
#include <fstream>
#include <list>
#include <stack>
#include <vector>
#include <numeric>
#include <set>
#include <math.h>
#include "NetworkAnalysis.h"
#include "scheduler.h"
#include "random.h"
#include "error.h"
#include "Statistics.h"
#include "three_column_files.h"
#include <gsl/gsl_roots.h>
#include "props_t.h"
#include "standardize.h"
#include "polyfit.h"
#include "s_t.h"
#include "topology_generator.h"
#include "Moran.h"
#include "MoranI.h"
#include "MoranIpp.h"

// these are included to force linking
#include "linpack_eigen.h"
#include "binomial_dist.h"
#include "niche_width_finder.h"

//#define NOLINKS
//#define DELETE_NEW

#ifdef DEBUGGING
#undef DEBUGGING
#endif

#ifdef DEBUGGING
static int counter=0;
#endif

using namespace std;

char * cfg_file_name="MoranF.cfg";

const double ln10=M_LN10; //==log(10.0);

char * in_file_name="/home/axel/Webs/St. Martin.web";

static int ndim=3;
static double C0=0.1;
static double beta_r=0.1;
static double beta_c=1;
static double p_f=-1; // used only if set to value >=0
static double p_v=-1; // used only if set to value >=0
static double relaxation_slowdown=1/0.2;
static double std_mutation_factor=0.001;
static double q=1;
static double saturation=0;
static Topology_Generator * web_generator=0;
static int cumulative_stats = 0;
static int number_of_repetitions=1000;
static int print_all_webs=0;
static int random_seed=1;
static int test_mode=0; // 1 = swap timing, 2 = braching prob, 
                        // 8 = save a web an exit
static int save_movie=0;
static int n_images=0;  // how many finished food-webs images to save
static int exact=0;
static int correct_correlations=0;
static int raw_target_S=0;

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGSTRING(in_file_name),
    CFGDOUBLE(C0),
    CFGDOUBLE(lambda),
    CFGDOUBLE(beta_r),
    CFGDOUBLE(beta_c),
    CFGDOUBLE(p_f),
    CFGDOUBLE(p_v),
    CFGDOUBLE(relaxation_slowdown),
    CFGDOUBLE(std_mutation_factor),
    CFGDOUBLE(q),
    CFGDOUBLE(saturation),
    CFGINT(ndim),
    CFGINT(cumulative_stats),
    CFGINT(number_of_repetitions),
    CFGINT(print_all_webs),
    CFGINT(random_seed),
    CFGINT(test_mode),
    CFGINT(save_movie),
    CFGINT(n_images),
    CFGINT(exact),
    CFGINT(correct_correlations),
    CFGINT(raw_target_S),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);

// dependent parameters and constants:
static double target_C;
static double target_S; // number of species
static const double max_s=1;
static const double min_s=0;
static double CC0; //1-C0
static double relaxation_time;
static const double rm=1;
static bool parameters_are_set=false;

//////////////////////////////////////////////////////////
// types and classes:


class bitrow {
  typedef unsigned int word;
  static const word wsize; //size of word in bits
  unsigned int _size;  //number of words in a row
  unsigned int _rowlength(){
    return _size*(sizeof(word)/sizeof(char));  //number of bytes in a row
  }
  unsigned int _bits;
  word * _array;
  void _copy(const bitrow & other){
    memcpy(_array,other._array, _rowlength() );
  }
  void set_size(unsigned int s){
    _bits=s;
    _size=((s+wsize-1)/wsize);
  }
public:
  bitrow(unsigned int s){
    set_size(s);
    _array = new word[_size];
    memset(_array, 0, _rowlength() );
  };
  bitrow(const bitrow & other):_size(other._size),
			       _bits(other._bits)
  {
    _array = new word[_size];
    _copy(other);
  };
  ~bitrow(){
    delete[] _array;
  }
  bitrow & operator=(const bitrow & other){
    if (this != &other){
      if(other._size!=_size){
	delete[] _array;
	_size=other._size;
	_bits=other._bits;
	_array = new word[_size];
      }
      _copy(other);
    }
    return *this;
  }
  void clear(){
    memset(_array, 0, _rowlength());
  }    
  void set(unsigned int i){
    ASSERT( i < _size*wsize );
    _array[i/wsize]|=(1<<(i%wsize));
  }
  void swap(unsigned int i){
    ASSERT( i < _size*wsize );
    _array[i/wsize]^=(1<<(i%wsize));
  }
  void unset(unsigned int i){
    ASSERT( i < _size*wsize );
    _array[i/wsize]&=-1^(1<<(i%wsize));
  }
  // not needed:
//   void set_to(bool b,unsigned int i){
//     ASSERT( i < _size*wsize );
//     if(b)
//       set(i);
//     else
//       unset(i);
//   }
  bool test(unsigned int i)
  {
    ASSERT( i < _size*wsize );
    return  _array[i/wsize] & (1<<(i%wsize));
  }
  void bit_copy(unsigned int i,unsigned int j){
    ASSERT( i < _size*wsize );
    ASSERT( j < _size*wsize );
    if(test(i))
      set(j);
    else
      unset(j);
  }
  void set_with_probability(double p0)
  {// here we expect p0 to be close to zero
    double p=1-p0;
    ASSERT(0 <= p && p < 1);
    unsigned int n=0;
    double ccrand=unirand();
    double cc=p;
    do{
      word w=0;
      for(unsigned int b=wsize ; b>0 ; b--){
	w<<=1;
	if(cc > ccrand){
	  cc*=p;
	}else{
	  cc=p;
	  ccrand=unirand();
	  w|=1;
	}
      }
      _array[n++]=w;
    }while( n < _size);
  }
};


class bitarray : public vector<bitrow>{ 
  //we define a bitarray from scratch for speed
  unsigned int _msize;
public:
  bitarray(unsigned int s=1):vector<bitrow>(s,bitrow(s)),_msize(s){
    for(unsigned int i=0;i<s;i++){
      (*this)[i].clear();
    }
  }
};

const bitrow::word bitrow::wsize = 8*sizeof(bitrow::word);

class robust_stack_t : public stack<unsigned int>{
public:
  robust_stack_t(unsigned int n){
    for(int i=n-1;i>=0;i--){
      push(i);
    }
  }
  unsigned int rpop(){ //Robust, Returning, Removing pop.
    if(empty()) FATAL_ERROR("No more indices");
    unsigned int t=top();
    pop();
    return t;
  }
};

// class rare_true_generator {
//   double random;
//   double threshold;
//   const double p;
//   rare_true_generator():p(0){};
// public:
//   rare_true_generator(double pr):
//     random(unirand()),
//     threshold(1-pr),
//     p(1-pr){};
//   bool operator()(){
//     if(random<threshold){
//       threshold*=p;
//       return false;
//     }else{
//       threshold=p;
//       random=unirand();
//       return true;
//     }
//   }
// };


  
void test_random_int(){
  for(int i=0;i<50;i++){
    unsigned int r=random_int();
    for(int b=0;b<props_t_::wbits;b++){
      cout << (r&1) ;
      r>>=1;
    }
    cout << endl;
  }
  exit(1);
}


static Interaction_Matrix im,exact_im;

typedef list<s_t> slist_tt;
typedef slist_tt::iterator siter;
class slist_t : public slist_tt {
public:
  siter smallest_potential_predator(siter & i,
				    double L=lambda){
    siter k=i;
    if(L > 0){
      double target_s=i->s+L;
      siter b=begin();
      while(k!=b && k->s < target_s) k--;
      if(k->s < target_s) k++;
    }else if(L < 0){
      double target_s=i->s+L;
      siter e=end();
      while(k!=e && k->s > target_s) k++;
    }
    return k;
  }
  siter smallest_too_large_to_prey(siter & i){
    return smallest_potential_predator(i,-lambda);
  }
};


#define forall_species(s,l) for((s)=(l).begin();(s)!=(l).end();(s)++)

typedef NetworkAnalysis::Interaction interaction;
const interaction none=NetworkAnalysis::none;

typedef sequence<int> histogram;
typedef sequence<double> distribution;

static int s_sum;

static slist_t slist;

///////////////////////////////////////////////////////////
// helper inlines:

#ifndef DEBUGGING
inline siter speciate(slist_t &slist,
		     siter i) __attribute__ ((always_inline));
inline 
#endif
siter speciate(slist_t &slist,
	      siter i){
  
  s_sum+=1;

  double delta;
  double new_s;
  do{
    delta=gaussian(0,std_mutation_factor);
    new_s=i->s+delta;
    // no-flux boundary conditions:
    if(new_s>max_s) 
      new_s=(max_s+max_s)-new_s; 
    if(new_s<min_s) new_s=(min_s+min_s)-new_s;
  }while(new_s>max_s || new_s<min_s);

#if 1 // this is faster:
  slist.push_front(*i); //new species
  // old species is i;
  siter j=slist.begin();
#else // this makes nicer movies:
  siter j=slist.insert(i,*i);
#endif

  j->s=new_s;
#ifndef NOLINKS
  j->mutate();
#endif
  return j;
}

#ifndef DEBUGGING
inline void extinguish(slist_t &slist,
		       siter i) __attribute__ ((always_inline));
inline 
#endif
void extinguish(slist_t &slist,
		       siter i){

  s_sum-=1;
  
#ifndef NOLINKS
#endif

  slist.erase(i);
}

// graph display *(*(_List_node<s_t> *) i._M_node )._M_data.predators._array
// graph display *(*(_List_node<s_t> *) j._M_node )._M_data.predators._array

#ifndef DEBUGGING
inline siter invade(slist_t &slist) __attribute__ ((always_inline));
inline 
#endif
siter invade(slist_t &slist){
  s_sum+=1;
  double new_s=min_s+unirand();
  slist.push_front(s_t(new_s)); // new species
  return slist.begin();
}



///////////////////////////////////////////////////////
// other helpers:



Interaction_Matrix get_im(){

  if(exact) return exact_im;

  int s=slist.size();
  Interaction_Matrix rim(s);
  simple_vector<bool> survivors(s,false);
  siter i,j,end=slist.end();
  int ii=s;
  for(i=slist.begin();i!=end;i++){
    ii--;
    if(i->is_plant)
      survivors[ii]=true;//plants always survive
    int jj=s;
    for(j=slist.begin();j!=end;j++){
      jj--;
      if(j->eats(*i)){
	rim[jj][ii]=NetworkAnalysis::eats;
	// This rule is too simple! There may be loops of animals that
	// eat each other, but do not get any input:
	survivors[jj]=true;//animals that find food survive
      }else{
	rim[jj][ii]=NetworkAnalysis::none;
      }
    }
  }
  if(plant_fraction>=0)
    return rim.select(survivors);
  else
    return rim;
}

double get_C0(slist_t & slist){
  int s=slist.size();
  average_meter avC0; 
  siter i,j,end=slist.end();
  int ii=s;
  for(i=slist.begin();i!=end;i++){
    ii--;
    int jj=s;
    for(j=slist.begin();j!=end;j++){
      jj--;
      avC0.sample(i->eats(*j) ? 1 : 0);
    }
  }
  return avC0.readout();
}

Interaction_Matrix real_im(Interaction_Matrix & im, 
			   slist_t & slist){
  Interaction_Matrix gim=get_im();
  return standardize(gim);
}

histogram index(int n){
  histogram h;
  for(int i=n-1;i>=0;i--)
    h[i]=i;
  return h;
}

histogram s_hist(slist_t & slist){
  histogram h;
  siter s;
  forall_species(s,slist){
    h[int(10*(s->s))]++;
  }
  return h;
}


// compute the machine precision:
double machine_precision(){
  double precision=1.0;
  double one=1.0;
  while((one+precision)!=one) precision/=2;
  return 2*precision;
}

void read_arguments(int argc,char *argv[]){
  int c;
  
  while(1){
    c=getopt(argc, argv, "hv");
    if (c == -1)
      break;
    
    switch(c){
    case 'v':
      printf("%s\n",version.c_str());
      exit(0);
      break;
    case 'h':
    default:
        fprintf(stdout,"usage: %s [config_file_name]\n",argv[0]);
        exit(c=='h'?0:1);
      }
  }
  argc-=optind-1;
  argv+=optind-1;
  
  if(argc>1) cfg_file_name=argv[1];
}

typedef unsigned long int loop_counter_t;

/// member randint(i) generates random numbers < i
class randint{
public:
  unsigned operator()(unsigned i){
    return unsigned(i*unirand());
  }
};

void make_new_exact_web(int S0); // see below

void make_new_web(double rS0){

  if(!exact){
    FATAL_ERROR("only exact webs are implemented");
  }
  exact_im=web_generator->draw_sample(rS0);
  cout << "+" ;
  cout.flush();
  if(save_movie && parameters_are_set) exit(0);
}

average_meter measure_S_at_precision(const double precision,const double S0, 
				     average_meter S=average_meter()){
  ASSERT(precision > 0);

  while(S.n()<50){
    make_new_web(S0);
    Interaction_Matrix im=get_im();
    S.sample(standardize(im).Number_of_Species_S());
  }
  while(S.error() > precision){
    //REPORT(S.error() / precision );
    make_new_web(S0);
    Interaction_Matrix im=get_im();
    S.sample(standardize(im).Number_of_Species_S());
  }
  cout << endl;
  return S;
}

#define SQUARED(X) ((X)*(X))

average_meter measure_S_next_step(double S0,double S1, double D){

  average_meter S;

  for(int i=0;i<50;i++){
    make_new_web(S0);
    Interaction_Matrix im=get_im();
    S.sample(standardize(im).Number_of_Species_S());
  }
  while(true){
    // stop if hit:
    if(fabs(S.readout()-S1)+S.error()<D) break;
    // stop if apparently not hit:
    if(S.error()<fabs(S.readout()-S1)/3) break;
    make_new_web(S0);
    Interaction_Matrix im=get_im();
    S.sample(standardize(im).Number_of_Species_S());
  }
  cout << endl;
  return S;
}

double var_of_estimate(double x1,
		       double x2,
		       average_meter Y1,
		       average_meter Y2,
		       double y0){
  return SQUARED(x1-x2)*
    (Y2.error_var()*SQUARED(Y1.readout()-y0)+
     Y1.error_var()*SQUARED(Y2.readout()-y0) )/
    SQUARED(SQUARED(Y1.readout()-Y2.readout()))+
    //extra penalty for distant points (nonlinearities or so):
    SQUARED(SQUARED(2*(x1-x2)/(x1+x2)))*SQUARED((x1+x2)/2);
}
    
void report_this(sequence<double> & S01,
		 sequence<average_meter> &S1,
		 int n){
  cout << "S01[" << n << "] = " 
       << S01[n] << endl;
  cout << "S1[" << n << "] = " 
       << S1[n].readout() << " +/- " << S1[n].error() << endl;
}

struct f_S_differencer_par {
  fitted_function * f;
  double S;
};

double f_S_differencer(double S0,void * vpar){
  f_S_differencer_par * par = (f_S_differencer_par *) vpar;
  return (*(par->f))(S0) - (par->S);
}

#define ADJUST_MAX_KAPPA_R(X) do{if(S01.last()>max_S0) max_S0=S01.last();}while(0)

double get_S0(){
  if(raw_target_S){
    return raw_target_S;
  }
  REPORT(target_S);
  const double D=0.08*target_S;

  //initials:
  sequence<double> S01;
  sequence<average_meter> S1;
  S01[0]=target_S;
  S01[1]=1.5*S01[0];
  double max_S0=S01[1];
  double precision=target_S/5.0;
  S1[0]=measure_S_at_precision(precision,S01[0]);
  report_this(S01,S1,0);
  S1[1]=measure_S_at_precision(precision,S01[1]);
  report_this(S01,S1,1);
  
  int rep=0;
  bool last_fit=false;
  double result=0;

  do{
    fitted_function f=fitted_function(S01,S1);
    //guarantee a zero of fitted function
    while(f(max_S0) <= target_S  || f(0) >= target_S){
      S01[S1.size()]=2*max_S0;
      S1[S1.size()]=measure_S_at_precision(precision,S01.last());
      report_this(S01,S1,S1.size()-1);
      ADJUST_MAX_KAPPA_R();
      f=fitted_function(S01,S1);
      while(f(0) >= target_S){
	S01[S1.size()]=max_S0/pow(2.0,double(S1.size()));
	S1[S1.size()]=measure_S_at_precision(precision,S01.last());
	report_this(S01,S1,S1.size()-1);
	ADJUST_MAX_KAPPA_R();
	f=fitted_function(S01,S1);
      };
    };

    REPORT(f(0));
    REPORT(f(max_S0));
    
    ////find root f(S0)==target_S
    // set up gsl root finder:
    double lower,upper;
    const gsl_root_fsolver_type * T =
      //gsl_root_fsolver_bisection;
      gsl_root_fsolver_brent;
    gsl_root_fsolver * s =
      gsl_root_fsolver_alloc (T);
    gsl_function F;
    F.function = &f_S_differencer;
    f_S_differencer_par par;
    par.f=&f;
    par.S=target_S;
    F.params = &par;
    if(int error=gsl_root_fsolver_set (s, &F, 0, max_S0)){
      REPORT(error);
      FATAL_ERROR("initializing gsl root solver");
    }
    cout << "using " << gsl_root_fsolver_name(s) << 
      " algorithm for root finding" << endl;
    
    // find root:
    int root_finding_iterations=0;
    do{
      if(gsl_root_fsolver_iterate(s))
	FATAL_ERROR("problem computing swapping prob");
      upper=gsl_root_fsolver_x_upper (s);
      lower=gsl_root_fsolver_x_lower (s);
      REPORT(upper);
      REPORT(lower);
      if(++root_finding_iterations>1000)
	FATAL_ERROR("too many iterations computing next S0");
    }while(upper-lower> 0.000001);
    REPORT(root_finding_iterations);
    REPORT(f(gsl_root_fsolver_root(s)));
    if(last_fit){
      result=gsl_root_fsolver_root(s);
      gsl_root_fsolver_free(s);
      break;
    }
    S01[S1.size()]=gsl_root_fsolver_root(s);
    ADJUST_MAX_KAPPA_R();
    average_meter Snew=measure_S_next_step(S01[S1.size()],target_S,D);
    S1[S1.size()]=Snew;
    report_this(S01,S1,S1.size()-1);

    if(fabs(Snew.readout()-target_S)+Snew.error()<D ||
       (test_mode&2) )
      last_fit=true;

  }while((++rep < 100));
  cout << "last estimate:" << endl;
  cout << "S = " << S1.last().readout() << " +/- " 
       << S1.last().error() << endl;
  if(rep>=100) 
    FATAL_ERROR("could not determine correct S0");
  
  return result;
}

void exiter(int i){
  exit(i);
}

void signal_handling(){
  signal(SIGUSR1,&exiter);
}

/////////////////////////////////////////////
// main

int main(int argc,char *argv[]){
  read_arguments(argc,argv);
  cout << source << endl;
  cout << version << endl;
  system("/bin/date");
  signal_handling();
  //test_BinomialDist();
  //test_random_int();

  REPORT(cfg_file_name);

  if (cfgParse(cfg_file_name, full_cfg_list(), CFG_SIMPLE) == -1)
    FATAL_ERROR("error reading parameter file");

  REPORT(random_seed);
  set_random_seed(random_seed);

  if(test_mode&1){
    props_t::measure_swapping_methods();
    exit(1);
  }
  if(test_mode&2){
    // try all kinds of swapping methods:
    const int n_reps=int(1e5);
    int i;
    double p;
    props_t props;
    props_t::swapping_probability prob;
    for(p=0.6,i=0;p>0.0001;p*=(1/1.3),i++){
      prob=p;
      for(int j=n_reps;j>0;j--)
	props.swap_a_fraction(prob);
    }
  }

  if(test_mode&4){
    // try all kinds of swapping methods:
    props_t::test_modes();
  }

  REPORT(in_file_name);

  string file_name=in_file_name;
  Interaction_Matrix raw_empirical=load_three_column_file(file_name);
  Interaction_Matrix empirical=standardize(raw_empirical);

  {// save empirical web:
    distribution distr;
    if(cumulative_stats)
      distr=
	distribution(empirical.cumulative_prey_hist())/
	empirical.Number_of_Species_S();
    else
      distr=distribution(empirical.prey_hist())/empirical.Number_of_Species_S();
    std::ofstream XY("XY.dat") ;
    XY << distribution(index(distr.size())).format("%#8.3g ") +
      distr.format("%#8.3g ")+
      "";
    
    if(cumulative_stats)
      distr=
	distribution(empirical.cumulative_predator_hist())/
	empirical.Number_of_Species_S();
    else
      distr=distribution(empirical.predator_hist())/empirical.Number_of_Species_S();
    std::ofstream XR("XR.dat") ;
    XR << distribution(index(distr.size())).format("%#8.3g ") +
      distr.format("%#8.3g ")+
      "";
    empirical.random_shuffle().tsort().pgm_write("imgX.pgm");
  }

  distribution eprops=empirical.props();
  target_S=eprops[Interaction_Matrix::pS];
  target_C=eprops[Interaction_Matrix::pC];
  REPORT(target_S);
  REPORT(target_C);

  simple_vector<bool> select(Interaction_Matrix::pend,false);
  simple_vector<bool> fix(Interaction_Matrix::pend,false);
//   enum {pS,pC,pZ,pT,pI,pB,pGenSD,pVulSD,pMxSim,
// 	pChnLg, pChnSD,pChnNo,pLoop,pCannib,pOmniv,
// 	poChnLg, poChnSD,poChnNo,poLoop,poOmniv,
// 	pDdiet,
// 	pend} property_t;
  select[Interaction_Matrix::pS]=true;
  //  fix[Interaction_Matrix::pS]=true;
  select[Interaction_Matrix::pC]=true;
  //  fix[Interaction_Matrix::pC]=true;
  select[Interaction_Matrix::pT]=true;
  if(plant_fraction>=0){
    select[Interaction_Matrix::pB]=true;
  }
  select[Interaction_Matrix::pGenSD]=true;
  select[Interaction_Matrix::pVulSD]=true;
  select[Interaction_Matrix::pMxSim]=true;
  select[Interaction_Matrix::pCannib]=true;
  select[Interaction_Matrix::poLoop]=true;
  select[Interaction_Matrix::poOmniv]=true;
  select[Interaction_Matrix::poChnLg]=true;
  select[Interaction_Matrix::poChnSD]=true;
  select[Interaction_Matrix::poChnNo]=true;
  select[Interaction_Matrix::pfDdiet]=true;
  select[Interaction_Matrix::pClust]=true;

  select[Interaction_Matrix::pI]=true;
  select[Interaction_Matrix::pB]=true;
  select[Interaction_Matrix::pZ]=true;

  static chi_square_meter multivar(select);

  ALWAYS_ASSERT(number_of_repetitions>=0);
  ALWAYS_ASSERT(lambda>=0);
  ALWAYS_ASSERT(C0>0);
  ALWAYS_ASSERT(C0<1);

  s_t::set_C0(C0);
  double true_C0=s_t::get_true_C0();
  REPORT(true_C0);
  if(p_v>=0){
    ALWAYS_ASSERT(p_v>=0);
    ALWAYS_ASSERT(p_v<=1);
    s_t::set_swap_armor(p_v);
  }else{
    FATAL_ERROR("setting beta_c won't work anymore");
  }
  if(p_f>=0){
    ALWAYS_ASSERT(p_f>=0);
    ALWAYS_ASSERT(p_f<=1);
    s_t::set_swap_weapons(p_f);
  }else{
    FATAL_ERROR("setting beta_r won't work anymore");
  }

//   web_generator = new Moran(std_mutation_factor*std_mutation_factor,
// 			    C0,p_f,p_v,lambda,q,correct_correlations);
#if 1
  web_generator = new MoranI(ndim,std_mutation_factor*std_mutation_factor,
			     C0,lambda,saturation,q,correct_correlations);
#else
   web_generator = 
     new MoranIpp(ndim,std_mutation_factor*std_mutation_factor,
 		 C0,lambda,saturation,q,correct_correlations,
 		 0.05/*T*/,1000/*D_f*/,1/*over variability*/,0.001/*sensitivity*/,
 		 0.5/*granularity*/);
#endif

  const double S0=get_S0();
  REPORT(S0);
  const double predation_area_fraction=
    (0.5+lambda-lambda*lambda/2);
  REPORT(predation_area_fraction);
  //C0=C/predation_area_fraction; // density of links in upper triangle
  CC0=1-C0;
  REPORT(C0);
  const double predation_bias_oversize=lambda;
  REPORT(predation_bias_oversize);
  //// prepare for progressive refinement of time scale:
  const double effective_invasion_time=S0;
  REPORT(effective_invasion_time);
  relaxation_time=0;//relaxation_slowdown*(effective_invasion_time+effective_clade_lifetime);
  
  int loop_counter=0;

  parameters_are_set=true;
  for(int repetitions_left=number_of_repetitions;
      repetitions_left > 0;/* decrement only of hit target S and C */){
    loop_counter++;
    if(loop_counter > 40*number_of_repetitions ) 
      FATAL_ERROR("not enought hits");

    
    make_new_web(S0); cout << endl;
    ///////////////////////////////////////////////
    // End of algorithm.  From here on evaluation.

#ifndef NOLINKS
    static sequence<int> N;
    static sequence<distribution> Ydist; // sum of distributions
    static sequence<distribution> Ydist2; // sum of squares of distributions
    static sequence<distribution> Rdist; // sum of distributions
    static sequence<distribution> Rdist2; // sum of squares of distributions
    static int NSdist=0; // number of distributions added;
    static average_meter Smeter,Cmeter,rawSmeter,rawCmeter,empC0;
    static distribution Sdist; // sum of squares of distributions
    static distribution Sdist2; // sum of squares of distributions
    static sequence<average_meter> avs;
    //REPORT(Z);
    //REPORT(relative_deltaZ);
    //REPORT(deltaZ);
    im = get_im();
    if(test_mode & 8){
      save_three_column_file(im,"sample.web");
      exit(0);
    }
    im.Print();
    int rawS=im.Number_of_Species_S();
    REPORT(rawS);
    rawSmeter.sample(rawS);
    double rawC;
    if(rawS){
      if(rawS>1)
	empC0.sample(get_C0(slist));
      rawC=im.Number_of_Links_L()/(double(rawS)*rawS);
      rawCmeter.sample(rawC);
    }else{
      rawC=C0*predation_area_fraction;
    }
    REPORT(rawC);
    cout << "mean rawS: " << rawSmeter.readout() << " +/- " << rawSmeter.error() << " theory: " << S0 << endl;
    cout << "mean rawC: " << rawCmeter.readout() << " +/- " << rawCmeter.error() << " theory: " << C0*predation_area_fraction << endl;
    //    cout << "mean empC0: " << empC0.readout() << " +/- " << empC0.error() << " theory: " << C0 << endl;
    
    Interaction_Matrix rim=standardize(im);
    //Interaction_Matrix rim=im;
    if(true || 0<rim.Number_of_Species_S()){
      Interaction_Matrix & dist_im = rim;
      double actualS = rim.Number_of_Species_S();
      REPORT(actualS/target_S);
      double actualL = 
	double(rim.Number_of_Links_L());
      //REPORT(actualL);
      double actualC =
	actualL/(actualS*actualS);
      REPORT(actualC/target_C);
      distribution dist;
      if(!raw_target_S && fabs(actualS-target_S)/target_S > 0.3){
 	cout << "-" << endl;
 	goto skip;
      }
      cout << "accepted" << endl;
      repetitions_left--;
      //int Zbin=int(10*actualS*actualC/(target_S*target_C)+0.5);
      int Zbin=int(actualS+0.5);
      if(Zbin==0) goto skip;
      N[Zbin]++;
      N[0]++;
      cout << "N[" << Zbin << "] " << N[Zbin] << endl;
      if(cumulative_stats)
	dist=
	  distribution(dist_im.cumulative_prey_hist())/
	  dist_im.Number_of_Species_S();
      else
	dist=distribution(dist_im.prey_hist())/dist_im.Number_of_Species_S();
      Ydist[Zbin]+=dist;
      Ydist2[Zbin]+=dist*dist;
      Ydist[0]+=dist;
      Ydist2[0]+=dist;
      if(cumulative_stats)
	dist=
	  distribution(dist_im.cumulative_predator_hist())/
	  dist_im.Number_of_Species_S();
      else
	dist=distribution(dist_im.predator_hist())/dist_im.Number_of_Species_S();
      Rdist[Zbin]+=dist;
      Rdist2[Zbin]+=dist*dist;
      Rdist[0]+=dist;
      Rdist2[0]+=dist;
 
      Smeter.sample(rim.Number_of_Species_S());
      Cmeter.sample(rim.connectance_C());

      distribution p=rim.props();
      for(unsigned int i=0;i<p.size();i++){
	if(!finite(p[i])){
	  REPORT(p);
	  WARNING("non-finite prop " << 
		  Interaction_Matrix::prop_names()[i] << 
		  " encountered. skipping" );
	  goto skip;
	}
      }
      avs+=p;
      multivar.sample(p);
      dist=distribution(s_hist(slist));
      cout << "Scurr " << dist << endl;
      cout << "Sdist " << accumulate(Sdist.begin(),Sdist.end(),0.0)/NSdist 
	   << endl;
      Sdist+=dist;
      Sdist2+=dist*dist;
      NSdist++;
      //Smeter.sample(accumulate(dist.begin(),dist.end(),0.0));
    }else{
      // all species extinct:
      Smeter.sample(0);
    }
#else //defined(NOLINKS)
    static int NSdist=0; // number of distributions added;
    static distribution Sdist; // sum of squares of distributions
    static distribution Sdist2; // sum of squares of distributions
    static average_meter Smeter;
    if(true/*sample_Sdist.due_now(t)*/){
      distribution dist;
      dist=distribution(s_hist(slist));
      cout << "Scurr " << dist << endl;
      cout << "Sdist " << accumulate(Sdist.begin(),Sdist.end(),0.0)/NSdist 
	   << endl;
      dist=distribution(s_hist(slist));// /slist.size();
      Sdist+=dist;
      Sdist2+=dist*dist;
      NSdist++;
      Smeter.sample(accumulate(dist.begin(),dist.end(),0.0));
    }
#endif
    if(true/*analysis.due_now(t)*/){
      int rawS=slist.size();
      cout << endl;
      REPORT(repetitions_left);
      REPORT(rawS);
      REPORT(NSdist);
#ifndef NOLINKS
      REPORT(N);
      Interaction_Matrix rim=real_im(im,slist);
      double actualZ = 
	double(rim.Number_of_Links_L())/
	rim.Number_of_Species_S();
      cout << "actualS " << rim.Number_of_Species_S() << endl;
      cout << "actualZ " << actualZ << endl;
      cout << "actualC " << actualZ/rim.Number_of_Species_S() 
	   << " =?= " << target_C << " ?" << endl;

      static sequence<int> previous_N;
      for(unsigned int Zbin=0; Zbin<Ydist.size(); Zbin++){
	if(N[Zbin] > previous_N[Zbin]){
	  Interaction_Matrix & dist_im = im;
	  distribution Ymean=Ydist[Zbin]/double(N[Zbin]);
	  distribution Ystd=sqrt(double(N[Zbin])/double(N[Zbin]-1))*
	    Map(::sqrt,Ydist2[Zbin]/double(N[Zbin])-Ymean*Ymean);
	  distribution Yerror=Ystd/sqrt(double(N[Zbin]));
	  double meanZ;
	  distribution meanZbasis;
	  if(cumulative_stats)
	    meanZbasis=Ymean;
	  else
	    meanZbasis=Ymean.reverse_cumulative_sum();
	  meanZ=accumulate(meanZbasis.begin(),meanZbasis.end(),0.0);
	  distribution Ytheory;
	  if(cumulative_stats)
	    Ytheory = dist_im.theoretical_cumulative_prey_dist(Yerror.size(),meanZ);
	  else
	    Ytheory = dist_im.theoretical_prey_dist(Yerror.size(),meanZ);
	  std::ofstream EY(("EY"+format("%02i",Zbin)+".dat").c_str()) ;
	  std::ofstream TY(("TY"+format("%02i",Zbin)+".dat").c_str()) ;
	  EY << 
	    distribution(index(Yerror.size())).format("%#8.3g ") +
	    Ymean.format("%#8.3g ")+
	    (Ymean+Yerror).format("%#8.3g ")+
	    (Ymean-Yerror).format("%#8.3g ")+
	    (Ymean+Ystd).format("%#8.3g ")+
	    (Ymean-Ystd).format("%#8.3g ")+
	    "";
	  TY <<
	    distribution(index(Yerror.size())).format("%#8.3g ") +
	    Ytheory.format("%#8.3g ");

// 	  EY <<
// 	    (distribution(index(Yerror.size()))/meanZ).format("%#8.3g ") +
// 	    (Ymean*meanZ).format("%#8.3g ")+
// 	    ((Ymean+Yerror)*meanZ).format("%#8.3g ")+
// 	    ((Ymean-Yerror)*meanZ).format("%#8.3g ")+
// 	    ((Ymean+Ystd)*meanZ).format("%#8.3g ")+
// 	    ((Ymean-Ystd)*meanZ).format("%#8.3g ")+
// 	    "";
// 	  TY <<
// 	    (distribution(index(Yerror.size()))/meanZ).format("%#8.3g ") +
// 	    (Ytheory*meanZ).format("%#8.3g ");
	  
	  
	  distribution Rmean=Rdist[Zbin]/double(N[Zbin]);
	  distribution Rstd=sqrt(double(N[Zbin])/double(N[Zbin]-1))*
	    Map(::sqrt,Rdist2[Zbin]/double(N[Zbin])-Rmean*Rmean);
	  distribution Rerror=Rstd/sqrt(double(N[Zbin]));
	  distribution Rtheory;
	  if(cumulative_stats)
	    Rtheory = 
	      dist_im.theoretical_cumulative_predator_dist(Rerror.size(),meanZ);
	  else
	    Rtheory = dist_im.theoretical_predator_dist(Rerror.size(),meanZ);
	  std::ofstream ER(("ER"+format("%02i",Zbin)+".dat").c_str()) ;
	  std::ofstream TR(("TR"+format("%02i",Zbin)+".dat").c_str()) ;
	  ER << 
	    distribution(index(Rerror.size())).format("%#8.3g ") +
	    Rmean.format("%#8.3g ")+
	    (Rmean+Rerror).format("%#8.3g ")+
	    (Rmean-Rerror).format("%#8.3g ")+
	    (Rmean+Rstd).format("%#8.3g ")+
	    (Rmean-Rstd).format("%#8.3g ")+
	    "";
	  TR <<
	    distribution(index(Rerror.size())).format("%#8.3g ") +
	    Rtheory.format("%#8.3g ");
// 	  ER <<
// 	    (distribution(index(Rerror.size()))/meanZ).format("%#8.3g ") +
// 	    (Rmean*meanZ).format("%#8.3g ")+
// 	    ((Rmean+Rerror)*meanZ).format("%#8.3g ")+
// 	    ((Rmean-Rerror)*meanZ).format("%#8.3g ")+
// 	    ((Rmean+Rstd)*meanZ).format("%#8.3g ")+
// 	    ((Rmean-Rstd)*meanZ).format("%#8.3g ")+
// 	    "";
// 	  std::ofstream TR(("TR"+format("%02i",Zbin)+".dat").c_str()) ;
// 	  TR <<
// 	    (distribution(index(Rerror.size()))/meanZ).format("%#8.3g ") +
// 	    (Rtheory*meanZ).format("%#8.3g ");
	} // N[Zbin]!=0
      } // for Zbin
      previous_N=N;
#endif
      
      distribution Smean=Sdist/double(NSdist);
      distribution Sstd=sqrt(double(NSdist)/double(NSdist-1))*
	Map(::sqrt,Sdist2/double(NSdist)-Smean*Smean);
      distribution Serror=Sstd/sqrt(double(NSdist));

      cout << "Smean " << Smean << endl;
      cout << "Sstd  " << Sstd/sqrt(double(NSdist)) << endl;
      cout << "Sav = " << Smeter.readout() << " +/- " 
	   << Smeter.error() 
	   << endl;
#ifndef NOLINKS
      cout << "Cav = " << Cmeter.readout() << " +/- " 
	   << Cmeter.error() 
	   << " (theory: "
	   << "not available"
	   << ")"
	   << endl;
#endif
      

      std::ofstream ES("ES.dat") ;
      ES <<
	distribution(index(Smean.size())).format("%#8.3g ") +
	(Smean).format("%#8.3g ") +
	(Serror).format("%#8.3g");

#ifndef NOLINKS

      if(print_all_webs ||
	 repetitions_left < 4 || 
	 repetitions_left > number_of_repetitions-4 ){
	// sorting takes a lot of time, so we won't always do this:
	rim.random_shuffle().tsort().Print();
      }
      if(n_images>=0){
	static char img_name[20];
	sprintf(img_name,"img%d.pgm",n_images--);
	rim.random_shuffle().tsort().pgm_write(img_name);
      }
      //       rim.prey_hist(cout,"AYH");
      //       rim.predator_hist(cout,"ARH");
      //       rim.cumulative_prey_hist(cout,"AYCH");
      //       rim.cumulative_predator_hist(cout,"ARCH");
      sequence<string> names=Interaction_Matrix::prop_names();
      for(unsigned int j=0;j<avs.size();j++){
	average_meter av(avs[j]);
	cout << names[j] << " " 
	     << eprops[j] << " in "
	     << av.readout()-av.std() 
	     << " -- "
	     << av.readout()+av.std() 
	     << " ( "
	     << av.readout()
	     << " +/- " 
	     << av.error()
	     << " )? "
	     << av.sample_var()
	     << endl;
      }
#endif
    } // analysis
#ifndef NOLINKS
  skip:;
#endif
  } // repetitions

  sequence<string> names=Interaction_Matrix::prop_names();
  multivar.estimate(eprops,fix).save("model");
  save(eprops,select,fix,"data");
  simple_vector<weighted_average_meter> xavs=
    multivar.mean_and_var(eprops,fix);
  for(unsigned int j=0;j<xavs.size();j++){
    weighted_average_meter av(xavs[j]);
    cout << names[j] << " " 
	 << eprops[j] << " in "
	 << av.readout()-av.sample_std() 
	 << " -- "
	 << av.readout()+av.sample_std() 
	 << " ( "
	 << av.readout()
	 << " ) "
	 << endl;
  }
  cout << endl;

  double log_det_cov;
  double chi2=multivar.chi_square(&log_det_cov, eprops,fix);
  cout << "chi^2 = " << chi2 << endl;
  simple_vector<double> dev=multivar.deviation(eprops,fix);
  // redeclared sequence<string> names=Interaction_Matrix::prop_names();

  for(unsigned int j=0;j<xavs.size();j++){
    weighted_average_meter av(xavs[j]);
    cout << names[j] 
	 << "&" 
	 << av.readout()
	 << "&"
	 << av.sample_std() 
	 << "&"
	 << endl;
  }
  cout << endl;

  cout << "deviations:" << endl;
  for(unsigned int i=0;i<dev.size();i++){
    if(!isnan(dev[i])){
      cout << names[i] << " " << dev[i] << endl;
    }
  }
  cout << "for 'space':" << endl;
  cout << chi2 + log_det_cov << "  " << chi2 << endl;
} // main
