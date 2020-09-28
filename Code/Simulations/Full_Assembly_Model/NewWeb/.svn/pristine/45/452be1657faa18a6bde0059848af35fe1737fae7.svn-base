// $Id:$
#include <string>
const std::string version("$Revision$");
const std::string source("$Source$ $Date$");

#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <complex>
typedef std::complex<double> cmplx;
#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


// This would be required only for n != 3/4 
// // Using AEAE_v1_0 from http://www.cpc.cs.qub.ac.uk/, must
// // acknowledge the paper Fast computation of the Gauss hypergeometric
// // function with all its parameters complex with application to the
// // Pöschl-Teller-Ginocchio potential wave functions by N. Michel and
// // M.V. Stoitsov

// using namespace std;
// #define SIGN(a) (((a) < 0) ? (-1) : (1))
// #include "complex_functions.H"
// #include "hyp_2F1.cpp"



#include <fftw3.h>

#include <plotter.h>

#include <time.h>

#include "ODE.h"

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
#include "polyfit.h"
#include "NewMatrix.h"
#endif
int max_num_threads=1;

using namespace std;

// Enable algebra mixing complex and int:

inline cmplx operator*(int i, cmplx z){
  return z*double(i);
}
inline cmplx operator*(cmplx z,int i){
  return z*double(i);
}
inline cmplx operator+(int i, cmplx z){
  return z+double(i);
}
inline cmplx operator+(cmplx z,int i){
  return z+double(i);
}
inline cmplx operator/(int i, cmplx z){
  return double(i)/z;
}
inline cmplx operator/(cmplx z,int i){
  return z/double(i);
}
inline cmplx operator-(int i, cmplx z){
  return double(i)-z;
}
inline cmplx operator-(cmplx z,int i){
  return z-double(i);
}


// A.K. Hui, B.H. Armstrong, A.A. Wray, Rapid computation of the Voigt
// and complex error functions, J. Quant. Spectrosc. Radiat. Transfer
// 19 (1978) 509–516, doi:10.1016/0022-4073(78)90019-5.
//
// This is supposed to have RELATIVE accuracy of 6 digits all over the
// complex plane.
class cef_implementation_hui_et_al_1978 {
public:
  cef_implementation_hui_et_al_1978(int dummy){};
  cmplx operator()(cmplx z0) const {
    bool negate=(z0.imag()<0);

    cmplx z;
    if(negate){
      z=cmplx(-z0.imag(),z0.real());
    }else{
      z=cmplx(z0.imag(),-z0.real());
    }

    const double 
      a0 = 122.607931777104326, b0 = 122.607931773875350,
      a1 = 214.382388694706425, b1 = 352.730625110963558,
      a2 = 181.928533092181549, b2 = 457.334478783897737,
      a3 = 93.155580458138441, b3 = 348.703917719495792,
      a4 = 30.180142196210589, b4 = 170.354001821091472,
      a5 = 5.912626209773153, b5 = 53.992906912940207,
      a6 = 0.564189583562615, b6 = 10.479857114260399,
      a7 = 0, b7 = 1;
    cmplx numerator_sum=a0,denominator_sum=b0;
    cmplx z_power=z;
    numerator_sum+=z_power*a1;
    denominator_sum+=z_power*b1;
    z_power*=z;
    numerator_sum+=z_power*a2;
    denominator_sum+=z_power*b2;
    z_power*=z;
    numerator_sum+=z_power*a3;
    denominator_sum+=z_power*b3;
    z_power*=z;
    numerator_sum+=z_power*a4;
    denominator_sum+=z_power*b4;
    z_power*=z;
    numerator_sum+=z_power*a5;
    denominator_sum+=z_power*b5;
    z_power*=z;
    numerator_sum+=z_power*a6;
    denominator_sum+=z_power*b6;
    z_power*=z;
    numerator_sum+=z_power*a7;
    denominator_sum+=z_power*b7;

    if(negate){
      return -numerator_sum/denominator_sum+double(2)*exp(z*z);
    }else{
      return numerator_sum/denominator_sum;
    }
  }
};

class cef_implementation_weideman_1994{
  // implements an approximation of cef(z)=exp(-z^2) erfc(-i z) 

  // See Weideman 1994.

  int N,M,M2;
  double L;
  vector<double> a;
public:
  cef_implementation_weideman_1994(int N);
  cmplx operator()(cmplx z) const;
};

cef_implementation_weideman_1994::cef_implementation_weideman_1994(int n):
  N(n),M(2*N),M2(2*M),L(sqrt(N/sqrt(double(2)))),a(N)
{
  fftw_complex *in, *out;
  fftw_plan p;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M2);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * M2);
  p = fftw_plan_dft_1d(M2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  for(int k=-M+1;k<=M-1;k++){
    const double theta=k*M_PI/M;
    const double t=L*tan(theta/2);
    const double f= exp(-t*t)*(L*L+t*t);
    in[(k+M2)%M2][0]=f;
    in[(k+M2)%M2][1]=0;
  }
  in[M][0]=0;
  in[M][1]=0;
  fftw_execute(p);
  for(int i=1; i<=N; i++){
    //Order of polynomial coefficients inverted as compared to
    //Weideman 1994 !
    a[i-1]=out[i][0]/M2;
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}
    
cmplx cef_implementation_weideman_1994::operator()(cmplx z) const{
  cmplx denom=double(1)/(L-cmplx(0,1)*z);
  cmplx Z=(L+cmplx(0,1)*z)*denom;
  cmplx p=0;
  cmplx powZ=1;
  for(int i=0;i<N;i++){
    p+=powZ*a[i];
    powZ*=Z;
  }
  return (p+p)*denom*denom+denom/sqrt(M_PI);
}

#define CEF_TEST(x,y) REPORT(cef1(cmplx(x,y))); REPORT(cef2(cmplx(x,y))); REPORT(exp(y*y)*erfc(y));

typedef cef_implementation_hui_et_al_1978 cef_implementation;

class erf_implementation : private cef_implementation{
  cmplx cef(cmplx z) const {
    // cout << "cef[" << z.real() << "+ I * " << z.imag() 
    // 	 << "] ==" << cef_implementation::operator()(z) << endl;
    return cef_implementation::operator()(z);
  }
public:
  erf_implementation(int N):cef_implementation(N){};
  cmplx operator()(cmplx z) const;
};

cmplx erf_implementation::operator()(cmplx z) const {
  return double(1)-exp(-z*z)*cef(cmplx(0,1)*z);
}

class erfc_implementation : private cef_implementation{
  cmplx cef(cmplx z) const {
    return cef_implementation::operator()(z);
  }
public:
  erfc_implementation(int N):cef_implementation(N){};
  cmplx operator()(cmplx z) const;
};

cmplx erfc_implementation::operator()(cmplx z) const {
    // cout << "Erfc[" << z.real() << "+ I * " << z.imag() 
    // 	 << "] ==" << exp(-z*z)*cef(cmplx(0,1)*z) << endl;
  return exp(-z*z)*cef(cmplx(0,1)*z);
}

const int BITMAPWIDTH=900, BITMAPHEIGHT=300;

// adjustable parameters:
int balance_case=2;
double alpha=0.6;
double h=85;
const double n=3/double(4); // fixed
double k=10;
double beta=100;
double w=1;
double search_gamma=1;
double gamma_boost=1;
int gamma_as_pp_boost=false;
double q=0.8;
double eta4=pow(4,1/double(4));
double sigma=1*log(10);
double mesh_size=0;
double absoluteF=0;
double lambda=2+q-n;
double relativeF=0;
double rho=1;
double wD=w;
double bRange=2;  
double x=1;//Scalefactor for transtition to eutrophication
int max_zeros=0;
int doSimulate=true;
double pulse_perturbation=0;
int compute_linear_mode=0;
int noSpecies=0;
int communitySpectrum=0;
int maturation_schedule=1;
double x0=0; // (offspring size) / (maturation size)
double stoppingTime=1e14;
char * show_zero_movie=0;
double upper_cutoff=0; 
int nonlinear_version=0; //0: infinite size axis, 1: using upper
			 //cutoff, 2: LV model, 3: fully_nonlinear
int phyto_spectrum=0; // we use g, m, yr units following Hartvig et al.
double phyto_u_min=log(1e-15);  
double phyto_u_max=log(1e-8);   

double phyto_zoo_transition_width=log(1e1);
double phyto_damping_rate=10;
double phyto_boost=1;

#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGINT(balance_case),
    CFGDOUBLE(alpha),
    CFGDOUBLE(h),
    CFGDOUBLE(k),
    CFGDOUBLE(beta),
    CFGDOUBLE(w),
    CFGDOUBLE(search_gamma),
    CFGDOUBLE(gamma_boost),
    CFGINT(gamma_as_pp_boost),
    CFGDOUBLE(q),
    CFGDOUBLE(eta4),
    CFGDOUBLE(sigma),
    CFGDOUBLE(mesh_size),
    CFGDOUBLE(absoluteF),
    CFGDOUBLE(relativeF),
    CFGDOUBLE(lambda),
    CFGDOUBLE(rho),
    CFGDOUBLE(x),
    CFGDOUBLE(wD),
    CFGDOUBLE(bRange),
    CFGDOUBLE(sigma),
    CFGINT(max_zeros),
    CFGINT(doSimulate),
    CFGINT(maturation_schedule),
    CFGDOUBLE(x0),
    CFGINT(noSpecies),
    CFGINT(compute_linear_mode),
    CFGINT(communitySpectrum),
    CFGDOUBLE(pulse_perturbation),
    CFGDOUBLE(stoppingTime),
    CFGDOUBLE(upper_cutoff),
    CFGINT(nonlinear_version),
    CFGINT(phyto_spectrum),
    CFGDOUBLE(phyto_u_min),
    CFGDOUBLE(phyto_u_max),
    CFGDOUBLE(phyto_zoo_transition_width),
    CFGDOUBLE(phyto_damping_rate),
    CFGDOUBLE(phyto_boost),
    CFGINT(show_zero_movie),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);

//dependent constants:
double I10,I20,phiTilde,NTilde,NTilde_crit,muTilde,f0;
double maturation_time,tildeBtot,time_unit;
double lower_tail_exponent,upper_tail_exponent;
double F;

char * gif_image_file=0;
char * output_file=0;
char * MERP_output_file=0;
char * initial_state_file=0;

// defined later:
void write_poles(ostream & plane);

class somePlotter{
public:
  somePlotter(){};
  virtual 
  void prepare(double min_u,double max_u, double time)=0;
  virtual
  void point(double x,double y)=0;
  virtual 
  double height()=0;
};

template< class PLOTTER = XPlotter >
class myPlotter : public PLOTTER, public somePlotter {
  bool _active;
  double _nominalHeight;
  FILE *_outfile;
public:
  myPlotter(bool active=true);
  myPlotter(const char * filename);
  ~myPlotter();
  void prepare(double min_u,double max_u, double time);
  void point(double x,double y);
  double height(){return _nominalHeight;};
};

template< class PLOTTER >
 myPlotter< PLOTTER >::myPlotter(bool active): 
   _active(active), PLOTTER(cin,cout,cerr) 
 {
   if(_active){
     if (PLOTTER::openpl () < 0)          // open Plotter
       FATAL_ERROR("could not open plotter");
   }
 }

template< class PLOTTER >
 myPlotter< PLOTTER >::myPlotter(const char * filename): 
   _active(filename),
   _outfile(filename ? fopen(filename,"w") : 0),
   //PLOTTER(stdin,_outfile,stderr) 
   PLOTTER(stdin,stderr,stdout) 
 {
   if(_active){
     if(_outfile == NULL){
       FATAL_ERROR("Could not open " << filename << " for writing");
     }
     if (PLOTTER::openpl () < 0)          // open Plotter
       FATAL_ERROR("could not open plotter");
   }
 }

template< class PLOTTER >
 myPlotter< PLOTTER >::~myPlotter()
 {
   if(_outfile){
     fclose(_outfile);
   };
 }

template< class PLOTTER >
void myPlotter< PLOTTER >::prepare(double min_u,double max_u, double time){
  if(!_active) return;
  _nominalHeight=(max_u-min_u)/M_LN10 * (double(BITMAPHEIGHT)/BITMAPWIDTH);
  PLOTTER::fspace(min_u/M_LN10,-_nominalHeight/2,
		  max_u/M_LN10,+_nominalHeight/2);
  PLOTTER::erase();
  PLOTTER::color(0,0,0);
  
  char time_string[100];
  PLOTTER::fontname("HersheySans");
  PLOTTER::fontsize(2);
  sprintf(time_string,"%.3f",time/time_unit);
  PLOTTER::fmove(max_u/log(10.0),_nominalHeight/2);
  PLOTTER::alabel ('r','t',time_string);

  if(mesh_size){
    PLOTTER::
      fline(log(mesh_size)/M_LN10,-_nominalHeight,
	    log(mesh_size)/M_LN10,+_nominalHeight/2);
    PLOTTER::endpath();
  }else{
    PLOTTER::
      fline(+2*sigma/M_LN10,-_nominalHeight,+2*sigma/M_LN10,_nominalHeight/2);
    PLOTTER::endpath();
    PLOTTER::
      fline(-2*sigma/M_LN10,-_nominalHeight,-2*sigma/M_LN10,_nominalHeight/2);
    PLOTTER::endpath();
  }

  const double tick_dist=1;
  const double tick_length=1;

  double lower_xtick = -int(-min_u/M_LN10/tick_dist) * tick_dist;
  double lower_ytick = -int(_nominalHeight/2/tick_dist) * tick_dist;
  
  for(double tick=lower_xtick;tick < max_u/M_LN10;tick+=tick_dist){
    PLOTTER::
      fline(tick,-_nominalHeight/2,tick,-_nominalHeight/2+tick_length);
    PLOTTER::endpath();
    PLOTTER::
      fline(tick,+_nominalHeight/2,tick,+_nominalHeight/2-tick_length);
    PLOTTER::endpath();
  }
  for(double tick=lower_ytick;tick < max_u/M_LN10;tick+=tick_dist){
    PLOTTER::
      fline(min_u/M_LN10,tick,min_u/M_LN10+tick_length,tick);
    PLOTTER::endpath();
    PLOTTER::
      fline(max_u/M_LN10,tick,max_u/M_LN10-tick_length,tick);
    PLOTTER::endpath();
  }
  
  PLOTTER::fmove(min_u/log(10.0),0);
}

template< class PLOTTER >
void myPlotter< PLOTTER >::point(double x,double y){
  if(!_active) return;
  PLOTTER::fcont(x/M_LN10,y/M_LN10);
}    

class fft_field {
  int _N;
  cmplx * _work;
  fftw_plan _fftw_forward_plan;
  fftw_plan _fftw_backward_plan;
public:
  fft_field(unsigned int N);
  ~fft_field();
  void fft_forward(){fftw_execute(_fftw_forward_plan);};
  void fft_backward(){fftw_execute(_fftw_backward_plan);};
  cmplx & operator[](unsigned int i) {
    return _work[i];
  }
  const cmplx & operator[](unsigned int i) const{
    return _work[i];
  }
  operator double * () const{
    return (double *)_work;
  }
  double & re(unsigned int i) const{
    return ((double *)_work)[i];
  }
};

fft_field::fft_field(unsigned int N):_N(N){
  _work = (cmplx*) fftw_malloc(sizeof(fftw_complex) * (_N+1));
  _fftw_forward_plan=
    fftw_plan_dft_r2c_1d(2*_N, (double *)_work, (fftw_complex *)_work, 
			 FFTW_ESTIMATE);
  _fftw_backward_plan=
    fftw_plan_dft_c2r_1d(2*_N, (fftw_complex *)_work, (double *)_work, 
			 FFTW_ESTIMATE);
}

fft_field::~fft_field(){
  fftw_destroy_plan(_fftw_forward_plan);
  fftw_destroy_plan(_fftw_backward_plan);
  fftw_free(_work);
}

class SizeSpectrum : public ODE_dynamical_object
{
  int _N;
  double _delta_u;
  double _delta_xi;
  double _min_u;
  vector<double> exp_minus_u_over_4,fishing,damper,equilibriumB;
  vector<double> phyto_zoo_mixer;
  vector<double> equilibrium_food_availability,h_mn,mn1;
  vector<double> search_gamma_mq,search_gamma_mq1;
  vector<cmplx> normalizing_kernel;
  vector<cmplx> normalizing_beta_hat;
  vector<cmplx> normalizing_exclusionHat,precomputed_sHat;
  vector<cmplx> normalizing_predation_window;
  vector<double> _b;
  fft_field _work,_work_X,_work_f,_work_mu,_work_cc;
public:
  double u(int i) const {return _delta_u*i+_min_u;}
  int index_of(double u) const {return int((u-_min_u)/_delta_u + 0.5);}
  double centered_u(int i) const{
    return _delta_u*(i < _N ? i : i-2*_N );
  }
  double xi(int i) const {
    return _delta_xi*i;
  }
  double min_u() const {return u(0);}
  double max_u() const {return u(_N-1);}
private:
  void initialize();
  void initialize_dynamics();
public:
  void to_xmgr(ostream &os, const double *x, double exponent=(-lambda+2)) const;
public:
  SizeSpectrum(int nPoints=0,double delta_u=1);
  SizeSpectrum(const SizeSpectrum & other);
  ~SizeSpectrum();
  virtual int dynamics(ODE_vector const & state, 
		       ODE_vector & time_derivative);
  virtual void precondition(ODE_vector const & state,
  			    ODE_vector const & in,
  			    ODE_vector & out,
  			    realtype gamma,
  			    bool left_rather_than_right);
  virtual bool has_preconditioner(){return true;};
  virtual void write_state_to(ODE_vector & state) const;
  virtual void read_state_from(const ODE_vector & state);
  virtual int number_of_variables() const;
  virtual void line_print(ODE_vector const & state,std::ostream &co);

  double individual_trophic_levels();
  void to_xmgr(ostream& os) const;
  void to_xmgr(const char * filename) const;
  void save_state(const char * filename) const;
  void load_state(const char * filename);
  void MERP_write_state(ostream& os);
  void MERP_report_size_classes(ostream & os=std::cout);
  void plot(somePlotter & p, const vector<double> & field);
  void plot(somePlotter & p);
  double plot_individuals(somePlotter & p); //returns system biomass
  vector<double> individual_spectrum(); 
  vector<double> species_spectrum(); 
  vector<double> l2_normalize(const vector<double> & b) const; 
  void mirror_the_kernel();
  void save_beta_tilde(const char * filename);
  void save_PPMR_window(const char * filename);
  void save_PPMR_mechanism(cmplx xi, const char * filename);
  void no_more_fishing();
  void no_more_damping();
  // extremely bad style: put this here because currently _b is
  // encapsulated, so adjoint_eigenfunction_extractor has not access:
  void adjoint_correct(){
    for(int i=_N;i-->0;){
      _b[i]*=exp(u(i)/4);
    }
  }
  double effective_density_dependence(const ODE_state & state);
};


// use to find linear mode with growth rate -1
class eigenfunction_extractor : public SizeSpectrum
{
  int _fixed_component;  // Fix component with this index to 1 to
 			 // ensure non-trivial normalization.
  double _eigenvalue;
public:
  eigenfunction_extractor(int nPoints=0,double delta_u=1,double eig_val=-1);
  virtual int dynamics(ODE_vector const & state, 
		       ODE_vector & time_derivative);
  virtual void read_state_from(const ODE_vector & state);
  eigenfunction_extractor compute_residuum();
};

// use to find linear mode with growth rate -1
class adjoint_eigenfunction_extractor : public eigenfunction_extractor
{
public:
  adjoint_eigenfunction_extractor(int nPoints=0,double delta_u=1,double eig_val=-1);
};

adjoint_eigenfunction_extractor::
adjoint_eigenfunction_extractor(int nPoints,double delta_u,double eig_val):
  eigenfunction_extractor(nPoints,delta_u,eig_val){
  mirror_the_kernel();
}

void eigenfunction_extractor::
read_state_from(const ODE_vector & state){
  ODE_vector modified_state(state);
  modified_state[_fixed_component]=1;
  SizeSpectrum::read_state_from(modified_state);
}

eigenfunction_extractor::
eigenfunction_extractor(int nPoints,double delta_u,double eig_val):
  _eigenvalue(eig_val),
  SizeSpectrum(nPoints,delta_u)
{
  no_more_fishing(); // set inhomogeneity to zero

  // determine the index of the _fixed_component:
  for(int i=0;i<number_of_variables();i++){
    if(u(i)*u(i-1)<=0){
      _fixed_component=i;
      break;
    }
  }
}

int eigenfunction_extractor::
dynamics(ODE_vector const & const_state, ODE_vector & time_derivative){

  // Remove the const.  The change we apply is only temporary
  ODE_vector & state=const_cast<ODE_vector &>(const_state);
  double hold=state[_fixed_component];
  state[_fixed_component]=1;

  // Set up eigenvalue equation for growth rate -1:
  if(SizeSpectrum::dynamics(state,time_derivative)){
    return 1;
  }
  for(int i=number_of_variables();i-->0;){
    time_derivative[i]-=state[i]*_eigenvalue;
  }

  // Tell solver that fixed component has value 1.
  int last=0;
  time_derivative[last]+=(1-hold)*exp(-(1-n)*u(last));
  
  state[_fixed_component]=hold;

  return 0;
}


eigenfunction_extractor eigenfunction_extractor::
compute_residuum(){
  ODE_vector state(number_of_variables());
  ODE_vector residuum(number_of_variables());

  write_state_to(state);
  dynamics(state,residuum);
  REPORT(residuum[_fixed_component]);
  residuum[0]=0;
  
  eigenfunction_extractor resi(*this);
  
  resi.SizeSpectrum::read_state_from(residuum);
  return resi;
}

vector<double> SizeSpectrum::l2_normalize(const vector<double> & b) const{
  
  double sum2=0;
  for(int i=b.size();i-->0;){
    sum2+=b[i]*b[i];
  }
  double l2_norm=sqrt(sum2*(u(1)-u(0)));

  vector<double> normalized(b.size());
  for(int i=number_of_variables();i-->0;){
    normalized[i]=b[i]*(1/l2_norm);
  }

  return normalized;
}

void SizeSpectrum::precondition(ODE_vector const & state,
		  ODE_vector const & in,
		  ODE_vector & out,
		  realtype gamma,
		  bool left_rather_than_right){
  if(!left_rather_than_right){
    for(int i=_N;i-->0;){
      out[i]=in[i]/(1+gamma*exp_minus_u_over_4[i]);
    }
  }else{
    for(int i=_N;i-->0;){
      out[i]=in[i];
    }
  }
}

void SizeSpectrum::initialize(){
  initialize_dynamics();
}

SizeSpectrum::SizeSpectrum(int nPoints,double delta_u) : 
  _N(nPoints), 
  _delta_u(delta_u), 
  _delta_xi(2*M_PI/(2*_N*delta_u)),
  exp_minus_u_over_4(_N),
  fishing(_N),
  damper(_N),
  phyto_zoo_mixer(phyto_spectrum ? _N : 0),
  normalizing_kernel(_N+1),
  normalizing_beta_hat(_N+1),
  _b(_N),
  _work(_N),
  _work_mu(phyto_spectrum || nonlinear_version >= 3 ? _N : 0),
  _work_cc(phyto_spectrum || nonlinear_version >= 3 ? _N : 0),
  _work_f(nonlinear_version >= 3 ? _N : 0),
  _work_X(nonlinear_version >= 3 ? _N : 0)
{
  initialize();
}

SizeSpectrum::SizeSpectrum(const SizeSpectrum & other) :
  _N(other._N), 
  _delta_u(other._delta_u), 
  _delta_xi(other._delta_xi),
  exp_minus_u_over_4(_N),
  fishing(_N),
  damper(_N),
  phyto_zoo_mixer(phyto_spectrum ? _N : 0),
  normalizing_kernel(_N+1),
  normalizing_beta_hat(_N+1),
  _b(other._b),
  _work(_N),
  _work_mu(phyto_spectrum || nonlinear_version >= 3 ? _N : 0),
  _work_cc(phyto_spectrum || nonlinear_version >= 3 ? _N : 0),
  _work_f(nonlinear_version >= 3 ? _N : 0),
  _work_X(nonlinear_version >= 3 ? _N : 0)
{
  initialize();
}

SizeSpectrum::~SizeSpectrum(){
}

int SizeSpectrum::number_of_variables() const{
  return _N;
}

void SizeSpectrum::write_state_to(ODE_vector & state) const{
  for(int i=_N;i-->0;){
    state[i]=_b[i];
  }
}

void SizeSpectrum::read_state_from(const ODE_vector & state){
  for(int i=_N;i-->0;){
    _b[i]=state[i];
  }
}  

cmplx sHat(cmplx zeta){
  return sqrt(2*M_PI) * w * exp(-((w*w)/2)*(zeta*zeta)) * 
    pow(beta,cmplx(0,-1)*zeta);
}

cmplx kappa0Hat_being_eaten(cmplx zeta)
{
  if(balance_case==1){
    return 
      h/phiTilde
      *(-conj(sHat(conj(zeta)))
	+
	((1/I10)*
	 sHat(zeta-cmplx(0,3-lambda-n))*
	 conj(sHat(conj(zeta)))
	 )
	);
  }else{
    return 
      -
      x*search_gamma*h/
      (x*search_gamma*phiTilde+h)*
      conj(sHat(conj(zeta)))
      +
      x*search_gamma*x*search_gamma*NTilde*h/
      (x*search_gamma*phiTilde+h)/
      (x*search_gamma*phiTilde+h)*
      sHat(zeta-cmplx(0,3-lambda-n))*
      conj(sHat(conj(zeta)));
  }
}

cmplx kappa0Hat_eating(cmplx zeta)
{
  if(balance_case==1){
    return double(0);
  }else{
    return 
      x*alpha*h*h*search_gamma/
      (x*search_gamma*phiTilde+h)/
      (x*search_gamma*phiTilde+h)*
      sHat(zeta-cmplx(0,3-lambda-n));
  }
}
  

cmplx kappa0Hat(cmplx zeta){
  return kappa0Hat_eating(zeta)+kappa0Hat_being_eaten(zeta);
}

// (betaHat[zeta] /. E^x_ -> exp[x] //. aa_^b_ -> pow[aa, b] // Evaluate // 
//                     Hold) /. pow[x_, 2] :> x*x /. pow[x_, 1/2] -> sqrt[x] /. 
//             pow[x_, -2] :> 1/(x*x) /. pow[x_, -1] :> 1/x /. 
//         pow[x_, 3] :> x*x*x /. pow[x_, 4] :> (x*x)*(x*x) /. 
//     Complex[xx_, yy_] :> cmplx[xx, yy] // CForm

cmplx betaHat(cmplx zeta)
{ 
  switch(maturation_schedule){
  case 1:
    return 
      3*pow(eta4,cmplx(0,-4)*zeta)/
      (3 + cmplx(0,-25)*zeta - 
       70*(zeta*zeta) + 
       cmplx(0,80)*
       (zeta*zeta*zeta) + 
       32*(zeta*zeta*(zeta*zeta))
       );
  case 2:
    return 
      (cmplx(0,5)+4*zeta)/
      (cmplx(0,5)+20*zeta);
  case 3:
    ALWAYS_ASSERT(n==double(3)/4);
    const double eta1=1/eta4;
    const double eta2=eta1*eta1;
    const double eta3=eta2*eta1;
    const double eta=eta2*eta2;
#define SQUARE(X) (X)*(X)
#define CUBE(X)  ((X)*(X))*(X)
#define FOURTH(X) ((X)*(X))*((X)*(X))
#define Power(X,Y) pow(X,Y)
    if(x0 == 0){
      return 
	(3*Power(eta,zeta*cmplx(0,1)) + 2*eta2*(5 + zeta*cmplx(0,-4))*(zeta + cmplx(0,1))*
	 (4*zeta + cmplx(0,3)) + Power(eta,1.25)*(zeta + cmplx(0,1))*
	 (2*zeta + cmplx(0,1))*(4*zeta + cmplx(0,3))*cmplx(0,4) + 
	 3*eta*(3 + zeta*cmplx(0,-4))*(2*zeta + cmplx(0,1))*(4*zeta + cmplx(0,5)) + 
	 eta3*(zeta + cmplx(0,1))*(2*zeta + cmplx(0,1))*(4*zeta + cmplx(0,5))*
	 cmplx(0,12))/
	(Power(-1 + eta1,4)*(1 + 4*eta1)*(zeta + cmplx(0,1))*(2*zeta + cmplx(0,1))*
	 (4*zeta + cmplx(0,1))*(4*zeta + cmplx(0,3)));
    }else{// x0 > 0
      return
	-((Power(x0,zeta*cmplx(0,-1))*
	   (Power(x0,zeta*cmplx(0,1))*
	    (3*Power(eta,zeta*cmplx(0,1)) + 
	     2*eta2*(5 + zeta*cmplx(0,-4))*(zeta + cmplx(0,1))*
	     (4*zeta + cmplx(0,3)) + 
	     Power(eta,1.25)*(zeta + cmplx(0,1))*(2*zeta + cmplx(0,1))*
	     (4*zeta + cmplx(0,3))*cmplx(0,4) + 
	     3*eta*(3 + zeta*cmplx(0,-4))*(2*zeta + cmplx(0,1))*
	     (4*zeta + cmplx(0,5)) + 
	     eta3*(zeta + cmplx(0,1))*(2*zeta + cmplx(0,1))*(4*zeta + cmplx(0,5))*
	     cmplx(0,12)) + 4*eta1*Power(x0,0.25)*(3 + zeta*cmplx(0,-4))*
	    (zeta + cmplx(0,1))*(2*zeta + cmplx(0,1))*FOURTH(-1 + eta1)))/
	  (Power(-1 + eta1,4)*(-1 + 4*eta1*(-1 + Power(x0,0.25)))*(zeta + cmplx(0,1))*
	   (2*zeta + cmplx(0,1))*(4*zeta + cmplx(0,1))*(4*zeta + cmplx(0,3))));
    }
  }
#undef SQUARE
#undef CUBE
#undef FOURTH
#undef Power
}

cmplx PhiHat(cmplx zeta){
  cmplx sz=sigma*zeta;
  return -F*sqrt(2*M_PI)*sigma*exp(-sz*sz/double(2));
}

cmplx exclusionHat(cmplx z){
  //  z-=cmplx(0,lambda-2);
  if(wD != 0){
    return 
      (exp(-(wD*z)*(wD*z)/double(2))-exp((wD*(2-lambda))*(wD*(2-lambda))/double(2)))*rho;
  }else{
    FATAL_ERROR("wD==0 not currently implemented");
    //return -D*z*z;
  }
}

cmplx kernel_function(cmplx z){
  if(noSpecies){
    return kappa0Hat(z+cmplx(0,3-lambda-n))+exclusionHat(z)/(conj(betaHat(conj(z+cmplx(0,(3-lambda-n)))))*betaHat(z));
  }else{
    cmplx kappaHat=
      kappa0Hat(z+cmplx(0,3-lambda-n)) * 
      conj(betaHat(conj(z+cmplx(0,(3-lambda-n)))));
    return 
      (tildeBtot*kappaHat * 
       (communitySpectrum ? betaHat(z+cmplx(0,(1-n))) : betaHat(z) ) + 
       exclusionHat(z) );
  }
}

vector<double> poles(){
  vector<double> p;
  if(!noSpecies){
    p.push_back(-(1-n) - ( communitySpectrum ? 1-n : 0 ));
    if(not (balance_case == 1)){
      p.push_back((1-n)-(3-lambda-n));
    }
    switch(maturation_schedule){
    case 1:
      for(double i=1;i>(1-n)*3/2;i-=(1-n)){
	p.push_back(-i - ( communitySpectrum ? 1-n : 0 ) );
	p.push_back( i-(3-lambda-n) );
      }
      break;
    case 2:
    case 3:
      break;
    }
  }
  return p;
}  

void write_poles(ostream & plane){
  vector<double> p=poles();
  for(int i=0;i<p.size();i++){
    plane << 0 << " " << p[i] << endl;
  }
}

double eutro_lambda(){
  return 2-log(beta)/(w*w)+
    sqrt((-1+n)*(-1+n)*w*w*w*w+log(beta)*(2*(-2+n)*w*w+log(beta))-
	 2*w*w*log((alpha*h-k)/(beta*h)) )/
    (w*w);
}

void set_dependent_constants(){
  double I10_eutro=pow(beta,-2+lambda)*
    exp((-2+lambda)*(-2+lambda)*w*w/double(2))*
    sqrt(2*M_PI)*w;
  if(balance_case==2) {
    lambda=2 + q - n; 
    REPORT(lambda);
  }
  I10=pow(beta,-2+lambda)*
    exp((-2+lambda)*(-2+lambda)*w*w/double(2))*
    sqrt(2*M_PI)*w;
  I20=pow(beta,-1+n)*
    exp((-1+n)*(-1+n)*w*w/double(2))*
    sqrt(2*M_PI)*w;
  if(balance_case==1){
    k=h*(alpha-(I20/I10));
    REPORT(k);
  }else if(k<0){
    REPORT(I20/I10_eutro);
    double k_suggest=h*(alpha-(I20/I10_eutro));
    REPORT(k_suggest);
    k=h*(alpha-(I20/I10_eutro));
    REPORT(k);
  }    
  NTilde=k/
    (search_gamma*((alpha-k/h)*I10-I20));
  if(NTilde<0 || balance_case==1){
    WARNING("NTilde SET TO UNITY!");
    NTilde=1;
  }
  phiTilde=
    NTilde*I10;
  NTilde_crit=
    ((h/search_gamma/sHat(cmplx(0,-1+q)))*
     (sHat(cmplx(0,-1+q))-alpha*sHat(double(0)))/
     (sHat(double(0))-sHat(cmplx(0,lambda-2)))
     ).real();
  tildeBtot=NTilde/betaHat(cmplx(0,lambda-2)).real();
  REPORT(tildeBtot);
  const double tildePhi=NTilde*I10;
  REPORT(tildePhi);
  REPORT(search_gamma*tildePhi);
  REPORT(h);
  f0=(balance_case==1 ? 1 :
      search_gamma*tildePhi/(search_gamma*tildePhi+h) );
  REPORT(f0);
  const double tildeG0=alpha*f0*h-k;
  REPORT(tildeG0);
  if(balance_case == 2){
    muTilde = search_gamma*h*NTilde/(search_gamma*tildePhi+h)*sHat(cmplx(0,-1+n)).real();
  }else{//balance_case == 1
    muTilde = (h*NTilde/tildePhi)*sHat(cmplx(0,-1+n)).real();
  }
  maturation_time=1/(tildeG0*(1-n));
  REPORT(maturation_time);
  REPORT(pow(1000,0.25)*maturation_time);
  time_unit=1;//maturation_time*tildeBtot;
  F=relativeF/maturation_time;
  if(noSpecies){
    lower_tail_exponent=HUGE_VAL;
    upper_tail_exponent=-HUGE_VAL;
  }else{ 
    lower_tail_exponent=(1-n) + (communitySpectrum ? (1-n) : 0);
    upper_tail_exponent=2-lambda;
    if(balance_case==1){
      // exclude shadowed pole
      upper_tail_exponent-=(1-n);
    }
  }
  double kens_a=pow(beta,2*n-q-1)*exp((2*n*(q-1)-q*q+1)*w*w/2)/alpha;
  REPORT(kens_a);
  double kens_a_simple=pow(beta,2*n-q-1)/alpha;
  REPORT(kens_a_simple);
}


void SizeSpectrum::initialize_dynamics(){
  WARNING("initialising dynamics");
  double lower_cutoff=(u(0)+1*log(beta));
  if(!phyto_spectrum){
    _min_u = -1e-9 + log(1e-10) - log(beta);// -_N*_delta_u/3;
    lower_cutoff = ( _min_u + 1*log(beta) );
  }else{
    _min_u = phyto_u_min;
    lower_cutoff = phyto_u_max;
  }
  REPORT(lower_cutoff/M_LN10);

  for(int i=_N;i-->0;){
    exp_minus_u_over_4[i]=exp(-u(i)/4);

    if(mesh_size){
      fishing[i]=(u(i)>log(mesh_size) ? -absoluteF : 0);
      //use with "poisoning.cfg"
      //fishing[i]=(u(i)>log(1e-9) && u(i)<log(1e-5) ? -absoluteF : 0);
    }else{
      fishing[i]=-F*exp(-0.5*u(i)*u(i)/(sigma*sigma) + (2-lambda)*u(i));
    }

    if(!phyto_spectrum){
      double ex=exp(-10.6*(u(i)-lower_cutoff));
      damper[i]=-8*ex/(ex+1);
    }else{
      double ex=exp(-(u(i)-lower_cutoff)/phyto_zoo_transition_width);
      phyto_zoo_mixer[i]=ex/(ex+1);
      //cout << exp(u(i)) << " " << phyto_zoo_mixer[i] << endl;
      damper[i]=-phyto_damping_rate*exp(-u(i)/4);
    }
  }
  
  for(int i=_N+1;i-->0;){
    normalizing_kernel[i]=kernel_function(xi(i))*(1/double(2*_N));
    normalizing_beta_hat[i]=betaHat(xi(i))*(1/double(2*_N));
  }

  if(false){
    // convolute fishing pressure with reverse species size structure:
    double * re_work=(double *)_work;
    for(int i=2*_N;i-->_N;){
      re_work[i]=0;
    }
    for(int i=_N;i-->0;){
      re_work[i]=fishing[i];
    }
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      _work[i]*=conj(normalizing_beta_hat[i]);
    }
    _work.fft_backward();
    for(int i=_N;i-->0;){
      fishing[i]=re_work[i];
    }
  }
    

  if(communitySpectrum){
    // convolute fishing pressure with species size structure:
    double * re_work=(double *)_work;
    for(int i=2*_N;i-->_N;){
      re_work[i]=0;
    }
    for(int i=_N;i-->0;){
      re_work[i]=fishing[i];
    }
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      _work[i]*=normalizing_beta_hat[i];
    }
    _work.fft_backward();
    for(int i=_N;i-->0;){
      fishing[i]=re_work[i];
    }
  }
  if(nonlinear_version >= 2){
    equilibriumB.resize(2*_N);
    for(int i=2*_N;i-->0;){
      int ii=( i<3*_N/2 ? i :  i-2*_N);
      equilibriumB[i] = tildeBtot * exp((2-lambda) * u(ii));
    }
    if(nonlinear_version >= 3){
      normalizing_exclusionHat.resize(_N+1);
      precomputed_sHat.resize(_N+1);
      normalizing_predation_window.resize(_N+1);
      for(int i=_N+1;i-->0;){
	normalizing_exclusionHat[i] = exclusionHat(xi(i)) * (1/double(2*_N));
	precomputed_sHat[i] = sHat(xi(i));
	normalizing_predation_window[i] = 
	  -conj(sHat(cmplx(xi(i),n-1))) * (1/double(2*_N)); 
      }
      equilibrium_food_availability.resize(2*_N);
      search_gamma_mq.resize(2*_N);
      search_gamma_mq1.resize(2*_N);
      h_mn.resize(2*_N);
      mn1.resize(2*_N);
      for(int i=2*_N;i-->0;){
	int ii=( i<3*_N/2 ? i :  i-2*_N);
	equilibrium_food_availability[i] = 
	  sHat(cmplx(0,lambda-2)).real()*NTilde*exp(u(ii)*(-lambda+2));
	search_gamma_mq[i]=gamma_boost*search_gamma*exp(u(ii)*q);
	search_gamma_mq1[i]=gamma_boost*search_gamma*exp(u(ii)*(q-1));
	h_mn[i]=h*exp(u(ii)*n);
	mn1[i]=exp(u(ii)*(n-1));  // FIXME: same as exp_minus_u_over_4
      }
    }
    if(gamma_as_pp_boost){
      for(int i=_N;i-->0;){
	_b[i]=-log(gamma_boost)/equilibriumB[i];
      }
    }
  }
}

void SizeSpectrum::mirror_the_kernel(){ 
  for(int i=_N+1;i-->0;){
    normalizing_kernel[i]=conj(normalizing_kernel[i]);
  }
}

int SizeSpectrum::dynamics(ODE_vector const & state, 
	     ODE_vector & time_derivative){
  double * re_work=(double *)_work;
  for(int i=2*_N;i-->_N;){
    re_work[i]=0;
  }
  if(nonlinear_version >= 2){
    for(int i=2*_N;i-->0;){
      int ii=( i<3*_N/2 ? i :  i-2*_N);
      if(ii<0){
	re_work[i] = 0;
      }else if(ii >= index_of(upper_cutoff)){
	re_work[i] = -equilibriumB[i];
      }else{
	re_work[i] = (exp(state[i]/equilibriumB[i])-1)*equilibriumB[i];
      }
    }
    // ofstream pl("tmp.dat");
    // to_xmgr(pl,&re_work[0]);FATAL_ERROR("check tmp.dat!");
  }else{
    for(int i=_N;i-->0;){
      re_work[i]=state[i];
    }
  }
  //***********************
  //for(int i=2*_N;i-->0;){re_work[i]=0;}re_work[index_of(0)]=0.01;
  //***********************
  if( nonlinear_version <= 2 ){
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      _work[i]*=normalizing_kernel[i];
    }
    _work.fft_backward();
  }else{ //fully_nonlinear
    //ofstream pl("tmp.dat");
    ALWAYS_ASSERT(balance_case==2);
    vector<double> phyto_imbalance;
    if(phyto_spectrum){
      double * cc = (double *) _work_cc;
      // Compute consumer community size spectrum:
      for(int i=2*_N;i-->0;){
	cc[i]=0;
      }
      for(int i=index_of(upper_cutoff);i-->phyto_imbalance.size();){
	cc[i]=re_work[i]+equilibriumB[i];
      }
      _work_cc.fft_forward();
      for(int i=_N+1;i-->0;){
	_work_cc[i] *= normalizing_beta_hat[i];
      }
      _work_cc.fft_backward();
      // take out phytoplankton species size specrum.
      phyto_imbalance.resize(index_of(phyto_u_max)+1);
      for(int i=phyto_imbalance.size();i-->0;){
	phyto_imbalance[i]=re_work[i];
	re_work[i]=-equilibriumB[i];
      }
    }    
    // Compute X
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      _work_X[i] = _work[i] * normalizing_exclusionHat[i];
    }
    _work_X.fft_backward();
    double * X = (double *) _work_X;
    // X computed.
    // Compute consumer cHat(u)...
    for(int i=_N+1;i-->0;){
      _work[i] *= normalizing_beta_hat[i];
    }
    if(phyto_spectrum){
      // Add phytoplankton availability
      _work.fft_backward();
      const double betafac = betaHat(cmplx(0,lambda-2)).real();
      for(int i=phyto_imbalance.size();i-->0;){
	re_work[i]+=(equilibriumB[i]+phyto_imbalance[i])*betafac;
      }
      for(int i=2*_N;i-->0;){
	re_work[i] *= (1/double(2*_N));
      }
      _work.fft_forward();
    }
    // Compute off-equilibrium food availability in _work_f
    for(int i=_N+1;i-->0;){
      _work_f[i]=_work[i]*precomputed_sHat[i];
    }
    _work_f.fft_backward();
    // Compute full food availability in _work_f
    for(int i=2*_N;i-->0;){
      _work_f.re(i) += equilibrium_food_availability[i];
    }    
    // Compute feeding level in _work_f
    double * f = (double *)_work_f;
    for(int i=2*_N;i-->0;){
      f[i] =
	(search_gamma_mq[i]*_work_f.re(i)) / 
	(search_gamma_mq[i]*_work_f.re(i) + h_mn[i]);
    }
    //@@ to_xmgr(pl,&f[0],0);FATAL_ERROR("check tmp.dat!");
    if(phyto_spectrum){
      // compute predator activity (/ mn1) in _work
      for(int i=2*_N;i-->0;){
	re_work[i] = (1-f[i]) * search_gamma_mq1[i] * 
	  _work_cc.re(i) / mn1[i];
      }
    }else{
      // _work still contains cHat
      _work.fft_backward();
      // _work now contains c
      const double * const c = (const double *) _work;
      // compute predator activity (/ mn1) in _work
      const double betafac = betaHat(cmplx(0,lambda-2)).real();
      for(int i=2*_N;i-->0;){
	re_work[i] = (1-f[i]) * search_gamma_mq1[i] * 
	  (c[i] + betafac * equilibriumB[i]) / mn1[i];
      }
    }
    //to_xmgr(pl,&re_work[0],n-1);FATAL_ERROR("check tmp.dat!");
    // // Clean upper work range:
    // for(int i=3*_N/2;i-->_N;){
    //   re_work[i]=0;
    // }
    // to_xmgr(pl,&re_work[0],n-1);FATAL_ERROR("check tmp.dat!");
    // Compute predation mortality (* -1) in _work:
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      //_work[i]*=-conj(sHat(cmplx(xi(i),0   ))) * (1/double(2*_N)); 
      _work_mu[i] = _work[i]*normalizing_predation_window[i]; 
      // ... includes normalization
    }
    _work_mu.fft_backward();
    // ... _work_mu now contains predation mortalities 
    // Clean upper re_work range: !! FIXME, SEE TWO BLOCKS BELOW !!
    for(int i=3*_N/2;i-->_N;){
      re_work[i]=0;
    }
    // and compute individual-level proportional biomass gains in re_work:
    for(int i=2*_N;i-->0;){
      re_work[i] = (_work_mu.re(i)+(alpha * f[i]*h - k))*mn1[i];
    }
    //to_xmgr(pl,&re_work[0],0);FATAL_ERROR("check tmp.dat!");
    // Clean upper work range:
    for(int i=3*_N/2;i<3*_N/2+(_N/4);++i){
      re_work[i]=0;
    }
    //to_xmgr(pl,&re_work[0],0);FATAL_ERROR("check tmp.dat!");
    // Compute population-level linear growth rates (of consumers):
    _work.fft_forward();
    for(int i=_N+1;i-->0;){
      _work[i]*=conj(betaHat(cmplx(xi(i),0)))  * (1/double(2*_N)); 
    }
    _work.fft_backward();
    //to_xmgr(pl,&re_work[0],0);FATAL_ERROR("check tmp.dat!");
    for(int i=_N;i-->0;){
      re_work[i] = re_work[i]*equilibriumB[i]/(mn1[i]) + X[i];
    }
    if(phyto_spectrum){// overwrite result for the phytoplankton part
		       // of the spectrum
      
      for(int i=phyto_imbalance.size();i-->0;){
	re_work[i]=
	  -phyto_damping_rate * (phyto_imbalance[i]/equilibriumB[i]+1-phyto_boost)*(1/phyto_boost) // logistic growth
	  + (_work_mu.re(i)+muTilde);
	// REPORT(phyto_damping_rate * phyto_imbalance[i]/equilibriumB[i]);
	// REPORT(_work_mu.re(i)+muTilde);
      }
      // WARNING("++++++++++++++++++++++++++");
    }    
  }
  //******************
  // ofstream pl("nkernel.dat");
  // to_xmgr(pl,&re_work[0],0);FATAL_ERROR("check nkernel.dat!");
  //******************
  if(phyto_spectrum){
    if(nonlinear_version <= 2){
      for(int i=_N;i-->0;){
	time_derivative[i]=
	  exp_minus_u_over_4[i]*(1-phyto_zoo_mixer[i])*re_work[i]+
	  phyto_zoo_mixer[i]*damper[i]*(state[i]+0)+
	  +fishing[i];
      }
    }else{ // nonlinear_version > 2
      for(int i=_N;i-->0;){
	time_derivative[i]=
	  exp_minus_u_over_4[i]*re_work[i]+fishing[i];
      }
    }
  }else{
    for(int i=_N;i-->0;){
      time_derivative[i]=
	exp_minus_u_over_4[i]*(re_work[i]+damper[i]*state[i])
	+fishing[i];
    }
  }
#if 0 // truncate spectrum from above
  for(int i=_N;i-->0;){
    if(u(i)>+2*sigma)
      time_derivative[i]=0;
  }
#endif
  if(nonlinear_version >= 1){
    for(int i=_N;i-->index_of(upper_cutoff);){
      time_derivative[i]=0;
    }
  }  
  // ofstream pl("plots.dat");
  // to_xmgr(pl,&time_derivative[0],0);FATAL_ERROR("check plots.dat!");
  return 0;
}

double SizeSpectrum::individual_trophic_levels(){
  ALWAYS_ASSERT(nonlinear_version == 3);
  ALWAYS_ASSERT(phyto_spectrum == 0);
  double * re_work=(double *)_work;
  double * X_re_work=(double *)_work_X;
  for(int i=2*_N;i-->_N;){
    re_work[i]=0;
  }
  for(int i=2*_N;i-->0;){
    int ii=( i<3*_N/2 ? i :  i-2*_N);
    if(ii<0){
      re_work[i] = 0;
    }else if(ii >= index_of(upper_cutoff)){
      re_work[i] = 0;
    }else{
      re_work[i] = exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    }
  }
  // re_work contains species size spectrum
  //ofstream pl("tmp.dat");
  ALWAYS_ASSERT(balance_case==2);
  vector<double> phyto_imbalance;
  _work.fft_forward();
  for(int i=_N+1;i-->0;){
    _work[i] *= normalizing_beta_hat[i];
  }
  // _work contains FT of individual size spectrum
  // Compute food availability in _work_f
  for(int i=_N+1;i-->0;){
    _work_f[i]=_work[i]*precomputed_sHat[i];
  }
  _work_f.fft_backward();

  // clear _work_X
  for(int i=2*_N;i-->0;){
    X_re_work[i]=1;
  }
  // compute individual spectrum
  _work.fft_backward();
  // now iteratively, compute TL in _work_X
  double old_max_level,max_level=0;
  do{
    for(int i=2*_N;i-->0;){
      X_re_work[i]*=_work.re(i);
    }
    _work_X.fft_forward();
    for(int i=_N+1;i-->0;){
      _work_X[i] *= precomputed_sHat[i] * (1/double(2*_N)); 
    }
    _work_X.fft_backward();
    old_max_level=max_level;
    max_level=0;
    for(int i=2*_N;i-->_N;){
      X_re_work[i]=0;
    }
    const int producer_start=index_of(min_u()+log(beta));
    for(int i=_N;i-->producer_start;){
      X_re_work[i]=1 + X_re_work[i]/_work_f.re(i);
      if(X_re_work[i] > max_level){
	max_level=X_re_work[i];
      }
    }
    for(int i=producer_start;i-->0;){
      X_re_work[i]=1;
    }
    REPORT(max_level);
  }while(abs(old_max_level-max_level)>0.0001);
  ofstream TL("TL.dat");
  to_xmgr(TL,X_re_work,0);
  return max_level;
}

double SizeSpectrum::
effective_density_dependence(const ODE_state & state){
  ODE_vector derivative(state.size());
  dynamics(state,derivative);
  const int i0=index_of(0);
  return -(derivative[i0]-fishing[i0])/state[i0];
}

void SizeSpectrum::no_more_fishing(){
  for(int i=2*_N;i-->0;){
    fishing[i]=0;
  }
}

void SizeSpectrum::no_more_damping(){
  for(int i=2*_N;i-->0;){
    damper[i]=0;
  }
}

void SizeSpectrum::save_PPMR_window(const char * filename){
  ofstream PPMR(filename);
  double * re_work=(double *)_work;

  for(int i=_N+1;i-->0;){
    cmplx z=xi(i);
    _work[i]=kappa0Hat_being_eaten(z+cmplx(0,3-lambda-n)) * 
      conj(betaHat(conj(z) + cmplx(0,-(3-lambda-n)))) * betaHat(z)
      /double(2*_N);
  }
  _work.fft_backward();
  for(int i=2*_N;i-->0;){
    PPMR << -centered_u(i)/M_LN10 << " " << re_work[i] << endl;
  }
  PPMR << endl;

  for(int i=_N+1;i-->0;){
    cmplx z=xi(i);
    _work[i]=kappa0Hat_eating(z+cmplx(0,3-lambda-n)) * 
      conj(betaHat(conj(z) + cmplx(0,-(3-lambda-n)))) * betaHat(z)
      /double(2*_N);
  }
  _work.fft_backward();
  for(int i=2*_N;i-->0;){
    PPMR << -centered_u(i)/M_LN10 << " " << re_work[i] << endl;
  }
  PPMR << endl;

  for(int i=_N+1;i-->0;){
    cmplx z=xi(i);
    _work[i]=kappa0Hat(z+cmplx(0,3-lambda-n)) * 
      conj(betaHat(conj(z) + cmplx(0,-(3-lambda-n)))) * betaHat(z)
      /double(2*_N);
  }
  _work.fft_backward();
  for(int i=2*_N;i-->0;){
    PPMR << -centered_u(i)/M_LN10 << " " << re_work[i] << endl;
  }
  PPMR << endl;

  for(int i=_N+1;i-->0;){
    _work[i]=kappa0Hat(xi(i))/double(2*_N);
  }
  _work.fft_backward();
  for(int i=2*_N;i-->0;){
    PPMR << -centered_u(i)/M_LN10 << " " << re_work[i] << endl;
  }
  PPMR << endl;

  for(int i=_N+1;i-->0;){
    cmplx z=xi(i);
    _work[i]=(
	      kappa0Hat(z+cmplx(0,3-lambda-n)) * 
	      conj(betaHat(conj(z+cmplx(0,(3-lambda-n)))))*betaHat(z)
	      +
	      exp(-(wD*(z-cmplx(0,lambda-2)))*(wD*(z-cmplx(0,lambda-2)))/double(2))
	      )
      /double(2*_N);
    _work[i]=normalizing_kernel[i];
  }
  _work.fft_backward();
  for(int i=2*_N;i-->0;){
    PPMR << -centered_u(i)/M_LN10 << " " << re_work[i] << endl;
  }
  PPMR << endl;
  double initial_density_dependence = re_work[0];
  REPORT(initial_density_dependence);
  REPORT(initial_density_dependence*time_unit);
}

void SizeSpectrum::save_PPMR_mechanism(cmplx zeta, const char * filename){
  ofstream PPMR(filename);
  double * re_work=(double *)_work;
  
  for(int i=_N+1;i-->0;){
    cmplx z=xi(i);
    _work[i]=//kappa0Hat(z+cmplx(0,3-lambda-n))/
      kernel_function(z)/
      //(betaHat(z)*conj(betaHat(conj(z+cmplx(0,(3-lambda-n))))))/
      double(2*_N);
  }
  _work.fft_backward();
  for(int ii=_N;ii-->-_N;){
    int i=(ii+2*_N)%(2*_N);
    PPMR << -centered_u(i) << " " << re_work[i] << endl;
  }
  PPMR << endl;
  for(int ii=_N;ii-->-_N;){
    int i=(ii+2*_N)%(2*_N);
    PPMR << -centered_u(i) << " " 
	 << (re_work[i]*exp(-cmplx(0,1)*zeta*centered_u(i))).real() << endl;
  }
  PPMR << endl;
  for(int ii=_N;ii-->-_N;){
    int i=(ii+2*_N)%(2*_N);
    PPMR << -centered_u(i) << " " 
	 << (re_work[i]*exp(-cmplx(0,1)*zeta*centered_u(i))).imag() << endl;
  }
  PPMR << endl;
}

void SizeSpectrum::save_beta_tilde(const char * filename){
  ofstream BETA(filename); 
  double * re_work=(double *)_work;

  for(int i=2*_N;i-->1;){
    re_work[i]=0;
  }
  re_work[0]=1/_delta_u;
  _work.fft_forward();
  for(int i=_N+1;i-->0;){
    _work[i]*=normalizing_beta_hat[i];
  }
  _work.fft_backward();
  for(int ii=_N;ii-->-_N;){
    int i=(ii+2*_N)%(2*_N);
    BETA << centered_u(i) << " " 
	 << re_work[i] << endl;
  }
  BETA << endl;
}

void SizeSpectrum::line_print(ODE_vector const & state,std::ostream &co){
  for(int i=0;i<_N;i++){
    co << state[(i+_N)%_N] << " ";
  }
}


void SizeSpectrum::save_state(const char * filename) const{
  ofstream os(filename);
  for(int i=0;i<number_of_variables();i++){
    os << u(i) << " " << _b[i] << endl;
  }
}

void SizeSpectrum::load_state(const char * filename){
  ifstream is(filename);
  double uu;
  for(int i=0;i<number_of_variables();i++){
    is >> uu;
    is >> _b[i];
    if(gamma_as_pp_boost){
      _b[i]-=log(gamma_boost)*equilibriumB[i];
    }
  }
}

void SizeSpectrum::to_xmgr(const char * filename) const{
  ofstream os(filename);
  to_xmgr(os);
}

void SizeSpectrum::to_xmgr(ostream &os, const double *x, double exponent) const{
  for(int ii=-_N/2;ii<3*_N/2;ii++){
    int i= (ii+2*_N)%(2*_N);
    os << u(ii)/M_LN10 << " " 
       << x[i]*exp(-exponent*u(ii)) << endl;
  }
  os << endl;
}

void SizeSpectrum::to_xmgr(ostream &os) const{
  to_xmgr(os,&_b[0]);
}

void SizeSpectrum::plot(somePlotter & p, const vector<double> & field){
  p.prepare(min_u(),max_u(),current_time);
  for(int i=0;i<_N;i++){
    int ii=(i+_N)%_N;
    if(field[ii]> -5*p.height() && !isnan(field[ii])){
      p.point(u(ii),field[ii]);
    }else{
      p.point(u(ii),-5*p.height());
    }
  }
}

void SizeSpectrum::plot(somePlotter & p){
  vector<double> scaled_b(_b.size());
  transform(_b.begin(), _b.end(), scaled_b.begin(),
  	    std::bind1st(std::multiplies<double>(),p.height()/bRange*M_LN10/2));
  plot(p,scaled_b);
}

void SizeSpectrum::MERP_write_state(ostream& os) {

  // convolute b with beta
  double * re_work=(double *)_work;
  for(int i=2*_N;i-->_N;){
    re_work[i]=0;
  }
  if(nonlinear_version >= 2){
    for(int i=_N;i-->index_of(upper_cutoff);){
      re_work[i]=0;
    }
    for(int i=index_of(upper_cutoff);i-->0;){
      re_work[i]=exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    }
  }else{
    for(int i=_N;i-->0;){
      re_work[i]=_b[i];
    }
  }
#if 0 // if 1 output community size spectrum
  _work.fft_forward();
  for(int i=_N+1;i-->0;){
    _work[i]*=normalizing_beta_hat[i];
  }
  _work.fft_backward();
#endif

  vector<double> c(_N);

  if(nonlinear_version >= 2){
    int i=0;
    for(; i<_N; ++i){
      if(gamma_as_pp_boost){
	c[i]=re_work[i]*gamma_boost;
      }else{	
	c[i]=re_work[i];
      }
      if(isnan(c[i]) || c[i]<0) break;
    }
    for(; i<_N; ++i){
      c[i] = 0;
    }
  }else{
    for(int i=_N;i-->0;){
      c[i] = re_work[i];
    }
  }

  os << current_time;

  int n_size_classes=(upper_cutoff-min_u())/M_LN10;
  for(int i=0;i<n_size_classes;i++){
    double Bsum=0;
    for(int j=index_of(min_u()+M_LN10*i);j<index_of(min_u()+M_LN10*(i+1));++j){
      Bsum+=c[j];
    }
    os << "," << Bsum;
  }
  os << endl;    
}

void SizeSpectrum::MERP_report_size_classes(ostream & os){
  int n_size_classes=(upper_cutoff-min_u())/M_LN10;
  for(int i=0;i<n_size_classes;i++){
    os << (i+1) << " " 
       << pow(eta4,4)*exp(u(index_of(min_u()+M_LN10*(i+1))))*1e-3 << endl;
  }
}

double SizeSpectrum::plot_individuals(somePlotter & p){

  double LSI;

  vector<double> c=individual_spectrum();
  
  plot(p,c);

  if(nonlinear_version >= 2){
    // compute LSI
    double delta_S=0,delta_L=0;
    double B_S=0,B_L=0;
    for(int i=index_of(-log(1000));i<index_of(-0*log(2));i++){
      B_S+=exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    }
    for(int i=index_of(-log(2));i<index_of(upper_cutoff);i++){
      B_L+=exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    } 
    LSI=B_L/(B_L+B_S);

  }
  
  return LSI;
}

vector<double> SizeSpectrum::individual_spectrum(){

  // convolute b with beta
  double * re_work=(double *)_work;
  for(int i=2*_N;i-->_N;){
    re_work[i]=0;
  }
  if(nonlinear_version >= 2){
    for(int i=_N;i-->index_of(upper_cutoff);){
      re_work[i]=0;
    }
    for(int i=index_of(upper_cutoff);i-->0;){
      re_work[i]=exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    }
  }else{
    for(int i=_N;i-->0;){
      re_work[i]=_b[i];
    }
  }
  _work.fft_forward();
  for(int i=_N+1;i-->0;){
    _work[i]*=normalizing_beta_hat[i];
  }
  _work.fft_backward();

  vector<double> c(_N);

  if(nonlinear_version >= 2){
    int i=0;
    for(; i<_N; ++i){
      if(gamma_as_pp_boost){
	c[i]=log(re_work[i]*(gamma_boost/NTilde));
      }else{	
	c[i]=log(re_work[i]/NTilde);
      }
      if(isnan(c[i])){
	c[i-1]= - HUGE_VAL;
	c[i]= - HUGE_VAL;
      }
    }
    for(; i<_N; ++i){
      c[i] = -HUGE_VAL;
    }
  }else{
    for(int i=_N;i-->0;){
      c[i] = re_work[i];
    }
  }
  return c;
}

vector<double> SizeSpectrum::species_spectrum(){

  double * re_work=(double *)_work;
  for(int i=2*_N;i-->_N;){
    re_work[i]=0;
  }
  if(nonlinear_version >= 2){
    for(int i=_N;i-->index_of(upper_cutoff);){
      re_work[i]=0;
    }
    for(int i=index_of(upper_cutoff);i-->0;){
      re_work[i]=exp(_b[i]/equilibriumB[i])*equilibriumB[i];
    }
  }else{
    for(int i=_N;i-->0;){
      re_work[i]=_b[i];
    }
  }

  vector<double> c(_N);

  if(nonlinear_version >= 2){
    int i=0;
    for(; i<_N; ++i){
      if(gamma_as_pp_boost){
	c[i]=log(re_work[i]*(gamma_boost/NTilde));
      }else{	
	c[i]=log(re_work[i]/NTilde);
      }
      if(isnan(c[i])) break;
    }
    for(; i<_N; ++i){
      c[i] = -HUGE_VAL;
    }
  }else{
    for(int i=_N;i-->0;){
      c[i] = re_work[i];
    }
  }
  return c;
}

int number_of_CPUs(){
  int number_of_cpus;
#ifdef _SC_NPROCESSORS_ONLN
  number_of_cpus = sysconf(_SC_NPROCESSORS_ONLN); 
#elif defined(HW_NCPU)
  {
    int mib[2];
    size_t len;

    mib[0] = CTL_HW;
    mib[1] = HW_NCPU;
    len = sizeof(number_of_cpus);
    sysctl(mib, 2, &number_of_cpus, &len, NULL, 0);
  } 
#else
  WARNING("Could not detect number of CPUs.");
  number_of_cpus = 1;
#endif
  REPORT(number_of_cpus);
  return number_of_cpus;
}

typedef cmplx analytic(cmplx); 

class analytic_root_finder_failed {};

cmplx analytic_root(analytic f,cmplx start,double epsilon=1e-10){
  // Implements secant root finder in complex plane, assuming f
  // analytic.
  cmplx z1=start,z2=start+cmplx(0,log(beta)/10);
  cmplx f1=f(z1),f2=f(z2);
  int iterations_left=30;
  while(abs(f1)>epsilon && iterations_left-- > 0  && !isnan(f1.real())){
    cmplx z3=f1*(z1-z2)/(f2-f1)+z1;
    if(abs(z3-z1) < abs(z3-z2)){
      z2=z3;
      f2=f(z2);
    } else {
      z1=z3;
      f1=f(z1);
    }
  }
  if(iterations_left <= 0 or isnan(f1.real())){
    //WARNING("Root finder did not converge");
    throw analytic_root_finder_failed();
  }
  return z1;
}

cmplx analytic_derivative(analytic f,cmplx z,double epsilon=1e-8){
  // Numerical derivative of f at z:
  return (f(z+epsilon)-f(z-epsilon))/(2*epsilon);
}

cmplx analytic_second_derivative(analytic f,cmplx z,double epsilon=1e-8){
  // Numerical derivative of f at z:
  return (f(z+epsilon)+f(z-epsilon)-double(2)*f(z))/(epsilon*epsilon);
}

class TheoreticalSizeSpectrum {
  static const int n_zeros;
  vector<cmplx> zero;
  vector<cmplx> v;
  vector<cmplx> D;
  double umin,umax;
  erf_implementation cerf;
  erfc_implementation cerfc;
  void comprehensive_zero_search();
  int zeros_found;
public:
  TheoreticalSizeSpectrum(double umin, double umax);
  double operator() (double u,double t) const;
  void plot(somePlotter & p, double t);
  void plot_zeros(XPlotter & p,double value=0);
  void to_xmgr(ostream& os, double t) const;
  void to_xmgr(const char * filename, double t) const;
  cmplx get_zero(int i){
    return zero[i];
  }
  void save_zeros(const char * filename);

};

const int TheoreticalSizeSpectrum::n_zeros=3;

class continue_root_search{};

class real_smaller_in {
  const vector<cmplx> & z;
public:
  real_smaller_in(const vector<cmplx> & v):z(v){};
  bool operator()(int a,int b) const{
    return z[a].real()<=z[b].real();
  }
};

class abs_smaller_in {
  const vector<cmplx> & z;
public:
  abs_smaller_in(const vector<cmplx> & v):z(v){};
  bool operator()(int a,int b) const{
    return abs(z[a])<=abs(z[b]);
  }
};

void TheoreticalSizeSpectrum::comprehensive_zero_search(){
  const double unit=M_PI/log(beta);
  const double raster=0.5/8;
  
  zeros_found=0;
  const double epsilon=1e-12;
  
  for(double x=0; x<8.1*unit; x+=raster*unit){
    for(double y=-4*unit; y<4.1*unit; y+=raster*unit){
      try{
	cmplx z=analytic_root(kernel_function,cmplx(x,y),epsilon);
	// keep only roots in positive half plane
	if(z.real()<0) z=-conj(z);
	cmplx dkdz=analytic_derivative(kernel_function,z);
	// make sure imaginary roots are immaginary
	if(fabs(z.real()) < abs(10*epsilon/dkdz)){
	  z=cmplx(0,z.imag());
	}
	if(z.imag() < -(1-n) || z.imag() > lambda-2){
	  //throw continue_root_search();
	}
	// // drop roots that are TOO unstable
	// if(dkdz.imag()*z.imag()>0 && 
	//    2*abs(dkdz.imag()) < abs(dkdz.real())){
	//   throw continue_root_search();
	// }
	// make sure root was not found, yet
	for(int i=zeros_found;i-->0;){
	  if(abs(z-zero[i]) < abs(10*epsilon/dkdz)){
	    throw continue_root_search();
	  }
	}
	
	int i=zeros_found;
	zeros_found++;
	zero.resize(zeros_found);
	v.resize(zeros_found);
	D.resize(zeros_found);
	zero[i]=z;
	v[i]=cmplx(0,1)*dkdz*(1-n);
	D[i]=-analytic_second_derivative(kernel_function,z);
      }catch(analytic_root_finder_failed){
      }catch(continue_root_search){
      }
    }
  }
  vector<int> index(zeros_found);
  vector<cmplx> sorted_zeros(zeros_found);
  for(int i=zeros_found;i-->0;){
    index[i]=i;
  }
  sort(index.begin(),index.end(),abs_smaller_in(zero));
  REPORT(zeros_found);
  REPORT(M_PI/log(beta));
  ofstream plane("complex_plane.dat");
  ofstream arrows("arrows.txt");
  if(max_zeros && zeros_found>max_zeros)
    zeros_found=max_zeros;
  for(int ii=0;ii<zeros_found;ii++){
    int i=index[ii];
    REPORT(zero[i]);
    REPORT(v[i]);
    REPORT(D[i]);
    cmplx z=zero[i];
    sorted_zeros[ii]=z;
    cmplx v_est=(1-n)*
      abs(betaHat(conj(unit+cmplx(0,(3-lambda-n))))*
    	  betaHat(unit)*
	  (
	   abs(analytic_derivative(kappa0Hat_eating,unit+cmplx(0,3-lambda-n)))+
	   abs(analytic_derivative(kappa0Hat_being_eaten,unit+cmplx(0,3-lambda-n)))
	   ));
    //REPORT(v_est);
    // REPORT(abs(betaHat(conj(z+cmplx(0,(3-lambda-n))))*
    // 	       betaHat(z)));
    // REPORT(abs(betaHat(conj(unit+cmplx(0,(3-lambda-n))))*
    // 	       betaHat(unit)));
    // REPORT(abs(analytic_derivative(kappa0Hat,z+cmplx(0,3-lambda-n))));
    // REPORT(abs(analytic_derivative(kappa0Hat_eating,unit+cmplx(0,3-lambda-n))));
    // REPORT(abs(analytic_derivative(kappa0Hat_being_eaten,unit+cmplx(0,3-lambda-n))));
    // REPORT(abs(analytic_derivative(kappa0Hat_eating,unit+cmplx(0,3-lambda-n)))+(abs(analytic_derivative(kappa0Hat_being_eaten,unit+cmplx(0,3-lambda-n)))));
    if(v[i].real()*zero[i].imag()<0)
      cout << "^^^^^^^^^^^^^^^^^^^^^^^" << endl;
    plane << zero[i].real() << " " << zero[i].imag() << endl;
    cmplx p1=zero[i],p2=zero[i]+0.3*v[i]/cmplx(0,abs(v[i]));
    
    if(abs(p1.real())<=3 && abs(p1.imag())<=3){
      arrows << "@with line" << endl;
      arrows << "@    line on" << endl;
      arrows << "@    line loctype world" << endl;
      arrows << "@    line g0" << endl;
      arrows << "@    line " 
	   << p1.real() << ", " 
	   << p1.imag() << ", " 
	   << p2.real() << ", " 
	   << p2.imag() << endl;
      arrows << "@    line linewidth 1.0" << endl;
      arrows << "@    line linestyle 1" << endl;
      arrows << "@    line color 1" << endl;
      arrows << "@    line arrow 2" << endl;
      arrows << "@    line arrow type 1" << endl;
      arrows << "@    line arrow length 1.000000" << endl;
      arrows << "@    line arrow layout 1.000000, 0.500000" << endl;
      arrows << "@line def" << endl;
      if(p1.real()>0){
	arrows << "@with line" << endl;
	arrows << "@    line on" << endl;
	arrows << "@    line loctype world" << endl;
	arrows << "@    line g0" << endl;
	arrows << "@    line " 
	     << -p1.real() << ", " 
	     << p1.imag() << ", " 
	     << -p2.real() << ", " 
	     << p2.imag() << endl;
	arrows << "@    line linewidth 1.0" << endl;
	arrows << "@    line linestyle 1" << endl;
	arrows << "@    line color 1" << endl;
	arrows << "@    line arrow 2" << endl;
	arrows << "@    line arrow type 1" << endl;
	arrows << "@    line arrow length 1.000000" << endl;
	arrows << "@    line arrow layout 1.000000, 0.500000" << endl;
	arrows << "@line def" << endl;
      }
    }
    if(zero[i].real()>0){
      plane << -zero[i].real() << " " << zero[i].imag() << endl;
    }
  }
  plane << endl;

  //find smallerst distance between zeros:
  int imin=1,jmin=0;
  double dmin=abs(sorted_zeros[imin]-sorted_zeros[jmin]);
  for(int i=zeros_found;i-->0;){
    for(int j=i;j-->0;){
      if(abs(sorted_zeros[i]-sorted_zeros[j])<dmin){
	imin=i;
	jmin=j;
	dmin=abs(sorted_zeros[imin]-sorted_zeros[jmin]);
      }
    }
  }
  REPORT(imin);
  REPORT(jmin);
  REPORT(dmin);

  write_poles(plane);
}


void TheoreticalSizeSpectrum::plot_zeros(XPlotter & p,double value){
  const double range=1.5*2;
  const double radius=0.012*range;
  p.fspace(-range,-range,range,range);
  p.erase();
  p.color(0,0,0);

  char value_string[100];
  p.fontname("HersheySans");
  p.ffontsize(8*radius);
  sprintf(value_string,"%+.3g",value);
  p.fmove(-range+radius,range-radius);
  p.alabel ('l','t',value_string);

  p.pencolor(0xcfff,0xcfff,0xcfff);
  p.fline(0,-range,0,range);
  p.endpath();
  p.fline(-range,0,range,0);
  p.endpath();
  p.linemod("shortdashed"); 
  p.fline(-range,lambda-2,range,lambda-2);
  p.linemod("solid"); 
  p.endpath();
  p.pencolor(0,0,0);
  for(int i=0;i<zeros_found;i++){
    cmplx direction=(cmplx(0,-1)*v[i])/abs(v[i]);
    cmplx head=zero[i]+3*radius*direction;
    p.fcircle(zero[i].real(),zero[i].imag(),radius);
    p.fline(zero[i].real(),zero[i].imag(),head.real(),head.imag());
    p.endpath();
    p.fcircle(-zero[i].real(),zero[i].imag(),radius);
    p.fline(-zero[i].real(),zero[i].imag(),-head.real(),head.imag());
    p.endpath();
  }
  vector<double> pp=poles();
  p.fillcolor (0,0,0);
  p.filltype(1);
  for(int i=0;i<pp.size();i++){
    p.fcircle(0,pp[i],radius);
  }   
  p.filltype(0);

}

void TheoreticalSizeSpectrum::save_zeros(const char * filename){
  ofstream file(filename);

  for(int i=0;i<zeros_found;i++){
    file << zero[i].real() << " " << zero[i].imag() << endl;
  }
  file << endl;
  vector<double> pp=poles();
  for(int i=0;i<pp.size();i++){
    file << 0 << " " << pp[i] << endl;
  }   
}


TheoreticalSizeSpectrum::TheoreticalSizeSpectrum(double min, double max):
  zero(n_zeros),v(n_zeros),umin(min),umax(max),cerf(40),cerfc(40) {
#if 0
  for(int i=0;i<n_zeros;++i){
    REPORT(i);
    zero[i]=analytic_root(kernel_function,M_PI/log(beta)*i);
    REPORT(zero[i]);
    v[i]=analytic_derivative(kernel_function,zero[i])/double(4)*cmplx(0,1);
    REPORT(v[i]);
  }
  if(n_zeros >0){
    zero[0]=cmplx(0,zero[0].imag());
    v[0]=cmplx(v[0].real(),0);
  }
  zeros_found=n_zeros;
#else
  comprehensive_zero_search();
#endif
}


double TheoreticalSizeSpectrum::operator()(double u, double t) const
{
  if(t<=0) return 0;

  double b=0;
  double z=exp(u/4);
  const double cutoff=8*sigma;

  for(int i=0;i<zeros_found;i++){
    // work in progress...
    int s=(v[i].real()>0 ? + 1 : -1);
    if(zero[i].real()!=0) s*=2;

    // determine upper and lower bounds of integral:
    double z1=z-v[i].real()*t;
    if(zero[i].imag()>-upper_tail_exponent && v[i].real() < 0){
      z1=0;
      if(v[i].real() < 0)
	s*=-1;
    }
    if(zero[i].imag()<-lower_tail_exponent && v[i].real() > 0){
      z1=HUGE_VAL;
      if(v[i].real() > 0)
	s*=-1;
    }
    double u1,u2;
    if(s>0){
      if(z1>0){
	u1=max<double>(-cutoff,4*log(z1));
      }else{
	u1=-cutoff;
      }
      u2=min<double>(cutoff,u);
    }else{
      u2=min<double>(cutoff,4*log(z1));
      u1=max<double>(-cutoff,u);
    }
    if(u2>u1){
      // analytically compute integrals
      if(false){// THIS HAS THE WRONG ALLOMETRIC WEIGHT FACTOR FOR FISHING
	cmplx A=cmplx(0,1)*(sigma*zero[i])*sigma;
	double B=M_SQRT2*sigma;
	cmplx erf_diff;
	if(u1+A.real()>0){
	  erf_diff=-(cerfc((u1+A)/B)-cerfc((u2+A)/B));
	}else{
	  erf_diff=(cerfc(-(u1+A)/B)-cerfc(-(u2+A)/B));
	}
	b=b+(double(s)*exp(cmplx(0,1)*u*zero[i]-(sigma*zero[i])*(sigma*zero[i])/double(2))/
	     (double(8)*v[i])*erf_diff).real();
      }else{// This one is correct
	cmplx n_xi=cmplx(0,1)*zero[i]-(+3-lambda-n);
	cmplx A=(sigma*n_xi)*sigma;
	double B=M_SQRT2*sigma;
	cmplx erf_diff;
	if((u1+A).real()>0){
	  erf_diff=-(cerfc((u1+A)/B)-cerfc((u2+A)/B));
	}else{
	  erf_diff=(cerfc(-(u1+A)/B)-cerfc(-(u2+A)/B));
	}
	double contribution=(1-n)*F*s*sigma*sqrt(M_PI/2)*
	  (exp(cmplx(0,1)*u*zero[i]+(sigma*n_xi)*(sigma*n_xi)/
	       double(2))*
	   erf_diff/v[i]).real();
	if(contribution != 2*contribution){
	  b=b+contribution;
	}	  
      }	
    }
  }
  return b;
}

void TheoreticalSizeSpectrum::to_xmgr(const char * filename, double t) const{
  ofstream os(filename);
  to_xmgr(os,t);
}

void TheoreticalSizeSpectrum::to_xmgr(ostream & os, double t) const{
  const double du=0.05;
  
  for(double u=umin; u<=umax; u+=du){
    os << u/M_LN10 << " " 
       << this->operator()(u,t)/(exp((2-lambda)*u)*tildeBtot) << endl;
  }
}

void TheoreticalSizeSpectrum::
plot(somePlotter & p, double t){
  const double du=0.05;
  const double scaling = (p.height()/bRange*M_LN10/2);
  
  p.prepare(umin,umax,t);
  for(double u=umin; u<=umax; u+=du){
    p.point(u,this->operator()(u,t) * scaling);
  }
}



bool has_suffix(const char * suff,const char * filename){
  return 0==strcmp(suff,filename+strlen(filename)-strlen(suff));
}

void read_configuration(int argc,char** argv){
  int c;
  const char * formatstring = "vg:W:M:Z:";

  while(argc>1){
    if(has_suffix(".cfg",argv[1])){
      string in_file_name=argv[1];
      read_parameters_from_file(in_file_name);
      argv++;
      argc--;
    }else if(has_suffix(".sss",argv[1])){
      initial_state_file=argv[1];
      argv++;
      argc--;
    }else{
      c=getopt(argc, argv, formatstring);
      if (c == -1)
	break;
      
      switch(c){
      case 'Z':
	show_zero_movie=optarg;
	break;
      case 'W':
	output_file=optarg;
	break;
      case 'M':
	MERP_output_file=optarg;
	break;
      case 'v':
	printf("Main code version: %s\n",version.c_str());
	printf("Copyright (C) 2010-2014 Axel G. Rossberg\n");
	printf("License GPL 3\n");
	exit(0);
	break;
      case 'g':
	gif_image_file=optarg;
	if(!has_suffix(".gif",gif_image_file)){
	  FATAL_ERROR(" \"" << optarg << "\" is not a .gif file");
	}
	break;
      default:
	// print_usage(argv,formatstring);
	// exit(c=='h'?0:1);
	FATAL_ERROR("unknown command-line parameters");
      }
      argc-=optind-1;
      argv+=optind-1;
      optind=1;
    }
  }
}

void kernelHat_toxmgr()
{
  ofstream kernelHat("kernelHat.dat");
  // for(double z=-10;z<10;z+=0.01){
  //   kernelHat 
  //     << z << " " 
  //     << kernel_function(z).real() << " "
  //     << kernel_function(z).imag() << endl;
  // }
  // kernelHat << endl;

  for(double z=-10;z<10;z+=0.01){
    kernelHat 
      << z << " " 
      << abs(kernel_function(z)) << endl;
  }
  kernelHat << endl;

  for(double z=-10;z<10;z+=0.01){
    kernelHat 
      << z << " " 
      << abs(kappa0Hat(z+cmplx(0,3-lambda-n))) << endl;
  }
  kernelHat << endl;

  for(double z=-10;z<10;z+=0.01){
    kernelHat 
      << z << " " 
      << abs(conj(betaHat(conj(z+cmplx(0,(3-lambda-n)))))) << endl;
  }
  kernelHat << endl;

  for(double z=-10;z<10;z+=0.01){
    kernelHat 
      << z << " " 
      << abs(betaHat(z)) << endl;
  }
  kernelHat << endl;

  for(double z=-10;z<10;z+=0.01){
    kernelHat 
      << z << " " 
      << abs(exclusionHat(z)) << endl;
  }
  kernelHat << endl;
}

void kernelHat2()
{
  ofstream kernelHat("kernelHat2.dat");

  for(double y=-3;y<3;y+=0.01){
    kernelHat 
      << y << " " 
      << kernel_function(cmplx(M_PI/log(beta),y)).real() << " "
      << kernel_function(cmplx(M_PI/log(beta),y)).imag() << " "
      << endl;
  }
  kernelHat << endl;
}

void kernelHat3()
{
  ofstream kernelHat("kernelHat3.dat");

  for(double y=-1.5;y<1.5;y+=0.01){
    kernelHat 
      << y << " " 
      << kernel_function(cmplx(0,y)).real() << " "
      << endl;
  }
  kernelHat << endl;
}

void zero_movie(){
  Plotter::parampl ("PAGESIZE", const_cast<char*>("letter"));
  Plotter::parampl ("USE_DOUBLE_BUFFERING",const_cast<char*>("yes"));
  Plotter::parampl ("BITMAPSIZE",const_cast<char*>("600x600"));
  
  XPlotter plotter(cin,cout,cerr);
  if (plotter.openpl () < 0)          // open Plotter
    FATAL_ERROR("could not open plotter");

  double x=get_cfg_parameter(show_zero_movie);
  cout << show_zero_movie << " = " << x << endl;
  char new_x_str[33];
  for(double scale=0.1;scale<10;scale*=1.05){
    double new_x=scale*x;
    snprintf(new_x_str, 33, "%25g", new_x);
    //    cout << show_zero_movie << " = " << new_x_str << endl; 
    set_cfg_parameter(show_zero_movie,new_x_str);
    set_dependent_constants();
    TheoreticalSizeSpectrum b_theo(-10,10);
    timespec t = {0,0.02 * 1e9};
    nanosleep(&t,&t);
    b_theo.plot_zeros(plotter,atof(new_x_str));
  }
  exit(0);
}


int main(int argc,char** argv){
  ///// Initialization
  signal_handling();
  read_configuration(argc,argv);

  if(show_zero_movie){
    zero_movie();
  }

  //  fftw_init_threads();
  //  fftw_plan_with_nthreads(number_of_CPUs());
  set_dependent_constants();
  set_cfg_parameter("MAX_STEP_SIZE","0.01");
  set_cfg_parameter("MAX_INTEGRATOR_STEPS","3200000");
  set_cfg_parameter("DEFAULT_ABSOLUTE_TOLERANCE","0.0001");

  // set Plotter parameter
  Plotter::parampl ("PAGESIZE", const_cast<char*>("letter"));
  Plotter::parampl ("USE_DOUBLE_BUFFERING",const_cast<char*>("yes"));
  ostringstream BITMAPSIZE;
  BITMAPSIZE << BITMAPWIDTH << "x" << BITMAPHEIGHT;
  Plotter::parampl ("BITMAPSIZE", const_cast<char*>(BITMAPSIZE.str().c_str()));
  Plotter::parampl ("GIF_DELAY", const_cast<char*>("5"));

  timespec tsleep = {0,0.1 * 1e9},trem;
  myPlotter<> plotter;
  nanosleep(&tsleep,&trem);
  myPlotter<> tplotter(true);
  nanosleep(&tsleep,&trem);
  myPlotter<> iplotter(!communitySpectrum && !noSpecies);
  myPlotter<GIFPlotter> gif_iplotter(gif_image_file);
  Plotter::parampl ("PAGESIZE", const_cast<char*>("letter"));
  Plotter::parampl ("USE_DOUBLE_BUFFERING",const_cast<char*>("yes"));
  Plotter::parampl ("BITMAPSIZE",const_cast<char*>("350x350"));
  nanosleep(&tsleep,&trem);
  myPlotter<> cplotter(true);

  REPORT(NTilde);
  REPORT(NTilde_crit);
  REPORT(kappa0Hat(cmplx(0,3-lambda-n)));

  // For whatever reason, eigenfunction_extractor leaves memory in a
  // corrupt state after exiting, so we have to separate the two
  // operations here:
  if(compute_linear_mode==1){//compute linear mode for growth rate -1
    eigenfunction_extractor b0(256*2*4,0.5/4,-1);
    {
      fixed_point_analyzer b0_finder(&b0,-200);
      b0_finder.snap_to_fixed_point();
    } 
    FATAL_ERROR("CODE BROKEN IN LINES FOLLOWING THIS ONE");
    // b0.l2_normalize(b0).to_xmgr("b0.dat");
    // b0.individual_spectrum().l2_normalize().to_xmgr("n0.dat");
  //   b0.compute_residuum().to_xmgr("b0_residuum.dat");
  //   exit(1);
  // }else if(compute_linear_mode==2){
  //   adjoint_eigenfunction_extractor c0(256*2*4,0.5/4,-1);
  //   {
  //     fixed_point_analyzer c0_finder(&c0,-200);
  //     c0_finder.snap_to_fixed_point();
  //   }
  //   c0.adjoint_correct();
  //   c0.l2_normalize().to_xmgr("c0.dat");
  //   c0.compute_residuum().to_xmgr("c0_residuum.dat");
  //   exit(1);
  }
  
  kernelHat_toxmgr();
  kernelHat2();
  kernelHat3();

  SizeSpectrum b(256*2,log(pow(10,0.1))/2);

  if(initial_state_file){
    b.load_state(initial_state_file);
  }


  b.save_PPMR_window("PPMR_window.dat");

  b.plot(plotter);
  b.plot(plotter);
  b.plot_individuals(iplotter);
  {
    double lambda_hold=lambda;
    int balance_case_hold=balance_case;
    REPORT(eutro_lambda());
    if(gamma_boost>1){
      lambda=eutro_lambda();
      balance_case=1;
    }
    set_dependent_constants();
    TheoreticalSizeSpectrum b_theo(b.min_u(),b.max_u());
    b.save_PPMR_mechanism(b_theo.get_zero(0),"zero0.dat");
    b.save_PPMR_mechanism(b_theo.get_zero(1),"zero1.dat");
    b.save_PPMR_mechanism(b_theo.get_zero(2),"zero2.dat");
    b.save_beta_tilde("beta_tilde.dat");
    b_theo.plot_zeros(cplotter,0);
    b_theo.plot_zeros(cplotter,0);
    b_theo.save_zeros("poles-zero.dat");
    lambda=lambda_hold;
    balance_case=balance_case_hold;
    set_dependent_constants();
  }
  TheoreticalSizeSpectrum b_theo(b.min_u(),b.max_u());
  b_theo.plot(tplotter,0);
  b_theo.plot(tplotter,0);
  //return 0;

  if(MERP_output_file){
    b.MERP_report_size_classes();
  }

  sleep(2);

  ofstream total_biomass_file("total_biomass.dat");
  ///// Main Work:
  double t_final=stoppingTime*time_unit;
  int plot_skip=100;
  int plot_step=0;
  double total_biomass=0;
  int stage=1;
  double t=0;
  int step_counter=0; // ode steps
  double next_gif_image_t = 0.001;
  if(!doSimulate){
    t=0.1;
  }else{
    b.current_time=t;
  }
  while(t<t_final){
    {
      ODE_state state(&b);
      while(t<t_final){	
	// if(stage<5 && t>1000000){
	//   stage++;
	//   set_cfg_parameter("MAX_STEP_SIZE","10.0");
	//   break;
	// }
	if(!doSimulate){}
	else if(stage==7 && t>10000000){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","30.0");
	  break;
	}
	else if(stage==6 && t>1000000){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","10.0");
	  break;
	}
	else if(stage==5 && t>100000){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","3.0");
	  break;
	}
	else if(stage==4 && t>10000){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","1.0");
	  break;
	}
	else if((maturation_schedule < 2) && stage==3 && t>1000){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","0.3");
	  plot_skip=1000;
	  break;
	}  
	else if(stage==2 && t>100){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","0.1");
	  plot_skip=1000;
	  break;
	}  
	else if(stage==1 && t>10){
	  stage++;
	  set_cfg_parameter("MAX_STEP_SIZE","0.03");
	  plot_skip=100;
	  //	  break;
	}
	else if(pulse_perturbation && t>pulse_perturbation*time_unit){
	  b.no_more_fishing();
	  //	  break;
	}

	if(t==0 || plot_step++>plot_skip){
	  REPORT(plot_step);
	  REPORT(t);
	  plot_step=0;
	  b.read_state_from(state);
	  b.plot(plotter);
	  b_theo.plot(tplotter,t);
	  total_biomass=b.plot_individuals(iplotter);
	  total_biomass_file << t << " " << total_biomass << endl;
	  total_biomass_file.flush();
	  // timespec sleeptime = {0,0.2e9}, reminder;
	  // nanosleep(&sleeptime,&reminder);
	}
	//cout << "\r" << t << " " << total_biomass << "               ";
	static double previous_digits=-10;
	double digits=floor(2*log10(t/time_unit))/2;
	if(digits!=previous_digits){
	  previous_digits=digits;
	  ostringstream filename;
	  filename << "change_" << digits << ".dat";
	  ofstream os(filename.str().c_str());
	  b.read_state_from(state);
	  b.to_xmgr(os);
	  os << endl;
	  //b_theo.to_xmgr(os,t);
	}
	if(t==0 || (gif_image_file && t > next_gif_image_t)){
	  b.read_state_from(state);
	  b.plot_individuals(gif_iplotter); //need to define a gif plotter!
	  next_gif_image_t*=1.1;
	  if(MERP_output_file){
	    static ofstream MERP_os(MERP_output_file);
	    b.MERP_write_state(MERP_os);
	    {
	      static int frame_number=0;
	      ostringstream filename;
	      filename  << "tmp/Bframe" << setfill('0') << setw(5) 
			<< frame_number << ".dat";
	      REPORT(frame_number);
	      ofstream frame(filename.str().c_str());
	      vector<double> c=b.individual_spectrum();
	      for(int i=0;i<b.number_of_variables();++i){
		frame << exp(b.u(i))*1e-3 << " " << exp(c[i]) << endl;
	      }
	      frame_number++;
	    }
	  }
	}
	if(doSimulate){
	  if(state.integrate_one_step(t+1000)){
	    state.diagnosis();
	    break;
	  }
	  step_counter++;
	  t=b.current_time;
	  //effective relaxation rate:
	  // double r = b.effective_density_dependence(state);
	  // cout << "RRR " << t/time_unit 
	  //      << " " << r*time_unit << endl;
	}else{
	  t*=1.01;
	}
      }
    } 
  }
  REPORT(step_counter);

  b.individual_trophic_levels();

  if(output_file){
    b.save_state(output_file);
  }

  b.to_xmgr("b_final.dat");
  ofstream B_final("B_final.dat");
  vector<double> c=b.individual_spectrum();
  for(int i=0;i<b.number_of_variables();++i){
    if(exp(c[i]) > 1e-20)
      B_final << b.u(i)/M_LN10 << " " << c[i]/M_LN10 << endl;
  }
  ofstream A_final("A_final.dat");
  vector<double> a=b.species_spectrum();
  for(int i=0;i<b.index_of(b.max_u());++i){
    A_final << b.u(i)/M_LN10 << " " << exp(a[i]) << endl;
  }
  // //  b.individual_spectrum().to_xmgr("b_individuals.dat");
  // b_theo.to_xmgr("b_theory.dat",t);

  ODE_vector state(b.number_of_variables());
  ODE_vector time_derivative(b.number_of_variables());
  b.write_state_to(state);
  b.dynamics(state,time_derivative);
  ofstream fitness("fitness");
  for(int i=0;i<b.index_of(b.max_u());++i){
    fitness << b.u(i)/M_LN10 << " " << time_derivative[i]/(k*exp(b.u(i)*(n-1))) << endl;
  }

  exit(0);
}
