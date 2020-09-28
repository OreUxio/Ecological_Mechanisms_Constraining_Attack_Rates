// $Id: Statistics.cc 2477 2016-10-23 11:38:50Z axel $

#include "Statistics.h"
#include "error.h"
#include "NewMatrix.h"
#include <math.h>
#include <fstream>

#ifndef NAN
static const double my_nan=sqrt(-1.0);
#define NAN my_nan
#endif
#ifndef INFINITY
static const double my_infinity=1.0/0.0;
#define INFINITY my_infinity
#endif

#ifdef ON_SX5FSV
// there seem to be two versions of iostreams, one in the std
// namespace, the other not, so we have to fight ourselves through
// this:
#define IOSTD 
#else
#define IOSTD std::
#endif

void set_format(std::ostream & os){
  os.width(15);
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(10);
}  

Histogram_Estimator::Histogram_Estimator(double min, double max, int n)
{
  if(n>100)
    the_histogram=std::simple_vector<int>(int(sqrt(double(n))));
  else
    the_histogram=std::simple_vector<int>(10);
  the_lowest_x_included=min;
  the_bin_width=(max-min)/the_number_of_bins();
}

Histogram_Estimator::Histogram_Estimator(double min, double width):
  the_histogram(),
  the_lowest_x_included(min),
  the_bin_width(width)
{
}

Histogram_Estimator::~Histogram_Estimator(){
}

void Histogram_Estimator::sample(double x)
{
  double normalized=(x-the_lowest_x_included)/the_bin_width;
  N++;
  if(normalized < 0) return;
  if(int(normalized) >= the_number_of_bins()){
    the_histogram.resize(int(normalized)+1);
  }
  the_histogram[int(normalized)]++;
  return;
}

double Histogram_Estimator::estimate_density(double x)
{
  double normalized=(x-the_lowest_x_included)/the_number_of_bins();
  if(N>0){
    if(normalized < 0 || int(normalized) >= the_number_of_bins())
      return 0;
    else
      return the_histogram[int(normalized)]/N/the_bin_width;
  }
  else
    FATAL_ERROR("attempt to calculate histrogam without taking samples");
  return 0;
}
  
  
void Histogram_Estimator::save(const char * name){
  std::ofstream os(name);
  set_format(os);
  for(int i=0;i<the_number_of_bins();i++){
    os << (i+0.5)*the_bin_width << " " 
       << the_histogram[i] << " "
       << sqrt(the_histogram[i]) << std::endl;
  }
}
  

void weighted_average_meter::sample(double x,double weight){
    sum+=x*weight;
    square_sum+=x*x*weight;
    samples+=weight;
}

double weighted_average_meter::readout() const {
  return sum/samples;
}

double weighted_average_meter::sample_var() const {
  return std::max<double>(0,(square_sum-sum*sum/samples))/(samples);
}

double weighted_average_meter::sample_std() const {
    return sqrt(sample_var());
}

double average_meter::std() const{
    return sqrt(sample_var()*samples/(samples-1));
}

double average_meter::var() const{
    return sample_var()*samples/(samples-1);
}

double average_meter::error() const{
    return std()/sqrt(double(samples));
}

double average_meter::error_var() const{
    return var()/double(samples);
}

void average_meter::sample(double x){
  weighted_average_meter::sample(x,1);
}
void average_meter::sample(double x,double weight){
  FATAL_ERROR("This is not a weighted average meter!");
}
int average_meter::n(){
  return int(0.5+samples);
}

std::ostream & operator<<(std::ostream &stream, 
			  const average_meter & av){
  stream << av.readout() << " +/- " << av.error();
  return stream;
}

// multinormal_distribution

multinormal_distribution::multinormal_distribution(NewVector &m, NewMatrix &c):mean(m),cov(c){
  ALWAYS_ASSERT(m.size()==c.SIZE1());
  ALWAYS_ASSERT(c.SIZE1()==c.SIZE2());
};

std::simple_vector<double> multinormal_distribution::main_axis() const {
  int n=cov.SIZE2();
  NewVector evals;
  NewMatrix m=cov;
  NewEigen(m,evals); // get eigensystem

  std::simple_vector<double> axis(n);
  for(int i=0;i<n;i++){
    axis[i]=m(i,n-1);
  }
  return axis;
}

std::simple_vector<double> multinormal_distribution::raw_main_axis() const {
  int n=cov.SIZE2();
  NewVector evals;
  NewMatrix m=cov+mean.t()*mean;
  NewEigen(m,evals); // get eigensystem
  REPORT(evals);

  std::simple_vector<double> axis(n);
  for(int i=0;i<n;i++){
    axis[i]=m(i,n-1);
  }
  return axis;
}

void multinormal_distribution::save(const char * name) const{
  std::ofstream os(name);
  set_format(os);

  os << mean.size() << std::endl;
  os << mean;
  os << cov;

  // as a service, include also correlation:
  NewMatrix corr_corrector(cov.SIZE1(),cov.SIZE2());
  corr_corrector = arma::diagmat(corr_corrector);
  
  for(int i=cov.SIZE1();i-->0;){
    if(cov(i,i))
      corr_corrector(i,i)=1/sqrt(cov(i,i));
  }
  os << prod(NewMatrix(prod(corr_corrector,cov)),corr_corrector);
};

void multinormal_distribution::load(const char * name){
  std::ifstream os(name);
  int size;
  os >> size;
  mean=NewVector(size);
  cov=NewMatrix(size,size);
  for(int i=0;i<size;i++){
    os >> mean[i];
  }
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      os >> cov(i,j);
    }
  }
  
  // drop correlation matrix
  double dummy;
  for(int i=0;i<size;i++){
    for(int j=0;j<size;j++){
      os >> dummy;
    }
  }
};

double multinormal_distribution::var_of_sum(){
  double sum=0;
  for(int i=cov.SIZE1();i-->0;){
    for(int j=cov.SIZE2();j-->0;){
      sum+=cov(i,j);
    }
  }
  return sum;
}

void save(std::simple_vector<double> const & data,
	  std::simple_vector<bool> const & select,
	  std::simple_vector<bool> const & fix,
	  const char * name){
  int size=0;
  for(unsigned int i=0;i<fix.size();i++) 
    if(select[i] && !fix[i]) size++;

  std::ofstream os(name);
  set_format(os);

  os << size << std::endl;;

  for(unsigned int i=0;i<fix.size();i++)
    if(select[i] && !fix[i]) os << data[i] << std::endl;;
}

void load(NewVector & data,
	  const char * name){
  int size;
  std::ifstream os(name);
  
  os >> size;
  
  data=NewVector(size);

  for(int i=0;i<size;i++)
    os >> data[i];
}

//chi_square_meter:

void chi_square_meter::_initialize_selection_given(){
  int count = 0;
  for(unsigned int i=0;i<_selection.size();i++)
    if(_selection[i]) count++;
  _size=count;
  _sum=NewVector(count,arma::fill::zeros);
  _square_sum=NewMatrix(count,count,arma::fill::zeros);
}


#if 0 // dead code:

static double norm(const HepMatrix & m){
  double sum=0;
  for(int i=0;i<m.num_row();i++)
    for(int j=0;j<m.num_col();j++)
      sum+=fabs(m[i][j]);
  return sum/(m.num_col()*m.num_row());
}
    

static HepMatrix sqrt2(const HepMatrix & m0){
  HepMatrix 
    previous(m0.num_row(),m0.num_col()),
    previous2(m0.num_row(),m0.num_col()),
    m=m0;
  m/=sqrt(norm(m0));
  double acc=1e-8*sqrt(norm(m0));
  int err;
  int count=0;
  do{
    previous2=previous;
    previous=m;
    if(norm(previous-previous2)<acc) break;
    m=0.5*(m+m0*m.inverse(err));
    count++;
  }while(norm(previous-m)>acc && ! err && count < 300);
  if(err) WARNING("matrix for sqrt is irregular");
  return m;
}

#endif

static NewMatrix sqrt(const NewMatrix & m0){
  int n=m0.SIZE2();
  NewVector evals;
  NewMatrix m=m0,evec;
  NewEigen(m,evals); // get eigensystem

  NewMatrix Lambda(n,n);
  Lambda = arma::diagmat(Lambda);
  for(int i=n;i-->0;){
    if(evals(i)<0){
      WARNING("Negative eigenvalue, cannot take square root");
    }else{
      Lambda(i,i) = sqrt(evals(i));
    }
  }

  FATAL_ERROR("This part of code needs to be checked, do we do the transposition of matrix m correctly?? Use e.g. chi_square_meter::deviation to test.");
  NewMatrix sqrt_mat = trans(prod(m,Lambda));
  sqrt_mat = prod(m,sqrt_mat);

  return sqrt_mat;
}

void chi_square_meter::sample(std::simple_vector<double> data){
  NewVector selected_data(_size,arma::fill::zeros);
  int j=0;
  for(unsigned int i=0;i<_selection.size() && i<data.size();i++){
    if(_selection[i]) selected_data[j++]=data[i];
  }
  _sum+=selected_data;
  _square_sum+=outer_prod(selected_data,selected_data);
  _n_samples++;
}
  
std::simple_vector<double> chi_square_meter::mean() const {
  std::simple_vector<double> result(_selection.size());
  int j=0;
  for(unsigned int i=0;i<_selection.size();i++){
    if(_selection[i]) 
      result[i]=_sum[j++]/_n_samples;
    else
      result[i]=NAN;
  }
  return result;
}

inline double dot(const NewVector & v1,const 
		  NewVector & v2){
  return inner_prod(v1,v2);
}

double chi_square_meter::chi_square_old(std::simple_vector<double> data) const {
  NewVector selected_data(_size);
  int j=0;
  for(unsigned int i=0;i<_selection.size() && i<data.size();i++){
    if(_selection[i]) selected_data[j++]=data[i];
  }
  //// Compute mean:
  NewVector mean=_sum/_n_samples;

  //// Compute Covariance Matrix:
  NewMatrix cov=
			 _square_sum/_n_samples-outer_prod(mean,mean);
  
  // correct for finite sample size (correct?):
  cov*=double(_n_samples)/(_n_samples-1);
  
  selected_data-=mean;
  int ierr;
  double result=dot(selected_data,solve(cov,selected_data));
  if(ierr){
    IOSTD cout << cov;
    FATAL_ERROR("covarance matrix is irregular");
  }
  return result;
}

void chi_square_meter::postfix(std::simple_vector<double> const & data,
			       std::simple_vector<bool> const & fixed0,
			       NewVector & v1, 
			       NewVector & mean_star, 
			       NewMatrix & cov11, 
			       NewMatrix & icov11) const {
  std::simple_vector<bool> fixed=fixed0;
  //fixed.resize(_selection.size(),false);
  fixed.resize(_selection.size()); //false is default element
  NewVector selected_data(_size);
  int j=0;
  for(unsigned int i=0;i<_selection.size() && i<data.size();i++){
    if(_selection[i]) selected_data[j++]=data[i];
  }
  //// Compute mean:
  NewVector mean=_sum/_n_samples;

  //// Compute Covariance Matrix:
  NewMatrix cov=
			 _square_sum/_n_samples-outer_prod(mean,mean);
  
  // correct for finite sample size (correct?):
  cov*=double(_n_samples)/(_n_samples-1);

  // regularize cov:
  double minc=INFINITY;
  for(int i=0;i<_size;i++){
    if(cov(i,i)!=0 && cov(i,i)<minc) minc=cov(i,i);
  }
  for(int i=0;i<_size;i++){
    if(cov(i,i)==0) cov(i,i)=minc/10;
  }
  
  NewMatrix icov=cov.i();

  //// now separate fixed and variable part:
  // count the number of fixed elements:
  int fcount=0;
  for(unsigned int i=0;i<fixed.size();i++){
    if(_selection[i] && fixed[i]) fcount++;
  }
  int vcount=_size-fcount;

  // compute sub-matrixes and vectors
  NewVector v0(fcount),m0(fcount),m1(vcount);
  v1=NewVector(vcount);
  NewMatrix icov10(vcount,fcount);
  icov11=NewMatrix(vcount,vcount);
  int k=0,k1=0; 
  for(unsigned int i=0;i<_selection.size();i++){
    if(_selection[i]){
      if(!fixed[i]){
	int l=0,l10=0,l11=0;
	for(unsigned int j=0;j<_selection.size();j++){
	  if(_selection[j]){
	    if(fixed[j]){
	      icov10(k1,l10)=icov(k,l);
	      v0[l10]=selected_data[l];
	      m0[l10]=mean[l];
	      l10++;
	    }else{
	      icov11(k1,l11)=icov(k,l);
	      v1[l11]=selected_data[l];
	      m1[l11]=mean[l];
	      l11++;
	    }
	    l++;
	  }
	}
	k1++;
      }
      k++;
    }
  }

  cov11=arma::inv(icov11);
      
  if(fcount>0){// armadillo does not like matrices with zero rows or cols
    mean_star=m1+prod(cov11,NewVector(prod(icov10,v0-m0)));
  }else{
    mean_star=m1;
  }
}

double chi_square_meter::chi_square(std::simple_vector<double> const & data,
				    std::simple_vector<bool> const & fixed) const {
  double dummy_log_det_cov;
  return chi_square(&dummy_log_det_cov,data,fixed);
}

double chi_square_meter::chi_square(double * log_det_cov,
				    std::simple_vector<double> const & data,
				    std::simple_vector<bool> const & fixed) const {
  NewVector v1,mean_star;
  NewMatrix cov11, icov11;
  //main calculations are done here:
  postfix(data,fixed,v1,mean_star,cov11,icov11);

  NewVector delta=v1-mean_star;
  
//   HepMatrix dd=cov11;
//   for(int i=0;i<dd.num_row();i++){
//     for(int j=0;j<dd.num_col();j++){
//       if(i!=j) dd[i][j]=0;
//     }
//   }
//   int ierr;
//   dd=sqrt(dd.inverse(ierr));
//   std::cout << "correlation matrix:" << std::endl;
//   IOSTD cout << dd*cov11*dd;

//   std::cout << "covariance matrix:" << std::endl;
//   IOSTD cout << cov11;

  NewVector ev;
  NewEigen(cov11,ev); // get eigensystem
  std::cout << std::endl << "eigenvalues: ";
  double log_det=0;
  double tr_cov2=0; //trace
  for(int i=0;i<cov11.SIZE2();i++){
    std::cout << ev[i] << " " ;
    if(!my_isnan(log_det)){
      if(ev[i]<0){
	log_det+=NAN;
      }else{
	log_det+=log(ev[i]);
      }
      tr_cov2+=ev[i]*ev[i];
    }
  }
  std::cout << std::endl << "log10(det): " << log_det/log(10.0) 
	    << ", TR(cov^2): " << tr_cov2 << std::endl;

  *log_det_cov=log_det; //secondary return value
  
  return dot(delta,prod(icov11,delta));
}

  
std::simple_vector<double> 
chi_square_meter::deviation(std::simple_vector<double> const & data,
			    std::simple_vector<bool> const & fixed) 
  const
{
  NewVector v1,mean_star;
  NewMatrix cov11, icov11;
  //main calculations are done here:
  postfix(data,fixed,v1,mean_star,cov11,icov11);

  NewVector delta=v1-mean_star;
  
  NewMatrix istderr_mat=sqrt(icov11);

  IOSTD cout << istderr_mat ;

  IOSTD cout << icov11-prod(istderr_mat,istderr_mat);

  NewVector dev=prod(istderr_mat,delta);
  
  std::simple_vector<double> result(_selection.size());
  int j=0;
  for(unsigned int i=0;i<_selection.size();i++){
    if(_selection[i] && !fixed[i])
      result[i]=dev[j++];
    else
      result[i]=NAN;
  }
  return result;
}


std::simple_vector<weighted_average_meter> chi_square_meter::
mean_and_var(std::simple_vector<double> const & data, 
	     std::simple_vector<bool> const & fixed) const {
  int s=_selection.size();
  std::simple_vector<weighted_average_meter> avs(s);
  NewVector v1,mean_star;
  NewMatrix cov11, icov11;
  //main calculations are done here:
  postfix(data,fixed,v1,mean_star,cov11,icov11);

  // fill 
  int j=0;
  for(unsigned int i=0;i<_selection.size();i++){
    if(_selection[i] && !fixed[i]){
      //std::cout << j << " " << mean_star[j] << std::endl ;
      avs[i]+=weighted_average_meter(mean_star[j],cov11(j,j));
      j++;
    }
  }
  
  return avs;
}

multinormal_distribution 
chi_square_meter::estimate(std::simple_vector<double> const & data, 
			   std::simple_vector<bool> const & fixed){
  NewVector v1,mean_star;
  NewMatrix cov11, icov11;
  //main calculations are done here:
  postfix(data,fixed,v1,mean_star,cov11,icov11);
  return multinormal_distribution(mean_star,cov11);
}

