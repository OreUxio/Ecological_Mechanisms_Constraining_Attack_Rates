// $Id:$ 

// A somewhat unrelated code I used to do the statistics for the
// "power" paper.

#include <fstream>
#include <algorithm>
#include <math.h>

#include "polyfit.h"
#include "random.h"
#include "error.h"
#include <sys/time.h>

#ifdef GET_INDIRECT_DEPENDENCIES
#include "Statistics.h"
//#include "linpack_eigen.h"
#include "NewMatrix.h"
#include "cfgList.h"
#endif

//#define ENVELOPES

using namespace std;

int main(int argc,char *argv[]){
  struct timeval t;
  if(gettimeofday(&t,0)) FATAL_ERROR("problem with system clock");
  srand(t.tv_usec); // a new seed every microsecond
  set_random_seed(rand());
  ifstream table( argc<2 ? "table.dat" : argv[1] );
  sequence<double> diversity, Sc, Sr, nu_mean, nu_sem, effort;
  sequence<double> lin_effort, inv_effort, no_effort, log_effort;
  sequence<average_meter> nu;
  string dummy;
  const double sem_fac=1;
  int n=0;
  int max_diversity=0;
  while(!table.eof()){
    table >> diversity[n];
    if(!diversity[n]){
      diversity.resize(n);
      break;
    }
    max_diversity=max<int>(diversity[n],max_diversity);
    table >> Sc[n] >> Sr[n] 
	  >> nu_mean[n] >>  nu_sem[n] >> lin_effort[n];
    getline(table,dummy);
    nu[n].sample(nu_mean[n]+sem_fac*nu_sem[n]);
    nu[n].sample(nu_mean[n]-sem_fac*nu_sem[n]);
    no_effort[n]=1;
    inv_effort[n]=1/lin_effort[n];
    log_effort[n]=log10(lin_effort[n]);
    REPORT(nu[n]);
    n++;
  }
  effort=inv_effort;
  double std_effort=1/(3e4);
  REPORT(diversity);
  REPORT(effort);

//   sequence<double> tmp=diversity;
//   diversity=effort;
//   effort=tmp;

  sequence< sequence<double> > X(3);

  X[0]=effort*0+1;
  X[1]=effort;
  X[2]=diversity;
  fitted_function nu_on_effort(effort,nu,2,false);
  fitted_function nu_on_lin_effort(lin_effort,nu,2,false);
  fitted_function nu_on_log_effort(log_effort,nu,2,false);
  fitted_function nu_on_no_effort(effort,nu,1,false);
  fitted_function nu_on_inv_effort(inv_effort,nu,2,false);
  fitted_function nu_on_diversity(diversity,nu,2,false);
  fitted_function nu_on_both(X,nu,false);

  sequence<average_meter> residual=nu;
  sequence<average_meter> residual0=nu;
  double SS=0,lin_SS=0,log_SS=0,no_SS=0,inv_SS=0;
  for(int i=0;i<n;i++){
    residual[i]-=nu_on_effort(effort[i]);
    residual0[i]-=nu_on_no_effort(effort[i]);
    SS+=residual[i]*residual[i]/(nu_sem[i]*nu_sem[i]);
    double res;
    res=nu[i]-nu_on_lin_effort(lin_effort[i]);
    lin_SS+=res*res/(nu_sem[i]*nu_sem[i]);
    res=nu[i]-nu_on_no_effort(effort[i]);
    no_SS+=res*res/(nu_sem[i]*nu_sem[i]);
    res=nu[i]-nu_on_inv_effort(inv_effort[i]);
    inv_SS+=res*res/(nu_sem[i]*nu_sem[i]);
    res=nu[i]-nu_on_log_effort(log_effort[i]);
    log_SS+=res*res/(nu_sem[i]*nu_sem[i]);

  }
  REPORT(SS);
  REPORT(lin_SS);
  REPORT(no_SS);
  REPORT(inv_SS);
  REPORT(log_SS);

  fitted_function residual_on_diversity(diversity,residual,2,false);
  fitted_function residual0_on_diversity(diversity,residual0,2,false);

  sequence<average_meter> log_corrected_nu(n);
  sequence<double> log_diversity(n);
  for(int i=n;i-->0;){
    log_corrected_nu[i].sample(log(nu_on_effort(std_effort)+residual_on_diversity(diversity[i]))+nu_sem[i]/(nu_on_effort(std_effort)+residual_on_diversity(diversity[i])));
    log_corrected_nu[i].sample(log(nu_on_effort(std_effort)+residual_on_diversity(diversity[i]))-nu_sem[i]/(nu_on_effort(std_effort)+residual_on_diversity(diversity[i])));
    log_diversity[i]=log(diversity[i]);
  }
  sequence< sequence<double> > Xlog=X;
  Xlog[2]=log_diversity;
  fitted_function power_law(Xlog,log_corrected_nu,false);
  
  cout << "nu_on_effort: " << nu_on_effort("f(effort)") << endl;
  cout << "nu_on_lin_effort: " << nu_on_lin_effort("f(effort)") << endl;
  cout << "nu_on_no_effort: " << nu_on_no_effort("f(effort)") << endl;
  cout << "nu_on_div: " << nu_on_diversity("div") << endl;
  cout << "residual_on_div: " << residual_on_diversity("div") << endl;
  cout << "residual0_on_div: " << residual0_on_diversity("div") << endl;

  ofstream corrected("corrected.dat");
  weighted_average_meter nu_av;
  weighted_average_meter diversity_av;
  for(int i=0;i<n;i++){
    double mean=nu[i]-nu_on_effort(effort[i])+nu_on_effort(std_effort);
    double error=nu_sem[i];
    corrected << diversity[i] << " "
	      << mean << " " 
	      << error << endl;
    nu_av.sample(mean,1/(error*error));
    diversity_av.sample(diversity[i],1/(error*error));
  }

  {// test here manual computation of beta, check if you understand the formula.
    double mean_x=diversity_av.readout();
    double mean_y=nu_av.readout();
    double Sxy=0,Sxx=0;
    for(int i=0;i<n;i++){
      double weight=1/(nu_sem[i]*nu_sem[i]);
      double dx=diversity[i]-mean_x;
      double mean=nu[i]-nu_on_effort(effort[i])+nu_on_effort(std_effort);
      double dy=mean-mean_y;
      Sxy+=dx*dy*weight;
      Sxx+=dx*dx*weight;
    }
    double beta=Sxy/Sxx;
    REPORT(beta);
    // compute CV propagator for errors in diversity:
    double propagator_sum=0;
    double dbeta_dx_hold=0;
    for(int i=0;i<n;i++){
      double weight=1/(nu_sem[i]*nu_sem[i]);
      double x=diversity[i];
      double dx=x-mean_x;
      double mean=nu[i]-nu_on_effort(effort[i])+nu_on_effort(std_effort);
      double dy=mean-mean_y;
      double dbeta_dx=(dy*weight-2*beta*dx*weight)/Sxx;
      dbeta_dx_hold+=dbeta_dx;
      if(i<n-1 && diversity[i+1]==diversity[i]){
	;//do nothing; continue accumulating dbeta_dx values
      }else{
	propagator_sum+=dbeta_dx_hold*dbeta_dx_hold*x*x;
	dbeta_dx_hold=0;
      }
    }
    cout << "additional error of beta = " << sqrt(propagator_sum)
	 << " * relative error of diversity " << endl;
  }
    

  sequence<average_meter> nu_null_model(nu);
  for(int i=0;i<n;i++){
    nu_null_model[i]-=nu_mean[i];
    REPORT(nu_null_model[i]);
  }


  double N=int(1e6);
  average_meter n0,ne0,ne1,nd0,nd1,rd0,rd1,pl0,pl1,pl2,ned0,ned1,ned2;
  average_meter pne1,pnd1_m,pnd1,prd1_m,prd1,pned1,pned2_m,pned2;
  average_meter ne0sig;
  average_meter pSS;
  average_meter pSS0;
  average_meter nu_av_null;
  sequence< sequence<double> > reg;
  sequence< double > exponent(N);
  for(int j=N;j-->0;){
    sequence<average_meter> nu_model(n);//(nu_null_model);
    typeof(diversity) perturbed_diversity=diversity;
    typeof(X) perturbed_X=X;
    for(int i=0;i<n;i++){
      double mean=gaussian(0,nu_sem[i]);
      nu_model[i].sample(mean+sem_fac*nu_sem[i]);
      nu_model[i].sample(mean-sem_fac*nu_sem[i]);
//       REPORT(nu_model[i]);
      //perturbed_diversity[i]*=exp(gaussian(0,0.3));
    }
    //perturbed_diversity[1]=perturbed_diversity[0];
//     REPORT(nu_model.size());
//     REPORT(effort.size());

    perturbed_X[2]=perturbed_diversity;

    fitted_function nu_model_on_effort(effort,nu_model,2,false);
    fitted_function nu_model_on_nothing(effort,nu_model,1,false);
    fitted_function nu_model_on_diversity(perturbed_diversity,nu_model,2,false);
    fitted_function nu_model_on_both(perturbed_X,nu_model,false);
    
    residual=nu_model;
    residual0=nu_model;
    double SS_model=0;
    double SS0_model=0;
    for(int i=0;i<n;i++){
      residual[i]-=nu_model_on_effort(effort[i]);
//       REPORT(residual[i]);
      SS_model+=residual[i]*residual[i]/(nu_sem[i]*nu_sem[i]);
      residual0[i]-=nu_model_on_nothing(effort[i]);
//       REPORT(residual[i]);
      SS0_model+=residual0[i]*residual0[i]/(nu_sem[i]*nu_sem[i]);
    }
    pSS.sample(SS_model>SS);
    pSS0.sample(SS0_model>no_SS);
    fitted_function residual_model_on_diversity(perturbed_diversity,residual,2,false);

//     cout << nu_model_on_effort("log10(effort)") << endl;
//     cout << nu_model_on_effort[0] << endl;
    weighted_average_meter nu_av_0;
    nu_av_0.sample(nu_model[1],1);
    nu_av_0.sample(nu_model[0],1);
    nu_av_0.sample(nu_model[6],1);
    //nu_av_null.sample(nu_av_0.readout());
    nu_av_null.sample(nu_model_on_effort(std_effort));

    ne0sig.sample(nu_model_on_effort[0]>0 ? 1 : -1);

    n0.sample(nu_model_on_nothing[0]);
    ne0.sample(nu_model_on_effort[0]);
    ne1.sample(nu_model_on_effort[1]);
    nd0.sample(nu_model_on_diversity[0]);
    nd1.sample(nu_model_on_diversity[1]);
    rd0.sample(residual_model_on_diversity[0]);
    rd1.sample(residual_model_on_diversity[1]);
    ned0.sample(nu_model_on_both[0]);
    ned1.sample(nu_model_on_both[1]);
    ned2.sample(nu_model_on_both[2]);

    pne1.sample(fabs(nu_model_on_effort[1])>fabs(nu_on_effort[1]));
    pnd1_m.sample(nu_model_on_diversity[1]<nu_on_diversity[1]);
    pnd1.sample(abs(nu_model_on_diversity[1])>fabs(nu_on_diversity[1]));
    prd1.sample(fabs(residual_model_on_diversity[1])>
		fabs(residual_on_diversity[1]));
    prd1_m.sample(residual_model_on_diversity[1]<
		  residual_on_diversity[1]);
    pned1.sample(fabs(nu_model_on_both[1])>
		 fabs(nu_on_both[1]) );
    pned2.sample(fabs(nu_model_on_both[2])>
		 fabs(nu_on_both[2]) );
    pned2_m.sample(nu_model_on_both[2]<
		   nu_on_both[2] );

    sequence<average_meter> log_corrected_nu_model(n);
    for(int i=n;i-->0;){
      sequence<double> xlog(2);
      xlog[0]=1;xlog[1]=effort[i];xlog[2]=log_diversity[i];
      double mean=log(exp(power_law(xlog))+
		      residual_model_on_diversity(perturbed_diversity[i]));
      double error=nu_sem[i]/mean;
      log_corrected_nu_model[i].sample(mean+error);
      log_corrected_nu_model[i].sample(mean-error);
    }
    // fitted_function power_law_model(Xlog,log_corrected_nu_model,false);
//     pl0.sample(power_law_model[0]);
//     pl1.sample(power_law_model[1]);
//     pl2.sample(power_law_model[2]);
//     exponent[j]=power_law_model[2];
    
#ifdef ENVELOPES
    int s=0;
    sequence<double> x(3);
    x[0]=1;
    x[1]=std_effort;//effort[n-1];
    for(double S=50;S<(max_diversity/100+1)*100;S+=25){
      //      reg[s][j]=exp(power_law_model(log(S)));//residual_model_on_diversity(S);
#if 1
      x[2]=S;
      reg[s][j]=nu_model_on_both(x)+nu_on_both(x);
#else
      x[2]=log(double(S));
      reg[s][j]=exp(power_law_model(x));
#endif
      s++;
    }
#endif
  } 
  REPORT(n0);REPORT(n0.std());
  REPORT(ne0sig);
  REPORT(ne0);REPORT(ne0.std());
  REPORT(ne1);REPORT(ne1.std());
  REPORT(pne1);
  REPORT(nd0);REPORT(nd0.std());
  REPORT(nd1);REPORT(nd1.std());
  REPORT(pnd1);
  REPORT(pnd1_m);
  REPORT(rd0);REPORT(rd0.std());
  REPORT(rd1);REPORT(rd1.std());
  REPORT(prd1);
  REPORT(prd1_m);
  REPORT(nu_on_both[0]);
  REPORT(nu_on_both[1]);
  REPORT(nu_on_both[2]);
  REPORT(ned0);REPORT(ned0.std());
  REPORT(ned1);REPORT(ned1.std());
  REPORT(ned2);REPORT(ned2.std());
  REPORT(pned1);
  REPORT(pned2);
  REPORT(pned2_m);
  REPORT(power_law[0]);REPORT(pl0);REPORT(pl0.std());
  REPORT(power_law[1]);REPORT(pl1);REPORT(pl1.std());
  REPORT(power_law[2]);REPORT(pl2);REPORT(pl2.std());
  sort(exponent.begin(),exponent.end());
  REPORT(exponent[int(0.975*N)]);
  REPORT(exponent[int(0.025*N)]);
  REPORT(SS);
  REPORT(pSS);
  REPORT(no_SS);
  REPORT(pSS0);
  REPORT(nu_av.readout());
  REPORT(nu_av.sample_std());
  REPORT(diversity_av.readout());
  REPORT(diversity_av.sample_std());
  REPORT(nu_av.sample_std()/diversity_av.sample_std());
  REPORT(nu_av.sample_std()/nu_av.readout());
  REPORT(diversity_av.sample_std()/diversity_av.readout());
  REPORT(nu_av.sample_std()/0.54);
  REPORT(nu_on_effort(std_effort));
  REPORT(nu_av_null);REPORT(nu_av_null.std());

  double elasticity=nu_on_both[2]*diversity_av.readout()/nu_av.readout();
  double elasticity_error=ned2.std()*diversity_av.readout()/nu_av.readout();
  REPORT(elasticity);
  REPORT(elasticity_error);

#ifdef ENVELOPES
  ofstream confidence("confidence.dat");
  int s=0;
  sequence<double> x(3);
  x[0]=1;
  x[1]=effort[n-1];
  for(double S=50;S<(max_diversity/100+1)*100;S+=25){
    x[2]=S;
    double offset=0;//nu_on_both(x);
    sort(reg[s].begin(),reg[s].end());
    confidence << S << " " 
	       << offset+reg[s][int(0.975*N)] << " " 
// 	       << offset+reg[s][int(0.95*N)] << " " 
// 	       << offset+reg[s][int(0.80*N)] << " " 
// 	       << offset+reg[s][int(0.5*N)] << " " 
// 	       << offset+reg[s][int(0.20*N)] << " " 
// 	       << offset+reg[s][int(0.05*N)] << " " 
	       << offset+reg[s][int(0.025*N)] << endl;
    s++;
  }
#endif
  WARNING("POWER-LAW estimates probably do not treat perturbed_diversity correctly");
}
