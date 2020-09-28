// -*- mode: c++ -*-
// $Id: MSY.cc 3 2005-12-01 07:13:32Z cvsrep $

#include "error.h"

#include "isres-variant.h" //overwrites the implementation of
			   //nlopt::GN_ISRES in nlopt.

#include "MSY.h"  
#include "Statistics.h"
#include "mass_integrator.h"
#include "snapshot.h"
#include <fstream>
#include <algorithm>
 
static my_evaluator_t eval_here;
static double catchability_cutoff = eval_here("1*gram");
double MSY_maxE=eval_here("5/year");
int MSY_penalize=false;
double MSY_decline_threshold=0.8;
double MSY_Threatened_threshold=eval_here("1e7*gram");
double MSY_rel_Threatened_threshold=0;
double MSY_productive_state_regularizer= 0e-10 ;
// isres_constraint_trickery depends on the fact that in our ISRES
// implementation the constraint is always evaluated just after the
// function call, so that we can save it from or evaluation of the
// function call and use it afterwards.  Of course, our trickery is not
// threat safe!
int isres_constraint_trickery=false;

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  {
    CFGDOUBLE(catchability_cutoff),
    CFGDOUBLE(MSY_maxE),
    CFGDOUBLE(MSY_decline_threshold),
    CFGDOUBLE(MSY_productive_state_regularizer),
    CFGDOUBLE(MSY_rel_Threatened_threshold),
    CFGINT(MSY_penalize),
    CFGINT(isres_constraint_trickery),
    {0, CFG_END, 0}
  };
static cfg_add dummy(cfg);

static sequence< double > 
fast_size_spectrum(const sequence<double> & B,
		   const sequence<double> & logM,
		   double bin_factor=100){
  
  double max_logM = *max_element(logM.begin(),logM.end());
  
  const double log_lowerM=log(catchability_cutoff);
  const double log_bin_factor=log(bin_factor);
  sequence< double > spectrum;
  
  for(int i=0;i<B.size();i++){
    if(logM[i]>=log_lowerM){
      spectrum[int((logM[i]-log_lowerM)/log_bin_factor)]+=
	B[i];
    }
  }
  return spectrum;
}

static sequence< double > 
fast_size_spectrum(const NewWeb & web, double bin_factor=100){

  sequence< double > B=web.get_biomass_B();
  sequence< double > logM=web.get_log_bodymass_M();

  // drop info on plants:
  B.resize(web.number_of_animals());
  logM.resize(web.number_of_animals());

  return fast_size_spectrum(B, logM, bin_factor);
}

//#define SPECTRAL_PENALIZATION

double MSY_t::worst_decline(const NewWeb & web){
#ifdef SPECTRAL_PENALIZATION // penalize based on size spectrum
  sequence< double > current_spectrum = 
    fast_size_spectrum(web);
  
  if(current_spectrum.size() < unperturbed_spectrum.size()){
    WARNING("spectrum truncated");
    return 1;
  }
  //REPORT(current_spectrum);// Trying to find out how worst_decline can be >1.
  //REPORT(unperturbed_spectrum);
  current_spectrum/=unperturbed_spectrum;
  //REPORT(current_spectrum);
  
  double worst_decline=1-
    *min_element(current_spectrum.begin(),current_spectrum.end());
#else  // penalize based on all species
  if(web.number_of_species() < this->web.number_of_species()){
    WARNING("species extinguished");
    return 1;
  }    
  double worst_decline=0;
  for(int i=web.number_of_species();i-->0;){
    worst_decline=
      std::max(worst_decline,
	       1 - 
	       web.s(i).biomass_abundance_B()/
	       this->web.s(i).biomass_abundance_B() );
  }
#endif
  REPORT(worst_decline);
  return worst_decline;
}

double MSY_t::sustainability_penalization(NewWeb & web){

  web.forbid_fishing();
  web.relax(eval("50*years"),NewWeb::species_set_t(),false);

  return simple_penalization(web,0.2/*=decline threshold*/);
}

double MSY_t::simple_penalization(NewWeb & web, double decline_threshold){

  const double decline=worst_decline(web);
  if(!isres_constraint_trickery){
    if(decline > decline_threshold){
      double penalization = exp(-10*(decline - decline_threshold));
      REPORT(penalization);
      return penalization;
    }
  }else{// with trickery
    saved_constraint = decline - decline_threshold;
    return 1; // no immediate penalization
  }
}

double equilibrium_yield(const std::vector<double> &x, 
			 std::vector<double> &grad, void* f_data){
  MSY_t & msy= *(MSY_t *)f_data;
  return msy(x,grad,f_data);
}

double isres_constraint(const std::vector<double> &x, 
			std::vector<double> &grad, void* f_data){
  MSY_t & msy= *(MSY_t *)f_data;
  REPORT(msy.saved_constraint);
  return msy.saved_constraint;
}  


MSY_t::MSY_t(const NewWeb & w,double maxE):
  web(w),max_E(maxE),best_Y_on_grid(0),
  unperturbed_spectrum(fast_size_spectrum(w))
{
};

void MSY_t::grid(const char * filename,int steps){
  std::ofstream os(filename);
  int depth=std::min(3,n())-1;
  grid(os,steps,depth,logE);
}
    
void MSY_t::grid(std::ostream & os,int steps,int depth,std::vector<double> logE){
  const char* sep[] = {" ","\n","\n"};
  std::vector<double> dummy;
  if(depth<0){
    double Y = operator()(logE,dummy,this);
    os << Y;
    os.flush();
    if(Y > best_Y_on_grid){
      best_Y_on_grid=Y;
      logE_on_grid=logE;
    }
  }else{
    for(int i=0;i<=steps;i++){
      logE[depth]=log(max_E)-log(2)*(steps-i);
      grid(os,steps,depth-1,logE);
      os << sep[depth];
    }
  }
}
    
void MSY_species::get_mortalities(const NewWeb & web, 
      			  std::vector<double> &E ){
  E.resize(n());
  for(int j=0;j<n();j++){
    E[j]=web.s[web.species_at(fished[j])].fishing_mortality();
  }
}

MSY_species::MSY_species(const NewWeb & w,std::vector<int> & f):
  MSY_t(w),fished(f) {
  if(fished.empty()){
    FATAL_ERROR("No species fished, cannot compute MSY.");
  }
}

void MSY_species::find_MSY(){
  // Create nlopt object: 
  //opt=nlopt::opt(nlopt::GN_DIRECT_L,n());
  opt=nlopt::opt(nlopt::GN_ISRES,n());
  opt.set_population(2*(n()+1));
  REPORT(opt.get_algorithm_name());

  double upper=log(max_E);
  double lower=log(1e-10*max_E);
  
  // Set magnitude of initial steps on log scale:
  double small_step=1;
    
  // Specify what to do:
  opt.set_max_objective(equilibrium_yield, (void*) this);
  opt.set_lower_bounds(lower);
  opt.set_upper_bounds(upper); //reasonable upper bound?
  //opt.set_maxtime(eval("5*hours/second"));
  opt.set_ftol_rel(-0.001);// negative value disables
  opt.set_xtol_rel(0.001);
  opt.set_initial_step(small_step);
  if(isres_constraint_trickery){
    opt.add_inequality_constraint(isres_constraint, (void*) this, 0.001);
  }
  
  // Find out if suitable initial fishing mortalities are set already:
  get_mortalities(web,E);
  double min_E_found=*std::min_element(E.begin(),E.end());

  REPORT(min_E_found);

  // Initialize optimization:
  if(min_E_found >= exp(lower)){
    logE.resize(n());
    for(int i=n();i-->0;){
      logE[i]=log(E[i]);
    }
  }else{
    logE=std::vector<double>(n(),lower);
  }

  // Optimize:
  try{
    result=opt.optimize(logE,msy);
    REPORT(result);
    REPORT(msy);
  }catch(std::runtime_error e){
    WARNING("runtime error exception while optimizing");
    E.resize(logE.size());
    for(int i=0;i<logE.size();++i){
      E[i]=exp(logE[i]);
    }
    return;
  }
  
  E.resize(logE.size());
  for(int i=0;i<logE.size();++i){
    E[i]=exp(logE[i]);
  }
  return;

  //////////////////////////////////////////////////////////
  std::cout << "*********************************************" << std::endl;
  std::cout << "***        Starting local search          ***" << std::endl;
  std::cout << "*********************************************" << std::endl;
  //////////////////////////////////////////////////////////

  // Optimize again with local search:
  opt = nlopt::opt(nlopt::LN_COBYLA,n());

  REPORT(opt.get_algorithm_name());

  // Specify what to do:
  opt.set_max_objective(equilibrium_yield, (void*) this);
  opt.set_lower_bounds(lower);
  opt.set_upper_bounds(upper); //reasonable upper bound?
  opt.set_ftol_rel(0.001);
  opt.set_xtol_rel(-1); // negative value disables
  opt.set_initial_step(0.1);

  // Optimize:
  try{
    result=opt.optimize(logE,msy);
    REPORT(result);
  }catch(std::runtime_error e){
    WARNING("runtime error exception while optimizing");
  }

  E.resize(logE.size());
  for(int i=0;i<logE.size();++i){
    E[i]=exp(logE[i]);
  }
}

void MSY_fleets::compute_catchabilities(){

  catchability.resize(catchability_function.size());

  for(int fleet=0;fleet<catchability.size();fleet++){

    REPORT(fleet);
    
    catchability[fleet].resize(web.species_at.size());

    for(int j=0;j<web.number_of_species();j++){ 
      my_evaluator_t eval_for_species;

      double Mmat = web.s[j].mean_bodymass_M();

      if(Mmat >= catchability_cutoff){
	
	eval_for_species.
	  set_variable("Mmat",Mmat);
	
	double c=
	  eval_for_species(catchability_function[fleet].c_str());

	REPORT(web.s[j].mean_bodymass_M());
	REPORT(c);
	
	catchability[fleet][web.assigned_column[j]]=c;
      }else{
	catchability[fleet][web.assigned_column[j]]=0;
      }
    }
  }
}


MSY_fleets::MSY_fleets(const NewWeb & w,std::vector<std::string> & cf):
  MSY_t(w),catchability_function(cf){
  compute_catchabilities();
}

void MSY_fleets::find_MSY(){
  // Create nlopt object: 
  //opt=nlopt::opt(nlopt::GN_DIRECT_L,n());
  opt=nlopt::opt(nlopt::LN_SBPLX,n());
  REPORT(opt.get_algorithm_name());

  double upper=log(max_E);
  double lower=log(1e-3*max_E);
  
  // Set magnitude of initial steps:
  double small_step=1;
    
  // Specify what to do:
  opt.set_max_objective(equilibrium_yield, (void*) this);
  opt.set_lower_bounds(lower);
  opt.set_upper_bounds(upper);
  opt.set_ftol_rel(0.1);
  opt.set_initial_step(small_step);

  // Initialize optimization:
  logE=std::vector<double>(n(),lower);

#if 1
  grid("Ygrid.dat",10);
  //exit(0);
#endif
  // Use experience on grid if there is:
  if(best_Y_on_grid > 0){
    logE=logE_on_grid;
  }

  // Optimize:
  try{
    result=opt.optimize(logE,msy);
    REPORT(msy);
  }catch(std::runtime_error e){
    WARNING("runtime error exception while optimizing");
    E.resize(logE.size());
    for(int i=0;i<logE.size();++i){
      E[i]=exp(logE[i]);
    }
    return;
  }

  
  //////////////////////////////////////////////////////////
  std::cout << "*********************************************" << std::endl;
  std::cout << "***        Starting local search          ***" << std::endl;
  std::cout << "*********************************************" << std::endl;
  //////////////////////////////////////////////////////////

  // Optimize again with local search:
  opt = nlopt::opt(nlopt::LN_COBYLA,n());

  REPORT(opt.get_algorithm_name());

  // Specify what to do:
  opt.set_max_objective(equilibrium_yield, (void*) this);
  opt.set_lower_bounds(lower); // Fishing mortalities must be positive!
  opt.set_upper_bounds(upper); //reasonable upper bound?
  opt.set_ftol_rel(0.01);
  opt.set_initial_step(0.1);

  // Optimize:
  try{
    result=opt.optimize(logE,msy);
  }catch(std::runtime_error e){
    WARNING("runtime error exception while optimizing");
  }
  
  E.resize(logE.size());
  for(int i=0;i<logE.size();++i){
    E[i]=exp(logE[i]);
  }
}

double MSY_t::operator()(const std::vector<double> &logE, 
		       std::vector<double> &grad, 
		       void* f_data, bool penalize){
  
  MSY_t & msy= *(MSY_t *)f_data;

  NewWeb web=msy.web;

  std::vector<double> E(logE.size());

  for(int i=0;i<logE.size();++i){
    E[i]=exp(logE[i]);
  }

  // Get yield:
  double yield=msy(web,E,grad);

  if(penalize){
    try{
      //yield*=sustainability_penalization(web);
      yield*=simple_penalization(web);
    }catch(terminal_condition t){
      std::cout << "caught " << t << std::endl;
      opt.force_stop();
    }
  }

  REPORT(yield);
  return yield;
}

void MSY_t::optimally_exploited_web(NewWeb & web){
  std::vector<double> dummy;
  std::cout << "E = ";
  for(int i=0;i<E.size();++i){
    std::cout << E[i] << " ";
  }
  std::cout << std::endl;
  set_mortalities(web,E);
    //operator()(web,E,dummy,(void *)this);
}


double MSY_t::operator()(NewWeb & web,
		       const std::vector<double> &E, 
		       std::vector<double> &grad){
  
  if(!grad.empty()){
    FATAL_ERROR("Can't use gradient methods for MSY yet.");
  }
    
  this->set_mortalities(web,E); // Virtual function.
    
  // Print mortalities:
  for(int i=0;i<n();i++){
    std::cout << E[i] << " ";
  }
  std::cout << std::endl;
  
  average_meter yield;
  try{
    // Relax:
    web.initialize_for_integration();
    web.relax(eval("200*years"),NewWeb::species_set_t(),false);
    
    // Get mean yield:
    double stretch=1;
    while( (yield.n() < 5 || yield.error()/yield > 0.01) && yield.n() < 100){
      web.relax(stretch*eval("-log(unirand)*years"),
		NewWeb::species_set_t(),false);
      yield.sample(web.biomass_yield());
      stretch=std::min(1.01*stretch,eval("10*years"));
    }
  }catch(terminal_condition t){
    std::cout << "caught " << t << std::endl;
    opt.force_stop();
  }
  if(exit_now){
    opt.force_stop();
  }

  // Print mortalities:
  for(int i=0;i<n();i++){
    std::cout << E[i] << " ";
  }
  std::cout << std::endl;
  
  return yield;
}

double MSY_species::worst_decline(const NewWeb & current){

  int S=fished.size();

  std::vector<double> rel_abundance(S);
  for(int i=S;i-->0;){
    int safi=current.species_at[fished[i]];
    if(safi == NewWeb::column_unused){
      return 1;
    }else{
      rel_abundance[i]=
	current.s(current.species_at[fished[i]]).biomass_abundance_B()/
	web.s(safi).biomass_abundance_B();
    }
  }
  
  return 
   1 - *std::min_element(rel_abundance.begin(),rel_abundance.end());
}

void MSY_species::set_mortalities(NewWeb & web, 
				    const std::vector<double> &E ){
  for(int j=0;j<n();j++){
    web.s[web.species_at[fished[j]]].
      set_fishing_mortality(E[j]);
  }
}

void MSY_fleets::set_mortalities(NewWeb & web, 
				    const std::vector<double> &E ){
  // Overwrite fishing mortalities
  for(int i=web.number_of_species();i-->0;){
    double F=0; // Fishing mortality
    for(int j=n();j-->0;){
      F += E[j] * catchability[j][web.assigned_column[i]];
    }
    web.s[i].set_fishing_mortality(F);
  }
}

void write_matrix(const NewMatrix & m,const char * file){
  std::ofstream f(file);
  for(int i=0;i<m.SIZE1();i++){
    for(int j=0;j<m.SIZE2();j++){
      f << m(i,j) << " ";
    }
    f << std::endl;
  }
}


void harvest_controller_t::
apply_fishing_mortalities(const ODE_vector & state, 
			  ODE_vector & time_derivative){

  NewVector F(fished.size());

  compute_fishing_mortalities(state,time_derivative,F);

  for(int i=fished.size();i-->0;){
    int safi=species_at[fished[i]];
    if(safi!=column_unused){
      time_derivative[safi]-=F[i];
    }
  }

  the_old_F = F; // perhaps there is a better context to put this line?
}

void harvest_controller_t::
set_fishing_mortalities(){

  const int S=number_of_variables();

  ODE_vector state(S),time_derivative(S);
  write_state_to(state);

  prepare_for_integration();
  if(NewWeb::dynamics(state,time_derivative)){
    FATAL_ERROR("problem evaluating dynamics");
  }

  NewVector F(fished.size());
  compute_fishing_mortalities(state,time_derivative,F);

  for(int i=fished.size();i-->0;){
    int safi=species_at[fished[i]];
    if(safi!=column_unused){
      s(safi).set_fishing_mortality(F[i]);
    }
  }
}


int harvest_controller_t::
dynamics(ODE_vector const & state, ODE_vector & time_derivative){

  int ret = NewWeb::dynamics(state,time_derivative);
  
  if(!constructing){
    this->apply_fishing_mortalities(state,time_derivative);
  }

  return ret;
}

void harvest_controller_t::recompute_LV_approximation(){

  int S=number_of_species();
  sequence<double> B=get_biomass_B();

  // Drop fished species that went extinct
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      for(int j=i; j<fished.size()-1;j++){
	the_old_F[j]=the_old_F[j+1];
      }
      the_old_F.resize(fished.size()-1);
      WARNING("erasing old_f " << i);
      fished.erase(fished.begin()+i);
    }
  }

  constructing=true;
  
  ODE_vector state(S),time_derivative(S);
  write_state_to(state);

  prepare_for_integration();
  if(NewWeb::dynamics(state,time_derivative)){
    FATAL_ERROR("problem evaluating dynamics");
  }

  std::cout << "Computing effective Jacobian, this may take a while." 
  	    << std::endl;
  std::cout.flush();
  // This J is for logarithmic abundances.  We need to keep this in mind!
  NewMatrix J=NewWeb::numerical_Jacobian();
  std::cout << "*"; std::cout.flush();
  
  // Split Bs and Fs into fish and unfished parts
  const int S_fished=fished.size();
  sequence<int>   is_fished(S,false);
  NewMatrix iBB_f(S_fished,S_fished);
  iBB_f = arma::diagmat(iBB_f);
  NewVector logB_dot_f(S_fished);
  the_production_rate.resize(S_fished);
  NewVector current_F(S_fished);
  the_old_B.resize(S_fished);
  for(int c=fished.size();c-->0;){
    int i=species_at[fished[c]];
    is_fished[i]=true;
    iBB_f(c,c)=1/B[i];
    logB_dot_f(c)=time_derivative[i];
    the_old_B(c)=B[i];
    current_F(c)=s(i).fishing_mortality();
  }

  const int S_unfished = S - S_fished;

  // index unfished species:
  std::vector<int> unfished(S_unfished);
   int i_u=S_unfished;
  for(int i=S;i-->0;){
    if(!is_fished[i]){
      unfished[--i_u]=i;
    }
  }

  // Now we have to split J into fished and unfished parts:
  NewMatrix J_ff(S_fished,S_fished);
  NewMatrix J_uu(S_unfished,S_unfished);
  NewMatrix J_uf(S_unfished,S_fished);
  NewMatrix J_fu(S_fished,S_unfished);

  for(int i=S_fished;i-->0;){
    for(int j=S_fished;j-->0;){
      J_ff(i,j) = J(species_at[fished[i]], species_at[fished[j]]);
    }
  }
  for(int i=S_fished;i-->0;){
    for(int j=S_unfished;j-->0;){
      J_fu(i,j) = J(species_at[fished[i]], unfished[j]);
    }
  }
  for(int i=S_unfished;i-->0;){
    for(int j=S_fished;j-->0;){
      J_uf(i,j) = J(unfished[i], species_at[fished[j]]);
    }
  }
  for(int i=S_unfished;i-->0;){
    for(int j=S_unfished;j-->0;){
      J_uu(i,j) = J(unfished[i], unfished[j]);
    }
  }
  std::cout << "*"; std::cout.flush();

  // Comput effective Jacobian for fished part:
  int err;
  NewMatrix Ji_uu=arma::inv(J_uu);
  if(err) FATAL_ERROR("Could not invert unfished Jacobian");
  //// Approximate unfished community as non-competing:
  //Ji_uu=ublas::diagonal_matrix<double>(Ji_uu);
  NewMatrix Jeff=
    -prod(NewMatrix(prod(J_fu,Ji_uu)),J_uf);
  std::cout << "*"; std::cout.flush();

  write_matrix(prod(Jeff,iBB_f),"Gindirect.dat");
  write_matrix(prod(J_ff,iBB_f),"Gdirect.dat");

  Jeff+=J_ff;

  //// Compare with single-species management:
  //Jeff=ublas::diagonal_matrix<double>(Jeff);

  // Now, compute effective interaction matrix for fished
  // community. The formula looks a bit odd because Jeff is for
  // logarithmic biomasses:
  the_interaction_matrix = prod(Jeff,iBB_f);
  std::cout << "*"; std::cout.flush();
  std::cout << "done." << std::endl;
  REPORT(current_F);

  REPORT(the_interaction_matrix);

  the_production_rate = logB_dot_f
    -prod(the_interaction_matrix,the_old_B) + current_F ;

  constructing=false;
}

harvest_controller_t::
harvest_controller_t(const MSY_species & msy) :
  MSY_species(msy),NewWeb(web),constructing(false){
  // Now compute interaction matrix:
  int S=number_of_species();

  forbid_fishing();
  initialize_for_integration();
  prepare_for_integration();

  recompute_strategy();
}

NewMatrix
harvest_controller_t::Gtilde_from_G(const NewMatrix & G){
  NewMatrix inv_G=arma::inv(G);
  NewMatrix Gtilde=
    arma::inv(arma::diagmat(inv_G));
  return Gtilde;
}

NewMatrix
harvest_controller_t::Gbar_from_G(const NewMatrix & G){
  return arma::diagmat(G);
}


transposed_interaction_controller_t::
transposed_interaction_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{
}

void transposed_interaction_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double bin_factor=100;
  const double log_bin_factor=log(bin_factor);
  const double log_lowerM=log(catchability_cutoff);

  sequence<double> B(fished.size());
  sequence<double> logB_dot(fished.size());
  sequence<double> logM(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B[i]=0;
      logB_dot[i]=0;
      logM[i]=log_lowerM;
    }else{
      B[i]=exp(state[safi]);
      logB_dot[i]=time_derivative[safi];
      logM[i]=s[safi].log_mean_bodymass_M();
    }
  }

  const NewVector vB(B);
  
  F=-prod(vB,the_interaction_matrix);
  
  if(MSY_penalize){
#if 1 // none of this works yet
#ifndef SPECTRAL_PENALIZATION
    FATAL_ERROR("For TIC only spectral penalization is implemented");
#else
    sequence< double > spectrum=
      fast_size_spectrum(B,logM,bin_factor);
    
    if(spectrum.size()!=unperturbed_spectrum.size()){
      FATAL_ERROR("size spectrum truncated");
    }

    // REPORT(spectrum);
    // REPORT(logB_dot);
    // REPORT(F);
    for(int i=spectrum.size();i-->0;){
      double badness=
	(unperturbed_spectrum[i] > 0
	 ?
	 (1-MSY_decline_threshold)-
	 spectrum[i]/unperturbed_spectrum[i]
	 :
	 -1 );
      // REPORT(badness);
      if(badness > 0){
	// Get total rate of decline in this size bin.
	weighted_average_meter rate_of_decline;

	const double bin_lower = log_lowerM+i*log_bin_factor;
	const double bin_upper = bin_lower + log_bin_factor;

	for(int i=fished.size();i-->0;){
	  if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	    rate_of_decline.sample(-logB_dot[i]+F[i],B[i]);
	  }
	}

	// REPORT(i);
	// REPORT(rate_of_decline.readout());

	if(rate_of_decline > 0){
	  const double correction=rate_of_decline*
	    (1+badness/(1-MSY_decline_threshold));
	  
	  
	  for(int i=fished.size();i-->0;){
	    if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	      F[i]-=correction;
	    }
	  }
	}
      }
    }
    
#endif
#endif
  }

  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream TIC_F("TIC_F.dat");
  TIC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    TIC_F << F[i] << " ";
  TIC_F << std::endl;
  
  double Y=inner_prod(F,vB);
  static std::ofstream TIC_yield("TIC_yield.dat");
  TIC_yield << current_time << " " << Y << std::endl;

  return;
}

productive_state_controller_t::
productive_state_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

#if 1 // variant with constraints:
void productive_state_controller_t::recompute_strategy(){

  recompute_LV_approximation();

  int S=the_interaction_matrix.SIZE1();
  the_target_B.resize(S);

  NewMatrix interaction_matrix = the_interaction_matrix;
  NewVector production_rate = the_production_rate;
  NewVector old_B = the_old_B;
  NewVector target_B= NewVector(S,0);
  std::vector<int> backmapper(S);

  for(int i=S;i-->0;) backmapper[i]=i;
  
  while(S > 0){// (usually exits through break command)
    // Repeat computation of strategy until no species violates constraint
    NewMatrix denominator=
      -trans(interaction_matrix)-interaction_matrix+
      MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

    denominator=arma::inv(denominator);
    
    NewVector numerator=
      production_rate + MSY_productive_state_regularizer*old_B;
    
    // This could be computed without a full matrix inversion!
    target_B = prod(denominator,numerator);
    
    std::vector<double> rel_target(S);
    for(int i=S;i-->0;){
      int safbi=species_at[fished[backmapper[i]]];
      rel_target[i]=target_B[i]/
	web.s(safbi).biomass_abundance_B();
    }

    int worst_i=std::min_element(rel_target.begin(),rel_target.end())-
      rel_target.begin();
    
    if(MSY_rel_Threatened_threshold){
      if(rel_target[worst_i] >= MSY_rel_Threatened_threshold){
	break;// exit while loop
      }
    }else{
      if(rel_target[worst_i] < 0){
	WARNING("targeting negative population biomass for " 
		<< backmapper[worst_i] );
      }
      break;// exit while loop
    }
    

    REPORT(target_B[worst_i]);
    REPORT(rel_target[worst_i]);// Pin biomass of worst_i to threshold
    WARNING("Pinning fished species " << backmapper[worst_i]);
    double pinnded_B = 
      MSY_rel_Threatened_threshold *
      web.s(species_at[fished[backmapper[worst_i]]]).biomass_abundance_B();
    the_target_B[backmapper[worst_i]] = pinnded_B;
      

    // Recompute LV approximation with worst_i removed:
    production_rate +=
      pinnded_B * COLUMN(interaction_matrix,worst_i);

    S=S-1;
    if(worst_i < S){
      ROW(interaction_matrix,worst_i)=ROW(interaction_matrix,S);
      COLUMN(interaction_matrix,worst_i)=COLUMN(interaction_matrix,S);
      interaction_matrix(worst_i,worst_i)=interaction_matrix(S,S);
      production_rate[worst_i]=production_rate[S];
      old_B[worst_i]=old_B[S];
      backmapper[worst_i]=backmapper[S];
    }
    interaction_matrix.resize(S,S);
    production_rate.resize(S);
    old_B.resize(S);
    backmapper.resize(S);
    
  }// end of while loop

  REPORT(old_B);
  REPORT(target_B);
  for(int i=S;i-->0;){
    the_target_B[backmapper[i]]=target_B[i];
  }

  REPORT(the_old_B);
  REPORT(the_target_B);
  
  return;
}
#else // variant without constraints
void productive_state_controller_t::recompute_strategy(){

  recompute_LV_approximation();

  const int S=the_interaction_matrix.SIZE1();
  
  NewMatrix denominator=
    -trans(the_interaction_matrix)-the_interaction_matrix+
    MSY_productive_state_regularizer*NewIdentityMatrix(S);
  
  int err=0;
  denominator=inverse(denominator,err);
  if(err) FATAL_ERROR("Could not invert denominator");
  
  NewVector numerator=
    the_production_rate + MSY_productive_state_regularizer*the_old_B;
  
  // This could be computed without a full matrix inversion!
  the_target_B = prod(denominator,numerator);

  for(int i=0;i<S;i++){
    //    if(the_target_B(i) < MSY_Threatened_threshold)
    //  the_target_B(i)= MSY_Threatened_threshold;
    int safi=species_at[fished[i]];
    double MSY_Thr_thr;
    if (safi!=column_unused){
     MSY_Thr_thr=web.s(safi).biomass_abundance_B()*MSY_rel_Threatened_threshold;
    }else{
     MSY_Thr_thr=0; 
    }  
    if(the_target_B(i) < MSY_Thr_thr)  the_target_B(i)= MSY_Thr_thr;
  }
  REPORT(the_old_B);
  REPORT(the_target_B);
  return;
}
#endif // variants with or without contraints

void productive_state_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double relaxation_rate=1/eval("1*year"); // not optimal in any sense!
  
  NewVector B(fished.size());
  NewVector logB_dot(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
      logB_dot(i)=0;
      F(i)=0;
    }else{
      B(i)=exp(state[safi]);
      logB_dot(i)=time_derivative[safi];
      F(i)=
	logB_dot(i) + relaxation_rate * (1-the_target_B(i)/B(i));
    }
  }
  
  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream PSC_F("PSC_F.dat");
  PSC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    PSC_F << F[i] << " ";
  PSC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream PSC_yield("PSC_yield.dat");
  PSC_yield << current_time << " " << Y << std::endl;
}


target_pressure_controller_t::
target_pressure_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void target_pressure_controller_t::recompute_strategy(){

  recompute_LV_approximation();

  const int S=the_interaction_matrix.SIZE1();
  
  NewMatrix denominator=
    -trans(the_interaction_matrix)-the_interaction_matrix+
    MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

  denominator=arma::inv(denominator);

  NewVector numerator=
    the_production_rate + MSY_productive_state_regularizer*the_old_B;

  // This could be computed without a full matrix inversion!
  NewVector target_B = prod(denominator,numerator);

  REPORT(target_B);
  REPORT(the_production_rate);

  // // Older, equivalent formulation
  // the_target_F = 
  //   the_production_rate + prod(the_interaction_matrix,target_B);

  the_target_F = 
    -prod(trans(the_interaction_matrix),target_B);

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);

  return;
}

void target_pressure_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }

  F=the_target_F;
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream TPC_F("TPC_F.dat");
  TPC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    TPC_F << F(i) << " ";
  TPC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream TPC_yield("TPC_yield.dat");
  TPC_yield << current_time << " " << Y << std::endl;

  return;
}


//----newest HCR : IPSC

individual_productive_state_controller_t::
individual_productive_state_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void individual_productive_state_controller_t::recompute_target_B(){

  int S=the_interaction_matrix.SIZE1();
  the_target_B.resize(S);
  
  NewMatrix interaction_matrix = the_interaction_matrix;
  NewVector production_rate = the_production_rate;
  NewVector old_B = the_old_B;
  NewVector target_B = NewVector(S,0);
  std::vector<int> backmapper(S);

  for(int i=S;i-->0;) backmapper[i]=i;

  while(S > 0){// (usually exits through break command)
    // Repeat computation of strategy until no species violates constraint
    NewMatrix Gtilde=
      Gtilde_from_G(interaction_matrix);
  
    NewMatrix denominator=
      -Gtilde-the_interaction_matrix+
      MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

    denominator=arma::inv(denominator);

    NewVector numerator=
      production_rate + MSY_productive_state_regularizer*old_B;
    
    // This could be computed without a full matrix inversion!
    target_B = prod(denominator,numerator);

    std::vector<double> rel_target(S);
    for(int i=S;i-->0;){
      int safbi=species_at[fished[backmapper[i]]];
      rel_target[i]=target_B[i]/
	web.s(safbi).biomass_abundance_B();
    }

    int worst_i=std::min_element(rel_target.begin(),rel_target.end())-
      rel_target.begin();
    
    if(MSY_rel_Threatened_threshold){
      if(rel_target[worst_i] >= MSY_rel_Threatened_threshold){
	break;// exit while loop
      }
    }else{
      if(rel_target[worst_i] < 0){
	WARNING("targeting negative population biomass for " 
		<< backmapper[worst_i] );
      }
      break;// exit while loop
    }
    REPORT(target_B[worst_i]);
    REPORT(rel_target[worst_i]);
    // Pin biomass of worst_i to threshold 
    WARNING("Pinning fished species " << backmapper[worst_i]);
    double pinnded_B = 
      MSY_rel_Threatened_threshold *
      web.s(species_at[fished[backmapper[worst_i]]]).
      biomass_abundance_B();
    the_target_B[backmapper[worst_i]] = pinnded_B;
    // Recompute LV approximation with worst_i removed:
    production_rate +=
      pinnded_B * COLUMN(interaction_matrix,worst_i);
    
    S=S-1;
    if(worst_i < S){
      ROW(interaction_matrix,worst_i)=ROW(interaction_matrix,S);
      COLUMN(interaction_matrix,worst_i)=COLUMN(interaction_matrix,S);
      interaction_matrix(worst_i,worst_i)=interaction_matrix(S,S);
      production_rate[worst_i]=production_rate[S];
      old_B[worst_i]=old_B[S];
      backmapper[worst_i]=backmapper[S];
    }
    interaction_matrix.resize(S,S);
    production_rate.resize(S);
    old_B.resize(S);
    backmapper.resize(S);
  }// end of while loop
  REPORT(old_B);
  REPORT(target_B);
  for(int i=S;i-->0;){
    the_target_B[backmapper[i]]=target_B[i];
  }
  
  REPORT(the_old_B);
  REPORT(the_target_B);
  return;
}

void individual_productive_state_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double relaxation_rate=1/eval("1*year");
  
  NewVector B(fished.size());
  NewVector logB_dot(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
      logB_dot(i)=0;
      F(i)=0;
    }else{
      B(i)=exp(state[safi]);
      logB_dot(i)=time_derivative[safi];
      F(i)=
	logB_dot(i) + relaxation_rate * (1-the_target_B(i)/B(i));
    }
  }
  
  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream IPSC_F("IPSC_F.dat");
  IPSC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    IPSC_F << F[i] << " ";
  IPSC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream IPSC_yield("IPSC_yield.dat");
  IPSC_yield << current_time << " " << Y << std::endl;
}

void soft_individual_productive_state_controller_t::
recompute_strategy(){

  individual_productive_state_controller_t::recompute_strategy();

  const int S=the_interaction_matrix.SIZE1();
  
  NewMatrix Gtilde=
    Gtilde_from_G(the_interaction_matrix);

  REPORT(the_old_F);
  if(the_old_F.size() != S){
    WARNING("the_old_F has wrong size, assuming zero");
    the_old_F.clear();
    the_old_F.resize(S);
  }

  NewVector estimated_logistic_growth_rate = 
    the_old_F - prod(Gtilde,the_old_B);
  the_relaxation_rate = 0.5 * estimated_logistic_growth_rate;

  for(int i=the_relaxation_rate.size();i-->0;){
    if(the_relaxation_rate(i)<0)
      the_relaxation_rate(i)*=-1;
  }

  REPORT(the_relaxation_rate);
  
  return;
}

void soft_individual_productive_state_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  NewVector logB_dot(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
      logB_dot(i)=0;
      F(i)=0;
    }else{
      B(i)=exp(state[safi]);
      logB_dot(i)=time_derivative[safi];
      F(i)=
	logB_dot(i) + the_relaxation_rate(i) * (B(i)/the_target_B(i)-1);
    }
  }
  
  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream SIPSC_F("SIPSC_F.dat");
  SIPSC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    SIPSC_F << F[i] << " ";
  SIPSC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream SIPSC_yield("SIPSC_yield.dat");
  SIPSC_yield << current_time << " " << Y << std::endl;
}

//---- ITPC, based on F=Gtilde * (G+Gtilde)^-1 * r

individual_target_pressure_controller_t::
individual_target_pressure_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void individual_target_pressure_controller_t::recompute_target_F(){

  const int S=the_interaction_matrix.SIZE1();
  
  REPORT(the_production_rate);

  NewMatrix Gtilde=
    Gtilde_from_G(the_interaction_matrix);
 
  NewMatrix denominator=
    -Gtilde-the_interaction_matrix+
    MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

  denominator=arma::inv(denominator);

  NewVector numerator=
    the_production_rate + MSY_productive_state_regularizer*the_old_B;

  // This could be computed without a full matrix inversion!
  NewVector target_B = prod(denominator,numerator);

  REPORT(target_B);
  REPORT(the_production_rate);

  // // Old, equivalent formulation:
  // the_target_F = 
  //   the_production_rate + prod(the_interaction_matrix,target_B);

  the_target_F = -prod(Gtilde,target_B);

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);
    
  
  return;
}

void individual_target_pressure_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }

  F=the_target_F;
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream ITPC_F("ITPC_F.dat");
  ITPC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    ITPC_F << F(i) << " ";
  ITPC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream ITPC_yield("ITPC_yield.dat");
  ITPC_yield << current_time << " " << Y << std::endl;

  return;
}



// individual TIC

individual_transposed_interaction_controller_t::
individual_transposed_interaction_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{
}

void individual_transposed_interaction_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double bin_factor=100;
  const double log_bin_factor=log(bin_factor);
  const double log_lowerM=log(catchability_cutoff);

  sequence<double> B(fished.size());
  sequence<double> logB_dot(fished.size());
  sequence<double> logM(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B[i]=0;
      logB_dot[i]=0;
      logM[i]=log_lowerM;
    }else{
      B[i]=exp(state[safi]);
      logB_dot[i]=time_derivative[safi];
      logM[i]=s[safi].log_mean_bodymass_M();
    }
  }

  const NewVector vB(B);    
    
  NewMatrix Gtilde=
    Gtilde_from_G(the_interaction_matrix);
  
  F=-prod(Gtilde,vB);
  
  if(MSY_penalize){
#if 1 // none of this works yet
#ifndef SPECTRAL_PENALIZATION
    FATAL_ERROR("For TIC only spectral penalization is implemented");
#else
    sequence< double > spectrum=
      fast_size_spectrum(B,logM,bin_factor);
    
    if(spectrum.size()!=unperturbed_spectrum.size()){
      FATAL_ERROR("size spectrum truncated");
    }

    // REPORT(spectrum);
    // REPORT(logB_dot);
    // REPORT(F);
    for(int i=spectrum.size();i-->0;){
      double badness=
	(unperturbed_spectrum[i] > 0
	 ?
	 (1-MSY_decline_threshold)-
	 spectrum[i]/unperturbed_spectrum[i]
	 :
	 -1 );
      // REPORT(badness);
      if(badness > 0){
	// Get total rate of decline in this size bin.
	weighted_average_meter rate_of_decline;

	const double bin_lower = log_lowerM+i*log_bin_factor;
	const double bin_upper = bin_lower + log_bin_factor;

	for(int i=fished.size();i-->0;){
	  if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	    rate_of_decline.sample(-logB_dot[i]+F[i],B[i]);
	  }
	}

	// REPORT(i);
	// REPORT(rate_of_decline.readout());

	if(rate_of_decline > 0){
	  const double correction=rate_of_decline*
	    (1+badness/(1-MSY_decline_threshold));
	  
	  
	  for(int i=fished.size();i-->0;){
	    if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	      F[i]-=correction;
	    }
	  }
	}
      }
    }
    
#endif
#endif
  }

  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream ITIC_F("ITIC_F.dat");
  ITIC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    ITIC_F << F[i] << " ";
  ITIC_F << std::endl;
  
  double Y=inner_prod(F,vB);
  static std::ofstream ITIC_yield("ITIC_yield.dat");
  ITIC_yield << current_time << " " << Y << std::endl;

  return;
}


// --- B-Nash based controllers

// coded the entire series by simply replacing 'tilde' with 'bar',
// '"I' with ''"IB', ' I' with '' IB' and 'individual_' 'individual_B_' in
// individual-based code above.

individual_B_productive_state_controller_t::
individual_B_productive_state_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void individual_B_productive_state_controller_t::recompute_target_B(){

  const int S=the_interaction_matrix.SIZE1();
  
  REPORT(the_production_rate);

  NewMatrix Gbar=
    Gbar_from_G(the_interaction_matrix);
  
  NewMatrix denominator=
    -Gbar-the_interaction_matrix+
    MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

  denominator=arma::inv(denominator);

  NewVector numerator=
    the_production_rate + MSY_productive_state_regularizer*the_old_B;

  // This could be computed without a full matrix inversion!
  the_target_B = prod(denominator,numerator);

  // !! Missing here: check if B falls below threshold or zero.
  // !! Keeping it as it is for now for consistency with previous
  // !! simulations.

  REPORT(the_old_B);
  REPORT(the_target_B);
  // for (int i=0;i<S;i++){
  //  int safi=species_at[fished[i]];  
  //  REPORT(i);REPORT(safi);REPORT(column_unused); 
  //  REPORT(web.s(safi).biomass_abundance_B());
  // }
  return;
}

void individual_B_productive_state_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double relaxation_rate=1/eval("1*year");
  
  NewVector B(fished.size());
  NewVector logB_dot(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
      logB_dot(i)=0;
      F(i)=0;
    }else{
      B(i)=exp(state[safi]);
      logB_dot(i)=time_derivative[safi];
      F(i)=
	logB_dot(i) + relaxation_rate * (1-the_target_B(i)/B(i));
    }
  }
  
  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream IBPSC_F("IBPSC_F.dat");
  IBPSC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    IBPSC_F << F[i] << " ";
  IBPSC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream IBPSC_yield("IBPSC_yield.dat");
  IBPSC_yield << current_time << " " << Y << std::endl;
}



//---- IBBTPC, based on F=Gbar * (G+Gbar)^-1 * r

individual_B_target_pressure_controller_t::
individual_B_target_pressure_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void individual_B_target_pressure_controller_t::recompute_target_F(){

  const int S=the_interaction_matrix.SIZE1();
  
  REPORT(the_production_rate);

  NewMatrix Gbar=
    Gbar_from_G(the_interaction_matrix);
 
  NewMatrix denominator=
    -Gbar-the_interaction_matrix+
    MSY_productive_state_regularizer*NewIdentityMatrix(S,S);

  denominator=arma::inv(denominator);

  NewVector numerator=
    the_production_rate + MSY_productive_state_regularizer*the_old_B;

  // This could be computed without a full matrix inversion!
  NewVector target_B = prod(denominator,numerator);

  REPORT(target_B);
  REPORT(the_production_rate);

  // // Old, equivalent formulation:
  // the_target_F = 
  //   the_production_rate + prod(the_interaction_matrix,target_B);

  the_target_F = -prod(Gbar,target_B);

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);
    
  
  return;
}

void individual_B_target_pressure_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }

  F=the_target_F;
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream IBTPC_F("IBTPC_F.dat");
  IBTPC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    IBTPC_F << F(i) << " ";
  IBTPC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream IBTPC_yield("IBTPC_yield.dat");
  IBTPC_yield << current_time << " " << Y << std::endl;

  return;
}



// individual TIC

individual_B_transposed_interaction_controller_t::
individual_B_transposed_interaction_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{
}

void individual_B_transposed_interaction_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  const double bin_factor=100;
  const double log_bin_factor=log(bin_factor);
  const double log_lowerM=log(catchability_cutoff);

  sequence<double> B(fished.size());
  sequence<double> logB_dot(fished.size());
  sequence<double> logM(fished.size());
  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B[i]=0;
      logB_dot[i]=0;
      logM[i]=log_lowerM;
    }else{
      B[i]=exp(state[safi]);
      logB_dot[i]=time_derivative[safi];
      logM[i]=s[safi].log_mean_bodymass_M();
    }
  }

  const NewVector(vB);
    
    
  NewMatrix Gbar=
    Gbar_from_G(the_interaction_matrix);
  
  F=-prod(Gbar,vB);
  
  if(MSY_penalize){
#if 1 // none of this works yet
#ifndef SPECTRAL_PENALIZATION
    FATAL_ERROR("For TIC only spectral penalization is implemented");
#else
    sequence< double > spectrum=
      fast_size_spectrum(B,logM,bin_factor);
    
    if(spectrum.size()!=unperturbed_spectrum.size()){
      FATAL_ERROR("size spectrum truncated");
    }

    // REPORT(spectrum);
    // REPORT(logB_dot);
    // REPORT(F);
    for(int i=spectrum.size();i-->0;){
      double badness=
	(unperturbed_spectrum[i] > 0
	 ?
	 (1-MSY_decline_threshold)-
	 spectrum[i]/unperturbed_spectrum[i]
	 :
	 -1 );
      // REPORT(badness);
      if(badness > 0){
	// Get total rate of decline in this size bin.
	weighted_average_meter rate_of_decline;

	const double bin_lower = log_lowerM+i*log_bin_factor;
	const double bin_upper = bin_lower + log_bin_factor;

	for(int i=fished.size();i-->0;){
	  if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	    rate_of_decline.sample(-logB_dot[i]+F[i],B[i]);
	  }
	}

	// REPORT(i);
	// REPORT(rate_of_decline.readout());

	if(rate_of_decline > 0){
	  const double correction=rate_of_decline*
	    (1+badness/(1-MSY_decline_threshold));
	  
	  
	  for(int i=fished.size();i-->0;){
	    if(logM[i] >= bin_lower and logM[i] < bin_upper and 
	     species_at[fished[i]] != column_unused){
	      F[i]-=correction;
	    }
	  }
	}
      }
    }
    
#endif
#endif
  }

  for(int i=fished.size();i-->0;){
    if(F(i)<0) F(i)=0;
  }

  static std::ofstream IBTIC_F("IBTIC_F.dat");
  IBTIC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    IBTIC_F << F[i] << " ";
  IBTIC_F << std::endl;
  
  double Y=inner_prod(F,vB);
  static std::ofstream IBTIC_yield("IBTIC_yield.dat");
  IBTIC_yield << current_time << " " << Y << std::endl;

  return;
}




//----growth rate controller

growth_rate_controller_t::
growth_rate_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void growth_rate_controller_t::recompute_target_F(){

  const int S=the_interaction_matrix.SIZE1();
  
  REPORT(the_production_rate);

  the_target_F = 0.5*the_production_rate;

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);
      
  return;
}

void growth_rate_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }

  F=the_target_F;
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream GRC_F("GRC_F.dat");
  GRC_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    GRC_F << F(i) << " ";
  GRC_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream GRC_yield("GRC_yield.dat");
  GRC_yield << current_time << " " << Y << std::endl;

  return;
}



//--- CFP controller

CFP_controller_t::
CFP_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void CFP_controller_t::recompute_target_F(){

  const int S=the_interaction_matrix.SIZE1();
  
  REPORT(the_production_rate);

  NewMatrix Gtilde=
    Gtilde_from_G(the_interaction_matrix);
 
  
  REPORT(the_old_F);
  if(the_old_F.size() != S){
    WARNING("the_old_F has wrong size, assuming zero");
    the_old_F.clear();
    the_old_F.resize(S);
  }
  
  // Rationale here:

  // - tildeG * B^2  is the non-linear term in the effective single-species LV model, where B is current biomass.

  // - the linear term is (r-F)*B, where F is current fishing mortality, and r is unknown.

  // - assuming equilibrium, 0= r - F + tildeG*B, so r = F - tildeG*B is the linear growth rate seen in single-species LV models.
  
  // - This would suggest to set Fmsy = r/2 = (F - tildeG*B)/2.
  
  // - CFP would then keep Fmsy fixed for each 50 year cycle in the model.

  NewVector estimated_logistic_growth_rate = 
    the_old_F - prod(Gtilde,the_old_B);
  the_target_F = 0.5 * estimated_logistic_growth_rate;

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);
    
  
  return;
}

void CFP_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }
  F=the_target_F;
  //REPORT(F);
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream CFP_F("CFP_F.dat");
  CFP_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    CFP_F << F(i) << " ";
  CFP_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream CFP_yield("CFP_yield.dat");
  CFP_yield << current_time << " " << Y << std::endl;

  return;
}

//--- DPF controller: data poor fishing, assumes an allometric
// relationship between size and optimal fishing effort.  Coupling
// constant is estimated a r/2.

DPF_controller_t::
DPF_controller_t(const MSY_species & msy):
  harvest_controller_t(msy)
{}

void DPF_controller_t::recompute_target_F(){

  const int S=the_interaction_matrix.SIZE1();
  
  NewMatrix Gtilde=
    Gtilde_from_G(the_interaction_matrix);
 
  REPORT(the_old_F);
  if(the_old_F.size() != S){
    WARNING("the_old_F has wrong size, assuming zero");
    the_old_F.clear();
    the_old_F.resize(S);
  }
  
  // Rationale here: see CFP_controller_t
  NewVector estimated_logistic_growth_rate = 
    the_old_F - prod(Gtilde,the_old_B);
  the_target_F = 0.5 * estimated_logistic_growth_rate;

  double exponent=get_cfg_parameter("respiration_rate_allometric_exponent");
  
  NewVector inv_scale_factor(S);
  for(int i=S;i-->0;){
    inv_scale_factor[i]=pow(s[web.species_at(fished[i])].bodymass(),
			    -exponent);
  }

  double mean_scaled_growth_rate=
    inner_prod(inv_scale_factor,estimated_logistic_growth_rate)/S;

  the_target_F.resize(S);
  for(int i=S;i-->0;){
    the_target_F[i]=(0.5*mean_scaled_growth_rate)/inv_scale_factor[i];
  }

  for(int i=0;i<S;i++)
    if(the_target_F(i) < 0)
      the_target_F(i)= 0;

  REPORT(the_target_F);
  
  return;
}

void DPF_controller_t::
compute_fishing_mortalities(const ODE_vector & state,
			    const ODE_vector & time_derivative,
			    NewVector & F){

  NewVector B(fished.size());
  const int S=the_interaction_matrix.SIZE1();

  for(int i=fished.size();i-->0;){
    const int safi=species_at[fished[i]]; // abbreviation
    if(safi==column_unused){
      B(i)=0;
    }else{
      B(i)=exp(state[safi]);
    }
  }
  F=the_target_F;
  //REPORT(F);
  
//  for(int i=0;i<S;i++)
//    if(B(i) < MSY_Threatened_threshold)
//      F(i)*=B(i)/MSY_Threatened_threshold;
  
  static std::ofstream DPF_F("DPF_F.dat");
  DPF_F << current_time << " ";
  for(int i=0;i<F.size();i++)
    DPF_F << F(i) << " ";
  DPF_F << std::endl;
  
  double Y=inner_prod(F,B);
  static std::ofstream DPF_yield("DPF_yield.dat");
  DPF_yield << current_time << " " << Y << std::endl;

  return;
}


