// -*- mode: c++ -*-
// $Id: snapshot.cc 2466 2016-05-01 23:27:44Z axel $

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>

#include "NewMatrix.h"

#include <math.h>

#include "snapshot.h"
#include "evaluate.h"
#include "matrix_transformers.h"
#include "xy_graph.h"
#include "Statistics.h"
#include "random.h"

using namespace std;

static my_evaluator_t eval_here;

const double gram=eval_here("1*gram"); //get a dimensionless gram
const double unit_mass=eval_here("1*kilogram"); //for output files
const double meter2=eval_here("1*meter^2");
const double unit_area=eval_here("1*meter^2"); //for output files
const double year=eval_here("1*year");


class level_comp{
  sequence<double> & level;
public:
  level_comp(sequence<double> & l):level(l){};
  bool operator()(int i,int j){
    return level[i]>level[j];
  }
};

Snapshot::Snapshot():S(0){};

void Snapshot::set_number_of_compartments(int s){
  ALWAYS_ASSERT(S==0);
  ALWAYS_ASSERT(s>0);
  S=s;
}
void Snapshot::set_biomasses(const sequence<double> & B){
  if(B.size()!=S){
    REPORT(B.size());
    REPORT(S);
    ALWAYS_ASSERT(B.size()==S);
  }
  ALWAYS_ASSERT(the_biomass_B.size()==0);
  the_biomass_B=B;
}
void Snapshot::set_bodymasses(const sequence<double> & M){
  ALWAYS_ASSERT(M.size()==S);
  ALWAYS_ASSERT(the_bodymass_M.size()==0);
  the_bodymass_M=M;
}
sequence< double > Snapshot::get_bodymasses(){
  return the_bodymass_M;
}
sequence< double > Snapshot::get_biomasses(){
  return the_biomass_B;
}
void Snapshot::set_flows(const link_strength_matrix & flow){
  ALWAYS_ASSERT(flow.size()==S);
  ALWAYS_ASSERT(the_flow.size()==0);
  the_flow=flow;
}
void Snapshot::set_links(const Interaction_Matrix & link){
  ALWAYS_ASSERT(the_link.size()==0);
  the_link=link;
  S=link.size();
}
void Snapshot::set_bottom(const sequence<int> & bottom){
  ALWAYS_ASSERT(bottom.size()==S);
  ALWAYS_ASSERT(the_bottom.size()==0);
  the_bottom=bottom;
}
void Snapshot::set_area(double area){
  the_area_per_compartment=area;
}
void Snapshot::adjust_links_given_flows(double thresh){
  ALWAYS_ASSERT(the_link.size()==0);
  compute_ifrac_if_required();
  the_link=threshold(the_ifrac,thresh);
}
void Snapshot::adjust_bottom_given_links(){
  ALWAYS_ASSERT(the_bottom.size()==0);
  WARNING("computing bottom from links");
  for(int i=S;i-->0;){
    the_bottom[i]=1;
    for(int j=S;j-->0;){
      if(the_link[i][j]==NetworkAnalysis::eats){
	the_bottom[i]=0;
	break;
      }
    }
  }
}
link_strength_matrix Snapshot::get_ifrac(){
  compute_ifrac_if_required();
  return the_ifrac;
}
const link_strength_matrix & Snapshot::ifrac(){
  compute_ifrac_if_required();
  return the_ifrac;
}
Interaction_Matrix Snapshot::get_im(){
  if(the_link.size()==0){
    adjust_links_given_flows(0);
  }
  return the_link;
}
int Snapshot::number_of_species(){
  return S;
}
int Snapshot::number_of_plants(){
  return sum(the_bottom);
}
int Snapshot::number_of_animals(){
  return S-number_of_plants();
}
int Snapshot::number_of_species_Mthreshold(double Mlowerthreshold, double Mupperthreshold){
  int S_Mthreshold=0;
  for(int i=S;i-->0;){
    if((the_bodymass_M[i]>Mlowerthreshold)&&(the_bodymass_M[i]<Mupperthreshold)){
      S_Mthreshold++;
    }
    //cout << "the_bodymass_M[i] = " << the_bodymass_M[i] << endl;  
  }
  return S_Mthreshold;
}
int Snapshot::number_of_species_TLthreshold(double TLlowerthreshold, double TLupperthreshold){
  int S_TLthreshold=0;
  compute_ifrac_if_required();
  compute_level_if_required();
  for(int i=S;i-->0;){
    if((the_level[i]>TLlowerthreshold)&&(the_level[i]<TLupperthreshold)){
      S_TLthreshold++;
    }
    //cout << "the_level[i] = " << the_level[i] << endl;  
  }
  return S_TLthreshold;
}
//const link_strength_matrix & get_intake_matrix();  //trash? 
link_strength_matrix Snapshot::intake_matrix(){
  FATAL_ERROR("TRASH?");
}              //trash?
void Snapshot::rank_abundance_plot(const char * filename,bool with_plants){
  std::ofstream os(filename);
  sequence<double> abund;
  for(int i=S;i-->0;){
    if(with_plants || ! the_bottom[i]){
      abund[abund.size()]=the_biomass_B[i]/the_bodymass_M[i];
    }
  }
  sort(abund.begin(),abund.end());
  for(int i=abund.size();i-->0;){
    os << abund.size()-i << " " << abund[i] << std::endl;
  }
}
void Snapshot::size_spectrum(const char * filename,
			     bool plants, bool animals){
  sequence<double> abundance;
  sequence<double> bodymass;
  double area=the_area_per_compartment/unit_area;
  
  int j=0;
  for(int i=S;i-->0;){
    if((the_bottom[i] && plants) || (!the_bottom[i] && animals)){
      abundance[j]=the_biomass_B[i]/the_bodymass_M[i]/
	area;
      bodymass[j]=the_bodymass_M[i]/unit_mass;
      j++;
    }
  }
  std::ofstream os(filename);
  os << log_spectrum(xy_graph(bodymass,abundance),10).log_xy();
}
void Snapshot::biomass_spectrum(const char * filename,
				bool plants, bool animals){
  sequence<double> biomass;
  sequence<double> bodymass;
  double area=the_area_per_compartment/unit_area;

  int j=0;
  for(int i=S;i-->0;){
    if((the_bottom[i] && plants) || (!the_bottom[i] && animals)){
      biomass[j]=the_biomass_B[i]/area/unit_mass;
      bodymass[j]=the_bodymass_M[i]/unit_mass;
      j++;
    }
  }
  std::ofstream os(filename);
  os << log_spectrum(xy_graph(bodymass,biomass),10).log_xy();
}

void Snapshot::trophic_level_structure(const char *filename){
  compute_level_if_required();
  std::ofstream os(filename);
  sequence<double> l=the_level;
  std::sort(l.begin(),l.end());
  for(int i=S;i-->0;){
    os << (S-i)/double(S) << " " << l[i] << endl;
  }
}
  

static sequence<double> sequence_double(const NewVector & v){
  // convert to sequence
  sequence<double> result(v.size());
  for(int i=v.size();i-->0;){
    result(i)=v[i];
  }
  
  return result;
}

static sequence<double> compute_level_from_ifrac_simple(const link_strength_matrix & ifrac){

  // prepare trophic level computation:
  int err;
  const int S=ifrac.size();
  NewVector ones(S);
  for(int i=S;i-->0;)
    ones[i]=1;

  WARNING("starting TL computation");  
#if 1
  // compute the trophic level vector:
  NewMatrix F=ifrac;
  NewMatrix D=NewIdentityMatrix(S,S)-F;
  NewVector level=arma::solve(D,ones);
#else
  NewMatrix F=ifrac;
  NewVector level=ones;
  NewVector level2;
  const int maxreps=20;
  int reps = 0;
  double diff = 0;
  do{
    level2=level;
    level=prod(F,level)+ones;
    diff=abs(level-level2);
    reps++;
  }while(diff>=0.0001 && reps<maxreps);
  if(reps == maxreps) WARNING("TL computation did not converge");
#endif
  WARNING("ending TL computation");


  return sequence_double(level);
}

static sequence<double> compute_level_from_ifrac(const link_strength_matrix & ifrac){

  // compute bottom from ifrac:
  sequence<int> bottom;
  const int S=ifrac.size();
  int n_animals=0;
  for(int i=S;i-->0;){
    bottom[i]=1;
    for(int j=S;j-->0;){
      if(ifrac[i][j]>0){
	bottom[i]=0;
	n_animals++;
	break;
      }
    }
  }


  const int Sa=n_animals;
  link_strength_matrix ifrac2;
  ifrac2.resize(Sa+1);
  
  int ii=0;
  for(int i=0;i<S;i++){
    if(!bottom[i]){
      int jj=0;
      double plant_fraction_sum=0;
      for(int j=0;j<S;j++){
	if(!bottom[j]){
	  ifrac2[ii][jj++]=ifrac[i][j];
	}else{
	  plant_fraction_sum+=ifrac[i][j];
	}
      }
      ifrac2[ii][Sa]=plant_fraction_sum;
      ii++;
    }
  }

  sequence<double> level2=compute_level_from_ifrac_simple(ifrac2);
  sequence<double> level(S);

  ii=0;
  for(int i=0;i<S;i++){
    if(!bottom[i]){
      level[i]=level2[ii++];
    }else{
      level[i]=1;
    }
  }
  
  return level;
}

void Snapshot::compute_ifrac_if_required(){
  if(the_ifrac.size()==0){
    if(the_flow.size()>0)
      the_ifrac=in_fraction(the_flow);
    else if(the_link.size()>0){
      WARNING("computing ifrac from links");
      the_ifrac=in_fraction(the_link);
    }
    else
      FATAL_ERROR("no data to compute ifrac");
  }
}

void Snapshot::compute_level_if_required(){
  if(the_level.size()==0){
    if(the_flow.size()){
      compute_ifrac_if_required();
      the_level=compute_level_from_ifrac(the_ifrac);
    }else{
      compute_level1_if_required();
      the_level=the_level1;
      WARNING("approximating h by h1");
    }
  }
}
void Snapshot::compute_level0_if_required(){
  if(the_level0.size()==0){
    ALWAYS_ASSERT(the_link.size()>0);
    link_strength_matrix ifrac0=in_fraction(the_link);
    the_level0=compute_level_from_ifrac(ifrac0);
  }
}
void Snapshot::compute_level1_if_required(){
  if(the_level1.size()==0){
    compute_level0_if_required();
    the_level1=
      0.5*(the_level0+the_link.shortest_path_level_vector());
    
  }
}
void Snapshot::compute_cStar_if_required(){
  if(the_cStar.size()==0){
    ALWAYS_ASSERT(the_link.size()>0);
    compute_ifrac_if_required();
    if(the_bottom.size()==0){
      adjust_bottom_given_links();
    }
    for(int i=number_of_species();i-->0;){
      if(the_bottom[i]){
	the_cStar[i]=0;
      }else{
	double c_sum=0;
	for(int j=number_of_species();j-->0;){
	  c_sum+=the_ifrac[i][j]*the_ifrac[i][j];
	}
	the_cStar[i]=c_sum;
      }
    }
  }
}


class rank_comp{
  sequence<double> & level;
  sequence<double> & prey_rank;
public:
  rank_comp(sequence<double> & l,sequence<double> & r):
    level(l),prey_rank(r){};
  bool operator()(int i,int j){
    if(int(level[i]+0.5) != int(level[j]+0.5)){
      return level[i]<level[j];
    }else{
     return prey_rank[i]<prey_rank[j];
    }
  }
};

void Snapshot::pajek_graph(const std::string & filename, double th){

  std::ofstream os(filename.c_str());
  compute_ifrac_if_required();
  compute_level_if_required();

  // Find permutation sorting species by trophic level.  The initial
  // "rectification" tries to find an ordering of species such that
  // the criss-cross in the plant-animal links in the food-web graph
  // is reduced.  However, this can take a bit of time.
#if 1  // Set this to zero if rectification() takes too much time.  
  permutation p=rectification();
#else
  permutation p(S);
  for(int i=S;i-->0;) p[i]=i;
#endif
  permutation rank=p.inverse();
  for(int r=2*max_level();r-->0;){
    sequence<double> mean_rank_of_prey(S);
    for(int i=number_of_animals();i-->0;){
      average_meter mean_rank;
      for(int j=S;j-->0;){
	if(the_ifrac[i][j]>=th && 
	   int(the_level[i]-0.5)==int(the_level[j]+0.5) )
	  mean_rank.sample(rank[j]);
      }
      if(mean_rank.n()) mean_rank_of_prey[i]=mean_rank;
    }
    stable_sort(p.begin(), 
		p.end(),
		rank_comp(the_level,mean_rank_of_prey) 
		);
  }


  // Compute how many species are at each trophic level
  sequence<int> species_at_level;
  for(int i=S;i-->0;){
    species_at_level[int(the_level[i]-0.5)]+=1;
  }
  REPORT(species_at_level);
  int max_l=max_level();

  // Compute starting positions for each level:
  double pos1=0;
  sequence<double> gap, pos(species_at_level.size(),pos1);
  gap=1.0/sequence<double>(species_at_level);
  pos+=0.5*gap;
  REPORT(gap);

  // Print output for Pajek
  // First the vertices, which are labelled with the original indices (dummyindex) for ease of reference
  os << "*Network " << endl;
  os << "*Vertices " << S << endl;

  for(int i=0;i<S;i++){
    int l=int(the_level[p[i]]+0.5);

    os << i+1 << " " << p[i]+1 << " ";
    os << pos[l-1] << " " 
       << 1-(the_level[p[i]]-0.8)/max_l << " " 
       << 0.5 << " ic ";
    switch(l){
    case 1:
      os << "Green"; break;
    case 2:
      os << "Yellow"; break;
    case 3:
      os << "Red"; break;
    case 4:
      os << "Blue"; break;
    case 5:
      os << "Black"; break;
    default:
      os << "White"; break;
    }
    pos[l-1]+=gap[l-1];
    os << " bc Black" << endl;
  } 

  // now the edges
  os << "*Arcslist " << endl;
  for(int i=0;i<S;i++){
    os << i+1;

    for(int k=0;k<S;k++){
      if(the_ifrac(p[i],p[k])>=th){ os << " " << k+1; }
    }
    os << endl;
  }
}
      
void Snapshot::species_table(const char * filename,
			     bool with_plants,
			     double th,
			     const sequence<string> & additional_headings,
			     const sequence<string> & additional_data){
  double area=the_area_per_compartment/unit_area;
  link_strength_matrix & intake = the_flow;

  compute_ifrac_if_required();
  compute_level_if_required();
  compute_level0_if_required();
  //  compute_level1_if_required();
  compute_cStar_if_required();

  // compute the mean trophic level and level sharpness:
  double max_level=0;
  average_meter mean_level;
  weighted_average_meter level_sharpness;
  for(int i=S;i-->0;){
    mean_level.sample(the_level[i]);
    if(the_level[i]>1.5)
      level_sharpness.sample(cos(2*M_PI*the_level[i]),1);
    if(the_level[i]>max_level) max_level=the_level[i];
  }
  cout << "max_level " << max_level << endl;
  cout << "mean_level " << mean_level.readout() << endl;
  cout << "level_sharpness " << level_sharpness.readout() << endl;
    
  // compute the trophic level vector based on topolgy:
  if(the_link.size()==0)
    adjust_links_given_flows(th);


  // compute generalities
  sequence<int> Gen;
  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      if(the_ifrac(j,i)>=th) Gen[i]++;
    }
  }
  
  // compute vulnerability
  sequence<int> Vul;
  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      if(the_ifrac(i,j)>=th) Vul[i]++;
    }
  }
  
  // compute prey log bodymass, log total biomass_B:
  NewVector log_bodymass(S), log_biomass_abundance_B(S), biomass(S);
  for(int n=S;n-->0;){
    log_bodymass[n]=log10(the_bodymass_M[n]/unit_mass);
    log_biomass_abundance_B[n]=log10(the_biomass_B[n]/unit_mass/area);
    biomass[n]=the_biomass_B[n]/unit_mass/area;
  }

  cout << "biomass_density " << sum(the_biomass_B)/unit_mass/area << endl;
  cout << "plant_biomass_density " << sum(the_biomass_B*the_bottom)/unit_mass/area << endl;
  cout << "animal_biomass_density " << sum(the_biomass_B*(1-the_bottom))/unit_mass/area << endl;

  NewMatrix F=the_ifrac;
  NewVector log_prey_bodymass=prod(F,log_bodymass);
  NewVector log_prey_biomass=prod(F,log_biomass_abundance_B);
  NewVector prey_biomass=prod(F,biomass);
  
  std::cout << "writing a species table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram; "
	    << "area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "meaning of colums is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": plant (1 if plant, 0 if animal)" << std::endl;
  std::cout << "column " << column++ << ": logM  (log10 mean bodymass)" << std::endl;
  std::cout << "column " << column++ << ": logB  (log10 total biomass density)" << std::endl;
  std::cout << "column " << column++ << ": logN  (log10 population)" << std::endl;
  std::cout << "column " << column++ << ": h     (trophic level)" << std::endl;
  std::cout << "column " << column++ << ": h0    (trophic level based on topolgy)" << std::endl;
  //  std::cout << "column " << column++ << ": h1    (short-weighted trophic level)" << std::endl;
  std::cout << "column " << column++ << ": cStar (specialization)" << std::endl;
  std::cout << "column " << column++ << ": logD  (log10 population density)" << std::endl;
  std::cout << "column " << column++ << ": Gen   (number of prey contributing " << th*100 << "% or more to diet)" << std::endl;
  std::cout << "column " << column++ << ": Vul   (number of predators it feeds to " << th*100 << "% or more)" << std::endl;
  std::cout << "column " << column++ << ": logMr (mean log resource bodymass)" << std::endl;
  std::cout << "column " << column++ << ": logBr (mean log resource total biomass density)" << std::endl;
  std::cout << "column " << column++ << ": logrB (log mean resource total biomass density)" << std::endl;
  std::cout << "column " << column++ << ": line  (line number)" << std::endl;
  //std::cout << "column " << column++ << ": pline (line number of main prey)" << std::endl;
  std::cout << "column " << column++ << ": recti (index in rectification)" << std::endl;
  for(int l=0;l<additional_headings.size();l++){
    std::cout << "column " << column++ << ": "
	      << additional_headings[l] << std::endl;
  }
  
  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */


  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(7);

  for(int n=0;n<S;n++){
    if(the_biomass_B[n]<the_bodymass_M[n]){
      WARNING("the_biomass_B[n]<the_bodymass_M[n] , i.e. " << the_biomass_B[n] << " < " << the_bodymass_M[n] );
    }
    //ALWAYS_ASSERT(the_biomass_B[n]>the_bodymass_M[n]);
    if(with_plants || !the_bottom[n]){
      os.width(1);
      os << the_bottom[n] << " ";
      os.width(9);
      os << log10(the_bodymass_M[n]/unit_mass) << " ";
      os.width(9);
      os << log10(biomass[n]) << " ";
      os.width(9);
      os << log10(the_biomass_B[n]/the_bodymass_M[n]) << " ";
      os.width(9);
      os << the_level[n] << " ";
      os.width(9);
      os << the_level0[n] << " ";
      os.width(9);
      //      os << the_level1[n] << " ";
      os << the_cStar[n] << " ";
      os.width(9);
      os << log10(the_biomass_B[n]/the_bodymass_M[n]/area) << " ";
      os.width(5);
      os << Gen[n] << " ";
      os.width(5);
      os << Vul[n] << " ";
      os.width(9);
      os << log_prey_bodymass[n] << " ";
      os.width(9);
      os << log_prey_biomass[n] << " ";
      os.width(9);
      os << log10(prey_biomass[n]) << " ";
      int maxi=0;
      for(int i=0;i<S;i++){
	if(F(n,i)>F(n,maxi))
	  maxi=i;
      }
      os << n+1 /*s[maxi].unique_id()*/ << " ";
      //os << maxi+1 /*s[maxi].unique_id()*/ << " ";
      os << (n<the_rectification.size()?the_rectification[n]:-1) << " ";
      os << additional_data[n];
      os << std::endl;
    }
  }
}

void Snapshot::species_table2(const char * filename,
			     bool with_plants,
			     double th,
			     const sequence<string> & additional_headings,
			     const sequence<string> & additional_data){
  double area=the_area_per_compartment/unit_area;
  
  // compute prey log bodymass, log total biomass_B:
  NewVector log_bodymass(S), log_biomass_abundance_B(S), biomass(S);
  for(int n=S;n-->0;){
    log_bodymass[n]=log10(the_bodymass_M[n]/unit_mass);
    log_biomass_abundance_B[n]=log10(the_biomass_B[n]/unit_mass/area);
    biomass[n]=the_biomass_B[n]/unit_mass/area;
  }
  
  std::cout << "writing a species table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram; "
	    << "area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "meaning of colums is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": logM  (log10 mean bodymass)" << std::endl;
  std::cout << "column " << column++ << ": logB  (log10 total biomass density)" << std::endl;
  for(int l=0;l<additional_headings.size();l++){
    std::cout << "column " << column++ << ": "
	      << additional_headings[l] << std::endl;
  }
  
  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */


  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(7);

  for(int n=0;n<S;n++){
    if(the_biomass_B[n]<the_bodymass_M[n]){
      WARNING("the_biomass_B[n]<the_bodymass_M[n] , i.e. " << the_biomass_B[n] << " < " << the_bodymass_M[n] );
    }
    //ALWAYS_ASSERT(the_biomass_B[n]>the_bodymass_M[n]);
    if(log10(the_bodymass_M[n]/unit_mass)>-3){
      if(with_plants){
	os.width(1);
	os << log10(the_bodymass_M[n]/unit_mass) << " ";
	os.width(9);
	os << log10(biomass[n]) << " ";
	os << additional_data[n];
	os << std::endl;
      }
    }
  }
}

void Snapshot::species_table3(const char * filename,
			     bool with_plants,
			     double th,
			     const sequence<string> & additional_headings,
			     const sequence<string> & additional_data){
  double area=the_area_per_compartment/unit_area;
  
  // compute prey log bodymass, log total biomass_B:
  NewVector log_bodymass(S), log_biomass_abundance_B(S), biomass(S);
  for(int n=S;n-->0;){
    log_bodymass[n]=log10(the_bodymass_M[n]/unit_mass);
    log_biomass_abundance_B[n]=log10(the_biomass_B[n]/unit_mass/area);
    biomass[n]=the_biomass_B[n]/unit_mass/area;
  }
  
  std::cout << "writing a species table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram; "
	    << "area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "meaning of colums is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": logM  (log10 mean bodymass)" << std::endl;
  std::cout << "column " << column++ << ": logB  (log10 total biomass density)" << std::endl;
  for(int l=0;l<additional_headings.size();l++){
    std::cout << "column " << column++ << ": "
	      << additional_headings[l] << std::endl;
  }
  
  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */


  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(7);

  for(int n=0;n<S;n++){
    if(the_biomass_B[n]<the_bodymass_M[n]){
      WARNING("the_biomass_B[n]<the_bodymass_M[n] , i.e. " << the_biomass_B[n] << " < " << the_bodymass_M[n] );
    }
    //ALWAYS_ASSERT(the_biomass_B[n]>the_bodymass_M[n]);
    if(log10(the_bodymass_M[n]/unit_mass)>-3){
      if(with_plants){
	os.width(1);
	os << log10(the_bodymass_M[n]/unit_mass) << " ";
	os.width(9);
	os << log10(biomass[n]) << " ";
	os << additional_data[n];
	os << std::endl;
      }
    }
  }
}

void Snapshot::link_table(const char * filename,
			  double th){
  double area=the_area_per_compartment/unit_area;
  
  if(the_flow.size()==0 && the_link.size()!=0)
    the_flow=the_link;

  compute_ifrac_if_required();
  compute_level_if_required();
  
  std::cout << "writing a link table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram;"
	    << " area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "the meaning of the colums is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": hc (consumer trophic level)" << std::endl;
  std::cout << "column " << column++ << ": hr (resource trophic level)" << std::endl;
  std::cout << "column " << column++ << ": Mc (log10 consumer body mass)" << std::endl;
  std::cout << "column " << column++ << ": Mr (log10 resource body mass)" << std::endl;
  std::cout << "column " << column++ << ": Bc (log10 consumer biomass density)" << std::endl;
  std::cout << "column " << column++ << ": Br (log10 resource biomass density)" << std::endl;
  std::cout << "column " << column++ << ": F  (biomass flow density [per year])" << std::endl;
  std::cout << "column " << column++ << ": f  (consumer diet fraction)" << std::endl;
  
  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */
  
  
  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(4);
  
  if(the_flow.size()==0) th=0;

  for(int n=0;n<S;n++){
    for(int m=0;m<S;m++){
      if(the_ifrac(m,n)>th){
	os.width(9);
	os << the_level[n] << " ";
	os.width(9);
	os << the_level[m] << " ";
	os.width(9);
	os << log10(the_bodymass_M[n]/unit_mass) << " ";
	os.width(9);
	os << log10(the_bodymass_M[m]/unit_mass) << " ";
	os.width(9);
	os << log10(the_biomass_B[n]/unit_mass/area) << " ";
	os.width(9);
	os << log10(the_biomass_B[m]/unit_mass/area) << " ";
	os.width(9);
	os << the_flow(m,n)/unit_mass/area/(1/year) << " ";
	os.width(9);
	os << the_ifrac(m,n) << " ";
	os << std::endl;
      }
    }
  }
}

void Snapshot::link_table2(const char * filename,
			  double th){
  double area=the_area_per_compartment/unit_area;
  
  if(the_flow.size()==0 && the_link.size()!=0)
    the_flow=the_link;

  compute_ifrac_if_required();
  compute_level_if_required();
  
  std::cout << "writing a link table to file \"" << filename << "\"" << std::endl;
  std::cout << "mass unit: " << unit_mass/gram << " gram;"
	    << " area unit: " << unit_area/meter2 << " meter^2" 
	    << std::endl;
  std::cout << "the meaning of the colums is " << std::endl;
  int column=1;
  std::cout << "column " << column++ << ": hc (consumer trophic level)" << std::endl;
  std::cout << "column " << column++ << ": hr (resource trophic level)" << std::endl;
  std::cout << "column " << column++ << ": Mc (log10 consumer body mass)" << std::endl;
  std::cout << "column " << column++ << ": Mr (log10 resource body mass)" << std::endl;
  std::cout << "column " << column++ << ": Bc (log10 consumer biomass density)" << std::endl;
  std::cout << "column " << column++ << ": Br (log10 resource biomass density)" << std::endl;
  std::cout << "column " << column++ << ": F  (biomass flow density [per year])" << std::endl;
  std::cout << "column " << column++ << ": f  (consumer diet fraction)" << std::endl;
  
  /* template:
     std::cout << "column " << column++ << ": ()" << std::endl;
  */
  
  
  std::ofstream os(filename);
  //os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(4);
  
  if(the_flow.size()==0) th=0;

  for(int n=0;n<S;n++){
    for(int m=0;m<S;m++){
      if(the_ifrac(m,n)>th){
	os.width(9);
	os << the_level[n] << " ";
	os.width(9);
	os << the_level[m] << " ";
	os.width(9);
	os << log10(the_bodymass_M[n]/unit_mass) << " ";
	os.width(9);
	os << log10(the_bodymass_M[m]/unit_mass) << " ";
	os.width(9);
	os << log10(the_biomass_B[n]/unit_mass/area) << " ";
	os.width(9);
	os << log10(the_biomass_B[m]/unit_mass/area) << " ";
	os.width(9);
	os << the_flow(m,n)/unit_mass/area/(1/year) << " ";
	os.width(9);
	os << the_ifrac(m,n) << " ";
	os << std::endl;
      }
    }
  }
}


bool Snapshot::is_plant(int i){
  if(the_bottom.size()==0){
    if(the_link.size()==0)
      adjust_links_given_flows(0);
    adjust_bottom_given_links();
  }
  return the_bottom[i];
}
 
  
void Snapshot::fast_stats(){
  // compute data similar to NewWeb::fast_stats()

  // print some stats that are easily computed
  cout << "n_animals " << number_of_animals() << endl;
  cout << "n_plants " << number_of_plants() << endl;

  if(S<=0){
    WARNING("no species");
  }

  double max_log_bodymass=log(the_bodymass_M[0]);
  average_meter log10_plant_bodymass;
  average_meter log10_animal_bodymass;
  double log10=log(10);
  double log10unitmass=log(unit_mass)/log10;
  for(int i=S;i-->0;){
    if(is_plant(i)){
      log10_plant_bodymass.sample(log(the_bodymass_M[i])/log10);
    }else{
      log10_animal_bodymass.sample(log(the_bodymass_M[i])/log10);
    }
    if(log(the_bodymass_M[i]) > max_log_bodymass){
      max_log_bodymass=log(the_bodymass_M[i]);
    }
  }
  cout << "max_log10_bodymass " << max_log_bodymass/log10-log10unitmass<< endl;
  cout << "log10_animal_bodymass " << log10_animal_bodymass.readout() -log10unitmass
       << " +/- " << log10_animal_bodymass.std() << endl;
  cout << "log10_plant_bodymass " << log10_plant_bodymass.readout() -log10unitmass
       << " +/- " << log10_plant_bodymass.std() << endl;
  weighted_average_meter log10_bodymass;
  log10_bodymass+=log10_plant_bodymass;
  log10_bodymass+=log10_animal_bodymass;
  cout << "log10_bodymass " << log10_bodymass.readout() -log10unitmass
       << " +/- " << log10_bodymass.sample_std() << endl;
}


double Snapshot::Omega(permutation p,double level_threshold){
  
  if(the_bottom.size()==0)
    adjust_bottom_given_links();

  const int S=the_ifrac.size();
  const double oS=1.0/S,oSa=1.0/(S-sum(the_bottom));

  // ATTENTION: p must be a permutation of size S!!

  double omega=0;
  int count=0;
  for(int con=S;con-->0;){
    if(the_level(p[con])>=level_threshold){
      count++;
      for(int res=S;res-->0;){
	double fac=(con*oSa-res*oS);
	omega+=the_ifrac(p[res],p[con])*fac*fac;
      }
    }
  }

  return omega/count;
}

double Snapshot::Delta_Omega(const permutation & p, int i, int j){
  //BUG ??: In a few cases, the values of Delta_Omega and differences
  //between omega values between exchanging i and j do not agree.  I
  //do not understand why.  In most cases, they seem to agree.

  
  const int S=the_ifrac.size();
  const double oS=1.0/S,oSa=1.0/(S-sum(the_bottom));

  // ATTENTION: p must be a permutation of size S!!

  double delta=0;
  for(int con=S;con-->0;){
    double faci=(con*oSa-i*oS),facj=(con*oSa-j*oS);
    delta+=
      (the_ifrac(p[i],p[con])-the_ifrac(p[j],p[con]))*
      (facj*facj-faci*faci);
  }

  for(int res=S;res-->0;){
    double faci=(i*oSa-res*oS),facj=(j*oSa-res*oS);
    delta+=
      (the_ifrac(p[res],p[i])-the_ifrac(p[res],p[j]))*
      (facj*facj-faci*faci);
  }

  return delta/S;
}
  

permutation Snapshot::rectification(){
  // Metropolis algorithm to minimize Omega by permutation of species.

  const int S=the_ifrac.size();
  permutation p(S);
  for(int i=S;i-->0;) p[i]=i;
  
  compute_ifrac_if_required();
  compute_level_if_required();

  stable_sort(p.begin(),p.end(),level_comp(the_level));

  int pure_herbivore_start=S,plant_start=S;
  for(int i=0;i<S;i++){
    if(the_level(p[i])<2.001){
      pure_herbivore_start=i;
      break;
    }
  }
  for(int i=pure_herbivore_start;i<S;i++){
    if(the_level(p[i])<1.5){
      plant_start=i;
      break;
    }
  }
  
  const int n_plants=S-plant_start;
  const int n_herbivores=plant_start-pure_herbivore_start;
    
  double omega=Omega(p);
  double cooldown_factor=1-1e-5;

  for(double T=1;T>1e-11;T*=cooldown_factor){
    int i=random_integer(n_plants+n_herbivores)+pure_herbivore_start;
    int j;
    if(i>=plant_start){
      j=random_integer(n_plants)+plant_start;
    }else{
      j=random_integer(n_herbivores)+pure_herbivore_start;
    }

//     REPORT(i);
//     REPORT(j);
    
    double delta=Delta_Omega(p,i,j);

    //double last_omega=Omega(p);
    
    if(delta<0 || exp(-delta/T)>unirand()){
      //accept step
      omega+=delta;
      {int hold=p[i];p[i]=p[j];p[j]=hold;}
    }      

    if(int(30+log10(T))!=int(30+log10(T*cooldown_factor*cooldown_factor))){
      // reporting
//       cerr << T << " " << i << " " << j << " " << p[i] << " " << p[j] << " " << delta 
// 	   << " " << Omega(p)-last_omega << "            \r";
      cerr << int(log10(T)+30.5)-30 //<< " " << sqrt(omega) 
	   << " " << sqrt(Omega(p)) << "            \r";
      cerr.flush();
    }
//     REPORT(T);
//     REPORT(omega);
  }
  
  //cout << "ranked_niche_width " << sqrt(Omega(p)) << endl;
  //  cerr << endl;
  the_rectification=p.inverse();
  return p;
}

double Snapshot::max_level(){
  compute_ifrac_if_required();
  compute_level_if_required();
  return *max_element(the_level.begin(),the_level.end());
}

double Snapshot::biomass_weighted_mean_level(){
  compute_ifrac_if_required();
  compute_level_if_required();
  return dot(the_level,the_biomass_B)/sum(the_biomass_B);
}

double Snapshot::mean_level(){
  compute_ifrac_if_required();
  compute_level_if_required();
  return sum(the_level)/the_level.size();
}

double Snapshot::Simpson_s_diversity_D(bool with_plants, bool with_animals){
  //const int S=the_ifrac.size();
  double B_square_sum=0;
  double B_sum=0;
  for(int i=S;i-->0;){
    if((the_bottom[i]?with_plants:with_animals)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_square_sum+=B*B;
    }
  }
  return B_square_sum/(B_sum*B_sum);
}

double Snapshot::Simpson_s_diversity_D_Mthreshold(double Mlowerthreshold, double Mupperthreshold){
  //const int S=the_ifrac.size();
  double B_square_sum=0;
  double B_sum=0;
  for(int i=S;i-->0;){
    if((the_bodymass_M[i]>Mlowerthreshold)&&(the_bodymass_M[i]<Mupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_square_sum+=B*B;
    }
  }
  return B_square_sum/(B_sum*B_sum);
}

double Snapshot::Simpson_s_diversity_D_TLthreshold(double TLlowerthreshold, double TLupperthreshold){
  const int S=the_ifrac.size();
  double B_square_sum=0;
  double B_sum=0;
  compute_ifrac_if_required();
  compute_level_if_required();
  for(int i=S;i-->0;){
    if((the_level[i]>TLlowerthreshold)&&(the_level[i]<TLupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_square_sum+=B*B;
    } 
  }
  return B_square_sum/(B_sum*B_sum);
}

double Snapshot::Shanon_diversity_H(bool with_plants, bool with_animals){
  //const int S=the_ifrac.size();
  double B_log_B_sum=0;
  double B_sum=0;
  for(int i=S;i-->0;){
    if((the_bottom[i]?with_plants:with_animals)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_log_B_sum+=B*log(B);
    }
  }
  return -B_log_B_sum/B_sum+log(B_sum);
}

double Snapshot::Shanon_diversity_H_Mthreshold(double Mlowerthreshold, double Mupperthreshold){
  //const int S=the_ifrac.size();
  double B_log_B_sum=0;
  double B_sum=0;
  for(int i=S;i-->0;){
    if((the_bodymass_M[i]>Mlowerthreshold)&&(the_bodymass_M[i]<Mupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_log_B_sum+=B*log(B);
    }
  }
  return -B_log_B_sum/B_sum+log(B_sum);
}
    
double Snapshot::Shanon_diversity_H_TLthreshold(double TLlowerthreshold, double TLupperthreshold){
  const int S=the_ifrac.size();
  double B_log_B_sum=0;
  double B_sum=0;
  compute_ifrac_if_required();
  compute_level_if_required();
  for(int i=S;i-->0;){
    if((the_level[i]>TLlowerthreshold)&&(the_level[i]<TLupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_log_B_sum+=B*log(B);
    }
  }
  return -B_log_B_sum/B_sum+log(B_sum);
}

double Snapshot::MTI_Mthreshold(double Mlowerthreshold, double Mupperthreshold){
  const int S=the_ifrac.size();
  double B_TL_sum=0;
  double B_sum=0;
  compute_ifrac_if_required();
  compute_level_if_required();
  for(int i=S;i-->0;){
    if((the_bodymass_M[i]>Mlowerthreshold)&&(the_bodymass_M[i]<Mupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_TL_sum+=B*the_level[i];
    }
  }
  return B_TL_sum/B_sum;
}

double Snapshot::MTI_TLthreshold(double TLlowerthreshold, double TLupperthreshold){
  const int S=the_ifrac.size();
  double B_TL_sum=0;
  double B_sum=0;
  compute_ifrac_if_required();
  compute_level_if_required();
  for(int i=S;i-->0;){
    if((the_level[i]>TLlowerthreshold)&&(the_level[i]<TLupperthreshold)){
      double B=the_biomass_B[i];
      B_sum+=B;
      B_TL_sum+=B*the_level[i];
    }
  }
  return B_TL_sum/B_sum;
}
