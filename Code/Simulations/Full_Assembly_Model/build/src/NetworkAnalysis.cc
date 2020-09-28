// $Id: NetworkAnalysis.cc 2466 2016-05-01 23:27:44Z axel $

#include "NetworkAnalysis.h"
#include "NetworkHelpers.h"
#include <iomanip>
#include <fstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include "random.h"
#include "error.h"
#include "Statistics.h"
#ifdef WITH_GSL
#include <gsl/gsl_sf.h>              // special functions
#endif

using namespace std;


int Interaction_Matrix::Number_of_Species_S() const
{
  return the_size;
}


int Interaction_Matrix::Number_of_Links_L()
{
  int count=0;
  for(int i=0;i<the_size;i++)
    for(int j=0;j<the_size;j++)
      if((*this)[i][j]!=NetworkAnalysis::none) count++;
  return count;
}

void Interaction_Matrix::label_species_by_index(){
#if defined(NO_POPULATION_DYNAMICS) || defined(NEW_POPULATION_DYNAMICS)
#else
  for(int i=0;i<the_size;i++){
    the_species[i].push_back(new Index_Species(i));
  }
#endif
}


#ifndef NO_POPULATION_DYNAMICS
static double total_biomass_sum(const Species_List & l){
  double s=0;
  for(Species_List::const_iterator i=l.begin();
      i!=l.end();
      i++){
    s+=(*i)->biomass_abundance_B();
  }
  return s;
}
#endif

Interaction_Matrix::Interaction_Matrix(const Interaction_Matrix& matrix):
  IMatrix(matrix),the_size(matrix.size()),
  the_species(matrix.the_species)
{
};

Interaction_Matrix const & Interaction_Matrix::operator=
(Interaction_Matrix const& matrix){
  //ASSERT(this->the_size == matrix.the_size); ??? What is this good for ???
  if(this != &matrix){
    the_size=matrix.the_size;
    the_species=matrix.the_species;
    this->IMatrix::operator=(matrix);
  }
  return (*this);
};


void Interaction_Matrix::Print(std::ostream &os) const
{
  os << std::setw(3);
  
  // horizontal scale:
  os << "WEB   ";
  for(int j=0;j<the_size;j++){
    os << (j/10)%10;
  }
  os << endl;
  os << "WEB   ";
  for(int j=0;j<the_size;j++){
    os << (j)%10;
  }
  os << endl;

  for(int i=0;i<the_size;i++){
    os << "WEB"<< ((i%100)<10?"0":"") << i%100<<" ";
    for(int j=0;j<the_size;j++){
      string s=".";
      if((*this)[j][i])
	s[0]='0'+(*this)[j][i];
      os << s ;
    }

#ifndef NO_POPULATION_DYNAMICS
    os << "  " << total_biomass_sum(the_species[i]);

    for(Species_List::const_iterator j=the_species[i].begin();
	j!=the_species[i].end();
	j++){
#if !defined(NEW_POPULATION_DYNAMICS)
      Species_with_Individuals * rs=dynamic_cast<Species_with_Individuals *>(*j);
#else
      NewSpecies * rs = *j;
#endif
      if(rs){
	os << "  " << rs->bodymass();
      }
      else
	os << "  " << "xx";
    }
#endif
    os << endl;
  }// for i
}


void Interaction_Matrix::PPrint(std::ostream &os) const
{
  os << std::setw(3);
  
  // horizontal scale:
//   for(int j=0;j<the_size;j++){
//     os << (j)%10;
//       }
//   os << endl;

  for(int i=0;i<the_size;i++){
//     os << i%10<<" ";
    os << "WIB ";
    for(int j=0;j<the_size;j++){
      string s=".";
      if((*this)[j][i])
	s[0]='X';
      os << s ;
    }

    os << endl;
  }// for i
  os << "WIB " << endl;

}


istream & operator>>(istream &is, Interaction_Matrix & im){
  int S;
  string l;

  is >> S;
  im=Interaction_Matrix(S);
  for(int i=0;i<S;i++){
    is >> l;
    if(l.size()!=S){
      FATAL_ERROR("input data has wrong format");
    }
    for(int j=S;j-->0;){
      im[j][i]=(l[j]=='.' ? 
		NetworkAnalysis::none : 
		NetworkAnalysis::eats );
    }
  }
  return is;
}


void Interaction_Matrix::tPrint(std::ostream &os) const
{
  os << the_size << endl;
  
  for(int i=0;i<the_size;i++){
    for(int j=0;j<the_size;j++){
      os << ((*this)[j][i] ? 'X' : '.') ;
    }
    os << endl;
  }
}

const string PGM_magic="P5"; 

void Interaction_Matrix::pgm_write(const char * filename)
{
  std::ofstream os(filename);
  
  os << PGM_magic << endl;
  os << the_size << endl; //width
  os << the_size << endl; //height
  os << 1 << endl; //maximum graylevel, i.e. black/white image

  for(int i=0;i<the_size;i++){
    for(int j=0;j<the_size;j++){
      int value=1;
      if((*this)[j][i])
	value=0;
      os << char(value);
    }
  }// for i
}


double Interaction_Matrix::connectance_C(){
  return double(Number_of_Links_L())/
    (double(Number_of_Species_S())*double(Number_of_Species_S()));
}
  
double Interaction_Matrix::links_per_species_Z(){
  return double(Number_of_Links_L())/
    (double(Number_of_Species_S()));
}

// Interaction_Matrix Interaction_Matrix::tsort(){
//   string tmp_file_name = "/tmp/NAtmp.dat";
//   FILE * TS = popen(("tsort - > "+tmp_file_name+" 2> /dev/null").c_str(),"w");
//   for(int i=0;i<the_size;i++){
//     fprintf(TS,"%i ",i);
//     fprintf(stdout,"%i ",i);
//     for(int j=0;j<the_size;j++){
//       if((*this)[i][j]){
// 	fprintf(TS,"%i ",j);
// 	fprintf(stdout,"%i ",j);
//       }
//     }
//     fprintf(TS,"\n",i);
//     fprintf(stdout,"\n",i);
//   }
//   pclose(TS);

//   TS=fopen(tmp_file_name.c_str(),"r");
//   int perm[the_size];
//   cout << "perm: ";
//   for(int i=0;i<the_size;i++){
//     fscanf(TS,"%i",&perm[i]);
//     cout << perm[i] << " " ;
//   }
//   cout << endl;
//   fclose(TS);

//   IMatrix pm(the_size,the_size);
//   for(int i=0;i<the_size;i++){
//     for(int j=0;j<the_size;j++){
//       //      pm[i][j]=(*this)[perm[i]][perm[j]];
//       pm[perm[i]][perm[j]]=(*this)[i][j];
//     }
//   }
//   Interaction_Matrix new_im(the_size);

//   for(int i=0;i<the_size;i++){
//     for(int j=0;j<the_size;j++){
//       (new_im)[i][j]=pm[i][j];
//     }
//   }
//   return new_im;
// }

Interaction_Matrix 
Interaction_Matrix::permute(permutation new_pos) const{
  ASSERT((int)new_pos.size()==the_size);
  Interaction_Matrix pm(the_size);
  for(int i=0;i<the_size;i++){
#ifndef NO_POPULATION_DYNAMICS
    pm.the_species[new_pos[i]]=the_species[i];
#endif
    for(int j=0;j<the_size;j++){
      pm[new_pos[i]][new_pos[j]]=(*this)[i][j];
    }
  }
  return pm;
}

Interaction_Matrix 
Interaction_Matrix::select(simple_vector<bool> &sel) const {
  ASSERT((int) sel.size()==the_size);

  int n_sels=count_if(sel.begin(),sel.end(),
		      bind2nd(equal_to<bool>(),true) );
  
  // make interaction matrix only for sels:
  Interaction_Matrix sel_im(n_sels);
  int i,j,k,l;
  for(i=0,j=0; i<the_size; i++)
    if(sel[i]){
#ifndef NO_POPULATION_DYNAMICS
      sel_im.the_species[j]=the_species[i];
#endif
      for(k=0,l=0; k<the_size; k++)
	if(sel[k]){
	  NetworkAnalysis::Interaction ii=
	    (*this)[i][k];
	  sel_im[j][l]=ii;
	  l++;
	}
      j++;
    }
  return sel_im;
}

Interaction_Matrix 
Interaction_Matrix::remove_lowest_level(){
  simple_vector<bool> sel(the_size);
  // make interaction matrix only for sels:
  for(int i=0; i<the_size; i++){
    for(int j=0; j<the_size ;j++){//simple_vectorized
      if( (*this)[i][j] != NetworkAnalysis::none ){
	sel[i]=true;
	break;
      }
    }
  }
  return select(sel);
}

#if 1 //new definition of lump_lowest_level:
Interaction_Matrix 
Interaction_Matrix::lump_lowest_level(){
  simple_vector<bool> lowest_level(the_size,true);
  simple_vector<bool> selection(the_size,true);
  //this sometimes yields the complete web as the lowest level!!:

  // make interaction matrix only for sels:
  for(int i=0; i<the_size; i++){
    for(int j=0; j<the_size ;j++){//simple_vectorized
      if( eats(i,j) ){
	lowest_level[i]=false;
	break;
      }
    }
  }

  int j=0; // the first species in lowest level
  while(j<the_size && !lowest_level[j]) j++;

  if(j==the_size){
    REPORT(the_size);
    PPrint();
    WARNING("there is no lowest level in food web");
    return *this;
  }

  Interaction_Matrix m2=*this;
  // lump together
  for(int i=0;i<the_size;i++){
    if(i>j && lowest_level[i]){
      selection[i]=false;
#ifndef NO_POPULATION_DYNAMICS
      m2.the_species[j].insert(m2.the_species[j].end(),
			       m2.the_species[i].begin(),
			       m2.the_species[i].end() );
#endif
      for(int k=0;k<the_size;k++){
	NetworkAnalysis::Interaction ii=m2[k][i];
	NetworkAnalysis::Interaction ij=m2[k][j];
	if(ii>ij) m2[k][j]=ii;
      }
    }
  }
  
  return m2.select(selection);
}
#else
Interaction_Matrix 
Interaction_Matrix::lump_lowest_level(){
  simple_vector<bool> lowest_level(the_size);
  simple_vector<bool> selection(the_size,true);
  //this sometimes yields the complete web as the lowest level!!:
  find_lowest_level(*this,lowest_level);

  int j=0; // the first species in lowest level
  while(!lowest_level[j]) j++;

  Interaction_Matrix m2=*this;
  // lump together
  for(int i=0;i<the_size;i++){
    if(i>j && lowest_level[i]){
      selection[i]=false;
#ifndef NO_POPULATION_DYNAMICS
      m2.the_species[j].insert(m2.the_species[j].end(),
			       m2.the_species[i].begin(),
			       m2.the_species[i].end() );
#endif
      for(int k=0;k<the_size;k++){
	NetworkAnalysis::Interaction ii=m2[k][i];
	NetworkAnalysis::Interaction ij=m2[k][j];
	if(ii>ij) m2[k][j]=ii;
      }
    }
  }
  
  return m2.select(selection);
}
#endif

int Interaction_Matrix::prop_Cy4() const{
  return Cy4(*this);
}

Interaction_Matrix Interaction_Matrix::randomize(){
  Interaction_Matrix im(size());
  int S=size();
  int L=Number_of_Links_L();
  for(int i=0;i<L;i++){
    int x,y;
    do{
      x=int(random_double(S));
      y=int(random_double(S));
    }while(im[x][y] || x==y);
    im[x][y]=NetworkAnalysis::eats;
  }
  return im;
}

Interaction_Matrix Interaction_Matrix::tsort(){
  //return forgiving_tsort(*this);
  graphs_stats_data stats=
    graph_analyze(*this,false/*don't do multiplicities*/);
  return stats.sorted;
}

Interaction_Matrix Interaction_Matrix::tsort2(){
  const int unset=-1;
  permutation perm(the_size);
  
  for(int i=0;i<the_size;i++){
    perm[i]=unset;
  }
  
  for(int i=0;i<the_size;i++){
    //find species for i-th position:

    //find an unset species j that is eaten by most n other unset k, keeping
    //n as small as possible:
    for(int n=0;n<2*the_size;n++){
      for(int j=0;j<the_size;j++){
 	if(perm[j]==unset){
 	  int count=0;
 	  for(int k=0;k<the_size;k++){
 	    if(k!=j && perm[k]==unset){
 	      //count+=(*this)[j][k]==NetworkAnalysis::eats;
 	      count+=((int)(*this)[k][j]==(int)NetworkAnalysis::eats);
 	    }
 	  }// for k
 	  if(count==n){
 	    //register j as i-the species:
 	    perm[j]=i;
 	    goto next_i;
 	  }// if
 	}// if
      }// for j
    }// for n
     //a species can have at most the_size-1 predators, therefore:
    FATAL_ERROR("This point should not be reached!");
  next_i:;
  }// for i

  return (*this).permute(perm);
}

#ifndef NO_POPULATION_DYNAMICS
double average_mass(Species_List l){
  double sum=0;
  for(Species_List::iterator s=l.begin();
      s!=l.end();
      s++){
#if defined(NEW_POPULATION_DYNAMICS)
    sum+=(*s)->bodymass();
#else
    if(Species_with_Individuals *si=
       dynamic_cast<Species_with_Individuals *>(*s) ){
      sum+=si->bodymass();
    }
    else
      WARNING("Species without mass!!");
#endif
  }
  return sum/l.size();
}
#endif

#ifndef NO_POPULATION_DYNAMICS
class higher_average_mass_in {
private:
  Interaction_Matrix const & the_matrix;
public:
  higher_average_mass_in(Interaction_Matrix const & m):the_matrix(m){};
  bool operator()(int s1, int s2){
    return 
      (average_mass(the_matrix.the_species[s1]) >=
       average_mass(the_matrix.the_species[s2])    );
  }
};

			 
Interaction_Matrix Interaction_Matrix::msort(){
  // sort by mass
  permutation perm(the_size);
  for(int i=0;i<the_size;i++){
    perm[i]=i;
  }
  stable_sort(perm.begin(),perm.end(),higher_average_mass_in(*this));
  permutation p2=perm.inverse();
  return this->permute(p2);
}
#endif
  

class random_ranking {
private:
  simple_vector<double> ranking;
public:
  random_ranking(int size):ranking(simple_vector<double>(size)){
    for(int i=0;i<size;i++){
      ranking[i]=random_double(1);
    }
  }
  bool operator()(int s1, int s2){
    return 
      (ranking[s1] >= ranking[s2] );
  }
};

			 
Interaction_Matrix Interaction_Matrix::random_shuffle(){
  // sort by mass
  permutation perm(the_size);
  for(int i=0;i<the_size;i++){
    perm[i]=i;
  }
  stable_sort(perm.begin(),perm.end(),random_ranking(the_size));
  return this->permute(perm);
}
  

class precedes_in {
private:
  Interaction_Matrix const & the_matrix;
public:
  precedes_in(Interaction_Matrix const & m):the_matrix(m){};
  bool operator()(int s1, int s2){
    int N=the_matrix.size();
    for(int i=0;i<N;i++){
      if(the_matrix[i][s1] < the_matrix[i][s2])
	return true;
      else if (the_matrix[i][s1] > the_matrix[i][s2])
	return false;
    }
    for(int i=0;i<N;i++){
      if(the_matrix[s1][i] < the_matrix[s2][i])
	return true;
      else if (the_matrix[s1][i] > the_matrix[s2][i])
	return false;
    }
    return false;
  }
};

class equal_in {
private:
  Interaction_Matrix const & the_matrix;
public:
  equal_in(Interaction_Matrix const & m):the_matrix(m){};
  bool operator()(int s1, int s2){
    int N=the_matrix.size();
    for(int i=0;i<N;i++){// simple_vectorized
      if(the_matrix[i][s1] != the_matrix[i][s2])
	return false;
    }
    for(int i=0;i<N;i++){// simple_vectorized
      if(the_matrix[s1][i] != the_matrix[s2][i])
	return false;
    }
    return true;
  }
};

Interaction_Matrix Interaction_Matrix::trophic(){
  // condense an Interaction_Matrix to a trophical food web

  permutation perm(the_size);
  for(int i=0;i<the_size;i++){
    perm[i]=i;
  }

  // sort so that species to be lumped together follow each other:
  //sort(perm.begin(),perm.end(),precedes_in(*this)); 
  //we prefere stable sorting:
  stable_sort(perm.begin(),perm.end(),precedes_in(*this));
  permutation p2=perm.inverse();
  Interaction_Matrix m2=this->permute(p2);

  // generate trophical species:
  equal_in eq=equal_in(m2);
  simple_vector<bool> selected(the_size,false);
  int j=0;
  for(int i=0;i<the_size;i++){
    if(i>j && eq(i,j)){
#ifndef NO_POPULATION_DYNAMICS
      m2.the_species[j].insert(m2.the_species[j].end(),
			       m2.the_species[i].begin(),
			       m2.the_species[i].end() );
#endif
    }else{
      selected[i]=true;
      j=i;
    }
  }
  
  return m2.select(selected);
}

void Interaction_Matrix::prey_hist(std::ostream & co, 
				   string prefix /*=""*/){
  histogram hist=prey_hist();
  distribution p=theoretical_prey_dist(hist.size()+1);
  double S=Number_of_Species_S();
  double L=Number_of_Links_L();
  double z=L/Number_of_Species_S();

  histogram ii;
  for(int i=p.size()-1;i>=0;i--)
    ii[i]=i;

  hist.resize(p.size());

  if(prefix.size()) 
    co << "prey hist:" << endl;
  sequence<int> h=0.5+S*p;
  co << prefix+" "+
    ii.format("%4g ")+
    hist.format("%5g ")+
    h.format("%5i ")+
    (distribution(hist)/S).format("%#8.3g ")+
    p.format("%#8.3g ")+
    (distribution(ii)/(2*z)).format("%#8.3g ");

  if(prefix.size()){
    co << prefix;
    co << "basal species fraction: " << hist[0]/double(the_size) 
       << " (" << log(1+2*z)/(2*z) << ")" << endl;
  }
}

void Interaction_Matrix::cumulative_prey_hist(std::ostream & co, 
					      string prefix /*=""*/){
  histogram hist=cumulative_prey_hist();
  distribution p=theoretical_cumulative_prey_dist(hist.size()+1);
  double S=Number_of_Species_S();
  double L=Number_of_Links_L();
  double z=L/S;

  histogram ii;
  for(int i=p.size()-1;i>=0;i--)
    ii[i]=i;

  hist.resize(p.size());

  if(prefix.size()) 
    co << "cumulative prey hist:" << endl;
  sequence<int> h=0.5+S*p;
  co << prefix+" "+
    ii.format("%4g ")+
    hist.format("%5g ")+
    h.format("%5i ")+
    (distribution(hist)/S).format("%#8.3g ")+
    p.format("%#8.3g ")+
    (distribution(ii)/(2*z)).format("%#8.3g ");

  if(prefix.size()){
    co << prefix;
    co << "basal species fraction: " << hist[0]/double(the_size) 
       << " (" << log(1+2*z)/(2*z) << ")" << endl;
  }
}


#if 0// dead code (def ON_SX5FSV)
static double gsl_sf_gamma_inc_P(double a, double x){
  gsl_sf_result res;
  gsl_sf_gamma_inc_P_e(a,x,&res);
  return res.val;
}
static double gsl_sf_expint_Ei(double x){
  gsl_sf_result res;
  gsl_sf_expint_Ei_e(x,&res);
  return res.val;
}
#endif


double P_pred(int m,double z){
#ifdef WITH_GSL
  double sum=0;
  double sum_old=-1;
  for(int mm=m;sum!=sum_old;mm++){
    sum_old=sum;
    sum+=gsl_sf_gamma_inc_P(mm+1,2*z);
  }
  return sum/(2*z);
#else
  return 0;
#endif
}
      


void Interaction_Matrix::predator_hist(std::ostream & co, 
				       string prefix /*=""*/){
  histogram hist=predator_hist();
  distribution p=theoretical_predator_dist(hist.size()+1);
  double S=Number_of_Species_S();
  double L=Number_of_Links_L();
  double z=L/Number_of_Species_S();

  histogram ii;
  for(int i=p.size()-1;i>=0;i--)
    ii[i]=i;

  hist.resize(p.size());

  if(prefix.size()) 
    co << "predator hist:" << endl;
  sequence<int> h=0.5+S*p;
  co << prefix+" "+
    ii.format("%4g ")+
    hist.format("%5g ")+
    h.format("%5i ")+
    (distribution(hist)/S).format("%#8.3g ")+
    p.format("%#8.3g ")+
    (distribution(ii)/(2*z)).format("%#8.3g ");

  if(prefix.size()){
    co << prefix;
    co << "top predator fraction: " << hist[0]/double(the_size) 
       << " (" << (1-exp(-2*z))/(2*z) << ")" << endl;
  }
}
void Interaction_Matrix::cumulative_predator_hist(std::ostream & co, 
						  string prefix /*=""*/){
  histogram hist=cumulative_predator_hist();
  distribution p=theoretical_cumulative_predator_dist(hist.size()+1);
  double S=Number_of_Species_S();
  double L=Number_of_Links_L();
  double z=L/Number_of_Species_S();

  histogram ii;
  for(int i=p.size()-1;i>=0;i--)
    ii[i]=i;

  hist.resize(p.size());

  if(prefix.size()) 
    co << "cumulative predator hist:" << endl;
  sequence<int> h=0.5+S*p;
  co << prefix+" "+
    ii.format("%4g ")+
    hist.format("%5g ")+
    h.format("%5i ")+
    (distribution(hist)/S).format("%#8.3g ")+
    p.format("%#8.3g ")+
    (distribution(ii)/(2*z)).format("%#8.3g ");

  if(prefix.size()){
    co << prefix;
    co << "top predator fraction: " << hist[0]/double(the_size) 
       << " (" << (1-exp(-2*z))/(2*z) << ")" << endl;
  }
}

Interaction_Matrix::histogram
Interaction_Matrix::prey_hist(){
  histogram hist;
  for(int i=0;i<the_size;i++){
    int count=0;
    for(int j=0;j<the_size;j++)
      if( (int)(*this)[i][j] == (int)NetworkAnalysis::eats ) count++;
    hist[count]++;
  }
  return hist;
}

Interaction_Matrix::histogram
Interaction_Matrix::predator_hist(){
  histogram hist;
  for(int i=0;i<the_size;i++){
    int count=0;
    for(int j=0;j<the_size;j++)
      if( (int)(*this)[j][i] == (int)NetworkAnalysis::eats ) count++;
    hist[count]++;
  }
  return hist;
}

Interaction_Matrix::histogram
Interaction_Matrix::cumulative_prey_hist(){
  return prey_hist().reverse_cumulative_sum();
}

Interaction_Matrix::histogram
Interaction_Matrix::cumulative_predator_hist(){
  return predator_hist().reverse_cumulative_sum();
}

Interaction_Matrix::distribution
Interaction_Matrix::theoretical_prey_dist(int n, double z){
#ifdef WITH_GSL
  distribution p;
  if(z<0) z=links_per_species_Z();
  if(z==0){
    WARNING("no links");
    return p;
  }
  for(int i=n-1;i>=0;i--){
    double k2z=double(i)/(2*z);
    double theo=(i?-1/(2*z)*gsl_sf_expint_Ei(-k2z): log(1+2*z)/(2*z) );
    p[i]=theo;
  }
  return p;
#else
  WARNING("theoretical_prey_dist needs gsl library");
  return 0;
#endif
}

Interaction_Matrix::distribution
Interaction_Matrix::theoretical_cumulative_prey_dist(int n, double z){
#ifdef WITH_GSL
  distribution p;
  if(z<0) z=links_per_species_Z();
  if(z==0){
    WARNING("no links");
    return p;
  }
  for(int i=n-1;i>=0;i--){
    double k2z=double(i)/(2*z);
    double theo=exp(-k2z)+(i?k2z*gsl_sf_expint_Ei(-k2z):0);
    p[i]=theo;
  }
  return p;
#else
  WARNING("theoretical_cumulative_prey_dist needs gsl library");
  return 0;
#endif
}

Interaction_Matrix::distribution
Interaction_Matrix::theoretical_predator_dist(int n, double z){
#ifdef WITH_GSL
  distribution p;
  if(z<0) z=links_per_species_Z();
  if(z==0){
    WARNING("no links");
    return p;
  }
  for(int m=n-1;m>=0;m--){
    double theo= gsl_sf_gamma_inc_P(m+1,2*z)/(2*z);
    p[m]=theo;
  }
  return p;
#else
  WARNING("theoretical_predator_dist needs gsl library");
  return 0;
#endif
}

Interaction_Matrix::distribution
Interaction_Matrix::theoretical_cumulative_predator_dist(int n, double z){
#ifdef WITH_GSL
  distribution p;
  if(z<0) z=links_per_species_Z();
  if(z==0){
    WARNING("no links");
    return p;
  }
  for(int m=n-1;m>=0;m--){
    double theo= P_pred(m,z);
    p[m]=theo;
  }
  return p;
#else
  WARNING("theoretical_cumulative_predator_dist needs gsl library");
  return 0;
#endif
}

double Interaction_Matrix::prop_T()//fraction of top predators
{
  histogram h=predator_hist();
  return prop_T(h);
}

double Interaction_Matrix::prop_I()//fraction of intermediate species
{
  histogram h1=predator_hist();
  histogram h2=prey_hist();
  return prop_I(h1,h2);
}

double Interaction_Matrix::prop_B()//fraction of bottom species
{
  histogram h=prey_hist();
  return prop_B(h);
}

double Interaction_Matrix::prop_GenSD()//std of prey count normalized to (L/S)
{
  histogram h=prey_hist();
  return prop_GenSD(h);
}

double Interaction_Matrix::prop_VulSD()//std of predator count
				       //normalized to (L/S)
{
  histogram h=predator_hist();
  return prop_VulSD(h);
}

double Interaction_Matrix::prop_T(histogram & h)//fraction of top predators
{
  return double(h[0])/size();
}

double Interaction_Matrix::prop_I(histogram & h1, histogram & h2)
  //fraction of intermediate species
{
  return 1-double(h1[0]+h2[0])/size();
}

double Interaction_Matrix::prop_B(histogram & h)//fraction of bottom species
{
  return double(h[0])/size();
}

double Interaction_Matrix::prop_GenSD(histogram & h)//std of prey count normalized to (L/S)
{
  weighted_average_meter av;
  for(unsigned int i=0;i<h.size();i++){
    av.sample(i,h[i]);
  }
  double S=Number_of_Species_S();
  // correct for an error made by Williams and Martinez
  return av.sample_std()/av.readout()*sqrt(S/(S-1));
}

double Interaction_Matrix::prop_VulSD(histogram & h)//std of predator count
  //normalized to (L/S)
{
  return prop_GenSD(h); // the method to calculate it is the same;
}


Interaction_Matrix::similarity_matrix_t Interaction_Matrix::similarity_s()
//fraction of prey and predators shared by a pair
{
  similarity_matrix_t s(size(),size());
  for(int i=0;i<size();i++){
    for(int j=0;j<size();j++){
      if(j<=i){ // we include i==j for testing
	// compute similarity:
// 	int linked_to_i=0;
// 	int linked_to_j=0;
// 	int linked_to_both=0;
// 	for(int k=0;k<size();k++){
// 	  if(true || (k!=i && k!=j)){
// 	    if(eats(i,k)){
// 	      linked_to_i++;
// 	      if(eats(j,k)){
// 		linked_to_j++;
// 		linked_to_both++;
// 	      }
// 	    }else if(eats(j,k)){
// 	      linked_to_j++;
// 	    }
// 	    if(eats(k,i)){
// 	      linked_to_i++;
// 	      if(eats(k,j)){
// 		linked_to_j++;
// 		linked_to_both++;
// 	      }
// 	    }else if(eats(k,j)){
// 	      linked_to_j++;
// 	    }
// 	  }
// 	}
// 	s[i][j]=
// 	  double(2*linked_to_both)/(linked_to_j+linked_to_i);

// 	int linked_to_i=0;
// 	int linked_to_j=0;
// 	int linked_to_both=0;
// 	for(int k=0;k<size();k++){
// 	  if(true || (k!=i && k!=j)){
// 	    if(eats(i,k) || eats(k,i)){
// 	      linked_to_i++;
// 	      if(eats(j,k) || eats(k,j)){
// 		linked_to_j++;
// 		linked_to_both++;
// 	      }
// 	    }else if(eats(j,k) || eats(k,j)){
// 	      linked_to_j++;
// 	    }
// 	  }
// 	}
// 	s[i][j]=
// 	  double(linked_to_both)/(linked_to_j+linked_to_i-linked_to_both);


	// this makes more sense, but is not what W & M used:
	average_meter ss;
	for(int k=0;k<size();k++){
	  if(true || (k!=i && k!= j)){
	    if(eats(i,k) || eats(j,k))
	      ss.sample((eats(i,k) && eats(j,k))?1:0);
	    if(eats(k,i) || eats(k,j))
	      ss.sample((eats(k,i) && eats(k,j))?1:0);
	  }
	}
	s(i,j)=ss.readout();
      }
    }
  }
  WARNING("Check code if symmetry is used correctly here!");
  s = symmatu(s);
  //cout << s << endl;
  return s;
}
double Interaction_Matrix::prop_MxSim()//web average over max s
{
  average_meter mxsim;
  similarity_matrix_t s=similarity_s();
  for(int i=0;i<size();i++){
    double mx=0;
    for(int j=0;j<size();j++)
      if(i!=j && s(i,j)>mx) mx=s(i,j);
    mxsim.sample(mx);
  }
  return mxsim.readout();
}
    
typedef Interaction_Matrix::histogram histogram;

histogram Interaction_Matrix::chain_hist(){
  //all loop-free directed chains starting from bottom
  return histogram(); //dummy
//   Interaction_Matrix im=(*this).tsort();
//   for(int i=0;i<im.size();i++){
//     for(int j=0;j<i+1;j++){
//       im[i][j]=NetworkAnalysis::none;
//     }
//   }
//  return loop_free_chain_hist(*this);
}

double Interaction_Matrix::prop_ChnLg(histogram & h){ //average chain length
  weighted_average_meter av;
  for(unsigned int i=0;i<h.size();i++){
    av.sample(i,h[i]);
  }
  return av.readout();
}

double Interaction_Matrix::prop_ChnSD(histogram & h){ //std chain length
  weighted_average_meter av;
  for(unsigned int i=0;i<h.size();i++){
    av.sample(i,h[i]);
  }
  return av.sample_std();
}

double Interaction_Matrix::prop_ChnNo(histogram & h){ //log # of chains
  return log10(accumulate(h.begin(),h.end(),0.0));
}

double Interaction_Matrix::prop_Cannib(){ //cannibal fraction
  int s=size();
  int count=0;
  for(int i=0;i<s;i++){
    if(eats(i,i)) count++;
  }
  return double(count)/s;
}

double Interaction_Matrix::prop_Loop(){ //fraction of species involved in loops
  return loop_species_fraction(*this);
}

double Interaction_Matrix::prop_Omniv(){ 
  // fraction of species that have food chains of different length
  return omnivore_fraction(*this);
}

double Interaction_Matrix::prop_Clust(){ // Clustering Coefficient
  const int s=size();
  average_meter C;
  for(int i=0;i<s;i++){
    average_meter localC;
    for(int j=s-1;j>=0;j--){
      if(connected(i,j) && j!=i){
	for(int k=j-1;k>=0;k--){
	  if(connected(i,k) && k!=i)
	    localC.sample(connected(j,k));
	}
      }
    }
    if(localC.n()>=1){
      C.sample(localC.readout());
    }
  }
  return C.readout();
}
  

//// helper class for prop_Ddiet
class diet_sorter_t {
  enum {before,on_line,after};
  typedef int position_to_species_line;
  static simple_vector<bool> _can_be_linearly_ordered;
  void assign_orderable(int order_class);
  void advance_positions_starting_with(int n,
				       simple_vector<position_to_species_line> & 
				       position,
				       int sharing_class,
				       int order_class);
public:
  static const int i_bit;
  static const int j_bit;
  static const int k_bit;
  static const int n_order_classes;
  void initialize__can_be_linearly_ordered();
  diet_sorter_t(){
    if(int(_can_be_linearly_ordered.size())!=n_order_classes){
      initialize__can_be_linearly_ordered();
    }
  }
  bool can_be_linearly_ordered(int order_class){
    return _can_be_linearly_ordered[order_class];
  }
};
const int diet_sorter_t::i_bit=1<<0;
const int diet_sorter_t::j_bit=1<<1;
const int diet_sorter_t::k_bit=1<<2;
const int diet_sorter_t::n_order_classes=(1<<((i_bit|j_bit|k_bit)+1))/*=256*/;
simple_vector<bool> diet_sorter_t::_can_be_linearly_ordered=simple_vector<bool>(0);

void diet_sorter_t::assign_orderable(int order_class){
  _can_be_linearly_ordered[order_class]=true;
}
void diet_sorter_t::initialize__can_be_linearly_ordered(){
  _can_be_linearly_ordered=simple_vector<bool>(n_order_classes,false);
  simple_vector<position_to_species_line> position(3,before);
  int sharing_class=0;
  int order_class=0;
  assign_orderable(order_class);
  // the rest is recursion:
  for(int m=0;m<3;m++){
    advance_positions_starting_with(m,position,sharing_class,order_class);
  }
  order_class|=(1<<sharing_class);
  assign_orderable(order_class);
  for(int m=0;m<3;m++){
    advance_positions_starting_with(m,position,sharing_class,order_class);
  }
}
  
void diet_sorter_t::advance_positions_starting_with
(int n,
 simple_vector<position_to_species_line> & position,
 int sharing_class,
 int order_class){
  if(position[n]==after) return; //stop recursion
  
  position[n]=position[n]+1;
  if(position[n]==after){
    sharing_class&=(i_bit|j_bit|k_bit)^(1<<n);
  }else{ //position[n]==on_line
    sharing_class|=(1<<n);
  }
  for(int m=0;m<3;m++){
    advance_positions_starting_with(m,position,sharing_class,order_class);
  }
  order_class|=(1<<sharing_class);
  assign_orderable(order_class);
  for(int m=0;m<3;m++){
    advance_positions_starting_with(m,position,sharing_class,order_class);
  }
  --position[n];
}

double Interaction_Matrix::prop_Ddiet(bool FrenchVariant){ 
  // the fraction of species triples with an irreducible dietary gap
  static diet_sorter_t sorter;
  int orderable_count=0,total_count=0;
  simple_vector<bool> possible(the_size,true);
  if(FrenchVariant){
    for(int i=0;i<the_size;i++){
      int link_count=0;
      for(int n=0;n<the_size && link_count <=2;n++){
	if(eats(i,n)) link_count++;
      }
      possible[i]=(link_count>=2);
    }
  }
  for(int i=the_size-1;i>=0;i--){
    if(possible[i]){
      for(int j=i-1;j>=0;j--){
	if(possible[j]){
	  for(int k=j-1;k>=0;k--){
	    if(possible[k]){
	      int order_class=0;
	      for(int n=the_size-1;n>=0;n--){
		int sharing_class=0;
		if(eats(i,n)) sharing_class|=diet_sorter_t::i_bit;
		if(eats(j,n)) sharing_class|=diet_sorter_t::j_bit;
		if(eats(k,n)) sharing_class|=diet_sorter_t::k_bit;
		order_class|=(1<<sharing_class);
	      }
	      if(sorter.can_be_linearly_ordered(order_class)){
		orderable_count++;
	      }else{
		//#define DDIET_CERTIFY
#ifdef DDIET_CERTIFY
// 		cout << the_size-i-1 << ", "<< 
// 		  the_size-j-1 << ", "<< 
// 		  the_size-k-1 << ", " ;
		cout << i << ", "<< 
		  j << ", "<< 
		  k << ", " ;
		int oc=order_class;
		for(int p=0;p<8;p++){
		  cout << ((oc<<=1)&256)/256 ;
		}	
		cout << endl;
#endif
	      }
	      
	      total_count++;
	    }
	  }
	}
      }
    }
  }
//   cout << total_count << " triples" << endl;
  return 1-orderable_count/double(total_count);
}



Interaction_Matrix::prop_vec_t Interaction_Matrix::props(){ // all properties.
  graphs_stats_data stats=graph_analyze(*this);
  Interaction_Matrix::prop_vec_t p;
  histogram h1=predator_hist();
  histogram h2=prey_hist();
  //histogram ch=chain_hist();

  p[pS]=size();
  p[pZ]=links_per_species_Z();
  p[pC]=connectance_C();
  p[pT]=prop_T(h1);
  p[pI]=prop_I(h1,h2);
  p[pB]=prop_B(h2);
  p[pGenSD]=prop_GenSD(h2);
  p[pVulSD]=prop_VulSD(h1);
  p[pMxSim]=prop_MxSim();
  p[pCannib]=prop_Cannib();
//   p[pChnLg]=prop_ChnLg(ch);
//   p[pChnSD]=prop_ChnSD(ch);
//   p[pChnNo]=prop_ChnNo(ch);
  p[pLoop]=prop_Loop();
  p[pOmniv]=prop_Omniv();
  p[poChnLg]=stats.oChnLg;
  p[poChnSD]=stats.oChnSD;
  p[poChnNo]=stats.oChnNo;
  p[poLoop]=stats.oLoop;
  p[poOmniv]=stats.oOmniv;
  p[pDdiet]=prop_Ddiet();
  p[pfDdiet]=prop_Ddiet(true);//French Variant
  p[pCy4]=prop_Cy4();
  p[pNest]=prop_Nest();
  p[pClust]=prop_Clust();
  p[pSStab]=structural_stability();
  return p;
}

sequence<const char *> Interaction_Matrix::prop_names(){ 
  // the names of these properties. 
  sequence<const char *> p;
  p[pS]="S";
  p[pZ]="Z";
  p[pC]="C";
  p[pT]="T";
  p[pI]="I";
  p[pB]="B";
  p[pGenSD]="GenSD";
  p[pVulSD]="VulSD";
  p[pMxSim]="MxSim";
  p[pChnLg]="ChnLg";
  p[pChnSD]="ChnSD";
  p[pChnNo]="ChnNo";
  p[pCannib]="Cannib";
  p[pLoop]="Loop";
  p[pOmniv]="Omniv";
  p[poChnLg]="aChnLg";
  p[poChnSD]="aChnSD";
  p[poChnNo]="aChnNo";
  p[poLoop]="aLoop";
  p[poOmniv]="aOmniv";
  p[pDdiet]="vDdiet";
  p[pfDdiet]="Ddiet";
  p[pCy4]="Cy4";
  p[pNest]="Nest";
  p[pClust]="Clust";
  p[pSStab]="SStab";
  return p;
} 

bool Interaction_Matrix::connected(){
  return is_connected(*this);
}

Interaction_Matrix Interaction_Matrix::largest_connected_subweb() const{
  return the_largest_connected_subweb(*this);
}

void Interaction_Matrix::two_column_write(const string & filename,
					  string delimiter/*=","*/)
{
  std::ofstream os(filename.c_str());
  for(int i=0;i<the_size;i++){
    for(int j=0;j<the_size;j++){
      if(eats(i,j))
	os << the_size-j << delimiter << the_size-i << endl;
    }
    
  }
}

Interaction_Matrix::distribution Interaction_Matrix::trophic_height_vector(){
  return ::trophic_height_vector(*this);
}

Interaction_Matrix::histogram 
Interaction_Matrix::shortest_path_level_vector(){
  int species_left=the_size;
  simple_vector<bool> lowest_level(the_size);
  histogram level(the_size);
  for(int i=the_size;i-->0;){
    bool lowest=true;
    for(int j=the_size;j-->0;){
      if(eats(i,j)){
	lowest=false;
	break;
      }
    }
    if(lowest){
      level[i]=1;
      species_left--;
    }
  }
  if(species_left==the_size) 
    FATAL_ERROR("no lowest level found!");
  int level_found=1;
  while(species_left){
    int last_species_left=species_left;
    for(int i=the_size;i-->0;){
      if(level[i]==level_found){
	for(int j=the_size;j-->0;){
	  if(!level[j] && eats(j,i)){
	    level[j]=level_found+1;
	    species_left--;
	  }
	}
      }
    }
    level_found++;
    if(last_species_left==species_left){
      WARNING("problem computing shortest_path_level_vector");
      break;
    }
  }
  return level;
}

double Interaction_Matrix::prop_Nest(){
  int violations=0;
  for(int i=the_size;i-->0;){
    for(int j=i-1;j-->0;){
      bool i_in_j=true,j_in_i=true,i_overlap_j=false;
      for(int k=the_size;k-->0;){
	if(eats(i,k) && !eats(j,k)) i_in_j=false;
	if(!eats(i,k) && eats(j,k)) j_in_i=false;
	if(eats(i,k) && eats(j,k)) i_overlap_j=true;
      }
      if(!(i_in_j || j_in_i || !i_overlap_j)){
	violations++;
      }
    }
  }
  return violations;
}
	

void Interaction_Matrix::dot_graph(const string & filename){
  std::ofstream os(filename.c_str());
  os << "strict digraph web {" << endl;
  os << "shape=circle;" << endl;
  os << "splines=false;" << endl;
  os << "size=\"7,11\";" << endl;
  os << "rankdir=BT;" << endl;
  os << "style=filled;" << endl;
  os << "height=1;" << endl;
  os << "width=1;" << endl;
  os << "fixedsize=true;" << endl;
  for(int i=0;i<the_size;i++){
    for(int j=0;j<the_size;j++){
      if(eats(i,j))
	os << "s" << j << " -> " << "s" << i << "[label = \"\"];" << endl;
    }
  }
  os << "}" << endl;
}

double Interaction_Matrix::structural_stability(){
  // "structural stability" according to J. Fox Oikos 115:97-109, 2006 
  NewMatrix plusminus(the_size,the_size);
  for(int i=the_size;i-->0;){
    for(int j=the_size;j-->0;){
      if(eats(i,j)){ 
	plusminus(i,j)=-1;
	if(eats(j,i)){
	  plusminus(j,i)=-1;
	}else{
	  plusminus(j,i)=1;
	}
      }
    }
  }
  
  int cannibs=0;
  for(int i=the_size;i-->0;){
    cannibs-=(int)plusminus(i,i);
  }

  NewMatrix save_plusminus=plusminus;

  average_meter stability;
  bool random_cannibs=true;//(cannibs<10);   // add cannibals
  
  for(int rep=(random_cannibs?100:1);rep-->0;){
    if(random_cannibs){
      int non_cannibs=the_size-cannibs;
      plusminus=save_plusminus;
      // add non_cannibs/2 cannibals
      for(int l=non_cannibs/2;l-->0;){
	//REPORT(non_cannibs);
	int nth=random_integer(non_cannibs); // = 0..non_cannibs-1
	//REPORT(nth);
	int n=-1,i=-1;
	while(n<nth){
	  i++;
	  while(plusminus(i,i)==-1) i++;
	  n++;
	}
	ALWAYS_ASSERT(i<the_size);
	//REPORT(i);
	plusminus(i,i)=-1;
	non_cannibs--;
      }
    }
    
//     PPrint();
//     for(int j=0;j<the_size;j++){
//       for(int i=0;i<the_size;i++){
// 	if(plusminus[i][j]==+1)
// 	  cout << '+';
// 	else if(plusminus[i][j]==-1)
// 	  cout << '-';
// 	else if(plusminus[i][j])
// 	  FATAL_ERROR("plusminus garbeled");
// 	else 
// 	  cout << '0';
// 	cout << ' ';
//       }
//       cout << endl;
//     }
//     cout << endl;

    NewVector evalr,evali;
    NewEigen(plusminus,evalr,evali);
    
    // find maximum real part:
    int maxind=0;
    for(int i=the_size;i-->0;){
      if(evalr[i]>evalr[maxind]){
      maxind=i;
      }
    }
    stability.sample(evalr[maxind]);
    if(rep==0){
      //REPORT(evalr);
      //REPORT(evali);
    }
  }
  //REPORT(stability.std());
  //REPORT(stability.error());
  return stability.readout();
}

// #include "consecutive_one.h"

// bool has_consecutive_ones(const CMatrix<NetworkAnalysis::Interaction> & m){
//   static consecutive_ones_tester cons_one;
//   return cons_one.test(m);
// }

#include "pqtree.h"

bool has_consecutive_ones(const CMatrix<NetworkAnalysis::Interaction> & m){
  set<int> U;
  for(int c=0;c<m.GetXSize();c++)
    U.insert(c);
  
  pqtree<int> T(U);

  try{
    for(int r=0;r<m.GetYSize();r++){
      // for each row ...
      set<int> s;
      for(int c=0;c<m.GetXSize();c++){
	// get the ones in each column:
	if(m[r][c]==NetworkAnalysis::eats){
	  s.insert(c);
	}
      }
      
      if(!T.reduce(s)){
	return false;
      }
    }
  }catch(AnalysisBug){
    Interaction_Matrix(m).PPrint();
    return false;
  }
  return true;
}


bool Interaction_Matrix::has_consecutive_ones(){
  return ::has_consecutive_ones(*this);
}

bool Interaction_Matrix::is_interval(){
  return niche_overlap_graph_analysis(*this) >= 2;
}

bool Interaction_Matrix::is_chordal(){
  return niche_overlap_graph_analysis(*this) >= 1;
}
