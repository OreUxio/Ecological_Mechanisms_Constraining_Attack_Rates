#include "error.h"
#include "NetworkAnalysis.h"
#include <iostream>
#include <fstream>
#include <string>

#ifdef GET_INDIRECT_DEPENDENCIES
#include "Statistics.h"
#include "linpack_eigen.h"
#include "random.h"
#endif

// Manage adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
{
  {0, CFG_END, 0}
};
static cfg_add dummy(cfg);


using namespace std;

Interaction_Matrix read_topological_web(istream & is){
  int S;
  is >> S;
  if(is.eof() || S<1) throw 1;
  Interaction_Matrix m(S);
  //REPORT(S);
  string line;
  for(int i=0;i<S;i++){
    is >> line;
    //REPORT(line);
    for(int j=0;j<S;j++){
      if(line[j]=='X'){
	m[j][i]=NetworkAnalysis::eats;
      }
    }
  }
  return m;
}


int main(int argc,char *argv[]){
  if(argc<2){
    FATAL_ERROR("no filename given");
  }
  ifstream is(argv[1]);

  int n=0;
  sequence<double> ss;

  try{
    while(!is.eof()){
      Interaction_Matrix m=read_topological_web(is);
      if(is.eof()) throw 1;
      ss[n++]=m.structural_stability();
      cout << ".";
      cout.flush();
    }
  }catch(int){
  }
  cout << endl;

  if(n==0) FATAL_ERROR("no data");

  sort(&ss[0],(&ss[0])+n);

  double alpha=0.05;
  double lower_n=(alpha/2)*n;
  double upper_n=(1-alpha/2)*n;


  double upper_limit,lower_limit;
  if(int(upper_n)==n-1){
    upper_limit=ss[n-1];
  }else{
    double resi=upper_n-int(upper_n);
    upper_limit=(1-resi)*ss[int(upper_n)]+resi*ss[int(upper_n)+1];
  }
  if(n==1){
    lower_limit=ss[0];
  }else{
    double resi=lower_n-int(lower_n);
    lower_limit=(1-resi)*ss[int(lower_n)]+resi*ss[int(lower_n)+1];
  }

  cout << lower_limit << "\t" << upper_limit << endl;

  
  if(argc>=3 && n>10){
    ofstream os(argv[2]);
    for(int i=0;i<n;i++){
      os << ss[i] << " " << i*(1/double(n)) << endl;
    }
  }

//   sequence<string> names=Interaction_Matrix::prop_names();
//   Interaction_Matrix::distribution pr=m.props();
  
//   for(unsigned int j=0;j<pr.size();j++){
//     std::cout << names[j] << " " 
// 	      << pr[j]
// 	      << endl;
//   }
}
