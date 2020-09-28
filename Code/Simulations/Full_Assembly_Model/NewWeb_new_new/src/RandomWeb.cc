#include "snapshot.h"
#include "random.h"
#include "matrix_transformers.h"

#ifdef GET_INDIRECT_DEPENDENCIES
#include "evaluate.h"
#include "Statistics.h"
#include "linpack_eigen.h"
#include "cfgList.h"
#include "xy_graph.h"
#include "polyfit.h"
#endif

const char * strength_spectrum_file="strength.dat";

int main(int argc,char *argv[]){

  int S=30;
  double alpha=1;
  double plant_fraction=2.0/3.0;

  if(argc>1){
    S=atoi(argv[1]);
  }
  if(argc>2){
    alpha=atof(argv[2]);
  }
  if(argc>3){
    plant_fraction=atof(argv[3]);
  }
  
  link_strength_matrix l;
  l.resize(S);

  int n_animals=int(S*(1-plant_fraction)+0.5);

  for(int i=n_animals;i-->0;){
    for(int j=S;j-->n_animals;){
      l(j,i)=exp(gaussian(0,alpha));
    }
    for(int j=i;j-->0;){
      l(j,i)=exp(gaussian(0,alpha));
    }
  }

  Snapshot s;
  s.set_number_of_compartments(S);
  s.set_flows(l);
  s.trophic_level_structure("levels.dat");
  strength_distribution(s.get_ifrac(),1e-3,strength_spectrum_file);
}
  

