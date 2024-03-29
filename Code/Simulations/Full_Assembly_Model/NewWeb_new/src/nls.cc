// $Id: nls.cc 307 2006-03-07 17:00:06Z cvsrep $


#include <iomanip>
#include "nls_web.h"
#include "error.h"
#include "random.h"
#include "Statistics.h"
#include <fstream>

#ifdef GET_INDIRECT_DEPENDENCIES
#include "linpack_eigen.h"
#endif

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);

using namespace std;

double LOG(double x){
  return log(x);
}

typedef sequence<double> vec;

/////////////////////////////////////////////
// main
int main(int argc,char *argv[]){
  set_random_seed(0);
  const int n_webs=argc-2;
  if(n_webs<1) FATAL_ERROR("please give at least one input web name");
  nls_web w[n_webs];
  for(int i=0;i<n_webs;i++){
    w[i]=nls_web(argv[i+1]);
  }
  nls_web web=w[0];
  for(int i=1;i<n_webs;i++){
    web+=w[i];
  }
  web.strength_distribution(0, argv[argc-1]);
  return 0;
}
