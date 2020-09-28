// $Id: IntervalTest.cc 1923 2010-11-04 16:36:21Z axel $
#include <string>
const std::string version("$Revision: 1923 $");
const std::string source("$Source: /home/cvsrep/CVS/NewWeb/src/IntervalTest.cc,v $ $Date: 2010-11-04 16:36:21 +0000 (Thu, 04 Nov 2010) $");

#include <signal.h>
#include "sequence.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "NetworkAnalysis.h"
#include "random.h"
#include "error.h"
#include "Statistics.h"
#include "standardize.h"
#include "topology_generator.h"
#include "random_pick_field.h"
#include "Moran.h"
#include "MoranI.h"
#include "MoranIpp.h"
#include "OUWeb.h"

// these are included to force linking
#include "linpack_eigen.h"
#include "binomial_dist.h"
#include "niche_width_finder.h"
#include "ODE.h"

//#define NOLINKS
//#define DELETE_NEW

#ifdef DEBUGGING
#undef DEBUGGING
#endif

#ifdef DEBUGGING
static int counter=0;
#endif

using namespace std;

const char * cfg_file_name="MoranF.cfg";

const double ln10=M_LN10; //==log(10.0);

const char * in_file_name="/home/axel/Webs/St. Martin.web";

static int ndim=3;
static double C0=0.1;
static double lambda=0;
static double p_f=-1; // used only if set to value >=0
static double p_v=-1; // used only if set to value >=0
static double std_mutation_factor=0.001;
static double q=1;
static double saturation=0;
static double rewire_fraction=0;
static Topology_Generator * web_generator=0;
static int cumulative_stats = 0;
static int number_of_repetitions=1000;
static int print_all_webs=0;
static int random_seed=1;
static int test_mode=0;
static int n_images=0;  // how many finished food-webs images to save
static int exact=0;
static int correct_correlations=0;
static int raw_target_S=0;
static const char* rec_file=0;
static const char* read_file=0;
static int adapt_foragers=0;

// adjustable parameters:
#include "cfgList.h"
static cfgStruct cfg[] = 
  { 
    CFGSTRING(in_file_name),
    CFGDOUBLE(C0),
    CFGDOUBLE(lambda),
    CFGDOUBLE(p_f),
    CFGDOUBLE(p_v),
    CFGDOUBLE(std_mutation_factor),
    CFGDOUBLE(q),
    CFGDOUBLE(saturation),
    CFGDOUBLE(rewire_fraction),
    CFGINT(ndim),
    CFGINT(cumulative_stats),
    CFGINT(number_of_repetitions),
    CFGINT(print_all_webs),
    CFGINT(random_seed),
    CFGINT(test_mode),
    CFGINT(n_images),
    CFGINT(exact),
    CFGINT(correct_correlations),
    CFGINT(raw_target_S),
    CFGINT(adapt_foragers),
    {0, CFG_END, 0}   /* no more parameters */
  };
static cfg_add dummy(cfg);


static Interaction_Matrix im,exact_im;
typedef sequence<int> histogram;
typedef sequence<double> distribution;


void read_arguments(int argc,char *argv[]){
  int c;
  
  while(1){
    c=getopt(argc, argv, "hvr:w:t");
    if (c == -1)
      break;
    
    switch(c){
    case 'w':
      rec_file=optarg;
      cout << "recording" << endl;
      break;
    case 'r':
      read_file=optarg;
      break;
    case 'v':
      printf("%s\n",version.c_str());
      exit(0);
      break;
    case 't':
      write_cfg_template("template.cfg");
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

/////////////////////////////////////////////
// rewire

// rewires a fraction frac of all trophic links in m, picked randomly
void rewire(double frac,Interaction_Matrix & m){
  ALWAYS_ASSERT(0 <= frac && frac <=1);
  typedef std::pair<int,int> ipair;
  const int S=m.size();
  int L=0;
  random_pick_field< ipair > linked;
  random_pick_field< ipair > not_linked;
  for(int i=S;i-->0;){
    for(int j=S;j-->0;){
      if(m[i][j]==NetworkAnalysis::eats){
	linked.insert(ipair(i,j));
	L++;
      }else{
	not_linked.insert(ipair(i,j));
      }
    }
  }
  const int to_rewire=int(L*frac+0.5);
  ALWAYS_ASSERT(to_rewire <= S*S-L);
  for(int i=to_rewire;i-->0;){
    random_pick_field< ipair >::position p,q;
    p=linked.random_pick(randint());
    q=not_linked.random_pick(randint());
    m[linked[p].first][linked[p].second]=NetworkAnalysis::none;
    m[not_linked[q].first][not_linked[q].second]=NetworkAnalysis::eats;
    linked.erase(p);
    not_linked.erase(q);
  }
}

histogram index(int n){
  histogram h;
  for(int i=n-1;i>=0;i--)
    h[i]=i;
  return h;
}


/////////////////////////////////////////////
// main

int main(int argc,char *argv[]){
  read_arguments(argc,argv);
  cout << source << endl;
  cout << version << endl;
  system("/bin/date");
  signal_handling();

  REPORT(cfg_file_name);

  if (cfgParse(cfg_file_name, full_cfg_list(), CFG_SIMPLE) == -1)
    FATAL_ERROR("error reading parameter file");

  REPORT(random_seed);
  set_random_seed(random_seed);

  ofstream rec;
  ifstream is;
  if(rec_file){
    cout << "opening " << rec_file << " for writing" << endl;
    rec.open(rec_file);
  }
  if(read_file){
    is.open(read_file);
  }else{
#if 1
//     web_generator = new Moran(std_mutation_factor*std_mutation_factor,
// 			      C0,p_f,p_v,lambda,q,correct_correlations);
    web_generator = new MoranI(ndim,std_mutation_factor*std_mutation_factor,
			       C0,lambda,saturation,q,correct_correlations,
			       adapt_foragers);
#elif 1
    web_generator = new OUWeb(ndim,std_mutation_factor*std_mutation_factor,
			      C0,lambda,saturation,q,correct_correlations);
#else
    web_generator = 
      new MoranIpp(ndim,std_mutation_factor*std_mutation_factor,
		   C0,lambda,saturation,q,correct_correlations,
		   0.05/*T*/,1000/*D_f*/,1/*over variability*/,0.001/*sensitivity*/,
		   0.5/*granularity*/);
#endif
  }
  
  average_meter av_Ddiet,av_C,av_Cy4,av_Nest,p_ones,p_interval,p_chordal;
  Interaction_Matrix im;
  int Ddiet_suggests_different=0;
  int nth_web=0;
  sequence<int> N;
  sequence<distribution> Ydist; // sum of distributions
  sequence<distribution> Ydist2; // sum of squares of distributions
  sequence<distribution> Rdist; // sum of distributions
  sequence<distribution> Rdist2; // sum of squares of distributions
  
  for(int repetitions_left=number_of_repetitions;repetitions_left-->0;){
    do{
      if(read_file){
	is >> im;
      }else{
	im=web_generator->draw_sample(raw_target_S);
      }
      nth_web++;
    }while(nth_web<0);
    
    if(rec_file)
      im.tPrint(rec);

    if(rewire_fraction){
      rewire(rewire_fraction,im);
    }
    
    Interaction_Matrix  dist_im = standardize(im);
    distribution dist; 
    const int Zbin=1;
    N[Zbin]++;
    N[0]++;
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
 
    int Cy4;
    double Ddiet;
    bool cons_ones,interval,chordal;
    av_Cy4.sample(Cy4=im.prop_Cy4());
//     REPORT(nth_web);
//     im.Print();
//     im.PPrint();
    av_C.sample(im.connectance_C());
    av_Nest.sample(im.prop_Nest());
    av_Ddiet.sample(Ddiet=im.prop_Ddiet(true));
//     try{
//       p_ones.sample(cons_ones=im.has_consecutive_ones());
//     }catch(AnalysisBug){
//       WARNING("Problem in consecutive ones test encountered, skipping this web:");
//       im.PPrint();
//     }
//     //     REPORT(Ddiet);
//     //     REPORT(cons_ones);
//     if(cons_ones && Ddiet>0){
//       static int counter=0;
//       im.PPrint();
//       REPORT(Ddiet);
//       REPORT(cons_ones);
//       FATAL_ERROR(++counter << " times cons-one inconsistency");
//     }
//     if(!cons_ones && Ddiet==0){
//       im.PPrint();
//       REPORT(Ddiet);
//       REPORT(cons_ones);
//       WARNING(++Ddiet_suggests_different << 
// 	      " times Ddiet suggests differently");
//     }
//     try{
//       p_interval.sample(interval=im.is_interval());
//       p_chordal.sample(chordal=im.is_chordal());
//     }catch(AnalysisBug){
//       WARNING("Problem in intervality test encountered, skipping this web:");
//       im.PPrint();
//     }
    
//     if(repetitions_left%(number_of_repetitions/10)==0){
//       cout << "start printing" << endl;
//       im./*tsort().*/PPrint();
//       cout << "finish printing" << endl;
//     }

//     ALWAYS_ASSERT((!chordal) || Cy4==0);
//     ALWAYS_ASSERT((!interval) || chordal);
//     ALWAYS_ASSERT((!cons_ones) || interval);
       
    
//     if(repetitions_left%(number_of_repetitions/10)==0){
//       REPORT(ndim);
//       REPORT(av_C);
//       REPORT(av_Ddiet);
//       REPORT(av_Cy4);
//       REPORT(av_Nest);
//       REPORT(p_ones);
//       REPORT(p_interval);
//       REPORT(p_chordal);
//     }
  }
  for(unsigned int Zbin=0; Zbin<Ydist.size(); Zbin++){
    double meanZ;
    distribution Ymean=Ydist[Zbin]/double(N[Zbin]);
    distribution Ystd=sqrt(double(N[Zbin])/double(N[Zbin]-1))*
      Map(::sqrt,Ydist2[Zbin]/double(N[Zbin])-Ymean*Ymean);
    distribution Yerror=Ystd/sqrt(double(N[Zbin]));
    distribution meanZbasis;
    if(cumulative_stats)
      meanZbasis=Ymean;
    else
      meanZbasis=Ymean.reverse_cumulative_sum();
    meanZ=accumulate(meanZbasis.begin(),meanZbasis.end(),0.0);
    distribution Ytheory;
    Interaction_Matrix dist_im=web_generator->draw_sample(raw_target_S);
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
  } // for Zbin
  REPORT(ndim);
  REPORT(av_C);
  REPORT(av_Ddiet);
  REPORT(av_Cy4);
  REPORT(av_Nest);
  REPORT(p_ones);
  REPORT(p_interval);
  REPORT(p_chordal);
  cout << ndim << " " 
       << av_Ddiet.readout() << " " 
       << av_Ddiet.error() << " "
       << av_Ddiet.std() << " "
       << av_Cy4.readout() << " "
       << av_Cy4.error() << " "
       << av_Cy4.std() << " "
       << av_Nest.readout() << " "
       << av_Nest.error() << " "
       << av_Nest.std() << " "
       << p_ones.readout() << " "
       << p_ones.error() << " "
       << p_ones.std() << " "
       << p_interval.readout() << " "
       << p_interval.error() << " "
       << p_interval.std() << " "
       << p_chordal.readout() << " "
       << p_chordal.error() << " "
       << p_chordal.std() << " "
       << Ddiet_suggests_different << " "
       << endl; 
} // main
