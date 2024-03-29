// -*- mode: c++ -*-
// $Id: perturb.cc 2464 2016-05-01 13:11:20Z axel $

#include "perturb.h"
#include "evaluate.h"

using namespace std;


// parse perturbation_type parameter 
bool parametrized(const char * perturbation_type,
		  const char * tag,
		  double & value){
  size_t len=strlen(tag);
  if(0==strncmp(perturbation_type,tag,len) && perturbation_type[len]=='='){
    value=eval(perturbation_type+len+1);
    return true;
  }else{
    return false;
  }
}

void perturb(NewWeb & web,const char * perturbation_type){
  if(0==strcmp(perturbation_type,"help")){
    cout << "known perturbation types:" << endl;
    cout << "'none': no change" << endl;
    cout << "'acid=x: extinguish x% of plants by acidification" << endl;
    cout << "'down=x': delete the largest x% of all species" << endl;
    cout << "'envi=x: change environment by x%" << endl;
    cout << "'fish=x': increase mortality of an abundand large species by x" << endl;
    cout << "'fishF=x': apply fishing mortality x to all animals >= 1g" << endl;
    cout << "'fishMrange=x': applies random fishing mortality to all consumers with M > x kg and M < fishMrangemax kg" << endl;
    cout << "'fishMrangeCS=x': applies fishing mortality to all consumers M > x kg and M < fishMrangemax kg according to specific values for 6 species" << endl;
    cout << "'fishMrangeconst=x': applies constant fishing mortality to all consumers with M > x kg and M < fishMrangemax kg" << endl;
    cout << "'fishMrangeline=x': applies fishing mortality to all consumers with M > x kg and M < fishMrangemax kg according to straight line" << endl;
    cout << "'fishTLrange=x': applies random fishing mortality to all consumers with TL > x and TL < fishTLrangemax" << endl;
    cout << "'fishTLrangeconst=x': applies constant fishing mortality to all consumers with TL > x and TL < fishTLrangemax" << endl;
    cout << "'fishTLrangeline=x': applies fishing mortality to all consumers with TL > x and TL < fishTLrangemax according to straight line" << endl;
    cout << "'fishall=x': applies random fishing mortality to all consumers with M > x kg" << endl;
    cout << "'mid10=x': delete 10 species closest to body mass x" << endl;
    cout << "'plants=x': delete x% of all plants" << endl;
    cout << "'rare=x: extinguish x% of (biomass) rare sepcies" << endl;
    cout << "'reduce=x: reduce all biomasses to x%" << endl;
    cout << "'deleterarefish=x': deletes rare fish species making up last x% of total fish biomass" << endl;
    cout << "'deleterarefish5per2=x': deletes rare fish species making up last 5% of total fish biomass, according to cut-off line with slope x" << endl;
    cout << "'deleterarefish10per2=x': deletes rare fish species making up last 10% of total fish biomass, according to cut-off line with slope x" << endl;
    cout << "'deletefishspeciesatrandom=x': deletes x randomly chosen fish species" << endl;
    cout << "'deletefishspeciesbysize=x': deletes x fish species by order of decreasing size" << endl;
    cout << "'deletefishspeciesbysizereverse=x': deletes x fish species by order of increasing size" << endl;
    cout << "'deletefishspeciesbybiomass=x': deletes x fish species by order of decreasing biomass" << endl;
    cout << "'deletefishspeciesbybiomassreverse=x': deletes x fish species by order of increasing biomass" << endl;
    cout << "'deletefishspeciesbyTL=x': deletes x fish species by order of decreasing TL" << endl;
    cout << "'deletefishspeciesbyTLreverse=x': deletes x fish species by order of increasing TL" << endl;
    cout << "'deletefishspeciesbyConn=x': deletes x fish species by order of decreasing Conn" << endl;
    cout << "'deletefishspeciesbyConnreverse=x': deletes x fish species by order of increasing Conn" << endl;
    cout << "'deleteallspecieswithbiomasslessthan=x': deletes species with biomass less than x" << endl;
    cout << "'deletefishspeciesexcept1=x': deletes all fish species except the one with column number x" << endl;
    cout << "'BottomUp=x': reduce net productivity of all producers by x%" << endl;
    // ... more types
    exit(0);
  }
  web.initialize_for_integration();
  if(0==strcmp(perturbation_type,"none")){
    return;
  }
  double x;
  if(parametrized(perturbation_type,"fish",x)){
    web.fish(x);
    return;
  }
  if(parametrized(perturbation_type,"fishF",x)){
    web.fishF(x);
    return;
  }
  if(parametrized(perturbation_type,"down",x)){
    web.delete_large_species_fraction(0.01*x);
    return;
  }
  if(parametrized(perturbation_type,"mid10",x)){
    web.delete_10_species_near_size(x);
    return;
  }
  if(parametrized(perturbation_type,"acid",x)){
    web.delete_plants_by_vulnerability_trait(0.01*x);
    return;
  }
  if(parametrized(perturbation_type,"rare",x)){
    web.delete_rare_species_fraction(0.01*x);
    return;
  }
  if(parametrized(perturbation_type,"plants",x)){
    web.delete_plant_fraction(0.01*x);
    return;
  }
  if(parametrized(perturbation_type,"envi",x)){
    web.model_environmenatal_change(0.01*x);
    return;
  }
  if(parametrized(perturbation_type,"reduce",x)){
    x*=0.01;
    REPORT(x);
    double mass_threshold=0;//eval("100*gram");
    REPORT(mass_threshold);
    for(int i=web.number_of_species();i-->0;){
      if(web.s(i).bodymass()>mass_threshold){
	web.s(i).
	  set_biomass_abundance_B(web.s(i).biomass_abundance_B()*x);
      }
    }
    return;
  }
  if(parametrized(perturbation_type,"fishall",x)){
    web.fishall(x);
    return;
  }
  if(parametrized(perturbation_type,"fishMrange",x)){
    web.fishMrange(x);
    return;
  }
  if(parametrized(perturbation_type,"fishTLrange",x)){
    web.fishTLrange(x);
    return;
  }
  if(parametrized(perturbation_type,"fishMrangeconst",x)){
    web.fishMrangeconst(x);
    return;
  }
  if(parametrized(perturbation_type,"fishTLrangeconst",x)){
    web.fishTLrangeconst(x);
    return;
  }
  if(parametrized(perturbation_type,"fishMrangeline",x)){
    web.fishMrangeline(x);
    return;
  }
  if(parametrized(perturbation_type,"fishTLrangeline",x)){
    web.fishTLrangeline(x);
    return;
  }
  if(parametrized(perturbation_type,"fishMrangeCS",x)){
    web.fishMrangeCS(x);
    return;
  }
  if(parametrized(perturbation_type,"deleterarefish",x)){
    web.delete_rare_fish_species(x);
    return;
  }
  if(parametrized(perturbation_type,"deleterarefish5per2",x)){
    web.delete_rare_fish_species_5per_2(x);
    return;
  }
  if(parametrized(perturbation_type,"deleterarefish10per2",x)){
    web.delete_rare_fish_species_10per_2(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesatrandom",x)){
    web.delete_fish_species_at_random(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbysize",x)){
    web.delete_fish_species_by_size(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbysizereverse",x)){
    web.delete_fish_species_by_size_reverse(x);
    return;
  }
  if(parametrized(perturbation_type,"custom",x)){
    web.custom_perturbation(x);
    return;
  }
  if(parametrized(perturbation_type,"deleteallspecieswithbiomasslessthan",x)){
    web.delete_all_species_with_biomass_less_than(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesexcept1",x)){
    web.delete_fish_species_except_1(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbybiomass",x)){
    web.delete_fish_species_by_biomass(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbybiomassreverse",x)){
    web.delete_fish_species_by_biomass_reverse(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbyTL",x)){
    web.delete_fish_species_by_TL(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbyTLreverse",x)){
    web.delete_fish_species_by_TL_reverse(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbyConn",x)){
    web.delete_fish_species_by_Conn(x);
    return;
  }
  if(parametrized(perturbation_type,"deletefishspeciesbyConnreverse",x)){
    web.delete_fish_species_by_Conn_reverse(x);
    return;
  }
  if(parametrized(perturbation_type,"BottomUp",x)){
    web.scale_primary_productivity(x/100+1);
  }
  FATAL_ERROR("unknown perturbation type");
}
