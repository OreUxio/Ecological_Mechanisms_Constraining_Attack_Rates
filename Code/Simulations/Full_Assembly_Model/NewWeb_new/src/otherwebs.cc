// -*- mode: c++ -*-
// $Id: otherwebs.cc 1945 2010-11-06 18:05:57Z axel $

#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <unistd.h>
#include <sys/types.h>
#include <fnmatch.h>
#include <string>
#include <limits.h> // to get the glibc version
#include "otherwebs.h"
#include "XMLStore.h"
#include "NewWeb.h"
#include "random.h"

const char * webdirlistfile="otherwebs.txt";

Otherwebs::Otherwebs(mode m):
  the_number_of_other_webs(0),
  the_total_number_of_animals(0),
  the_total_number_of_plants(0),
  the_total_animal_biomass(0),
  the_total_plant_biomass(0),
  the_mode(off),
  the_time_of_last_refresh(0)
{
  activate(m);
}

Otherwebs::~Otherwebs(){
}

int Otherwebs::activate(mode m){
  the_mode=m;
  if(m==off) return 0;
  // read names of directories containing other webs and check
  // consistency:
  errno=0;
  std::ifstream is(webdirlistfile);
  if(!is.good())
    FATAL_ERROR("problem opening " << webdirlistfile);
  struct stat statinfo,pwd_statinfo;
  if(stat(getenv("PWD"),&pwd_statinfo))
    SYS_ERROR();
  while(!is.eof()){
    int id; // id is currently not used!!
    is >> id;
    is >> the_dirs[the_number_of_other_webs];
    if(!is.good()) break;
    std::cerr << "this or other webdir " << the_dirs[the_number_of_other_webs] 
	      << std::endl;
    if(stat(the_dirs[the_number_of_other_webs].c_str(), &statinfo))
      SYS_ERROR();
    if(!S_ISDIR(statinfo.st_mode)){
      FATAL_ERROR(the_dirs[the_number_of_other_webs] << "is not a directory");
    }
    if(!((statinfo.st_mode & S_IXUSR)&&(statinfo.st_mode & S_IRUSR))){
      FATAL_ERROR("cannot read files in " << the_dirs[the_number_of_other_webs]);
    }      
    // exclude the current working directory
    if(statinfo.st_ino!=pwd_statinfo.st_ino){
      the_number_of_other_webs++;
    }else{
      this_dir=the_dirs[the_number_of_other_webs];
    }
  }
  the_dirs.resize(the_number_of_other_webs); //because current may be last

  // set the_other to empty webs:
  the_others=std::vector<NewWeb>(the_number_of_other_webs);
  the_file_stats=std::vector<struct stat>(the_number_of_other_webs);
  return the_number_of_other_webs;
}

#ifdef NO_CONST_WEBNAME_FILTER
static int webname_filter(struct dirent * dir_entry)
#else
static int webname_filter(const struct dirent * dir_entry)
#endif
{
  return !fnmatch("web*.xml*",dir_entry->d_name,0);
}


#if defined(__GLIBC__) && ((__GLIBC__ > 2) || (__GLIBC_MINOR__ >= 10))
static int index_sorter(const struct dirent **v1,
			const struct dirent **v2)
#else
static int index_sorter(const void * v1, const void * v2)
#endif
{
  const dirent * e1=*(dirent **)v1;
  const dirent * e2=*(dirent **)v2;
  return -atoi((e1->d_name)+3)+atoi((e2->d_name)+3);
}

static char * allocate_and_get_cwd(){
  size_t buffsize=256;
  char * buff=(char *)malloc(buffsize*sizeof(char));
  errno=0;
  while(!getcwd(buff,buffsize)){
    if(errno) SYS_ERROR();
    free(buff);
    buffsize+=256;
    buff=(char *)malloc(buffsize*sizeof(char));
  }
  return buff;
}


bool Otherwebs::get_newest_webs(){
  {
    // do nothing if we just got the newest webs:
    time_t current_time=time(0); //seconds accuracy
    if(current_time < the_time_of_last_refresh + 1)
      return false;
    the_time_of_last_refresh=current_time;
  }

  char * this_dir=allocate_and_get_cwd();
  for(int i=0;i<the_number_of_other_webs;i++){
    struct dirent **namelist;
    chdir(the_dirs[i].c_str());
    errno=0;
    int n=
      scandir(the_dirs[i].c_str(),&namelist,
	      webname_filter,index_sorter);
    if(n<0){
      REPORT(i);
      REPORT(the_dirs[i]);
      SYS_ERROR();
    }else if(n>0){
      // see if this is a new web:
      struct stat statinfo;
      stat(namelist[0]->d_name,&statinfo);
#define OTH_EQUAL_MEM(MEM) ((statinfo.MEM)==(the_file_stats[i].MEM))
      if(!( OTH_EQUAL_MEM(st_dev) && OTH_EQUAL_MEM(st_ino) && 
	    OTH_EQUAL_MEM(st_mtime) ) ){
#undef OTH_EQUAL_MEM
	the_file_stats[i]=statinfo;
	WARNING(this_dir << "reading " << the_dirs[i] << "/" 
		<< namelist[0]->d_name);
	// this is a new web, read it
	// wait for file to become readable:
	if(!((statinfo.st_mode) & S_IRUSR)){
	  WARNING("waiting for " << namelist[0]->d_name
		  << " to become readable");
	  while(!((statinfo.st_mode) & S_IRUSR)){
	    usleep((unsigned long)1e5);// wait 0.1 second
	    stat(namelist[0]->d_name,&statinfo);
	  }
	}
	read_web(namelist[0]->d_name,i);
      }

      // free scandir memory:
      while(n--) {
	free(namelist[n]);
      }
      free(namelist);
    }
  }
  // free getcwd memory:
  chdir(this_dir);
  free(this_dir);

  return true;
}    

bool Otherwebs::get_webs(int webnumber, int & exit_flag){
  if(webnumber<0){
    WARNING("negative web number");
    return false;
  }
  const char * webname_format="web%i.xml.bz2";
  char * this_dir=allocate_and_get_cwd();
  char * webname=(char *)
    malloc(1+snprintf(NULL,0,webname_format,webnumber)*sizeof(char));
  sprintf(webname,webname_format,webnumber);
  for(int i=0;i<the_number_of_other_webs;i++){
    chdir(the_dirs[i].c_str());
    struct stat statinfo;
    WARNING("reading " << the_dirs[i] << "/" << webname);
    errno=0;
    stat(webname,&statinfo);
    while(errno==ENOENT || !(statinfo.st_mode & S_IRUSR)) {
      // file does not exist or is not readable
      chdir(this_dir); // change to own diretory while waiting
      usleep((unsigned long)1e6);// wait 1 seconds
      chdir(the_dirs[i].c_str());
      errno=0;
      stat(webname,&statinfo);
      if(exit_flag) goto finish;
    }
    if(errno){
      SYS_ERROR();
    }
    the_file_stats[i]=statinfo;
    read_web(webname,i);
  }
 finish:
  // free getcwd memory:
  chdir(this_dir);
  free(webname);
  free(this_dir);
  return true;
}    

void Otherwebs::read_web(char * name, int i){
  the_total_number_of_animals-=
    the_others[i].number_of_animals();
  the_total_number_of_plants-=
    the_others[i].number_of_plants();
  the_total_animal_biomass-=
    the_others[i].animal_biomass();
  the_total_plant_biomass-=
    the_others[i].plant_biomass();

  the_others[i]=NewWeb(); // clear web
  XMLStore store(name);
  store.get(&the_others[i],"FoodWeb");

  the_total_number_of_animals+=
    the_others[i].number_of_animals();
  the_total_number_of_plants+=
    the_others[i].number_of_plants();
  the_total_animal_biomass+=
    the_others[i].animal_biomass();
  the_total_plant_biomass+=
    the_others[i].plant_biomass();
  if(the_total_plant_biomass<0){
    the_total_plant_biomass=0;
  }
  if(the_total_animal_biomass<0){
    the_total_animal_biomass=0;
  }
}


//     // get name of last file (we should do this with a call to
//     // scandir!)
//     const char * command_form="ls -1t %s/web*.xml*|head -1";
//     char * command=
//       malloc((strlen(command_form)*strlen(the_others[i])+1)*sizeof(char));
//     sprintf(command,command_form,the_others[i]);
//     last_webfile=popen(command,"r");
//     if(!feof(last_webfile)){
//       WARNING("directory " << the_others[i] << " contains no webs");
//     }else{
//       const int chunksize=256;
//       int webfilename_size=chunksize;
//       char * webfilename=malloc(webfilename_size);
//       do{
// 	read
//       }while(!feof(last_webfile));
    
int Otherwebs::exit_now_dummy=0;

const NewSpecies & Otherwebs::get_random_species(NewSpecies::taxon_t taxon){
  int specsum=
    (taxon==NewSpecies::plant?
     the_total_number_of_plants : the_total_number_of_animals );
  int specindex=random_integer(specsum);
  for(int i=the_number_of_other_webs;i-->0;){
    if(taxon==NewSpecies::plant){
      if(the_others[i].number_of_plants()>specindex){
	return the_others[i].
	  s[specindex+the_others[i].number_of_animals()];
      }else{
	specindex-=the_others[i].number_of_plants();
      }
    }else{
      if(the_others[i].number_of_animals()>specindex){
	return the_others[i].s[specindex];
      }else{
	specindex-=the_others[i].number_of_animals();
      }
    }
  }
  FATAL_ERROR("Could not find other species");
}
 
const NewSpecies & 
Otherwebs::get_random_species_by_biomass(NewSpecies::taxon_t taxon){
  double biomass_sum=
    (taxon==NewSpecies::plant ? 
     the_total_plant_biomass:
     the_total_animal_biomass);
  double specindex=biomass_sum*unirand();
  double cumulative_mass_sum=0;
  for(int i=the_number_of_other_webs;i-->0;){
    for(int j=the_others[i].number_of_species();j-->0;){
      if(the_others[i].s[j].taxon()==taxon){
	cumulative_mass_sum+=the_others[i].s[j].biomass_abundance_B();
	if(cumulative_mass_sum>specindex){
	  return the_others[i].s[j];
	}
      }
    }
  }
  WARNING("we have not found an appropriate species, returning anything");
  for(int i=the_number_of_other_webs;i-->0;){
    if(the_others[i].number_of_species()>0 and
       the_others[i].s[0].taxon()==taxon)
      return the_others[i].s[0];
  }
  FATAL_ERROR("cound not find a species");
}
  
