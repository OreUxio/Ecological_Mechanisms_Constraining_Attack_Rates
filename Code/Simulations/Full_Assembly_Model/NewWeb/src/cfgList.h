//$Id: cfgList.h 2451 2015-04-24 16:15:59Z axel $

#ifndef __CFGLIST_H__
#define  __CFGLIST_H__

/// \file High-level handling of configuration parameters.

#include "parsecfg.h"
#include <list>
#include <vector>
#include <string>

//using namespace std;

typedef cfgStruct *cfgStructP;

cfgStruct *full_cfg_list();

struct cfg_list_t;

struct cfg_list_t {
  cfgStructP current;
  cfg_list_t * next;
};

/// Used to map static configuration variables to a global configuration file.
/** Look for "Manage adjustable parameters:" to see how to use this. */
class cfg_add{
  //  static std::list<cfgStructP> cfg_list;
  static cfg_list_t * cfg_list;
  cfgStruct * the_s;
  friend cfgStruct *full_cfg_list();
 public:
  cfg_add(cfgStruct * s);
  ~cfg_add();
};

void read_parameters_from_file(const std::string in_file_name);
double get_cfg_parameter(const char * name);
int set_cfg_parameter(const char * name, const char * value);
void write_cfg_template(const char * filename);
void write_cfg_template(std::ostream & file);
void do_assignments(const std::vector< std::string > & assignments);

#define CFGDOUBLE(X) {#X, CFG_DOUBLE, &(X)}
#define CFGINT(X) {#X, CFG_INT, &(X)}
#define CFGBOOL(X) {#X, CFG_BOOL, &(X)}
#define CFGSTRING(X) {#X, CFG_STRING, &(X)}
#endif //  __CFGLIST_H__
