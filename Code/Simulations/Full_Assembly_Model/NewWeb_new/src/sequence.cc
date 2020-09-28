// $Id: sequence.cc 419 2006-04-20 08:30:18Z cvsrep $
#include "sequence.h"

#undef simple_

std::ostream & operator<<(std::ostream &stream, 
			  const sequence<std::string> &s){
  std::simple_vector<std::string>::const_iterator i=s.begin(),e=s.end();
  while(i!=e) {
    stream << *i << std::endl;
    i++;
  }
  return stream;
}

