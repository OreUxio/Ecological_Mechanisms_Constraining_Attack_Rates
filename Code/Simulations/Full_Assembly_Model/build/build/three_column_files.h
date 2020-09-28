// -*- c++ -*-
// $Id: three_column_files.h 838 2007-01-09 06:46:53Z cvsrep $

#ifndef __THREE_COLUMN_FILES__
#define __THREE_COLUMN_FILES__

#include "NetworkAnalysis.h"
#include <string>

Interaction_Matrix load_three_column_file(std::string & name);
void save_three_column_file(Interaction_Matrix & m,const char * name);

#endif // __THREE_COLUMN_FILES__
