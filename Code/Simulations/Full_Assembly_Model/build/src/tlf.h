// -*- mode: c++ -*-
// $Id: tlf.h 543 2006-05-25 08:37:03Z cvsrep $
#ifndef _TLF_H_
#define _TLF_H_

#include "link_strength_matrix.h"
#include "snapshot.h"


// This reads food-web data on Tuesday Lake as copy-pasted from
// AdvEcolRes paper.
void read_tlf_file(const char * file_name, 
		   Snapshot & data);

#endif // _TLF_H_
