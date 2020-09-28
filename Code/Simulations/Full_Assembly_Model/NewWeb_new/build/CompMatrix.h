// -*- mode: c++ -*-
// $Id: CompMatrix.h 3 2005-12-01 07:13:32Z cvsrep $
#ifndef _COMPMATRIX_H_
#define _COMPMATRIX_H_


#include "NewMatrix.h"

int CompMatrix(NewMatrix & Chat, NewMatrix & iChat,
	       NewMatrix & epsAT, NewMatrix & A, NewMatrix C);

int CompMatrixSchur(NewMatrix & Chat, NewMatrix & iChat,
		    NewMatrix & epsAT, NewMatrix & A, NewMatrix C);


#endif // _COMPMATRIX_H_
