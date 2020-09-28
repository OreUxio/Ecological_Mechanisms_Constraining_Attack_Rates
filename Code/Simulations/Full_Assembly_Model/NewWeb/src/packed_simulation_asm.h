// -*- mode: c++ -*-
// $Id: packed_simulation_asm.h 3 2005-12-01 07:13:32Z cvsrep $
#ifndef _PACKED_SIMULATION_ASM_H_
#define _PACKED_SIMULATION_ASM_H_

extern "C" {
  double asm_fast_sandwich_product(const double * v1,
				   const double * v2,
				   double factor,  
				   char * __restrict__ tsk
				   );
}

#endif // _PACKED_SIMULATION_ASM_H_

