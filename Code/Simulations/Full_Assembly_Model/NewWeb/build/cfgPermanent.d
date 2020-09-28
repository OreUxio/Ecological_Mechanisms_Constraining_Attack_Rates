cfgPermanent.o cfgPermanent.d1 : cfgPermanent.cc  cfgPermanent.h \
 remember.h  \
   error.h \
  cfgList.h parsecfg.h \
 evaluate.h CLHEP_Evaluator.h  \
  
cfgPermanent : cfgPermanent.o remember.o error.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
