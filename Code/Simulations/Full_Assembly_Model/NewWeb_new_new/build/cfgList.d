cfgList.o cfgList.d1 : cfgList.cc  \
  cfgList.h parsecfg.h evaluate.h \
 CLHEP_Evaluator.h  \
 error.h  \
 
cfgList : cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o error.o 
