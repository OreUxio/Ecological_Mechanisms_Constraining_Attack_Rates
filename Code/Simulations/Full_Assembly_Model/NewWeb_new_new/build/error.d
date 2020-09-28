error.o error.d1 : error.cc  error.h \
  cfgList.h parsecfg.h \
 evaluate.h CLHEP_Evaluator.h  \
  
error : error.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
