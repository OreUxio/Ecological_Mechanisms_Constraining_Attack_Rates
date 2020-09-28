cell.o cell.d1 : cell.cc   \
  cell.h \
  error.h  \
  NewSpecies.h \
 remember.h  \
  sequence.h  \
  simple_vector.h \
   md5.h \
 cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
 
cell : cell.o error.o NewSpecies.o remember.o sequence.o simple_vector.o md5.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
