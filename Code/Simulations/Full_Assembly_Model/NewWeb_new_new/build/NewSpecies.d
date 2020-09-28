NewSpecies.o NewSpecies.d1 : NewSpecies.cc  NewSpecies.h \
  remember.h \
   error.h \
  sequence.h \
  simple_vector.h \
   md5.h \
 random.h  \
  evaluate.h \
 CLHEP_Evaluator.h Statistics.h NewMatrix.h  \
  cfgList.h parsecfg.h \
 
NewSpecies : NewSpecies.o remember.o error.o sequence.o simple_vector.o md5.o random.o evaluate.o CLHEP_Evaluator.o Statistics.o NewMatrix.o cfgList.o parsecfg.o 
