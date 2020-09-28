props_t.o props_t.d1 : props_t.cc   \
 random.h  \
 sequence.h   \
 simple_vector.h error.h  \
 props_t.h cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
  Statistics.h NewMatrix.h \
 
props_t : random.o sequence.o simple_vector.o error.o props_t.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o Statistics.o NewMatrix.o 
