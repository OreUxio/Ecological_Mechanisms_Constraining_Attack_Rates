relax.o relax.d1 : relax.cc   \
 relax.h   \
 ODE.h   \
  error.h \
  sequence.h \
 simple_vector.h  \
 polyfit.h   \
 Statistics.h NewMatrix.h  \
  random.h \
  evaluate.h \
 CLHEP_Evaluator.h cfgList.h parsecfg.h  \
  
relax : relax.o ODE.o error.o sequence.o simple_vector.o polyfit.o Statistics.o NewMatrix.o random.o evaluate.o CLHEP_Evaluator.o cfgList.o parsecfg.o 
