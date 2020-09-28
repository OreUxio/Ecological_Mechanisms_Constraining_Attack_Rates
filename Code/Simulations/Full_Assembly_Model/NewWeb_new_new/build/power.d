power.o power.d1 : power.cc   \
  polyfit.h \
 Statistics.h NewMatrix.h  \
  simple_vector.h \
 error.h sequence.h  \
  random.h \
  cfgList.h parsecfg.h \
 evaluate.h CLHEP_Evaluator.h  \
  
power : polyfit.o Statistics.o NewMatrix.o simple_vector.o error.o sequence.o random.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
