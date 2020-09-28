period_cutter.o period_cutter.d1 : period_cutter.cc  \
 period_cutter.h sequence.h   \
  simple_vector.h error.h \
 cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
 
period_cutter : period_cutter.o sequence.o simple_vector.o error.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
