mass_integrator.o mass_integrator.d1 : mass_integrator.cc  \
  error.h  \
  mass_integrator.h \
 sequence.h simple_vector.h  \
 evaluate.h CLHEP_Evaluator.h cfgList.h parsecfg.h \
 
mass_integrator : error.o mass_integrator.o sequence.o simple_vector.o evaluate.o CLHEP_Evaluator.o cfgList.o parsecfg.o 
