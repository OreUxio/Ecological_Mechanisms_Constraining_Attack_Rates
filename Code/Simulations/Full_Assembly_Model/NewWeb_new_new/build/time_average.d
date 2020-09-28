time_average.o time_average.d1 : time_average.cc  time_average.h \
 sequence.h   \
  simple_vector.h error.h \
 NewWeb.h   \
  NetworkAnalysis.h random.h \
  CMatrix.h NewMatrix.h \
  NetworkHelpers.h \
 ODE.h  \
  NewSpecies.h remember.h \
  md5.h \
 link_strength_matrix.h snapshot.h otherwebs.h \
  Statistics.h SortedVector.h \
  vector_with_max.h \
 SortedMatrix.h relax.h ptr_vector.h CompMatrix.h period_cutter.h \
 Integrator.h cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h
time_average : time_average.o sequence.o simple_vector.o error.o NewWeb.o NetworkAnalysis.o random.o CMatrix.o NewMatrix.o NetworkHelpers.o ODE.o NewSpecies.o remember.o md5.o link_strength_matrix.o snapshot.o otherwebs.o Statistics.o SortedVector.o vector_with_max.o SortedMatrix.o relax.o ptr_vector.o CompMatrix.o period_cutter.o Integrator.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
