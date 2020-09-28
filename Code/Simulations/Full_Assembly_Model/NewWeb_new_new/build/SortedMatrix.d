SortedMatrix.o SortedMatrix.d1 : SortedMatrix.cc  SortedMatrix.h \
  vector_with_max.h \
  error.h \
 Statistics.h NewMatrix.h  \
  simple_vector.h \
 cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
 
SortedMatrix : SortedMatrix.o vector_with_max.o error.o Statistics.o NewMatrix.o simple_vector.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
