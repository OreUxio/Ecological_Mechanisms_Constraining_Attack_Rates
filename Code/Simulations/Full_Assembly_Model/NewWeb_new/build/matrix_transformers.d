matrix_transformers.o matrix_transformers.d1 : matrix_transformers.cc  \
  matrix_transformers.h \
 link_strength_matrix.h NetworkAnalysis.h error.h  \
  random.h  \
  CMatrix.h \
  sequence.h  \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h \
 polyfit.h   \
 Statistics.h evaluate.h CLHEP_Evaluator.h cfgList.h parsecfg.h \
 
matrix_transformers : matrix_transformers.o link_strength_matrix.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o polyfit.o Statistics.o evaluate.o CLHEP_Evaluator.o cfgList.o parsecfg.o 
