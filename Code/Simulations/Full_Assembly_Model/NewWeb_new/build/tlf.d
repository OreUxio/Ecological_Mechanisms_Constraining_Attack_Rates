tlf.o tlf.d1 : tlf.cc   \
  tlf.h \
 link_strength_matrix.h NetworkAnalysis.h error.h  \
  random.h \
 CMatrix.h  sequence.h \
 simple_vector.h  \
  NewMatrix.h  \
  NetworkHelpers.h \
 snapshot.h evaluate.h CLHEP_Evaluator.h
tlf : tlf.o link_strength_matrix.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o snapshot.o evaluate.o CLHEP_Evaluator.o 
