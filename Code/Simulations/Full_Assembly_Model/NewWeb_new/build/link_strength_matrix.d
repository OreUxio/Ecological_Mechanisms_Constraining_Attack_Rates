link_strength_matrix.o link_strength_matrix.d1 : link_strength_matrix.cc \
  link_strength_matrix.h NetworkAnalysis.h \
 error.h   \
  random.h \
 CMatrix.h  sequence.h  \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h
link_strength_matrix : link_strength_matrix.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o 
