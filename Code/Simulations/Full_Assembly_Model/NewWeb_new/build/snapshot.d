snapshot.o snapshot.d1 : snapshot.cc  \
 NewMatrix.h   \
  snapshot.h \
 sequence.h simple_vector.h error.h  \
 link_strength_matrix.h NetworkAnalysis.h random.h \
  CMatrix.h \
  NetworkHelpers.h evaluate.h CLHEP_Evaluator.h \
 matrix_transformers.h xy_graph.h Statistics.h
snapshot : NewMatrix.o snapshot.o sequence.o simple_vector.o error.o link_strength_matrix.o NetworkAnalysis.o random.o CMatrix.o NetworkHelpers.o evaluate.o CLHEP_Evaluator.o matrix_transformers.o xy_graph.o Statistics.o 
