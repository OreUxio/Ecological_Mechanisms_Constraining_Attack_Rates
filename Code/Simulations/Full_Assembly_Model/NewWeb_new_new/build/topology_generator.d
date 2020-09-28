topology_generator.o topology_generator.d1 : topology_generator.cc  \
 topology_generator.h NetworkAnalysis.h error.h  \
  random.h \
 CMatrix.h  sequence.h  \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h
topology_generator : topology_generator.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o 
