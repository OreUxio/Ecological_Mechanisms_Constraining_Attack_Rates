NetworkAnalysis.o NetworkAnalysis.d1 : NetworkAnalysis.cc  \
 NetworkAnalysis.h error.h   \
  random.h \
 CMatrix.h  sequence.h  \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h \
 Statistics.h   \
  pqtree.h  \
 pqnode.h setTemplates.h
NetworkAnalysis : NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o Statistics.o pqtree.o pqnode.o setTemplates.o 
