otherwebs.o otherwebs.d1 : otherwebs.cc  \
  otherwebs.h \
 NewSpecies.h   \
 remember.h  \
  error.h  \
  sequence.h \
 simple_vector.h  \
  md5.h \
  XMLStore.h \
  NewWeb.h \
  NetworkAnalysis.h random.h \
  CMatrix.h NewMatrix.h \
  NetworkHelpers.h \
 ODE.h  \
  link_strength_matrix.h snapshot.h \
 Statistics.h SortedVector.h \
  vector_with_max.h \
 SortedMatrix.h relax.h ptr_vector.h CompMatrix.h
otherwebs : otherwebs.o NewSpecies.o remember.o error.o sequence.o simple_vector.o md5.o XMLStore.o NewWeb.o NetworkAnalysis.o random.o CMatrix.o NewMatrix.o NetworkHelpers.o ODE.o link_strength_matrix.o snapshot.o Statistics.o SortedVector.o vector_with_max.o SortedMatrix.o relax.o ptr_vector.o CompMatrix.o 
