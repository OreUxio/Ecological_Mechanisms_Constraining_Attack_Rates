perturb.o perturb.d1 : perturb.cc  perturb.h NewWeb.h \
  NetworkAnalysis.h error.h \
 random.h  \
  CMatrix.h sequence.h \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h \
 ODE.h  \
  NewSpecies.h remember.h \
  md5.h \
 link_strength_matrix.h snapshot.h otherwebs.h \
  Statistics.h SortedVector.h \
  vector_with_max.h \
 SortedMatrix.h relax.h ptr_vector.h CompMatrix.h evaluate.h \
 CLHEP_Evaluator.h
perturb : perturb.o NewWeb.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o ODE.o NewSpecies.o remember.o md5.o link_strength_matrix.o snapshot.o otherwebs.o Statistics.o SortedVector.o vector_with_max.o SortedMatrix.o relax.o ptr_vector.o CompMatrix.o evaluate.o CLHEP_Evaluator.o 
