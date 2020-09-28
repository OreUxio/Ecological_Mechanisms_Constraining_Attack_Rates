CompMatrix.o CompMatrix.d1 : CompMatrix.cc  error.h \
  Statistics.h NewMatrix.h \
  simple_vector.h \
 CompMatrix.h
CompMatrix : error.o Statistics.o NewMatrix.o simple_vector.o CompMatrix.o 
