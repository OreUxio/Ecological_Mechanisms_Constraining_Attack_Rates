Statistics.o Statistics.d1 : Statistics.cc  Statistics.h \
 NewMatrix.h   \
  simple_vector.h \
 error.h   \
  
Statistics : Statistics.o NewMatrix.o simple_vector.o error.o 
