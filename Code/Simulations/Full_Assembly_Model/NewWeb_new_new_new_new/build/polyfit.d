polyfit.o polyfit.d1 : polyfit.cc  \
  polyfit.h \
 Statistics.h NewMatrix.h  \
  simple_vector.h \
 error.h sequence.h  \
 
polyfit : polyfit.o Statistics.o NewMatrix.o simple_vector.o error.o sequence.o 
