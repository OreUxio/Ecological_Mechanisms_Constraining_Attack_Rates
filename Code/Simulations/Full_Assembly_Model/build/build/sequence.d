sequence.o sequence.d1 : sequence.cc  sequence.h \
  simple_vector.h error.h \
  
sequence : sequence.o simple_vector.o error.o 
