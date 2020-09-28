index_set.o index_set.d1 : index_set.cc  index_set.h \
 sequence.h   \
  simple_vector.h error.h \
  
index_set : index_set.o sequence.o simple_vector.o error.o 
