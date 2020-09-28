memory_power.o memory_power.d1 : memory_power.cc  memory_power.h \
 sequence.h   \
  simple_vector.h error.h \
  
memory_power : memory_power.o sequence.o simple_vector.o error.o 
