nls_web.o nls_web.d1 : nls_web.cc  nls_web.h \
 sequence.h   \
  simple_vector.h error.h \
  
nls_web : nls_web.o sequence.o simple_vector.o error.o 
