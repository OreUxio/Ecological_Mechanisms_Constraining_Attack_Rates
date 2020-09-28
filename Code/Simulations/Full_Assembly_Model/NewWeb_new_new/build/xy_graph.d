xy_graph.o xy_graph.d1 : xy_graph.cc  \
  xy_graph.h sequence.h \
  simple_vector.h error.h \
  
xy_graph : xy_graph.o sequence.o simple_vector.o error.o 
