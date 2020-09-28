scheduler.o scheduler.d1 : scheduler.cc  scheduler.h \
  error.h  \
 
scheduler : scheduler.o error.o 
