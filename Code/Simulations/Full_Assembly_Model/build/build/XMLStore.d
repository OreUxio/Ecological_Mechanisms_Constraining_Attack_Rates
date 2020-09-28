XMLStore.o XMLStore.d1 : XMLStore.cc  \
  XMLStore.h \
  remember.h \
   error.h \
 
XMLStore : XMLStore.o remember.o error.o 
