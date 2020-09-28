evaluate.o evaluate.d1 : evaluate.cc  \
  evaluate.h CLHEP_Evaluator.h \
 error.h  \
 
evaluate : evaluate.o CLHEP_Evaluator.o error.o 
