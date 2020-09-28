s_t.o s_t.d1 : s_t.cc   \
   s_t.h \
 props_t.h binomial_dist.h error.h  \
   random.h \
 cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
 
s_t : s_t.o props_t.o binomial_dist.o error.o random.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
