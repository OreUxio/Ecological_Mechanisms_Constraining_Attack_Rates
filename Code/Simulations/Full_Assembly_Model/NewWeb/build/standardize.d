standardize.o standardize.d1 : standardize.cc  \
 NetworkAnalysis.h error.h   \
  random.h \
 CMatrix.h  sequence.h  \
  simple_vector.h \
 NewMatrix.h   \
  NetworkHelpers.h \
 cfgList.h parsecfg.h evaluate.h CLHEP_Evaluator.h \
 
standardize : NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o cfgList.o parsecfg.o evaluate.o CLHEP_Evaluator.o 
