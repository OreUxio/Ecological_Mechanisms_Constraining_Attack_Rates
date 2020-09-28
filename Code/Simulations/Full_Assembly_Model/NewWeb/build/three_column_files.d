three_column_files.o three_column_files.d1 : three_column_files.cc  \
  three_column_files.h \
 NetworkAnalysis.h error.h   \
  random.h \
 CMatrix.h  sequence.h \
 simple_vector.h  \
  NewMatrix.h  \
  NetworkHelpers.h
three_column_files : three_column_files.o NetworkAnalysis.o error.o random.o CMatrix.o sequence.o simple_vector.o NewMatrix.o NetworkHelpers.o 
