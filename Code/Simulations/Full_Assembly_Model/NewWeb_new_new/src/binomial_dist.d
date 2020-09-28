binomial_dist.o binomial_dist.d1 : binomial_dist.cc  \
 binomial_dist.h random.h  \
 error.h  \
 
binomial_dist : binomial_dist.o random.o error.o 
