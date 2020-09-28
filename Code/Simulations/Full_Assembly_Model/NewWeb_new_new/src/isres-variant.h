/* $Id:$ */
/* The file isres.c from nlopt 2.3, plus necessary header files to
 *  make this run, plus modifications by Axel G. Rossberg to adapt to
 *  the problem of finding MSY finding in the PDMM. */

/* Copyright (c) 2007-2012 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef NLOPT_UTIL_H
#define NLOPT_UTIL_H

#include <stdlib.h>
#include <math.h>

#include "nlopt.h"

/* workaround for Solaris + gcc 3.4.x bug (see configure.ac) */
#if defined(__GNUC__) && defined(REPLACEMENT_HUGE_VAL)
#  undef HUGE_VAL
#  define HUGE_VAL REPLACEMENT_HUGE_VAL
#endif

#ifndef HAVE_COPYSIGN
   /* not quite right for y == -0, but good enough for us */
#  define copysign(x, y) ((y) < 0 ? -fabs(x) : fabs(x))
#endif

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

int nlopt_isinf(double x);

/* re-entrant qsort */
extern void nlopt_qsort_r(void *base_, size_t nmemb, size_t size, void *thunk,
			  int (*compar)(void *, const void *, const void *));

/* seconds timer */
extern double nlopt_seconds(void);
extern unsigned long nlopt_time_seed(void);

/* pseudorandom number generation by Mersenne twister algorithm */
extern void nlopt_init_genrand(unsigned long s);
extern double nlopt_urand(double a, double b);
extern int nlopt_iurand(int n);
extern double nlopt_nrand(double mean, double stddev);

/* Sobol' low-discrepancy-sequence generation */
typedef struct nlopt_soboldata_s *nlopt_sobol;
extern nlopt_sobol nlopt_sobol_create(unsigned sdim);
extern void nlopt_sobol_destroy(nlopt_sobol s);
extern void nlopt_sobol_next01(nlopt_sobol s, double *x);
extern void nlopt_sobol_next(nlopt_sobol s, double *x,
			    const double *lb, const double *ub);
extern void nlopt_sobol_skip(nlopt_sobol s, unsigned n, double *x);

/* stopping criteria */
typedef struct {
     unsigned n;
     double minf_max;
     double ftol_rel;
     double ftol_abs;
     double xtol_rel;
     const double *xtol_abs;
     int nevals, maxeval;
     double maxtime, start;
     int *force_stop;
} nlopt_stopping;
extern int nlopt_stop_f(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_ftol(const nlopt_stopping *stop, double f, double oldf);
extern int nlopt_stop_x(const nlopt_stopping *stop, 
			const double *x, const double *oldx);
extern int nlopt_stop_dx(const nlopt_stopping *stop, 
			 const double *x, const double *dx);
extern int nlopt_stop_xs(const nlopt_stopping *stop, 
			 const double *xs, const double *oldxs,
			 const double *scale_min, const double *scale_max);
extern int nlopt_stop_evals(const nlopt_stopping *stop);
extern int nlopt_stop_time_(double start, double maxtime);
extern int nlopt_stop_time(const nlopt_stopping *stop);
extern int nlopt_stop_evalstime(const nlopt_stopping *stop);
extern int nlopt_stop_forced(const nlopt_stopping *stop);

/* for local optimizations, temporarily setting eval/time limits */
extern nlopt_result nlopt_optimize_limited(nlopt_opt opt, 
					   double *x, double *minf,
					   int maxevals, double maxtime);

/* data structure for nonlinear inequality or equality constraint
   (f <= 0 or f = 0, respectively).  tol (>= 0) is a tolerance
   that is used for stopping criteria -- the point is considered
   "feasible" for purposes of stopping if the constraint is violated
   by at most tol. */
typedef struct {
     unsigned m; /* dimensional of constraint: mf maps R^n -> R^m */
     nlopt_func f; /* one-dimensional constraint, requires m == 1 */
     nlopt_mfunc mf;
     nlopt_precond pre; /* preconditioner for f (NULL if none or if mf) */
     void *f_data;
     double *tol;
} nlopt_constraint;

extern unsigned nlopt_count_constraints(unsigned p, const nlopt_constraint *c);
extern unsigned nlopt_max_constraint_dim(unsigned p, const nlopt_constraint *c);
extern void nlopt_eval_constraint(double *result, double *grad,
				  const nlopt_constraint *c,
				  unsigned n, const double *x);

/* rescale.c: */
double *nlopt_compute_rescaling(unsigned n, const double *dx);
double *nlopt_new_rescaled(unsigned n, const double *s, const double *x);
void nlopt_rescale(unsigned n, const double *s, const double *x, double *xs);
void nlopt_unscale(unsigned n, const double *s, const double *x, double *xs);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
/* Copyright (c) 2007-2012 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#ifndef ISRES_H
#define ISRES_H

#include "nlopt.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

nlopt_result isres_minimize(int n, nlopt_func f, void *f_data,
			    int m, nlopt_constraint *fc, /* fc <= 0  */
			    int p, nlopt_constraint *h, /* h == 0 */
			    const double *lb, const double *ub, /* bounds */
			    double *x, /* in: initial guess, out: minimizer */
			    double *minf,
			    nlopt_stopping *stop,
			    int population); /* init. population */

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif


/* Copyright (c) 2010 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>


/* Copyright (c) 2010 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

/* #include "isres.h" done above */

/* Improved Stochastic Ranking Evolution Strategy (ISRES) algorithm
   for nonlinearly-constrained global optimization, based on
   the method described in:

   Thomas Philip Runarsson and Xin Yao, "Search biases in constrained
   evolutionary optimization," IEEE Trans. on Systems, Man, and Cybernetics
   Part C: Applications and Reviews, vol. 35 (no. 2), pp. 233-243 (2005).

   It is a refinement of an earlier method described in:

   Thomas P. Runarsson and Xin Yao, "Stochastic ranking for constrained
   evolutionary optimization," IEEE Trans. Evolutionary Computation,
   vol. 4 (no. 3), pp. 284-294 (2000).

   This is an independent implementation by S. G. Johnson (2009) based
   on the papers above.  Runarsson also has his own Matlab implemention
   available from his web page: http://www3.hi.is/~tpr
*/

static int key_compare(void *keys_, const void *a_, const void *b_)
{
     const double *keys = (const double *) keys_;
     const int *a = (const int *) a_;
     const int *b = (const int *) b_;
     return keys[*a] < keys[*b] ? -1 : (keys[*a] > keys[*b] ? +1 : 0);
}

static unsigned imax2(unsigned a, unsigned b) { return (a > b ? a : b); }

nlopt_result isres_minimize(int n, nlopt_func f, void *f_data,
			    int m, nlopt_constraint *fc, /* fc <= 0 */
			    int p, nlopt_constraint *h, /* h == 0 */
			    const double *lb, const double *ub, /* bounds */
			    double *x, /* in: initial guess, out: minimizer */
			    double *minf,
			    nlopt_stopping *stop,
			    int population) /* pop. size (= 0 for default) */
{
     const double ALPHA = 0.5; /* here hard-coded as "sqrt(..)" */
     const double GAMMA = 0.85; /* step-reduction factor from paper */
     const double PHI = 1.0/ALPHA; /* adjusted expected rate of convergence */
     const double PF = 0.45; /* fitness probability, from paper */
     const double SURVIVOR = 1.0/7.0; /* survivor fraction, from paper */
     int survivors;
     nlopt_result ret = NLOPT_SUCCESS;
     double *sigmas = 0, *xs; /* population-by-n arrays (row-major) */
     double *fval; /* population array of function vals */
     double *penalty; /* population array of penalty vals */
     double *x0;
     int *rank, *irank = 0;
     int k, i, j, c;
     int mp = m + p;
     double minf_penalty = HUGE_VAL, minf_gpenalty = HUGE_VAL;
     double taup, tau;
     double *results = 0; /* scratch space for mconstraint results */
     unsigned ires;
     int denovo = 1; /* indicates are we doing the first iteration */
     FILE* saved_state = NULL;



     *minf = HUGE_VAL;

     if (!population) population = 20 * (n + 1);
     if (population < 1) return NLOPT_INVALID_ARGS;
     survivors = ceil(population * SURVIVOR);

     taup = PHI / sqrt(2*n);
     tau = PHI / sqrt(2*sqrt(n));

     /* we don't handle unbounded search regions */
     for (j = 0; j < n; ++j) if (nlopt_isinf(lb[j]) || nlopt_isinf(ub[j]))
				  return NLOPT_INVALID_ARGS;

     ires = imax2(nlopt_max_constraint_dim(m, fc),
		  nlopt_max_constraint_dim(p, h));
     results = (double *) malloc(ires * sizeof(double));
     if (ires > 0 && !results) return NLOPT_OUT_OF_MEMORY;

     sigmas = (double*) malloc(sizeof(double) * (population*n*2
						 + population
						 + population
						 + n));
     if (!sigmas) { free(results); return NLOPT_OUT_OF_MEMORY; }
     xs = sigmas + population*n;
     fval = xs + population*n;
     penalty = fval + population;
     x0 = penalty + population;

     irank = (int*) malloc(sizeof(int) * 2 * population);
     if (!irank) { ret = NLOPT_OUT_OF_MEMORY; goto done; }
     rank = irank + population;

     saved_state=fopen("isres-state.bin","r");
     if(saved_state!=NULL){
       printf("!!!! READING OPTIMIZER STATE FROM isres-state.bin\n");
       fread(sigmas,sizeof(double),population*n,saved_state);
       fread(xs,sizeof(double),population*n,saved_state);
       if(ferror(saved_state)){
	 printf("!!!! READ OF SAVED STATE FAILED. EXITING\n");
	 exit(1);
       }
       if(fclose(saved_state)){
	 printf("!!!! READ OF SAVED STATE FAILED. EXITING\n");
	 exit(1);
       }
       printf("!!!! DONE\n");
     }else{
       for (k = 0; k < population; ++k) {
	 for (j = 0; j < n; ++j) {
	   sigmas[k*n+j] = (ub[j] - lb[j]) / sqrt(n);
	   xs[k*n+j] = nlopt_urand(lb[j], ub[j]);
	 }
       }
       memcpy(xs, x, sizeof(double) * n); /* use input x for xs_0 */
     }

     while (1) { /* each loop body = one generation */
	  
          /* First, save state just in case we are killed */
          if(!denovo){
	    saved_state=fopen("isres-state.bin","w");
	    if(saved_state!=NULL){
	      printf("!!!! WRITING OPTIMIZER STATE TO isres-state.bin\n");
	      int c=fwrite(sigmas,sizeof(double),population*n,saved_state);
	      c=fwrite(xs,sizeof(double),population*n,saved_state);
	      if(ferror(saved_state)){
		printf("!!!! WRITE OF SAVED STATE FAILED!\n");
	      }
	      fclose(saved_state);
	      printf("!!!! DONE\n");fflush(stdout);
	    }
	  }

	  int all_feasible = 1;

	  /* evaluate f and constraint violations for whole population */
	  for (k = 0; k < population; ++k) {
	       int feasible = 1;
	       double gpenalty;
	       stop->nevals++;
	       fval[k] = f(n, xs + k*n, NULL, f_data);
	       if (nlopt_stop_forced(stop)) { 
		    ret = NLOPT_FORCED_STOP; goto done; }
	       penalty[k] = 0;
	       for (c = 0; c < m; ++c) { /* inequality constraints */
		    nlopt_eval_constraint(results, NULL,
					  fc + c, n, xs + k*n);
		    if (nlopt_stop_forced(stop)) { 
			 ret = NLOPT_FORCED_STOP; goto done; }
		    for (ires = 0; ires < fc[c].m; ++ires) {
			 double gval = results[ires];
			 if (gval > fc[c].tol[ires]) feasible = 0;
			 if (gval < 0) gval = 0;
			 penalty[k] += gval*gval;
		    }
	       }
	       gpenalty = penalty[k];
	       for (c = m; c < mp; ++c) { /* equality constraints */
		    nlopt_eval_constraint(results, NULL,
					  h + (c-m), n, xs + k*n);
		    if (nlopt_stop_forced(stop)) { 
			 ret = NLOPT_FORCED_STOP; goto done; }
		    for (ires = 0; ires < h[c-m].m; ++ires) {
			 double hval = results[ires];
			 if (fabs(hval) > h[c-m].tol[ires]) feasible = 0;
			 penalty[k] += hval*hval;
		    }
	       }
	       if (penalty[k] > 0) all_feasible = 0;

	       /* convergence criteria (FIXME: improve?) */

	       /* FIXME: with equality constraints, how do
		  we decide which solution is the "best" so far?
		  ... need some total order on the solutions? */

	       if ((penalty[k] <= minf_penalty || feasible)
		   && (fval[k] <= *minf || minf_gpenalty > 0)
		   && ((feasible ? 0 : penalty[k]) != minf_penalty
		       || fval[k] != *minf)) {
		    if (fval[k] < stop->minf_max && feasible) 
			 ret = NLOPT_MINF_MAX_REACHED;
		    else if (!nlopt_isinf(*minf)) {
			 if (nlopt_stop_f(stop, fval[k], *minf)
			     && nlopt_stop_f(stop, feasible ? 0 : penalty[k], 
					     minf_penalty))
			      ret = NLOPT_FTOL_REACHED;
			 else if (nlopt_stop_x(stop, xs+k*n, x))
			      ret = NLOPT_XTOL_REACHED;
		    }
		    memcpy(x, xs+k*n, sizeof(double)*n);
		    *minf = fval[k];
		    minf_penalty = feasible ? 0 : penalty[k];
		    minf_gpenalty = feasible ? 0 : gpenalty;
		    if (ret != NLOPT_SUCCESS) goto done;
	       }

	       if (nlopt_stop_forced(stop)) ret = NLOPT_FORCED_STOP;
	       else if (nlopt_stop_evals(stop)) ret = NLOPT_MAXEVAL_REACHED;
	       else if (nlopt_stop_time(stop)) ret = NLOPT_MAXTIME_REACHED;
	       if (ret != NLOPT_SUCCESS) goto done;
	       
	  }

	  /* stop also if all x in population are close to each other */
	  for(k = 0; k < population ; ++k){
	    for(i = k+1; i < population ; ++i){
	      if( !nlopt_stop_x(stop, xs+k*n, xs+i*n) )
		goto not_internal_xtol;
	    }
	  }
	  ret = NLOPT_XTOL_REACHED;
	  goto done;
     not_internal_xtol:
	  
	  denovo = 0; /* no need to evaluate survivors again */

	  /* "selection" step: rank the population */
	  for (k = 0; k < population; ++k) irank[k] = k;
	  if (all_feasible) /* special case: rank by objective function */
	       nlopt_qsort_r(irank, population, sizeof(int), fval,key_compare);
	  else {
	       /* Runarsson & Yao's stochastic ranking of the population */
	       for (i = 0; i < population; ++i) {
		    int swapped = 0;
		    for (j = 0; j < population-1; ++j) {
			 double u = nlopt_urand(0,1);
			 if (u < PF || (penalty[irank[j]] == 0
					&& penalty[irank[j+1]] == 0)) {
			      if (fval[irank[j]] > fval[irank[j+1]]) {
				   int irankj = irank[j];
				   irank[j] = irank[j+1];
				   irank[j+1] = irankj;
				   swapped = 1;
			      }
			 }
			 else if (penalty[irank[j]] > penalty[irank[j+1]]) {
			      int irankj = irank[j];
			      irank[j] = irank[j+1];
			      irank[j+1] = irankj;
			      swapped = 1;
			 }
		    }
		    if (!swapped) break;
	       }
	  }
	  
	  /* compute rank from irank */
	  for(k = 0; k < population; ++k){
	    rank[irank[k]] = k;
	  }

	  /* OLD: evolve the population: differential evolution for
	          the best survivors, and standard mutation of the
	          best survivors for the rest: */
	  /* NOW: evolve the population: "survivors" actually survive,
	          differential evolution for the survivors gives some
	          offspring, the rest is generated through standard
	          mutation. */
	  for (k = survivors; k < population; ++k) { /* standard mutation */
	       double taup_rand = taup * nlopt_nrand(0,1);
	       int rk = irank[k], ri;
	       i = k % survivors;
	       ri = irank[i];
	       for (j = 0; j < n; ++j) {
		    double sigmamax = (ub[j] - lb[j]) / sqrt(n);
		    sigmas[rk*n+j] = sigmas[ri*n+j] 
			 * exp(taup_rand + tau*nlopt_nrand(0,1));
		    if (sigmas[rk*n+j] > sigmamax)
			 sigmas[rk*n+j] = sigmamax;
		    do {
			 xs[rk*n+j] = xs[ri*n+j] 
			      + sigmas[rk*n+j] * nlopt_nrand(0,1);
		    } while (xs[rk*n+j] < lb[j] || xs[rk*n+j] > ub[j]);
		    sigmas[rk*n+j] = sigmas[ri*n+j] * 
		      sqrt( sigmas[rk*n+j]/sigmas[ri*n+j] );
/* This gives a heritability of log(sigma) of h^2=0.5^2=0.25, which
 * appears reasonable in view of empirical values of h^2 around 0.3
 * in data such as in
 * http://www.genetics.org/content/130/1/195.full.pdf or
 * http://bio.fsu.edu/~dhoule/Publications/Hansen@@11h2ev.pdf (though
 * this is not heritability of evolvability). */
			}
	  }
	  memcpy(x0, xs+irank[0]*n, n * sizeof(double));
	  for (k = 0; k< survivors; ++k){/* differential variation */
	       double taup_rand = taup * nlopt_nrand(0,1);
	       int rk = irank[k], ri;
	       i = k;
	       ri = irank[i];
	       for (j = 0; j < n; ++j) {
		    double xi = xs[ri*n+j];
		    if (i+1 < survivors)
			 xs[rk*n+j] = 
			    xi + GAMMA * (x0[j] - xs[irank[i+1]*n+j]);
		    if (i+1 == survivors
			|| xs[rk*n+j] < lb[j] || xs[rk*n+j] > ub[j]) {
			 /* standard mutation for last survivor and
			    for any survivor components that are now
			    outside the bounds */
		         double sigmamax = (ub[j] - lb[j]) / sqrt(n);
			 double sigi = sigmas[ri*n+j]; 
			 sigmas[rk*n+j] = sigi 
			   * exp(taup_rand + tau*nlopt_nrand(0,1));
			 if (sigmas[rk*n+j] > sigmamax)
			      sigmas[rk*n+j] = sigmamax;
			 do {
			      xs[rk*n+j] = xi 
				   + sigmas[rk*n+j] * nlopt_nrand(0,1);
			 } while (xs[rk*n+j] < lb[j] || xs[rk*n+j] > ub[j]);
			 sigmas[rk*n+j] = sigi * 
			   sqrt( sigmas[rk*n+j]/sigi );
		    }
	       }
	  }
     }

done:
     if (irank) free(irank);
     if (sigmas) free(sigmas);
     if (results) free(results);
     return ret;
}
