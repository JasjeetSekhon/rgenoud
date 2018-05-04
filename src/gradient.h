/*

  RGENOUD

  Walter R. Mebane, Jr.
  Cornell University
  http://macht.arts.cornell.edu/wrm1
  <wrm1@macht.arts.cornell.edu>

  Jasjeet Singh Sekhon 
  UC Berkeley
  http://sekhon.polisci.berkeley.edu
  <sekhon@berkeley.edu>

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/gradient.h,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/

#define EPS 0.000000000000001   /* machine precision */

void gradient(SEXP fn, SEXP rho,
	      double *p, double *g, long nparms, short MinMax, short BoundaryEnforcement,
	      double **Domains);

double func4g(SEXP fn, SEXP rho,
	      double *X, long nvars, short MinMax, short BoundaryEnforcement, 
	      double **Domains);

void numgrad(SEXP fn, SEXP rho,
	     double *epsacc, double *optint,
	     int nparms, double *invals, double *grads, double *wrk,
	     double (*func)(SEXP fn, SEXP rho,
			     double *X, int nvars, short int MinMax), 
	     short MinMax);

void numgradc(SEXP fn, SEXP rho,
	      double *epsacc, double *optint,	      
	      int nparms, double *invals, double *grads, double *wrk,
	      double (*func)(SEXP fn, SEXP rho,
			     double *X, long nvars, short MinMax, 
			     short BoundaryEnforcement, double **Domains), 
	      short MinMax, short BoundaryEnforcement, double **Domains);

extern double *numopgc(double *epsacc, double *optint,
		int nparms, int nobs, double *invals, double *opg, double *wrk,
		int (*func)(double *X, double *outvec));

double **eaccuracy(SEXP fn, SEXP rho,
		   int nparms, int ndiffs, double h, double *invals,
		   double *wrk, 
		   double (*func)(SEXP fn, SEXP rho,
				  double *X, long nvars, short MinMax, 
				  short BoundaryEnforcement, double **Domains), 
		   short MinMax, short BoundaryEnforcement, double **Domains);

struct estints {
  int nparms;
  int *errors;  /* 0 == OK, >=1 == error */
  double
    *hf,   /* interval */
    *phi,  /* f' (first derivative) */
    *phic,  /* f' (first derivative, central-difference) */
    *phi2, /* f'' (second derivative) */
    *ef,   /* error bound */
    *hessian;  /* hessian matrix (lower triangle) */
};

struct estints *algfd(SEXP fn, SEXP rho,
		      int nparms, double *eps, double *invals, double *wrk,
		      double (*func)(SEXP fn, SEXP rho,
				     double *X, long nvars, short MinMax, 
				     short BoundaryEnforcement, double **Domains), 
		      short MinMax, short BoundaryEnforcement, double **Domains);

void fdestimates(SEXP fn, SEXP rho,
		 int parm, double fvalue, double *invals, double *wrk,
		 double eps, double h,
		 double *fplus, double *fminus,
		 double *phif, double *phib, double *phic, double *phi2,
		 double *cf, double *cb, double *c2,
		 double (*func)(SEXP fn, SEXP rho,
				double *X, long nvars, short MinMax, 
				short BoundaryEnforcement, double **Domains), 
		 int nparms, short MinMax, short BoundaryEnforcement, double **Domains);

struct estints *numhessian(struct estints *instruc, double *invals, double *wrk,
			   double (*func)(double *));

struct estints *numhessianc(SEXP fn, SEXP rho,
			    struct estints *instruc, double *invals, double *wrk,
			    double (*func)(SEXP fn, SEXP rho,
					   double *X, long nvars, short MinMax, 
					   short BoundaryEnforcement, double **Domains), 
			    short MinMax, short BoundaryEnforcement, double **Domains);

void estoptint(SEXP fn, SEXP rho,
	       double *epsacc, double *optint,
	       int nparms, int ndiffs, int pflag, double *invals,
	       double (*func)(SEXP fn, SEXP rho,
			      double *X, long nvars, short MinMax, short BoundaryEnforcement,
			      double **Domains), 
	       short MinMax, short BoundaryEnforcement, double **Domains);

void dohessians(SEXP fn, SEXP rho,
		double *epsacc, 
		int nparms, int nobs, int ndiffs, double *invals,
		double (*func)(SEXP fn, SEXP rho,
			       double *X, long nvars, short MinMax, 
			       short BoundaryEnforcment, double **Domains), 
		double (*funco)(double *, double *),
		short MinMax, short BoundaryEnforcement, double **Domains);
