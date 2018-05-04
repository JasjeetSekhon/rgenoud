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

  June 2, 2012

*/

#include "genoud.h"
#include "gradient.h"

double GammaLN(double xx)
{
	double cof[6] = {
						 76.18009173,
						-86.50532033,
						 24.01409822,
						 -1.231739516,
						  0.12085003E-2,
						 -0.536382E-5,
					};

	double stp  = 2.50662827465;
	double x, tmp, ser;
	int j;

	x = xx - 1.0;
	
	tmp = x + 5.5;
	
	tmp = (x + 0.5) * log(tmp) - tmp;

	ser = 1.0;
	
	for (j = 0; j < 6; j++)
	{
		x	+= 1.0;
		ser += cof[j] / x;
	}

	return tmp + log(stp * ser);
}

double VMgamma(double xx)
{
	#define PI2 3.141592654

	if		(xx > 0) return exp(GammaLN(xx));
	else if (xx < 0) return PI2 / exp(GammaLN(1 - xx)) / sin( PI2 * (1 - xx));
	else			 return 0;

}



/*
  func4g is required by the gradient code.  It takes what is in
  evaluate() and increments it down by 1
  
  NOTE: funco() (which not used here) does not increment evaluate() down by 1.
  The code internal to the function does this.
  */

double func4g(SEXP fn, SEXP rho,
	      double *X, long nvars, short MinMax, short BoundaryEnforcement, 
	      double **Domains)
{
  double fit;
  short BoundaryTrigger=0;
  long i;

  if (BoundaryEnforcement==2) {
    for (i=0; i<nvars; i++) {
      if (X[i] < Domains[(i+1)][1]) {
	BoundaryTrigger=1;
	break;
      }
      if (X[i] > Domains[(i+1)][3]){
	BoundaryTrigger=1;
	break;
      }
    }

    if (BoundaryTrigger > 0) {
      // *Status=-3;
      // min
      if (MinMax==0) return(-1*DOUBLEMAX);
      //max
      else return(DOUBLEMAX);
    }
  }

  if (MinMax==0) fit=evaluate(fn, rho, X-1, nvars, MinMax);
  else fit = -1*evaluate(fn, rho, X-1, nvars, MinMax);
  return(fit);
}

/* replace gradient() body with analytical gradient code, if desired.
   by default, numerical gradients with the intervals in *optint are used

   arguments (must be identical in any replacement code):
   *p, vector of parameters, p[0] is the first parameter;
   *g, vector that will hold the gradient, allocated by caller,
       g[0] is the derivative for the first parameter;
   nparms, the number of parameters, p[nparms-1] is the last parameter.
*/
void gradient(SEXP fn, SEXP rho,
	      double *p, double *g, long nparms, short MinMax, short BoundaryEnforcement,
	      double **Domains)
{

  double *wrk;
  int ndiffs;

  /* formally declared global in graddata.h */
  // double *epsacc, *optint, *ihessians;
  double *epsacc, *optint;

  optint = (double *) malloc(nparms*sizeof(double));
  epsacc = (double *) malloc(nparms*sizeof(double));
  wrk = (double *) malloc(nparms*sizeof(double));

  ndiffs = 9;  /* order of differences for num grad optimal interval calcs */

  estoptint(fn, rho, epsacc, optint, nparms, ndiffs, 2, p, func4g, 
	    MinMax, BoundaryEnforcement, Domains);
  
  /* numgradc:  numerical gradient, central-difference */
  numgradc(fn, rho, epsacc, optint, nparms, p, g, wrk, func4g, MinMax,
	   BoundaryEnforcement, Domains);

  free(wrk);
  free(epsacc);
  free(optint);

  return;
}

/* numerical gradient, forward-difference */
/* see Gill, Murray and Wright, 1981, "Practical Optimization," p. 342 */

/* invals, grads, wrk should point to double w[nparms+1] arrays;
   func is the function whose gradient is to be evaluated */
void numgrad(SEXP fn, SEXP rho,
	     double *epsacc, double *optint,
	     int nparms, double *invals, double *grads, double *wrk,
	     double (*func)(SEXP fn, SEXP rho, double *X, int nvars, short int MinMax), 
	     short MinMax)
{
  int i;
  //  double u, rf, fplus, fminus;
  double u, fplus, fminus;
  double epsuse;
  double duse;
  
  /* evaluate func at the input point */
  u = func(fn, rho, invals, nparms, MinMax);
  
  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];
  
  /* evaluate the gradient */
  for (i=0; i<nparms; i++) {
    epsuse = epsacc[i];
    duse = optint[i];
    wrk[i] += duse;
    grads[i] = (func(fn, rho, wrk, nparms, MinMax) - u) / duse ;

    /* check the gradient */
    if (2.0*epsuse / (duse*fabs(grads[i])) > 0.1) { /* switch to central-diff */
      duse = pow(duse, 2.0/3.0);  /* see GMW p 131 */
      wrk[i] = invals[i] + duse;
      fplus = func(fn, rho, wrk, nparms, MinMax);
      wrk[i] = invals[i] - duse;
      fminus = func(fn, rho, wrk, nparms, MinMax);
      grads[i] = (fplus-fminus) * 0.5 / duse ;
    }
    wrk[i] = invals[i];
  }
}


/* numerical gradient, central-difference */
/* see Gill, Murray and Wright, 1981, "Practical Optimization," p. 342 */

/* invals, grads, wrk should point to double w[nparms+1] arrays;
   func is the function whose gradient is to be evaluated */
void numgradc(SEXP fn, SEXP rho,
	      double *epsacc, double *optint,	      
	      int nparms, double *invals, double *grads, double *wrk,
	      double (*func)(SEXP fn, SEXP rho,
			     double *X, long nvars, short MinMax, 
			     short BoundaryEnforcement, double **Domains), 
	      short MinMax, short BoundaryEnforcement, double **Domains)
{
  int i;
  // double u, rf, fplus, fminus;
  double u, fplus, fminus;
  double epsuse;
  double duse;

  /* evaluate func at the input point */
  u = func(fn, rho, invals, nparms, MinMax, BoundaryEnforcement, Domains);

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  /* evaluate the gradient */
  for (i=0; i<nparms; i++) {

    epsuse = epsacc[i];
    duse = optint[i];
    duse = pow(duse, 2.0/3.0);  /* see GMW p 131 */
    wrk[i] = invals[i] + duse;
    fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    wrk[i] = invals[i] - duse;
    fminus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    grads[i] = (fplus-fminus) * 0.5 / duse ;
    wrk[i] = invals[i];
  }
}

#ifdef NEVERDEFINED
/* numerical outer-product of gradient, central-difference */
/* return pointer to full numerical OPG */
/*  uses central differences: */
/*   phi2 = (f(x+hi ei) - f(x-hi ei)) * (f(x+hj ej) - f(x-hj ej)) */
/* invals, grads, wrk should point to double w[nparms+1] arrays;
   func is the function whose gradient is to be evaluated */
double *numopgc(double *epsacc, double *optint,
		int nparms, int nobs, double *invals, double *opg, double *wrk,
		int (*func)(double *X, double *outvec))
{
  double *outvec, *outplus, *outminus, *outdiff;

  int i,j,k, idx;
  int ni = 0;
  double hi, ih;
  double phi2, dini;

  double *hfuse = optint;

  outvec = (double *) malloc(nobs*sizeof(double));
  outplus = (double *) malloc(nobs*nparms*sizeof(double));
  outminus = (double *) malloc(nobs*nparms*sizeof(double));
  outdiff = (double *) malloc(nobs*nparms*sizeof(double));

  /* allocate storage for the hessian */
  opg = (double *) calloc(nparms*nparms,sizeof(double));

  /* evaluate func at the input point */
  /* parameters, ooutvec=nobs */
  /* evaluates the function (i.e., the likelihood) at each observation */
  func(invals-1, outvec);

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  /* evaluate the gradient elements */
  for (i=0; i<nparms; i++) {
    hi = pow(hfuse[i], 2.0/3.0);
    wrk[i] = invals[i] + hi;
    func(wrk-1, outplus+i*nobs);
    wrk[i] = invals[i] - hi;
    func(wrk-1, outminus+i*nobs);
    wrk[i] = invals[i];

    ih = 0.5 / hi;
    for (k=0; k<nobs; k++) {
      idx = i*nobs+k;
      outdiff[idx] = (outplus[idx] - outminus[idx]) * ih ;
    }
  }
  for (i=0; i<nparms; i++) {
    idx = i*nparms;
    for (j=0; j<=i; j++) {
      phi2 = 0.0;
      for (k=0; k<nobs; k++) {
	phi2 += outdiff[i*nobs+k] * outdiff[j*nobs+k];
      }
      opg[i*nparms + j] = phi2;
      if (i!=j) opg[j*nparms + i] = phi2;
    }
  }
  free(outdiff);
  free(outminus);
  free(outplus);
  free(outvec);

  return opg;
}
#endif

/* estimate accuracy */
/* see Gill, Murray and Wright, 1981, "Practical Optimization," p. 337 */

double **eaccuracy(SEXP fn, SEXP rho,
		   int nparms, int ndiffs, double h, double *invals,
		   double *wrk, 
		   double (*func)(SEXP fn, SEXP rho,
				  double *X, long nvars, short MinMax, 
				  short BoundaryEnforcement, double **Domains), 
		   short MinMax, short BoundaryEnforcement, double **Domains)
{
  double **table;

  int i,j,k, idx;
  int nsteps = 1+2*ndiffs, nrows = nparms*nsteps, ncols = ndiffs+1;
  double u, huse, v;
  double scale = 2.0*pow(10.0,6.0);

  /* allocate storage for the table to differences to be returned */
  table = (double **) malloc(ncols*sizeof(double *));
  for (i=0; i<ncols; i++)
    table[i] = (double *) calloc(nrows, sizeof(double));

  /* evaluate func at the input point */
  u = func(fn, rho, invals, nparms, MinMax, BoundaryEnforcement, Domains);
  for (i=0; i<nparms; i++) table[0][i*nsteps] = u;

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  /* evaluate the offsets */
  for (i=0; i<nparms; i++) {
    /* make sure huse is sufficiently small */
    v = fabs(invals[i]);
    huse = h;
    if (v>EPS*scale) {
      while (huse > v/scale)
	huse *= 0.1;
    }
    for (j=1; j<nsteps; j++) {
      wrk[i] += huse;
      table[0][i*nsteps+j] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains) ;
    }
    wrk[i] = invals[i];
  }

  /* compute the differences */
  for (i=0; i<nparms; i++) {
    idx = i*nsteps;
    for (j=0; j<ndiffs; j++) {
      for (k=0; k<nsteps-j-1; k++) {
	table[j+1][idx+k] = table[j][idx+k+1] - table[j][idx+k];
      }
    }
  }
  return table;
}

/* estimate intervals for use in numerical gradients */
/* see Gill, Murray and Wright, 1981, "Practical Optimization," p. 343-344 */

struct estints *algfd(SEXP fn, SEXP rho,
		      int nparms, double *eps, double *invals, double *wrk,
		      double (*func)(SEXP fn, SEXP rho,
				     double *X, long nvars, short MinMax, 
				     short BoundaryEnforcement, double **Domains), 
		      short MinMax, short BoundaryEnforcement, double **Domains)
{
  struct estints *outstruc;

  // int i,j,k;
  int i, k;
  int K = 20;
  int errval = 0;
  double omega = 1.0, eta = 1.0;
  double u;
  double hf, hk, hbar, hs, hphi=0;
  double fplus, fminus, phi, phif, phib, phic, phi2, cf, cb, c2;
  // double fplus, fminus, phi, phic, phi2, c2;
  double ef, ebar;

  /* allocate structure to return */
  outstruc = (struct estints *) malloc(sizeof(struct estints));
  outstruc->errors = (int *) calloc(nparms, sizeof(int));
  outstruc->hf = (double *) calloc(nparms, sizeof(double));
  outstruc->phi = (double *) calloc(nparms, sizeof(double));
  outstruc->phic = (double *) calloc(nparms, sizeof(double));
  outstruc->phi2 = (double *) calloc(nparms, sizeof(double));
  outstruc->ef = (double *) calloc(nparms, sizeof(double));
  outstruc->nparms = nparms;

  /* evaluate func at the input point */
  u = func(fn, rho, invals, nparms, MinMax, BoundaryEnforcement, Domains);

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  for (i=0; i<nparms; i++) {
    // #ifdef NEVERDEFINED
  FD1:
    hbar = 2.0*(eta+fabs(invals[i]))* sqrt(eps[i]/(omega + fabs(u))) ;
    hk = 10.0 * hbar;
    k = 0;
    fdestimates(fn, rho, i, u, invals, wrk, eps[i], hk,
		&fplus, &fminus, &phif, &phib, &phic, &phi2, &cf, &cb, &c2,
		func, nparms, MinMax, BoundaryEnforcement, Domains);
    hs = -1.0;
  FD2:
    if ((cf>cb ? cf : cb) <= 0.1) hs = hk;
    if (0.001 <= c2 && c2 <= 0.1) {
      hphi = hk;
      goto FD5;
    }
    if (0.001 > c2) goto FD4;
  FD3:
    do {
      k++;
      hk *= 10.0;
      fdestimates(fn, rho, i, u, invals, wrk, eps[i], hk,
		  &fplus, &fminus, &phif, &phib, &phic, &phi2, &cf, &cb, &c2,
		  func, nparms, MinMax, BoundaryEnforcement, Domains);
      if (hs<0 && ((cf>cb ? cf : cb) <= 0.1)) hs = hk;
      if (c2 <= 0.1) {
	hphi = hk;
	goto FD5;
      }
    } while (k<K);
    if (k==K) goto FD6;
  FD4:
    do {
      k++;
      hk /= 10.0;
      fdestimates(fn, rho, i, u, invals, wrk, eps[i], hk,
		  &fplus, &fminus, &phif, &phib, &phic, &phi2, &cf, &cb, &c2,
		  func, nparms, MinMax, BoundaryEnforcement, Domains);
      if (c2 > 0.1) {
	hphi = hk * 10.0;
	goto FD5;
      }
      if ((cf>cb ? cf : cb) <= 0.1) hs = hk;
      if (0.001 <= c2 && c2 <= 0.1) {
	hphi = hk;
	goto FD5;
      }
    } while (k<K);
    if (k==K) goto FD6;
  FD5:
    hf = 2.0*sqrt(eps[i]/fabs(phi2));
    wrk[i] = invals[i] + hf;
    fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);
    phi = (fplus-u)/hf;
    wrk[i] = invals[i] + hphi;
    fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    wrk[i] = invals[i] - hphi;
    fminus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    wrk[i] = invals[i];
    phic = (fplus-fminus)/(2.0*hphi);

    ef = hf*fabs(phi2)*0.5 + 2.0*eps[i]/hf ;
    ebar = fabs(phi-phic);

    outstruc->hf[i] = hf;
    outstruc->phi[i] = phi;
    outstruc->phic[i] = phic;
    outstruc->phi2[i] = phi2;
    outstruc->ef[i] = ef;
    if ((ef>ebar ? ef : ebar) <= 0.5*fabs(phi)) {
      outstruc->errors[i] = 0;
    }
    else
      outstruc->errors[i] = 1;
    continue;
  FD6:
    if (hs<0) {
      hf = hbar;
      phi = phi2 = ef = 0.0;
      errval = 2;
    }
    else if (hs > 0 && c2 > 0.1) {
      hf = hs;
      wrk[i] = invals[i] + hf;
      fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[i] = invals[i];
      phi = (fplus-u)/hf;
      phi2 = 0.0;
      ef = 2.0*eps[i]/hf ;
      errval = 3;
    }
    else {
      hf = hk;
      wrk[i] = invals[i] + hf;
      fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      phi = (fplus-u)/hf;
      wrk[i] = invals[i] - hf;
      fminus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[i] = invals[i];
      phic = (fplus-fminus)/(2.0*hf);
      ef = hf*fabs(phi2)*0.5 + 2.0*eps[i]/hf ;
      errval = 4;
    }
    outstruc->hf[i] = hf;
    outstruc->phi[i] = phi;
    outstruc->phic[i] = phic;
    outstruc->phi2[i] = phi2;
    outstruc->ef[i] = ef;
    outstruc->errors[i] = errval;
    // #endif
  }
  return outstruc;
}

void fdestimates(SEXP fn, SEXP rho,
		 int parm, double fvalue, double *invals, double *wrk,
		 double eps, double h,
		 double *fplus, double *fminus,
		 double *phif, double *phib, double *phic, double *phi2,
		 double *cf, double *cb, double *c2,
		 double (*func)(SEXP fn, SEXP rho,
				double *X, long nvars, short MinMax, 
				short BoundaryEnforcement, double **Domains), 
		 int nparms, short MinMax, short BoundaryEnforcement, double **Domains)
{
  double ih = 1.0/h;

  wrk[parm] = invals[parm] + h;
  *fplus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

  wrk[parm] = invals[parm] - h;
  *fminus = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

  wrk[parm] = invals[parm];
  *phif = (*fplus-fvalue) * ih;
  *phib = (fvalue-*fminus) * ih;
  *phic = (*fplus-*fminus) * 0.5 * ih;
  *phi2 = (*phif-*phib) * ih;
  *cf = 2.0*eps*ih/fabs(*phif);
  *cb = 2.0*eps*ih/fabs(*phib);
  *c2 = 4.0*eps*ih*ih/fabs(*phi2);
}

/* put strict lower triangle of numerical hessian into instruc->hessian */
/* instruc should have been set by algfd */
struct estints *numhessian(struct estints *instruc, double *invals, double *wrk,
			   double (*func)(double *))
{
  int nparms;
  double *fplusi = NULL;

  int i,j;
  double hi, hj, ih, jh;
  double fvalue, fplus, phi2;

  nparms = instruc->nparms;
  fplusi = (double *) malloc(nparms*sizeof(double));


  /* allocate storage for the hessian */
  instruc->hessian = (double *) calloc((nparms*(nparms+1))/2,sizeof(double));

  /* evaluate func at the input point */
  fvalue = func(invals);

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  for (i=0; i<nparms; i++) {
    wrk[i] = invals[i] + instruc->hf[i];
    fplusi[i] = func(wrk);
    wrk[i] = invals[i];
  }
  for (i=1; i<nparms; i++) {
    hi = instruc->hf[i];
    ih = 1.0/hi;
    wrk[i] = invals[i] + hi;
    for (j=0; j<i; j++) {
      hj = instruc->hf[j];
      jh = 1.0/hj;
      wrk[j] = invals[j] + hj;
      fplus = func(wrk);
      wrk[j] = invals[j];
      phi2 = (fplus - fplusi[i] - fplusi[j] + fvalue) * ih * jh;
      instruc->hessian[(i*(i-1))/2 + j] = phi2;
    }
    wrk[i] = invals[i];
  }

  free(fplusi);
  return instruc;
}

/* put strict lower triangle of numerical hessian into instruc->hessian */
/*  uses central differences: */
/*   phi2 =
       f(x+hj ej+hi ei) - f(x+hj ej-hi ei) - f(x-hj ej+hi ei) + f(x-hj ej-hi ei) */
/* instruc should have been set by algfd */
struct estints *numhessianc(SEXP fn, SEXP rho,
			    struct estints *instruc, double *invals, double *wrk,
			    double (*func)(SEXP fn, SEXP rho,
					   double *X, long nvars, short MinMax, 
					   short BoundaryEnforcement, double **Domains), 
			    short MinMax, short BoundaryEnforcement, double **Domains)
{
  int nparms;
  int nelems;
  double *fplusi;
  double *fminusi;
  double *fppi;
  double *fpmi;
  double *fmmi;

  int i,j, idx;
  double hi, hj, ih, jh;
  // double fvalue, fplus, fminus, phi2;
  double fvalue, phi2;

  nparms = instruc->nparms;
  nelems = (nparms*(nparms-1))/2;
  fplusi = (double *) malloc(nparms*sizeof(double));
  fminusi = (double *) malloc(nparms*sizeof(double));
  fppi = (double *) malloc(nelems*sizeof(double));
  fpmi = (double *) malloc(nparms*nparms*sizeof(double));
  fmmi = (double *) malloc(nelems*sizeof(double));

  /* allocate storage for the hessian */
  instruc->hessian = (double *) calloc(nelems,sizeof(double));

  /* evaluate func at the input point */
  fvalue = func(fn, rho, invals, nparms, MinMax, BoundaryEnforcement, Domains);

  /* copy the parameter values for point at which to evaluate the gradient */
  for (i=0; i<nparms; i++) wrk[i] = invals[i];

  for (i=0; i<nparms; i++) {
    hi = pow(instruc->hf[i], 2.0/3.0);
    idx = (i*(i-1))/2;
    wrk[i] = invals[i] + 2.0*hi;
    fplusi[i] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    wrk[i] = invals[i] - 2.0*hi;
    fminusi[i] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

    for (j=0; j<i; j++) {
      hj = pow(instruc->hf[j], 2.0/3.0);
      wrk[i] = invals[i] + hi;
      wrk[j] = invals[j] + hj;
      fppi[idx + j] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[i] = invals[i] + hi;
      wrk[j] = invals[j] - hj;
      fpmi[i*nparms + j] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[i] = invals[i] - hi;
      wrk[j] = invals[j] + hj;
      fpmi[j*nparms + i] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[i] = invals[i] - hi;
      wrk[j] = invals[j] - hj;
      fmmi[idx + j] = func(fn, rho, wrk, nparms, MinMax, BoundaryEnforcement, Domains);

      wrk[j] = invals[j];
    }
    wrk[i] = invals[i];
  }
  for (i=0; i<nparms; i++) {
    hi = pow(instruc->hf[i], 2.0/3.0);
    ih = 1.0/hi;
    idx = (i*(i-1))/2;
    phi2 = (fplusi[i] - 2.0*fvalue + fminusi[i]) * ih * ih * 0.25;
    instruc->phi2[i] = phi2;
    for (j=0; j<i; j++) {
      hj = pow(instruc->hf[j], 2.0/3.0);
      jh = 1.0/hj;
      phi2 = (fppi[idx+j] - fpmi[j*nparms+i] - fpmi[i*nparms+j] + fmmi[idx+j])
	* ih * jh * 0.25;
      instruc->hessian[idx + j] = phi2;
    }
  }

  free(fmmi);
  free(fpmi);
  free(fppi);
  free(fplusi);

  return instruc;
}

void estoptint(SEXP fn, SEXP rho,
	       double *epsacc, double *optint,
	       int nparms, int ndiffs, int pflag, double *invals,
	       double (*func)(SEXP fn, SEXP rho,
			      double *X, long nvars, short MinMax, short BoundaryEnforcement,
			      double **Domains), 
	       short MinMax, short BoundaryEnforcement, double **Domains) 
{
  double *wrk;

  int i,j,k;
  // int nsteps=1+2*ndiffs, nrows=nparms*nsteps, ncols=ndiffs+1;
  int nsteps=1+2*ndiffs;
  double h, beta, dwrk;
  double **table;
  struct estints *estructure;

  wrk = (double *) malloc(nparms*(ndiffs+1)*sizeof(double));

  h = 0.0000002;

  table = eaccuracy(fn, rho, nparms, ndiffs, h, invals, wrk, func, MinMax,
		    BoundaryEnforcement, Domains);

  for (i=0; i<nparms*ndiffs; i++) wrk[i] = 0.0;

  for (i=0; i<nparms; i++) {
    for (j=0; j<ndiffs; j++) {
      for (k=1; k<=ndiffs; k++) {
	/* find largest difference of each order for each parm */
	dwrk = fabs(table[j+1][i*nsteps+k]);
	if (wrk[i*ndiffs+j] < dwrk) wrk[i*ndiffs+j] = dwrk;
      }
      beta = sqrt(VMgamma(1.0+2.0*(1.0+(double)j))/pow(VMgamma(2.0+(double)j),2.0));
      wrk[i*ndiffs+j] /= beta;
    }
  }
  /* put estimates for highest difference order into epsacc */
  for (i=0; i<nparms; i++) {
    dwrk = wrk[i*ndiffs+ndiffs-1];
    /* make sure epsacc values >= EPS */
    epsacc[i] = (dwrk > EPS ? dwrk : EPS) ; 
  }

#ifdef NEVERDEFINED
  Rprintf("accuracy estimates:\n");
  for (i=0; i<nparms; i++) {
    Rprintf("parm = %d\n", i+1);
    for (j=0; j<ndiffs; j++) {
      Rprintf(" %14.7e", wrk[i*ndiffs+j]);
    }
    Rprintf("\n");
  }

  Rprintf("difference table:\n");
  for (i=0; i<nparms; i++) {
    Rprintf("parm = %d\n", i+1);
    for (j=0; j<=ndiffs; j++) {
      for (k=0; k<=ndiffs; k++) {
	Rprintf(" %14.7e", table[j][i*nsteps+k]);
      }
      Rprintf("\n");
    }
  }
#endif

  estructure = algfd(fn, rho, nparms, epsacc, invals, wrk, func, MinMax, 
		     BoundaryEnforcement, Domains);

  if (pflag==1) {
    Rprintf("err   interval          f'                fc'               f''               errorbound\n");
    for (i=0; i<nparms; i++) {
      Rprintf(" %d  ", estructure->errors[i]);
      Rprintf(" %17.10e", estructure->hf[i]);
      Rprintf(" %17.10e", estructure->phi[i]);
      Rprintf(" %17.10e", estructure->phic[i]);
      Rprintf(" %17.10e", estructure->phi2[i]);
      Rprintf(" %17.10e", estructure->ef[i]);
      Rprintf("\n");
    }
  }

  /* put estimates for optimal interval into optint */
  for (i=0; i<nparms; i++) {
    optint[i] = estructure->hf[i];
  }
  
  free(table);
  free(wrk);
  /* free the estructure */
  free(estructure->errors);
  free(estructure->hf);
  free(estructure->phi);
  free(estructure->phic);
  free(estructure->phi2);
  free(estructure->ef);
  free(estructure);
}

void dohessians(SEXP fn, SEXP rho,
		double *epsacc, 
		int nparms, int nobs, int ndiffs, double *invals,
		double (*func)(SEXP fn, SEXP rho,
			       double *X, long nvars, short MinMax, 
			       short BoundaryEnforcment, double **Domains), 
		double (*funco)(double *, double *),
		short MinMax, short BoundaryEnforcement, double **Domains)
{
  double *wrk;

  int i,j;
  // int nsteps=1+2*ndiffs, nrows=nparms*nsteps, ncols=ndiffs+1;
  // int nsteps=1+2*ndiffs;
  struct estints *estructure;

  // double *opg;

  wrk = (double *) malloc(nparms*(ndiffs+1)*sizeof(double));

  estructure = algfd(fn, rho, nparms, epsacc, invals, wrk, func, MinMax, 
		     BoundaryEnforcement, Domains);

#ifdef NEVERDEFINED
  /* numerical hessian, using forward differences for off-diagonal elements */
  numhessian(estructure, invals, wrk, func);
  Rprintf("numerical hessian, forward differences:\n");
  for (i=0; i<nparms; i++) {
    for (j=0; j<nparms; j++) {
      if (i==j)
	Rprintf(" %19.12e", estructure->phi2[i] / 2.0);
      else if (i>j)
	Rprintf(" %19.12e", estructure->hessian[(i*(i-1))/2 + j] / 2.0);
      else if (j>i)
	Rprintf(" %19.12e", estructure->hessian[(j*(j-1))/2 + i] / 2.0);
    }
    Rprintf("\n");
  }
  /* fflush(stdout); */
#endif

  /* numerical hessian, using central differences */
  numhessianc(fn, rho, estructure, invals, wrk, func, MinMax, BoundaryEnforcement, 
	      Domains);

  Rprintf("numerical hessian, central differences:\n");
  for (i=0; i<nparms; i++) {
    for (j=0; j<nparms; j++) {
      if (i==j)
	Rprintf(" %19.12e", estructure->phi2[i] / 2.0);
      else if (i>j)
	Rprintf(" %19.12e", estructure->hessian[(i*(i-1))/2 + j] / 2.0);
      else if (j>i)
	Rprintf(" %19.12e", estructure->hessian[(j*(j-1))/2 + i] / 2.0);
    }
    Rprintf("\n");
  }
  /* fflush(stdout); */

#ifdef NEVERDEFINED
  /* numerical outer product of gradients, using central differences */
  if (funco != NULL) {
    opg = numopgc(nparms, nobs, invals, opg, wrk, funco);
    Rprintf("numerical outer product of gradients, central differences:\n");
    for (i=0; i<nparms; i++) {
      for (j=0; j<nparms; j++)
	Rprintf(" %19.12e", opg[i*nparms+j] / 2.0);
      Rprintf("\n");
    }
    /* fflush(stdout); */
  }
#endif 

  free(wrk);
}
