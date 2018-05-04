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

  June 3, 2012

*/

#include "genoud.h"

void add(double *in1, double *in2, double
         *out, int row, int col) 
{ int i, j, idx;

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
        idx=M(i,j,col);
        out[idx]=in1[idx]+in2[idx];
      }
  }
}


void copy(double *in, double *target, int row, int col)
{
  int i, j, idx;

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
        idx=M(i,j,col);
        target[idx]=in[idx];
      }
  }
}


void multi(double *in1, double *in2, double *out,
	   int row1, int col1, int row2, int col2, int outrowcol[2],
	   FILE *output)
{
  int oi, oj, i;

if (col1!=row2) {
  error("The matrices are not conformable for muliplication\n");
}

outrowcol[0]=row1;
outrowcol[1]=col2;

for (oi=0;oi<outrowcol[0];oi++) {
  for (oj=0;oj<outrowcol[1];oj++) {
    out[M(oi,oj,outrowcol[1])] = 0.0;
  }
}

for (oi=0;oi<outrowcol[0];oi++) {
   for (oj=0;oj<outrowcol[1];oj++) {
       for (i=0;i<col1;i++) {
                out[M(oi,oj,outrowcol[1])] += in1[M(oi,i,col1)]*in2[M(i,oj,col2)];

        }
    }
}

} /* end of multi.c */



void scalarmulti(double scalar, double *in1, double *out, int row, int col) 
{
  int i, j, idx;

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
        idx=M(i,j,col);
        out[idx]=scalar*in1[idx];
      }
  }
}

void scalarmultioffdiag(double scalar, double *in1, double *out, int row, int col) 
{
  int i, j, idx;

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
      idx=M(i,j,col);
      if (i==j)
	out[idx] = in1[idx];
      else
	out[idx] = scalar*in1[idx];
    }
  }
}

void subtract(double *in1, double *in2,
         double *out, int row, int col) { int i, j, idx;

  for (i=0;i<row;i++) {
    for (j=0;j<col;j++) {
        idx=M(i,j,col);
        out[idx]=in1[idx]-in2[idx];
      }
  }
}

double trace(double *a, int n)
{

int i;
double  out;
out=0.0;

 for (i=0;i<n;i++) {
     out += a[M(i,i,n)];
}

 return out;
}


void transpose(double *orig_matrix,
	       double *t_matrix,
	       int orig_rows, int orig_columns) 
{
  int i, j, idx, t_idx;
  
  
  for(i=0;i<orig_rows;i++) {
    for(j=0;j<orig_columns;j++) {
      idx = M(i,j,orig_columns);
      t_idx = M(j,i,orig_rows);
      t_matrix[t_idx]=orig_matrix[idx];
    }
  }
  
} /* end transpose */


void copy_matrix(MATRIX mat1, MATRIX mat2, int lr, int ur, int lc, int uc)
     /*
       int lr,ur,lc,uc;
       MATRIX mat1,mat2;
     */
{
  int i,j;

  for(i=lr; i<=ur; i++)
    for(j=lc; j<=uc; j++)
      mat2[i][j] = mat1[i][j];
}


/* this rounds a double to the closest integer */
int Iround(double in)
{
  double ip, frac;

  frac = modf(in, &ip);
  if (frac < 0.5) return( (int) ip);
  else return( (int) ip+ (int) 1.0);
}


/* This function computs some sample statistics */
void samplestats(double **obsdata, int numobsv, int novarsv, int weightflag, 
		 double *weightdata, FILE *output)
{
  double *mean, *var, *skew, *kur;
  double *rmu, *rvar, *rskew, *rkur;

  int i, j, maxnovars, nobs, nvars;
  double sum[4], x1, x2, wsum, iwsum;
  double dinobs;

  maxnovars = novarsv;

  mean = (double *) malloc((size_t) maxnovars*sizeof(double));
  var = (double *) malloc((size_t) maxnovars*sizeof(double));
  skew = (double *) malloc((size_t) maxnovars*sizeof(double));
  kur = (double *) malloc((size_t) maxnovars*sizeof(double));
  
  rmu = (double *) malloc((size_t) maxnovars*sizeof(double));
  rvar = (double *) malloc((size_t) maxnovars*sizeof(double));
  rskew = (double *) malloc((size_t) maxnovars*sizeof(double));
  rkur = (double *) malloc((size_t) maxnovars*sizeof(double));

  nobs = numobsv;
  nvars = novarsv;
  dinobs = 1.0 / nobs;

  if (weightflag==0) {
    for (j=0; j<nvars; j++) {
      sum[0] = 0.0;
      for (i=0; i<nobs; i++) {
	sum[0] += obsdata[i][j];
      }
      sum[0] *= dinobs;
      sum[1] = 0.0; sum[2] = 0.0; sum[3] = 0.0;
      for (i=0; i<nobs; i++) {
	x1 = obsdata[i][j] - sum[0];
	x2 = x1*x1;
	sum[1] += x2;
	x2 *= x1;
	sum[2] += x2;
	x2 *= x1;
	sum[3] += x2;
      }
      rmu[j] = sum[0];
      rvar[j] = sum[1] * dinobs;
      rskew[j] = sum[2] * dinobs;
      rkur[j] = sum[3] * dinobs;
    }
    
    for (j=0; j<nvars; j++) {
      mean[j] = rmu[j];
      x1 = rvar[j] ;
      var[j] = x1;
      x2 = 1.0 / (x1*x1) ;
      kur[j] = rkur[j] * x2;
      skew[j] = rskew[j] * sqrt(x2/x1);
      
      Rprintf("var %d:\n", j+1);
      Rprintf("sample mean = %f\n", mean[j]);
      Rprintf("sample var = %f\n", var[j]);
      Rprintf("sample skewness = %f\n", skew[j]);
      Rprintf("sample kurtosis = %f\n", kur[j]);
    }
  }
  else if (weightflag==1) {
    for (wsum = 0.0, i=0; i<nobs; i++) wsum += weightdata[i];
    iwsum = 1.0/wsum;
    
    for (j=0; j<nvars; j++) {
      sum[0] = 0.0;
      for (i=0; i<nobs; i++) {
	sum[0] += obsdata[i][j] * weightdata[i];
      }
      sum[0] *= iwsum;
      sum[1] = 0.0; sum[2] = 0.0; sum[3] = 0.0;
      for (i=0; i<nobs; i++) {
	x1 = obsdata[i][j] - sum[0];
	x2 = x1*x1;
	sum[1] += x2 * weightdata[i];
	x2 *= x1;
	sum[2] += x2 * weightdata[i];
	x2 *= x1;
	sum[3] += x2 * weightdata[i];
      }
      rmu[j] = sum[0];
      rvar[j] = sum[1] * iwsum;
      rskew[j] = sum[2] * iwsum;
      rkur[j] = sum[3] * iwsum;
    }
    
    for (j=0; j<nvars; j++) {
      mean[j] = rmu[j];
      x1 = rvar[j] ;
      var[j] = x1;
      x2 = 1.0 / (x1*x1) ;
      kur[j] = rkur[j] * x2;
	skew[j] = rskew[j] * sqrt(x2/x1);
	
	Rprintf("var %d:\n", j+1);
	Rprintf("sample mean = %f\n", mean[j]);
	Rprintf("sample var = %f\n", var[j]);
	Rprintf("sample skewness = %f\n", skew[j]);
	Rprintf("sample kurtosis = %f\n", kur[j]);
      }
  }

  /* free data and work storage */
  free(rkur);
  free(rskew);
  free(rvar);
  free(rmu);
  free(kur);
  free(skew);
  free(var);
  free(mean);
}


/* This function computs some sample statistics for the population matrix */
void populationstats(double **popdata, int numobsv, int novarsv, 
		     double *mean, double *var, double *skew, double *kur,
		     long *tobs)
{
  double *rvar, *rskew, *rkur;

  long i, j, maxnovars, nobs, nvars;
  double sum[4], x1, x2;
  double dinobs;

  maxnovars = novarsv+1;

  rvar = (double *) malloc(maxnovars*sizeof(double));
  rskew = (double *) malloc(maxnovars*sizeof(double));
  rkur = (double *) malloc(maxnovars*sizeof(double));

  nobs = numobsv;
  nvars = novarsv;
  dinobs = 1.0 / nobs;

  for (j=0; j<=nvars; j++) {
    tobs[j] = nobs;
    sum[0] = 0.0;
    for (i=1; i<=nobs; i++) {
      if (popdata[i][j] > DOUBLEMAX) {
	tobs[j]=tobs[j]-1;
	}
      if (popdata[i][j] < -1*DOUBLEMAX) {
	tobs[j]=tobs[j]-1;
      }
      else sum[0] += popdata[i][j];
    }
    dinobs = 1.0/tobs[j]; 
    sum[0] *= dinobs; 
    sum[1] = 0.0; sum[2] = 0.0; sum[3] = 0.0;
    for (i=1; i<=nobs; i++) {
      if ( (popdata[i][j] < DOUBLEMAX) && (popdata[i][j] > -1*DOUBLEMAX) )
	{
	  x1 = popdata[i][j] - sum[0];
	  x2 = x1*x1;
	  sum[1] += x2;
	  x2 *= x1;
	  sum[2] += x2;
	  x2 *= x1;
	  sum[3] += x2;
	} /* end of if */
    } /* end of i loop */
      mean[j] = sum[0];
      rvar[j] = sum[1] * dinobs;
      rskew[j] = sum[2] * dinobs;
      rkur[j] = sum[3] * dinobs;
  } /* end of j loop */
  
  for (j=0; j<=nvars; j++) {
    x1 = rvar[j] ;
    var[j] = x1;
    x2 = 1.0 / (x1*x1) ;
    kur[j] = rkur[j] * x2;
    skew[j] = rskew[j] * sqrt(x2/x1);
  }

  /* free data and work storage */
  free(rkur); 
  free(rskew);
  free(rvar);
}
