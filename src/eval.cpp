/*

  RGENOUD

  Walter R. Mebane, Jr.
  University of Michigan
  http://www-personal.umich.edu/~wmebane
  <wmebane@umich.edu>

  Jasjeet Singh Sekhon 
  UC Berkeley
  http://sekhon.berkeley.edu
  <sekhon@berkeley.edu>

  August 3, 2009

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "genoud.h"

extern "C" 
{
  double evaluate(SEXP fn, SEXP rho, double *X, long nvars, short int MinMax)
  {
    SEXP R_fcall, Rx;
    double fit;
    long i;
    int isFinite=0;
    
    PROTECT(Rx = allocVector(REALSXP, nvars));  
    
    for (i=0; i<nvars; i++)
      {
	REAL(Rx)[i] = X[(i+1)];
      }
    
    PROTECT(R_fcall = lang2(fn, R_NilValue));
    SETCADR(R_fcall, Rx);
    fit =  REAL(eval(R_fcall, rho))[0];
    UNPROTECT(2);
    
    isFinite = R_finite(fit);  
    if (!isFinite)
      {
	if (MinMax)
	  {
	    return(-1*DOUBLEMAX);
	  }
	else
	  {
	    return(DOUBLEMAX);
	  }
      }  
    
 
   return(fit);
  } /* double evaluate(SEXP fn, VECTOR X, long nvars, short int MinMax) */


  void EvaluateLexical(SEXP fn, SEXP rho,
		       double *X, long nvars, long lexical, short int MinMax, double *ret)
  {
    SEXP R_fcall, Rx, Rret;
    long i;
    int isFinite=0;
    
    PROTECT(Rx = allocVector(REALSXP, nvars));

    for (i=0; i<nvars; i++)
      {
	REAL(Rx)[i] = X[i+1];
      }

    PROTECT(R_fcall = lang2(fn, R_NilValue));
    SETCADR(R_fcall, Rx);
    Rret = eval(R_fcall, rho);

    for (i=0; i<lexical; i++)
      {
        ret[i] =  REAL(Rret)[i];
	isFinite = R_finite(ret[i]);  
	if (!isFinite)
	  {
	    if (MinMax)
	      {
		ret[i]=(-1*DOUBLEMAX);
	      }
	    else
	      {
		ret[i]=(DOUBLEMAX);
	      }
	  }  
      }
    UNPROTECT(2);
  } /*   void EvaluateLexical(SEXP fn, SEXP rho */
  

  void EvaluateTransform(SEXP fn, SEXP rho,
                       double *X, long nvars, long lexical, short int MinMax, double *ret)
  {
    SEXP R_fcall, Rx, Rret;
    long i;
    int isFinite=0;

    PROTECT(Rx = allocVector(REALSXP, nvars));

    for (i=0; i<nvars; i++)
      {
        REAL(Rx)[i] = X[i+1];
      }

    PROTECT(R_fcall = lang2(fn, R_NilValue));
    SETCADR(R_fcall, Rx);
    Rret = eval(R_fcall, rho);

    for (i=0; i<lexical; i++)
      {
        ret[i] =  REAL(Rret)[i];
        isFinite = R_finite(ret[i]);
        if (!isFinite)
          {
            if (MinMax)
              {
                ret[i]=(-1*DOUBLEMAX);
              }
            else
              {
                ret[i]=(DOUBLEMAX);
              }
          }
      }
    for(i=0; i<nvars; i++) // this is the only part that differs from EvaluateLexical
      {
        X[(i+1)] = REAL(Rret)[(i+lexical)];
      }
    UNPROTECT(2);
  } /*   void EvaluateTransform(SEXP fn, SEXP rho */

  void userGradientfn(SEXP fnGR, SEXP rho, double *parms, double *grad, long nvars)
  {
    SEXP Rparms, R_fcall, Rgrad;
    long i;
    
    PROTECT(Rparms = allocVector(REALSXP, nvars));    
    PROTECT(Rgrad  = allocVector(REALSXP, nvars));    
    
    for(i=0; i<nvars; i++)
      {
	REAL(Rparms)[i] = parms[i];
      }

    PROTECT(R_fcall = lang2(fnGR, R_NilValue));
    SETCADR(R_fcall, Rparms);    
    Rgrad = eval(R_fcall, rho);  
    
    for(i=0; i<nvars; i++)
      {
	grad[i] = REAL(Rgrad)[i];
      }
    
    UNPROTECT(3);    
  } /*   void userGradientfn(SEXP fnGR, SEXP rho, double *parms, double *grad, long nvars) */


  void RlexicalSort(SEXP fnLexicalSort, SEXP rho,
		    double **population, 
		    short int MinMax, long pop_size, long nvars, long lexical_end,
		    short int type)
  {
    SEXP parms, MAT, R_fcall, MATret;
    long i,j,k=0;

    /* MinMax: 0 min, 1 max */
    /* parms = (1) MinMax, (2) nvars, (3) lexical_end, (4) [nvars/or lexical sort] */
    /* using: #define M(ROW,COL,NCOLS) (((ROW)*(NCOLS))+(COL)) */
    
    PROTECT(MAT = allocMatrix(REALSXP, pop_size, lexical_end));
    PROTECT(parms = allocVector(REALSXP, 4));
    
    REAL(parms)[0] = MinMax;
    REAL(parms)[1] = nvars;
    REAL(parms)[2] = lexical_end;
    REAL(parms)[3] = type; /* 0=nvars, 1=lexical on obj function */

    for(j=0; j<lexical_end; j++)
      for (i=1; i<=pop_size; i++)
	{
	  {
	    REAL(MAT)[k] = population[i][j];
	    k++;
	  }
	}  
    
    PROTECT(R_fcall = lang3(fnLexicalSort, MAT, parms));
    SETCADR(R_fcall, parms);
    SETCADR(R_fcall, MAT);
    MATret = eval(R_fcall, rho);  

    k = 0;
    for(j=0; j<lexical_end; j++)
      for (i=1; i<=pop_size; i++)
	{
	  {
	    population[i][j] = REAL(MATret)[k];
	    k++;
	  }
	}  
    UNPROTECT(3);
  }

  long RmemoryMatrixEvaluate(SEXP fnMemoryMatrixEvaluate, SEXP rho,
			     double **Memory, double **population, 
			     short int MinMax, long pop_size, long UniqueCount,
			     long nvars, long lexical, long lexical_end)
  {
    SEXP parms, Rmemory, Rpopulation, R_fcall, Rret;
    long i,j,k;    

    /* MinMax: 0 min, 1 max */
    /* parms = (1) MinMax, (2) UniqueCount, (3) nvars, (4) lexical */
    
    PROTECT(Rmemory = allocMatrix(REALSXP, UniqueCount, lexical_end));
    PROTECT(Rpopulation = allocMatrix(REALSXP, pop_size, lexical_end));
    PROTECT(parms = allocVector(REALSXP, 3));

    REAL(parms)[0] = MinMax;
    REAL(parms)[1] = nvars;
    REAL(parms)[2] = lexical;

    if(UniqueCount > 1)
      {
	k=0;
	for(j=0; j<lexical_end; j++)
	  for (i=1; i<=UniqueCount; i++)
	    {
	      {
		REAL(Rmemory)[k] = Memory[i][j];
		k++;
	      }
	    }  	
      }

    k =0;
    for(j=0; j<lexical_end; j++)
      for (i=1; i<=pop_size; i++)
	{
	  {
	    REAL(Rpopulation)[k] = population[i][j];
	    k++;
	  }
	}  

    PROTECT(R_fcall = lang4(fnMemoryMatrixEvaluate, Rmemory, Rpopulation, parms));    
    SETCADR(R_fcall, parms);
    SETCADR(R_fcall, Rpopulation);
    SETCADR(R_fcall, Rmemory);
    Rret = eval(R_fcall, rho);      

    UniqueCount = (long) REAL(Rret)[0];
    k =1;
    for(j=0; j<lexical_end; j++)
      for (i=1; i<=UniqueCount; i++)
	{
	  {
	    Memory[i][j] = REAL(Rret)[k];
	    k++;
	  }
	}      

    for(j=0; j<lexical_end; j++)
      for (i=1; i<=pop_size; i++)
	{
	  {
	    population[i][j] = REAL(Rret)[k];
	    k++;
	  }
	}  

    UNPROTECT(4);
    return(UniqueCount);
  }
} /* end of extern "C" */

