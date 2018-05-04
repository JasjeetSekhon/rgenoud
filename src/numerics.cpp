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

double **JaMatrixAllocate(long nobs, long nvars)
{
  long i;
  double **M;

  M= (double **) malloc(nobs*sizeof(double *));
  for (i=0; i<nobs; i++) {
    M[i] = (double *) malloc(nvars*sizeof(double));
  }
  return(M);
}

void JaMatrixFree(double **M, long nobs)
{
  long i;

  if (M == NULL)
    return;
  else
    {
      for (i=0; i<nobs; i++) {
	free(M[i]);
      }
    }
  free(M);
}

short **JaShortMatrixAllocate(long nobs, long nvars)
{
  long i;
  short **M;

  M= (short **) malloc(nobs*sizeof(short *));
  for (i=0; i<nobs; i++) {
    M[i] = (short *) malloc(nvars*sizeof(short));
  }
  return(M);
}

void JaShortMatrixFree(double **M, long nobs)
{
  long i;

  if (M == NULL)
    return;
  else
    {
      for (i=0; i<nobs; i++) {
	free( (short *) M[i]);
      }
    }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   nrerror()                                    */
/*                                                                              */
/*           SYNOPSIS          :   void nrerror(error_text)                     */
/*                                                                              */
/*           DESCRIPTION       :   This function gives out an error message on  */
/*                                  to the standard output.                     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   Gvector() was vector()                       */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

#ifdef NEVERDEFINED
/* Dated GENOUD Memory run-time error */
void nrerror(char error_text[])
{
  error("Dated GENOUD Memory run-time error...\n %s.\n ...now exiting to system...\n",error_text);
}
#endif

/* Dated fault tolerant GENOUD Memory run-time error */
#ifdef NEVERDEFINED
void nrerror(char error_text[])
{
  warning("Dated GENOUD Memory run-time error...\n %s.\n I SHOULD exit the system, but I will not because of cheap fault tolerance.\n",error_text);
  return;
}
#endif

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   Gvector()                                     */
/*                                                                              */
/*           SYNOPSIS          :   double *Gvector(nl,nh)                         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a single dimensional   */
/*                                  double array after allocating memory from    */
/*                                  indices nl to nh                            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_org_in_eq(),                            */
/*                                 initialize_x2(),                             */
/*                                 oper1(),                                     */
/*                                 oper2(),                                     */
/*                                 oper3(),                                     */
/*                                 optimization(),                              */
/*                                 main().                                      */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



VECTOR Gvector(int nl, int nh)
{
        VECTOR v;

        if (nh <  nl)
          return(NULL);

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        /* if (!v) nrerror("allocation failure in Gvector()"); */
        return v-nl;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   ivector()                                    */
/*                                                                              */
/*           SYNOPSIS          :   int *vector(nl,nh)                           */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a single dimensional   */
/*                                  integer array after allocating memory from  */
/*                                  indices nl to nh                            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_probability(),                          */
/*                                 initialize_x2(),                             */
/*                                 main(),                                      */
/*                                 optimization(),                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



IVECTOR ivector(int nl, int nh)
{
        IVECTOR v;

        if (nh <  nl)
          return(NULL);

        v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
        /* if (!v) nrerror("allocation failure in ivector()"); */
        return v-nl;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   matrix()                                     */
/*                                                                              */
/*           SYNOPSIS          :   double *matrix(nrl,nrh,ncl,nch)               */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a two dimensional      */
/*                                  double array after allocating memory for the */
/*                                  rows from indices nrl to nrh, and for the   */
/*                                  columns from ncl to nch                     */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   det(),                                       */
/*                                 find_org_in_eq(),                            */
/*                                 initialize_x2(),                             */
/*                                 inverse(),                                   */
/*                                 main(),                                      */
/*                                 oper4(),                                     */
/*                                 oper5(),                                     */
/*                                 optimization(),                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



MATRIX matrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        MATRIX m;

        if (nrh <  nrl)
          return(NULL);
        if (nch <  ncl)
          return(NULL);

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        /* if (!m) nrerror("allocation failure 1 in matrix()"); */
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                /* if (!m[i]) nrerror("allocation failure 2 in matrix()"); */
                m[i] -= ncl;
        }
        return m;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   imatrix()                                    */
/*                                                                              */
/*           SYNOPSIS          :   int *imatrix(nrl,nrh,ncl,nch)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a two dimensional      */
/*                                  integer array after allocating memory for   */
/*                                  the rows from indices nrl to nrh, and for   */
/*                                  the columns from ncl to nch                 */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_probability(),                          */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



IMATRIX imatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        IMATRIX m;

        if (nrh <  nrl)
          return(NULL);
        if (nch <  ncl)
          return(NULL);

        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
        /* if (!m) nrerror("allocation failure 1 in imatrix()"); */
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
                /* if (!m[i]) nrerror("allocation failure 2 in imatrix()"); */
                m[i] -= ncl;
        }
        return m;
}

void free_vector(double *v, int nl)
{
     if (v == NULL)
      return;
     else
      free((double*) (v+nl));
}

void free_ivector(int *v, int nl)
{
     if (v == NULL)
      return;
     else
      free((unsigned int*) (v+nl));
}
void free_matrix(double **m, int nrl, int nrh, int ncl)
{
     int i;

     if (m == NULL)
      return;
     else
      {
        for(i=nrh;i>=nrl;i--) free((double*) (m[i]+ncl));
          free((double*) (m+nrl));
      }
}
void free_imatrix(int **m, int nrl, int nrh, int ncl)
{
     int i;

     if (m == NULL)
      return;
     else
      {
        for(i=nrh;i>=nrl;i--) free((unsigned int*) (m[i]+ncl));
          free((unsigned int*) (m+nrl));
      }
}

