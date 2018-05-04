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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/multiply.cpp,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/

#include "genoud.h"
/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   mmprod()                                     */
/*                                                                              */
/*           SYNOPSIS          :   void mmprod(m,nm,n,mul_cm,mul_am,mul_bm)     */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the product of two       */
/*                                  double matrices                              */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_org_in_eq(),                            */
/*                                 main().                                      */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void mmprod(int m, int nm, int n, MATRIX mul_cm, MATRIX mul_am, MATRIX mul_bm)
     /*
       int m,      row of mul_am matrix
       nm,         column of mul_am and row of mul_bm matrices
       n;          row of mul_bm matrix
       MATRIX mul_cm,   the final matrix
       mul_am,          the first matrix to be multiplied
       mul_bm;          the second matrix to be multiplied
     */
{
 int i,j,k;       /*counter variables, where i=m, j=nm, k=n*/

 for(i= 1; i<=m; i++)
  for(j = 1; j<=n; j++)
   {
     mul_cm[i][j] = 0.0;
     for (k= 1;k<nm+1;k++)
       mul_cm[i][j] = mul_am[i][k]*mul_bm[k][j]  + mul_cm[i][j];
    }
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   mvprod()                                     */
/*                                                                              */
/*           SYNOPSIS          :   void mvprod(m,nm,cm,am,bm)                   */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the product of a double   */
/*                                  vector and a double matrix                   */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   find_org_in_eq(),                            */
/*                                 main().                                      */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void mvprod(int m, int nm, VECTOR cm, MATRIX am, VECTOR bm)
     /*
       int m,     row of matrix am and length of vector cm
       nm;        column of matrix am
       VECTOR cm,  the final vector
       bm;         original vector
       MATRIX am;  original matrix
     */
{
int i,k;    /*counter variables, where i=m, k=nm*/

 for (i = 1; i<=m; i++)
  {
   cm[i]=0.0;
   for (k = 1; k<=nm; k++)
    cm[i] = cm[i] + am[i][k]*bm[k];
  }
}

