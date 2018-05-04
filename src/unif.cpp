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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/unif.cpp,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/


#include <stdio.h>
#include "unif.h"

void ruxorv (integer *iseed, int n, double *rvec, integer *aux)
{
/* the double precision random uniforms in rvec have only 31 random bits
     set *iseed!=0 to reinitialize the sequence, else use *iseed==0
   *iseed always equals 0 on return
   aux points to an integer vector of length > 1279+3
*/
  int i, ibloc, nw;
  static integer wrk[TLPDBSIZ];
  double fk = 1.0/((double) TLPMOD) ;

  /* if n==0, initialize the TLP seed vector in aux[] and return */
  if (n == 0) {
    /* initialize only if *iseed > 0 */
    if (*iseed != 0) tlpseq (iseed, 0, wrk, aux);
  } else {
    nw = n;
    ibloc = 0;
    while (nw > TLPDBSIZ) {
      tlpseq (iseed, TLPDBSIZ, wrk, aux);
      if (ibloc > 0) {
	for (i=0; i<TLPDBSIZ; i++) rvec[ibloc+i] = wrk[i] * fk ;
      } else {
	for (i=0; i<TLPDBSIZ; i++) rvec[i] = wrk[i] * fk ;
      }
      nw -= TLPDBSIZ;
      ibloc += TLPDBSIZ;
    } ;
    if (nw > 0) {
      tlpseq (iseed, nw, wrk, aux);
      if (ibloc > 0) {
	for (i=0; i<nw; i++) rvec[ibloc+i] = wrk[i] * fk;
      } else {
	for (i=0; i<nw; i++) rvec[i] = wrk[i] * fk;
      }
    }
  }
}

/* Tausworthe-Lewis-Payne generator,
    for primitive trinomial with k=K, q=Q.  Here K=1279, Q=216.
   see Bratley, Fox and Schrage (1983, 190) and Bright and Enison (1979)
*/
void tlpseq (integer *iseed, int n, integer *rvec, integer *aux)
{
  static integer k = 1279, q = 216, k0 = 89, q0 = 38,
    seed = (integer) 524287 ;
  static integer aux0[89+3];
  int i;
  integer seedw;

  /* aux is a vector of length > k+3 */
  /* initialize aux if nonzero *iseed is supplied or
     if aux[k+2] is not set to the magic value */
  if (*iseed != ZEROI || aux[k+2] != k) {
    if (*iseed < ZEROI) *iseed = -(*iseed);
    seedw = (*iseed != seed) ? (*iseed ^ seed) : seed ;
    *iseed = ZEROI;
    for (i=0; i<k0; i++) {
      tauint (&seedw) ;
      aux0[i] = seedw;
    }
    aux0[k0] = k0 - 1;
    aux0[k0+1] = aux0[k0] - q0;
    aux0[k0+2] = k0;
    /* initialize the k=K, q=Q generator using the k=89, q=38 generator */
    tlpcor (k0, k, aux, aux0);
    /* store initial recurrence location indexes and the magic value in aux */
    aux[k] = k - 1;
    aux[k+1] = aux[k] - q;
    aux[k+2] = k;
  }
  if (n > 0) tlpcor (k, n, rvec, aux);
}

/* Tausworthe-Lewis-Payne generator, core recursion.
   assumes aux is valid
   see Bratley, Fox and Schrage (1983, 190) and Bright and Enison (1979)
*/
void tlpcor (integer k, int n, integer *rvec, integer *aux)
{
  /* aux is a vector of length k+2 */
  int ii;
  integer i, j, inew;

  /* get the recursion location indexes out of AUX */
  i = aux[k];
  j = aux[k+1];
  for (ii=0; ii<n; ii++) {
    inew = aux[j] ^ aux[i] ;
    aux[i] = inew;
    rvec[ii] = inew;
    j = j==0 ? k-1 : j - 1;
    i = i==0 ? k-1 : i - 1;    
  }
  /* store the recursion location indexes in aux */
  aux[k] = i;
  aux[k+1] = j;
}

void tauint (integer *ix)
{
/* Tausworthe generator (integer part)
   requires full-word logical operations, 32-bit words, and
     two's complement arithmetic
   from Bratley, Fox and Schrage (1983, 203)
   a recommended initial value for IX is 524287
*/
  integer i, j;

  i = *ix;
  /* do the right shift */
  j = i/8192;
  /* exclusive-or */
  j ^= i;
  /* lower-order bits are done.  now do the rest */
  i = j;
  /* do the left shift */
  j <<= 18;
  /* exclusive-or, handling cases where the left shift makes j negative */
  *ix = j<0 ? i^(-j) : i^j ;
}
