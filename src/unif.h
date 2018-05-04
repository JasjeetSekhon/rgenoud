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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/unif.h,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/


/* Tausworthe-Lewis-Payne generator */

#ifndef DEFTLPGEN
typedef int integer;

#define TLPAUXSIZE 1300
#define TLPDBSIZ 2000
#define ZEROI ((integer)0)
#define TLPMOD 2147483647

void ruxorv (integer *, int, double *, integer *);
void tlpseq (integer *, int, integer *, integer *);
void tauint (integer *);
void tlpcor (integer, int, integer *, integer *);

#define DEFTLPGEN
#endif
