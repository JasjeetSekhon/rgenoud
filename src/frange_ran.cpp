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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/frange_ran.cpp,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/
#include <random>
#include "genoud.h"

/* replacements for random number functions in GENOCOP's frange_ran.c */

double frange_ran(double llim, double ulim)
     /*
       double ulim;
       double llim;
     */
{
  // Define the mersenne twister 19937 generator seeded with NewUnifSeed[ThreadNumber]
  extern std::mt19937 mt_engine_unif;
  
  // Create a uniform real distribution on [0,1]
  std::uniform_real_distribution<double> unif_dist(llim, ulim);
  
  double generated_num = unif_dist(mt_engine_unif);

  return(generated_num);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   irange_ran()                                 */
/*                                                                              */
/*           SYNOPSIS          :   int irange_ran(llim,ulim)                    */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a random integer       */
/*                                  between the llim and ulim.                  */
/*                                                                              */
/********************************************************************************/
int irange_ran(int llim, int ulim)
{
  extern std::mt19937 mt_engine_int;
  
  // Create a uniform int distribution on [llim,ulim]
  std::uniform_int_distribution<int> unif_int_dist(llim, ulim);

  int generated_int = unif_int_dist(mt_engine_int);
  
  return(generated_int);
}
