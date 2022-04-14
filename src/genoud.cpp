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
#include <random>
#include "genoud.h"

/* frange_ran integer definition */
int NewUnifSeed[MAXTHREADS];
int RandIntSeed[MAXTHREADS];
int ThreadNumber;
std::mt19937 mt_engine_int;
std::mt19937 mt_engine_unif;

extern double func4g(double *X);

void genoud(struct GND_IOstructure *Structure)
{

  extern int NewUnifSeed[MAXTHREADS];
  extern int RandIntSeed[MAXTHREADS];
  extern int ThreadNumber;
  extern std::mt19937 mt_engine_int;
  extern std::mt19937 mt_engine_unif;

  MATRIX 
    domains,      /*Matrix for Domains*/
    final_mat;    /*The final Domain*/

  VECTOR 
    LowerBounds,
    UpperBounds,
    X;            /*Initial values for the variables*/

  INDEX 
    fin;          /*Size of final matrix*/

  int 
    i,            /*Counter variable*/
    nvars;        /*Remaining variables after p-variables were eliminated*/

  time_t start_time,
         stop_time;
  double delta_time;
  long   hours, minutes, seconds;
  char   time_str[27];

  static short firsttime=1;

  /* FILE *output; */

  /* Lamarck Child Test Variables */
  // char LVMchar[1000];
  long LVMreturn;

/********************************************************************************/

  LVMreturn=0;
  time(&start_time);
  strcpy(time_str, ctime(&start_time));

  /* Fault Tolerant MinMax */
  if (Structure->MinMax > 0) Structure->MinMax=1;
  else Structure->MinMax=0;

  if (Structure->OutputType!=0) {
    error("output path/type must be the 'R' default");
  } /* else {
    output=stdout;
    } */

  if(Structure->PrintLevel>0)
    Rprintf("\n\n%s",time_str);

  ThreadNumber=Structure->ThreadNumber;
  if (ThreadNumber > MAXTHREADS) {
    error("No more than %d threads allowed\n\n", MAXTHREADS);
  }
  
  if (Structure->ProvideSeeds == 1) {
     //  Only toggle the instance number if we have threads!
    NewUnifSeed[ThreadNumber] = Structure->UnifSeed;
    RandIntSeed[ThreadNumber] = Structure->IntSeed; 
  } 
  else {
     /*If a Seed is NOT provided, use the base random number and run from that base!
       In other words, might as well the ThreadNumber equal to 0 
      */
    
    if (firsttime==1) {
      //NewUnifSeed[0] = BaseNewUnifSeed;
      //RandIntSeed[0] = BaseRandIntSeed;	
      firsttime=0;
    }
    ThreadNumber = 0;
 }
  
  // Now need to initialize the std::mt19937 objects
  mt_engine_int.seed(RandIntSeed[ThreadNumber]);
  mt_engine_unif.seed(NewUnifSeed[ThreadNumber]);
  
  
  fin.r =   Structure->nvars;            /*total number of inequalities + domains*/
  fin.c =   Structure->nvars+2;          /*x2 variables + lower limits + upper limits*/

  nvars = Structure->nvars;

  /*Allocating memory for all the vectors and matrices*/
  final_mat = matrix(1,fin.r,1,fin.c);
  domains = matrix(1,nvars,1,3);
  LowerBounds = Gvector(1, nvars);
  UpperBounds = Gvector(1, nvars);
  X = Gvector(1,nvars);

  /* SETUP DOMAINS */

  /* alter the domains when we are using integers because of the "open
     set" problem.  We only extend the UPPER domain bound */
  if (Structure->DataType==1) {
    for(i=0; i<Structure->nvars; i++)
      Structure->Domains[i][1] = Structure->Domains[i][1] + 0.99999;
  }
  
  for(i=1; i<=Structure->nvars; i++)
    {
      domains[i][1] = Structure->Domains[i-1][0];
      domains[i][2] =  (double) i;
      domains[i][3] = Structure->Domains[i-1][1];
    }
  
  for (i=1; i<=nvars; i++)
    {
      LowerBounds[i] = domains[i][1];
      UpperBounds[i] = domains[i][3];
    }

  /*Initialization*/
  if(Structure->PrintLevel>0)
    print_domains(domains,nvars,Structure->DataType);

  if (Structure->DataType==1) {
    JaIntegerOptimization(Structure, X, domains);
  }
  else {
    optimization(Structure, X, domains);
  }

  /* free memory */
  free_matrix(final_mat,1,fin.r,1);
  free_matrix(domains, 1, nvars, 1);
  free_vector(LowerBounds,1);
  free_vector(UpperBounds,1);
  free_vector(X,1);

  /* print final numbers and the time this has taken */
  if(Structure->PrintLevel>0)
    {
      Rprintf( "\n");
      Rprintf( "Solution Found Generation %d\n", Structure->oPeakGeneration);
      Rprintf( "Number of Generations Run %d\n", Structure->oGenerations);
    }
  time(&stop_time);

  strcpy(time_str, ctime(&stop_time));
  if(Structure->PrintLevel>0)
    Rprintf("\n%s",time_str);

  delta_time = difftime(stop_time, start_time);
  hours   = (int) (delta_time/3600);
  minutes = (int) (delta_time - (hours*3600))/60;
  seconds = (int) delta_time - (hours*3600) - (minutes*60);
  
  if(Structure->PrintLevel>0)
    {
      Rprintf("Total run time : %d hours %d minutes and %d seconds\n", 
	      hours, minutes, seconds);
      /* fflush(output); */
    }

  /* no longer allowed */
  /* if (Structure->OutputType==1) fclose(output); */
}

