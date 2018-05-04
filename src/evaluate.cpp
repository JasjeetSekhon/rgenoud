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

  June 27, 2013

*/

#include "genoud.h"
#include "gradient.h"

extern "C" 
{

  double genoud_optim(SEXP fn_optim, SEXP rho, double *X, long parameters);

  void RlexicalSort(SEXP fnLexicalSort, SEXP rho,
			   double **population, 
			   short int MinMax, long pop_size, long nvars, long lexical_end,
			   short int type);

  long RmemoryMatrixEvaluate(SEXP fnMemoryMatrixEvaluate, SEXP rho,
				    double **Memory, double **population, 
				    short int MinMax, long pop_size, long UniqueCount,
				    long nvars, long lexical, long lexical_end);

  void userGradientfn(SEXP fnGR, SEXP rho, double *parms, double *grad, long nvars);
}

long Gnvars[MAXINSTANCES];
struct GND_IOstructure *ExternStructure;

int JaIntegerCMP(double **a, double **b) 
{
  extern long Gnvars[MAXINSTANCES];
  extern struct GND_IOstructure *ExternStructure;


  long i = 0;
  long nvars;

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=nvars; i++) {
    if ( (int) a[0][i] != (int) b[0][i])
      break;
  }

  if ( (int) a[0][i] >  (int) b[0][i]) i = 1;
  else if ( (int) a[0][i] <  (int) b[0][i]) i = -1;

  return i;
} /* end of JaIntegerCMP */


int JaDoubleCMP(double **a, double **b) 
{
  extern long Gnvars[MAXINSTANCES];
  extern struct GND_IOstructure *ExternStructure;

  long i = 0;
  long nvars;

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=nvars; i++) {
    if ( a[0][i] != b[0][i])
      break;
  }

  if ( a[0][i] > b[0][i]) i = 1;
  else if ( a[0][i] < b[0][i]) i = -1;

  return i;
} /* end of JaCMP */


/* Cummulative probability on crossover */
/* Random probability on mutation       */
/* NO multiple hits per agent possible  */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   optimization()                               */
/*                                                                              */
/*           SYNOPSIS          :   double optimization(X,x1,x2,fin_mat,rc,tot_eq) */
/*                                                                              */
/*           DESCRIPTION       :   This procedure initializes the population    */
/*                                  with the values X passed from main, and     */
/*                                  evaluates them.  After assigning weight     */
/*                                  for each member of the populaiton, a group  */
/*                                  of them are chosen to reproduce and a group */
/*                                  is chosen to die.  Genetic operators are    */
/*                                  applied and a new generation is produced    */
/*                                  to replace the members that died.  This     */
/*                                  cycle continues for the number of times,    */
/*                                  user specifies in the input file            */
/*                                                                              */
/*           FUNCTIONS CALLED  :   assign_probab(),                             */
/*                                 evaluate(),                                  */
/*                                 find_cum_probab(),                           */
/*                                 find_live_die(),                             */
/*                                 find_parent(),                               */
/*                                 ivector(),                                   */
/*                                 matrix(),                                    */
/*                                 oper1(),                                     */
/*                                 oper2(),                                     */
/*                                 oper3(),                                     */
/*                                 oper4(),                                     */
/*                                 oper5(),                                     */
/*                                 oper6(),                                     */
/*                                 print_population(),                          */
/*                                 sort(),                                      */
/*                                 Gvector() was vector().                      */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void optimization(struct GND_IOstructure *Structure, VECTOR X, 
		    MATRIX domains)
{

  extern struct GND_IOstructure *ExternStructure;

  MATRIX 
    new_genera,   /*Temporary storage for the new generation*/
    population,   /*Population of x2 variables*/
    temp;

  VECTOR probab,       /*Probability of agents to die or live*/
         cum_probab,   /*Cumilative probability of agents*/
         t_vec;

  IVECTOR live;
  /* for oper4 */
  IVECTOR parents;



  long count_gener= 1; /*Counter to keep track of the number of generations*/
  unsigned long peak_cnt;

  int                     /*Total number of agents chosen to reproduce*/
    j1,
    j2,
    j3,
    j4,
    j5,
    j6,
    j7,
    j8,
    oper,
    ocnt,
    B,                     /*Parameter for the 3rd operator - nonuniform mutation*/
    STEP,                  /*Parameter for the 5th operator - simple arithmetical crossover*/
    first_live=0,            /*Index of the two parents for crossover parents*/
    second_live=0,
    first_die,             /*Index of the two parents for crossover death*/
    second_die,
    die_now;               /*index of agent to replace in current operation*/

  long i,j, s, k;
  /* for oper 4 */
  int p2use;


  double Q;                /*Probability of the best agent*/
  FLAG  same;
  double **Jnew; 

  double bfgsfit;

  double *grad, *evalX, *finalhessin, *bfgsoutX;

  int nochange_gen=0;

  double oldfitvalue=0;

  int IncreaseGenerations;
  short int GradientTrigger=0;
  short int BoundaryTrigger;
  long InstanceNumber;

  long nvars, MaxGenerations, WaitGenerations, count;
  long pop_size, P, P0, P1, P2, P3, P4, P5, P6, P7, P8;
  short int MinMax, GradientCheck, BoundaryEnforcement, UseBFGS, HardMaximumNumber=0;
  double SolutionTolerance, *Results, *Gradients;
  short PrintLevel, HardGenerationLimit;

  /* Old variables which may change when SetRunTimeParameters is run during a run! */
  long pop_size_old;

  /* Summary Statistics (mean, variance etc) */
  /* double popmean, popvar, popwrk, popstat; */

  /* Population Print population*/
  FILE *popout;
  long *tobs, nnull;
  double *mean, *var, *skew, *kur;

  /* Stuff for the Unique Stuff (how's that for an informative comment! */
  /* A big Matrix which remembers all of our past evaluations. It's
     maximum memory is set in genoud.h */
  extern long Gnvars[MAXINSTANCES];
  double **Memory;
  long MemorySize=0, UniqueCount=0, OldUniqueCount=0;

  /* fine two unique parents count */
  long SameCount, UniquePairs;

  ExternStructure=Structure;

  Results=Structure->oResults;
  Gradients=Structure->oGradients;

  /* Structure Done */
  SetRunTimeParameters(Structure, 1,
		       &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
		       &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
		       &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
		       &PrintLevel, &HardGenerationLimit);

  /*Space allocation for all the vectors and matrices involved*/
  long lexical_end = (Structure->Lexical-1)+nvars+2;
  /* population[][0] = fitness value (first)
     population[][1:nvars] = parameter values
     population[][nvars+1] = flag for fitting
     population[][(nvars+2):((Structure->Lexical-1)+nvars+2)] = other fitness for Lexical fitting
  */
  population    = JaMatrixAllocate(pop_size+2, lexical_end);
  new_genera    = JaMatrixAllocate(pop_size+2, lexical_end);

  /* reset population to get rid of odd things being passed to R */
  for(i=1; i<=pop_size; i++)
    {
      for(j=0; j<lexical_end; j++)
	{
	  population[i][j] = 0;
	}
    }

  VECTOR LexicalReturn;
  VECTOR oldfitvalueVEC;
  short int LexicalFitsImproving;
  /* Transform is logically similar to Lexical even if there is only one fit criterion */
  if(Structure->Lexical > 1 || Structure->Transform == 1)
    {
      LexicalReturn = (double *)  malloc(Structure->Lexical*sizeof(double));  
      oldfitvalueVEC = (double *)  malloc(Structure->Lexical*sizeof(double));  
    }

  temp       = matrix(0,nvars+1,0,nvars);
  probab     = Gvector(1,pop_size);
  t_vec      = Gvector(1,nvars);
  cum_probab = Gvector(1,pop_size);
  live       = ivector(1,pop_size);

  /*for oper4 Find max(2,nvars) parents for crossover operator 4*/
  p2use = nvars > 2 ? nvars : 2;
  parents    = ivector(1,p2use);

  Gnvars[Structure->InstanceNumber]=nvars;

  if (Structure->MemoryUsage==1)
    {
      if (HardGenerationLimit==0)
	MemorySize=3*(MaxGenerations+1)*pop_size+1+pop_size;
      else
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      
      Memory = JaMatrixAllocate(MemorySize, lexical_end);
    }

  grad = (double *) malloc((nvars)*sizeof(double));
  evalX = (double *) malloc((nvars)*sizeof(double));
  finalhessin = (double *) malloc(((nvars*nvars)+(nvars))*sizeof(double));
  bfgsoutX = (double *) malloc((nvars+1)*sizeof(double));

  /* populationstats variables */
  mean = (double *) malloc((nvars+1)*sizeof(double));
  var = (double *) malloc((nvars+1)*sizeof(double));
  skew = (double *) malloc((nvars+1)*sizeof(double));
  kur = (double *) malloc((nvars+1)*sizeof(double));
  tobs = (long *) malloc((nvars+1)*sizeof(long));

  Q=0.5;
  B=6;
  STEP=10;

  if(PrintLevel>0)
    {
      switch(MinMax) {
      case 0:
	Rprintf("Minimization Problem.\n");  
	break;
      case 1:
	Rprintf("Maximization Problem.\n");  
	break;
      }
    }

  /*
    if (PrintLevel>2) {
    Rprintf("Parameter B (hardcoded): %d\n", B); 
    Rprintf("Parameter Q (hardcoded): %f\n", Q);
    }
  */

  peak_cnt = 0;

  pop_size_old=0;
  if (Structure->ShareType == 1 || Structure->ShareType == 3) {

    if(PrintLevel>0)
      Rprintf( "Using old population file to initialize new population.\n");

    if((popout = fopen(Structure->ProjectPath, "r")) == NULL) {
      Rprintf("         Generating new population\n");
      warning("Unable to open the old project file: %s", Structure->ProjectPath);
    }
    else {
      pop_size_old=ReadPopulation(population, pop_size, nvars, popout, PrintLevel);
      fclose(popout);
      if (pop_size_old<2) {
	warning("The old population file appears to be from a different genoud specification.");
	pop_size_old=0;
      }
    }
    if (PrintLevel>1) {
      if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	warning("Unable to open the project file: %s", Structure->ProjectPath);
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
        /* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);

	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 0, nvars+1, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	free_ivector(parents, 1);

        if(Structure->Lexical > 1 || Structure->Transform == 1)
	  {
	    free(LexicalReturn);
	    free(oldfitvalueVEC);
	  }

	error("Fatal Error. See output for diagnostic information.");
      }
      fclose(popout);
    }
  } /* end of ShareType 0 */
  else {
    if (PrintLevel>1) {
      if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	Rprintf("Unable to open the project file: %s", 
		Structure->ProjectPath);

	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
        /* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);

	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 0, nvars+1, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	free_ivector(parents, 1);

        if(Structure->Lexical > 1 || Structure->Transform == 1)
	  {
	    free(LexicalReturn);
	    free(oldfitvalueVEC);
	  }

	error("Fatal Error. See output for diagnostic information.");
      }
      fclose(popout);
    }
  }

  /* The new initial value matrix: setting a new initial value for every individual */
  if (ExternStructure->nStartingValues > 0) 
    {
      /* Adjust old starting values (from ReadPopulation) so we have enough room for our 
	 starting.values */
      pop_size_old = pop_size_old-ExternStructure->nStartingValues-1;
      if (pop_size_old < 0)
	pop_size_old = 0;
      
      // seed the starting values until we run out of population or starting values!
      j = pop_size_old;
      
      for(s=0; s<ExternStructure->nStartingValues; s++) {
	j++;
	for(i=1; i<=nvars; i++) {
	  population[j][i] = ExternStructure->StartingValues[s][i-1];
	  population[j][nvars+1] = -1.0;
	}
      } // end of for loop
      pop_size_old = j;

      // randomly add on people if we still have population left over!
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of we have starting values!
  else 
    {
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of else

  if (Structure->MemoryUsage==1)
    {
      OldUniqueCount=UniqueCount;

      if (UniqueCount==0)
	UniqueCount = 1;

      UniqueCount = RmemoryMatrixEvaluate(Structure->fnMemoryMatrixEvaluate, Structure->rho,
					  Memory, population,
					  MinMax, pop_size, UniqueCount,
					  nvars, Structure->Lexical, lexical_end);

      if ( (UniqueCount+pop_size) >= MemorySize )
	{
	  Structure->MemoryUsage=0;
	  warning("Turned Off MemoryMatrix because memory usage was too great.");
	} /* end of if */
    } // end of Memory based evaluation
  else
    {
      for (i=1; i<=pop_size; i++) 
	{
	  if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0)
	    {
	      for(j=1; j<=nvars; j++)
		X[j] = population[i][j];
	      
              if (Structure->whichFUN == 1) // neither Lexical, nor Transform
		{
		  population[i][0] = evaluate(Structure->fn, Structure->rho, X, nvars, MinMax);
		} 
              else if(Structure->whichFUN == 2) // Lexical but not Transform
		{
		  EvaluateLexical(Structure->fn, Structure->rho, 
				  X, nvars, Structure->Lexical, MinMax, LexicalReturn);
		  
		  population[i][0] = LexicalReturn[0];
		  count = 0;
		  for(j=(nvars+2);j<lexical_end;j++)
		    {
		      count++;
		      population[i][j] = LexicalReturn[count];
		    }		      
                }  // else if
              else // Transform
                {
                  EvaluateTransform(Structure->fn, Structure->rho,
				    X, nvars, Structure->Lexical, MinMax, LexicalReturn);
		  
                  population[i][0] = LexicalReturn[0];
                  count = 0;
                  if(Structure->Lexical > 1) 
		    for(j=(nvars+2);j<lexical_end;j++)
		      {
			count++;
			population[i][j] = LexicalReturn[count];
		      }
                  if (BoundaryEnforcement==0) 
		    for(j=1; j<=nvars; j++)
		      {
			population[i][j] = X[j];
		      }
                  else 
		    for(j=1; j<=nvars; j++)
		      {
			if(X[j] < domains[j][1])
			  {
			    if(PrintLevel>0)
			      {
				Rprintf(
					"\nNOTE: Transformed individual below lower bound.\n");
				Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e \t Bound: %e\n\n",
					count_gener, j, X[j], domains[j][1]);
			      }
			    warning("Transformed individual below lower bound. Generation: %d; Parameter: %d; Value: %e; Bound: %e",
				    count_gener, j, X[j], domains[j][1]);
			  }
			if(X[j] > domains[j][3])
			  {
			    if(PrintLevel>0)
			      {
				Rprintf(
					"\nNOTE: Transformed individual above upper bound.\n");
				Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e \t Bound: %e\n\n",
					count_gener, j, X[j], domains[j][3]);
			      }
			    warning("Transformed individual above upper bound. Generation: %d; Parameter: %d; Value: %e; Bound: %e",
				    count_gener, j, X[j], domains[j][3]);
			  }
			population[i][j] = X[j]; // put the transformation in anyway
		      }//end for
                }//end Transform
	    }
	} //end of i loop
    } // end of default evaluation

  if(Structure->MemoryUsage!=1)
    {
      /*Sort the initial individuals based on their evaluation function*/
      if (Structure->Lexical < 2)
	{
	  sort(MinMax,population,pop_size,0);
	}
      else
	{
	  /* in eval.cpp because it is like the EvaluateLexical() function */
	  RlexicalSort(Structure->fnLexicalSort, Structure->rho,
		       population,
		       MinMax, pop_size, nvars, lexical_end, 1);
	}
    }

  peak_cnt = count_gener;

  /* since we are on generation 0 */
  oldfitvalue=population[1][0];
  if(Structure->Lexical  > 1)
    {
      oldfitvalueVEC[0]=population[1][0];
      k = 1;
      for (i=(nvars+2);i<lexical_end;i++)  {
	oldfitvalueVEC[k]=population[1][i];
	k++;  
      } /* for (i=(nvars+2);i<lexical_end;i++) */
    }

  /*
  if(PrintLevel>0)
    {
      Rprintf("\nThe best initial individual is:\n");
      print_vector(population[1],1,nvars,output);

      if (Structure->Lexical > 1)
	{
	  Rprintf("\nbest (lexical) fitness:\n");
	  Rprintf("%e  ", population[1][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[1][j]);
	    }		      
	  Rprintf("\n");
	} else {
	Rprintf("\nbest fitness: %e\n", population[1][0]);
      }
      Rprintf("\n");

      if (Structure->Lexical > 1)
	{      
	  Rprintf("The worst (lexical) fitness is:\n");
	  Rprintf("%e  ", population[pop_size][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[pop_size][j]);
	    }		   
	  Rprintf("\n");   	  
	} else {
	Rprintf("The worst fit is: %e\n", 
		population[pop_size][0]);
      }
      Rprintf("\n");
    }
  */

  if(PrintLevel==1)
    {
      if (Structure->Lexical > 1)
	{
	Rprintf("\n\nGeneration#\t    Solution Values (lexical)\n");	  
	Rprintf("\n%7d \t%e  ", 0, population[1][0]);
	for(j=(nvars+2);j<lexical_end;j++)
	  {
	    Rprintf("%e  ", population[1][j]);
	  }		      
	Rprintf("\n");	
	} else {
	Rprintf("\n\nGeneration#\t    Solution Value\n");
	Rprintf("\n%7d \t%e\n", 0, population[1][0]);
      }
    }

  /* compute and print mean and variance of population */
  if (PrintLevel>1) {
    Rprintf("GENERATION: 0 (initializing the population)\n");
    populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
    
    if(Structure->Lexical > 1)
      {
	Rprintf( "Lexical Fit..... %e  ", population[1][0]);
	for(j=(nvars+2);j<lexical_end;j++)
	  {
	    Rprintf("%e  ", population[1][j]);
	  }		      
	Rprintf("\n");	    
      }
    else
      {
	Rprintf( "Fitness value... %e\n", population[1][0]);
	Rprintf( "mean............ %e\n", mean[0]);
	Rprintf( "variance........ %e\n", var[0]);
	/*
	  Rprintf( "skewness........ %e\n", skew[i]);
	  Rprintf( "kurtosis........ %e\n", kur[i]);
	*/
      }

    nnull = pop_size-tobs[0];
    if(nnull > 0)
      Rprintf( "#null........... %d\n", nnull);
    if(Structure->MemoryUsage==1)
      Rprintf( "#unique......... %d, #Total UniqueCount: %d\n", 
	      UniqueCount-OldUniqueCount, UniqueCount);
    /* Rprintf( "tobs............ %d\n", tobs[i]); */
    
    for (i=1; i<=nvars; i++) {
      Rprintf( "var %d:\n", i);
      Rprintf( "best............ %e\n", population[1][i]);
      Rprintf( "mean............ %e\n", mean[i]);
      Rprintf( "variance........ %e\n", var[i]);
      /*
	Rprintf( "skewness........ %e\n", skew[i]);
	Rprintf( "kurtosis........ %e\n", kur[i]);
      */
      nnull = pop_size-tobs[i];
      if(nnull > 0)
	Rprintf( "#null........... %d\n", nnull);
      /* Rprintf( "tobs............ %d\n", tobs[i]); */
    }
  } /* end of printlevel if */
  
  /* if(PrintLevel>0)
    fflush(output);  */
      
  /* Print the population file */
  if ( PrintLevel == 1 ) {
    if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
      Rprintf("Unable to open the project file: %s", 
	      Structure->ProjectPath);

      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);

      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 0, nvars+1, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      free_ivector(parents, 1);

      if(Structure->Lexical > 1 || Structure->Transform == 1)
	{
	  free(LexicalReturn);
	  free(oldfitvalueVEC);
	}
			  
	error("Fatal Error. See output for diagnostic information.");
    }
    print_population((int) pop_size, (int) nvars, 0, (int) Structure->Lexical, population, popout);
    fclose(popout);
  } /* end of PrintLevel if */
  if ( PrintLevel>1 ) {
    if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
      Rprintf("Unable to open the project file: %s", 
	      Structure->ProjectPath);

      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);

      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 0, nvars+1, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      free_ivector(parents, 1);

      if(Structure->Lexical > 1 || Structure->Transform == 1)
	{
	  free(LexicalReturn);
	  free(oldfitvalueVEC);
	}
      
      error("Fatal Error. See output for diagnostic information.");
    }
    print_population((int) pop_size, (int) nvars, 0, (int) Structure->Lexical, population, popout);
    fflush(popout);
    fclose(popout);
  }

  /* Interrupt setup.  Let's print a nice message to recover the best
     solution so far if at least generation 0 has been run */
  if (PrintLevel > 0 & (strcmp(Structure->ProjectPath, "/dev/null")!=0))
    setVar(install("interrupted"), ScalarLogical(1), Structure->rho);

  /*Assigning probability of survival for each of the agent, with the*/
  /*probability provided by the user for the best agent*/
  assign_probab(probab,pop_size,Q); 

  /*Finding the cumulative probability of the agents*/
  find_cum_probab(cum_probab,probab,pop_size);

  /*Reproducing and evaluating for the total number of generations times*/
  do
    {

      /*Initializing the live vector*/
      for(j=1; j<=pop_size; j++)
        {
          live[j] = 0;
          for(i=0; i<lexical_end; i++)
            new_genera[j][i] = population[j][i];
	  new_genera[j][nvars+1]=0;
        }

      /*Finding the agents that will die and the agents that will reproduce*/
      find_live(cum_probab,live,pop_size,P);
      /* set die_now counter to start replacements with the worst agent.
         use of die_now is okay if the entire population (except the previous
         best) is to be replaced in each generation */
      die_now = pop_size;

      j1=j2=j3=j4=j5=j6=j7=j8=0;

      /* This was causing a difference to appear between MemoryMatrix and !MemoryMatrix runs see oper5 and oper7
      UniquePairs= UniqueCount-OldUniqueCount;
      UniquePairs= (int) (0.5*(UniquePairs*UniquePairs-UniquePairs));
      if ( MAX_OPER_UNIQUE_TRY < UniquePairs)
	UniquePairs = MAX_OPER_UNIQUE_TRY;
      */
      UniquePairs = MAX_OPER_UNIQUE_TRY;

      /* main operator loop */
      while(j1+j2+j3+j4+j4+j5+j5+j6+j7+j7+j8 < P)
        {
          oper = irange_ran(1,8);
          switch (oper)
            {
              case 1:
		/* JS Description: Uniform Mutation */
                     /*Applying the first operator, uniform mutation*/
                    if (j1 < P1)
                      {
			/*Find one parent for mutation operator 1*/
			first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}

			new_genera[die_now][nvars+1] = 1.0;
			for(i=1; i<=nvars; i++)
			  t_vec[i] = population[first_live][i];
			for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
			  oper1(t_vec,domains,nvars);
			for(i=1; i<=nvars; i++)
			  new_genera[die_now][i] = t_vec[i];
			die_now--;
			j1++;
		      }
                    break;
              case 2:
		/* JS Description: Boundary Mutation */
                    /*Applying the second operator, boundary mutation*/
                    if (j2 < P2)
                      {
                        /*Find one parent for mutation operator 2*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}

                        new_genera[die_now][nvars+1] = 2.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
			oper2(t_vec,domains,nvars);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j2++;
                      }
                    break;
              case 3:
		/* JS Description: Non-uniform Mutation */
                    /*Applying the third operator, non-uniform mutation*/
                    if (j3 < P3)
                      {
                        /*Find one parent for mutation operator 3*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}

                        new_genera[die_now][nvars+1] = 3.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
			for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
			  oper3(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j3++;
                      }
                    break;

              case 4:
                    /*Applying the fourth operator, GENOUD Polytope Crossover */
                    if (j4 < (int) P4)
                      {
                        /*Find max(2,nvars) parents for crossover operator 4*/
			for (i=1; i<p2use; i++) {
			  parents[i] = find_parent(live,pop_size);
			  live[parents[i]]++;  /* no decr. first p2use-1 parents */
			}
			parents[p2use] = find_parent(live,pop_size);
			/* check that agents to replace are in range */
			if (die_now < 2) {
			  error("No agents to be replaced\n");
			}
			new_genera[die_now][nvars+1]  = 4.0;
			for(j=1; j<=p2use; j++)
			  for(i=1; i<=nvars; i++)
			    temp[j][i] = population[parents[j]][i];
			oper4(temp,p2use,nvars);
			for(i=1; i<=nvars; i++)
			  new_genera[die_now][i]  = temp[1][i];
			die_now--;
                        j4++;
                      }		
                    break;
              case 5:
		/* JS Description: Simple Crossover
		   Applying the fifth operator, simple arithmetical crossover*/
                    if (j5 < (int) P5/2)
                      {
                        /*Find two distinct parents for crossover operator 5*/
                        same = TRUE;
			SameCount=0;
			while (same==TRUE) {
			  SameCount++;
			  
			  first_live  = find_parent(live,pop_size);
			  second_live = find_parent(live,pop_size);

			  if (SameCount >= (UniquePairs) ) 
			    break;

			  for(i=1; i<=nvars; i++)
			    {
			      if (population[first_live][i] != population[second_live][i])
				{
				  same = FALSE;
				  break;
				}
			    }
			} /* end of while same==TRUE loop */
			/* check that agents to replace are in range */
			if (die_now < 3) {
			  error("Not enough agents to be replaced\n");
			}
			live[first_live]--;
			live[second_live]--;
			first_die   = die_now-- ;
			second_die  = die_now-- ;
			new_genera[first_die][nvars+1]  = 5.0;
			new_genera[second_die][nvars+1] = 5.0;
                        if (!same)
                          {
                            for(i=1; i<=nvars; i++)
                              {
                                temp[1][i] = population[first_live][i];
                                temp[2][i] = population[second_live][i];
                              }
                            oper5(temp[1],temp[2],STEP,domains,nvars);
                            for(i=1; i<=nvars; i++)
                              {
                                new_genera[first_die][i]  = temp[1][i];
                                new_genera[second_die][i] = temp[2][i];
                              }
                          }
			else {
			  /* copy agent chosen twice into two new indivs */
			  for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
			  }
			}
                        j5++;
                      }
                    break;
              case 6:
		/* JS Description: Whole Non-uniform Mutation */
                    /*Applying the sixth operator, whole non-uniform mutation*/
                    if (j6 < P6)
                      {
                        /*Find one parent for mutation operator 6*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}

                        new_genera[die_now][nvars+1] = 6.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        oper6(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j6++;
                      }
                    break;
              case 7:
		/* JS Description: Heuristic Crossover */
                    /*Applying the seventh operator*/
                    if (j7 < (int) P7/2)
                      {
                        /*Find two distinct parents for operator 7*/
                        same = TRUE;
			SameCount=0;
			while (same==TRUE) {
			  SameCount++;

			  first_live  = find_parent(live,pop_size);
			  second_live = find_parent(live,pop_size);
			  
			  if (SameCount >= (UniquePairs) ) 
			    break;

			  for(i=1; i<=nvars; i++)
			    {
			      if (population[first_live][i] != population[second_live][i])
				{
				  same = FALSE;
				  break;
				}
			    }
			} /* end of while same==TRUE loop */
			/* check that agents to replace are in range */
			if (die_now < 3) {
			  error("Not enough agents to be replaced\n");
			}
			live[first_live]--;
			live[second_live]--;
			first_die   = die_now-- ;
			second_die  = die_now-- ;
			new_genera[first_die][nvars+1]  = 7.0;
			new_genera[second_die][nvars+1] = 7.0;
                        if (!same) {
			  if (first_live < second_live)
			    /* first agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[first_live][i];
			      temp[1][i] = population[second_live][i];
			    }
			  else
			    /* second agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[second_live][i];
			      temp[1][i] = population[first_live][i];
			    }
			  oper7(temp[1],temp[2],domains,nvars);
			  for(i=1; i<=nvars; i++)
			    new_genera[first_die][i]  = temp[1][i];
			  if (first_live < second_live)
			    /* first agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[first_live][i];
			      temp[1][i] = population[second_live][i];
			    }
			  else
			    /* second agent is better agent */
			    for(i=1; i<=nvars; i++) {
			      temp[2][i] = population[second_live][i];
			      temp[1][i] = population[first_live][i];
			    }
			  oper7(temp[1],temp[2],domains,nvars);
			  for(i=1; i<=nvars; i++)
			    new_genera[second_die][i]  = temp[1][i];
			}
			else {
			  /* copy agent chosen twice into two new indivs */
			  for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
			  }
			}
                        j7++;
                      }
              case 8:
		/* JS Description: Local-Minimum Crossover */
                     /*Applying the eighth operator, homotopy (BFGS) */
		if (j8 < P8 & (Structure->BFGSburnin >= 0) & (count_gener > Structure->BFGSburnin))
                      {
                        /*Find one parent for BFGS operator 1*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}

                        new_genera[die_now][nvars+1] = 8.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        oper8(Structure->fn_optim, Structure->rho, t_vec, domains, SolutionTolerance, 
			      nvars, BoundaryEnforcement, PrintLevel,
			      Structure->P9mix);

			for(i=1; i<=nvars; i++)
			  new_genera[die_now][i] = t_vec[i];
			die_now--;
			j8++;
                      }
                    break;
            }
        }

      /*Replace the population with the new generation */
      Jnew = new_genera;
      new_genera = population;
      population = Jnew;

      if (Structure->MemoryUsage==1)
	{
	  OldUniqueCount=UniqueCount;

	  UniqueCount = RmemoryMatrixEvaluate(Structure->fnMemoryMatrixEvaluate, Structure->rho,
					      Memory, population,
					      MinMax, pop_size, UniqueCount,
					      nvars, Structure->Lexical, lexical_end);	  

	  if ( (UniqueCount+pop_size) >= MemorySize )
	    {
	      Structure->MemoryUsage=0;
	      warning("Turned Off MemoryMatrix because memory usage was too great.");
	    } /* end of if */
	} // end of MemoryUsage==1
      else
	{
	  for (i=1; i<=pop_size; i++) 
	    {
	      if (population[i][nvars+1]!=0)
		{
		  for(j=1; j<=nvars; j++)
		    X[j] = population[i][j];
		  
                  if (Structure->whichFUN == 1) // neither Lexical, nor Transform
		    {
		      population[i][0] = evaluate(Structure->fn, Structure->rho, X, nvars, MinMax);
		    } 
                  else if(Structure->whichFUN == 2) // Lexical but not Transform
		    {
		      EvaluateLexical(Structure->fn, Structure->rho, 
				      X, nvars, Structure->Lexical, MinMax, LexicalReturn);

		      population[i][0] = LexicalReturn[0];
		      count = 0;
		      for(j=(nvars+2);j<lexical_end;j++)
			{
			  count++;
			  population[i][j] = LexicalReturn[count];
			}		      			  
		    }
                  else // Transform
                    {
                      EvaluateTransform(Structure->fn, Structure->rho,
                                      X, nvars, Structure->Lexical, MinMax, LexicalReturn);

                      population[i][0] = LexicalReturn[0];
                      count = 0;
                      if(Structure->Lexical > 1) for(j=(nvars+2);j<lexical_end;j++)
                        {
                          count++;
                          population[i][j] = LexicalReturn[count];
                        }
                      for(j=1; j<=nvars; j++)
                        {
                          population[i][j] = X[j];
                        }
                    }
		}
	    } //end of i loop	  
	} //end of default evaluation scheme

      if(Structure->MemoryUsage!=1)
	{
	  /*Sort the new population based on their evaluation function*/
	  if (Structure->Lexical < 2)
	    {
	      sort(MinMax,population,pop_size,0);
	    }
	  else
	    {
	      /* in eval.cpp because it is like the EvaluateLexical() function */
	      RlexicalSort(Structure->fnLexicalSort, Structure->rho,
			   population,
			   MinMax, pop_size, nvars, lexical_end, 1);
	    }
	}

      /* apply the bfgs to the best individual */
      if (UseBFGS != 0 & (Structure->BFGSburnin >= 0) & (count_gener > Structure->BFGSburnin)) {
	for (i=1; i<=nvars; i++)
	  {
	    bfgsoutX[i-1]=population[1][i];
	  }

	bfgsfit = genoud_optim(Structure->fn_optim, Structure->rho, bfgsoutX, nvars);

        if (Structure->whichFUN == 1) // neither Lexical, nor Transform
	  {
	    switch(MinMax) {
	    case 0:
	      if (population[1][0] > bfgsfit) /* minimize */
		{
		  /* is the BFGS individual in the bounds? */
		  BoundaryTrigger=0; /* outside of bounds ? */
		  for (i=0; i<nvars; i++) {
		    j = i+1;
		    if (bfgsoutX[i] < domains[j][1]) {
		      BoundaryTrigger=1;
		      if (PrintLevel>1)
			{
			  Rprintf(
				  "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
			  Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
				  count_gener, i+1, bfgsoutX[i]);
			  Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
			}
		      warning("BFGS hit on best individual produced Out of Boundary individual.");
		    }
		    if (bfgsoutX[i] > domains[j][3]) {
		      BoundaryTrigger=1;
		      if (PrintLevel>1)
			{
			  Rprintf(
				  "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
			  Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n", 
				  count_gener, i+1, bfgsoutX[i]);
			  Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
			}
		      warning("BFGS hit on best individual produced Out of Boundary individual.");
		    }
		  } /* end for loop */
		  
		  /* if we use out of bounds individuals then proceed */
		  /* 0=anything goes, 1: regular; 2: no trespassing! */
		  if (BoundaryEnforcement==0) {
		    for(i=1;i<=nvars;i++) 
		      {
			population[1][i]=bfgsoutX[i-1];
		      }
		    population[1][0]=bfgsfit;
		  }
		  else if (BoundaryTrigger==0) {
		    for(i=1;i<=nvars;i++) 
		      {
			population[1][i]=bfgsoutX[i-1];
		      }
		    population[1][0]=bfgsfit;
		  }
		} /* end if (population[1][0] > bfgs) */
	    case 1:
	      if (population[1][0] < bfgsfit) /* maximize */
		{
		  /* is the BFGS individual in the bounds? */
		  BoundaryTrigger=0; /* outside of bounds ? */
		  for (i=0; i<nvars; i++) {
		    j = i+1;
		    if (bfgsoutX[i] < domains[j][1]) {
		      BoundaryTrigger=1;
		      if (PrintLevel>1)
			{
			  Rprintf(
				  "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
			  Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
				  count_gener, i+1, bfgsoutX[i]);
			}
		      warning("BFGS hit on best individual produced Out of Boundary individual.");
		    }
		    if (bfgsoutX[i] > domains[j][3]) {
		      BoundaryTrigger=1;
		      if (PrintLevel>1)
			{
			  Rprintf(
				  "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
			  Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
				  count_gener, i+1, bfgsoutX[i]);
			}
		      warning("BFGS hit on best individual produced Out of Boundary individual.");
		    }
		  } /* end for loop */
	      
		  /* if we we use out of bounds individuals then proceed */
		  /* 0=anything goes, 1: regular; 2: no trespassing! */
		  if (BoundaryEnforcement==0) {
		    for(i=1;i<=nvars;i++) 
		      {
			population[1][i]=bfgsoutX[i-1];
		      }
		    population[1][0]=bfgsfit;
		  }
		  else if (BoundaryTrigger==0) {
		    for(i=1;i<=nvars;i++) 
		      {
			population[1][i]=bfgsoutX[i-1];
		      }
		    population[1][0]=bfgsfit;
		  }
		} /* end if (population[1][0] < bfgsfit) */
	    } /* end switch */
	  }/*end of NOT lexical bfgs hit */
        else if(Structure->whichFUN == 2) // Lexical but not Transform
	  {
	    /* is the BFGS individual in the bounds? */
	    BoundaryTrigger=0; /* outside of bounds ? */
	    for (i=0; i<nvars; i++) {
	      j = i+1;
	      if (bfgsoutX[i] < domains[j][1]) {
		BoundaryTrigger=1;
		if (PrintLevel>1)
		  {
		    Rprintf(
			    "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
		    Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n\n", 
			    count_gener, i+1, bfgsoutX[i]);
		    Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
		  }
		warning("BFGS hit on best individual produced Out of Boundary individual.");
	      }
	      if (bfgsoutX[i] > domains[j][3]) {
		BoundaryTrigger=1;
		if (PrintLevel>1)
		  {
		    Rprintf(
			    "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
		    Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n", 
			    count_gener, i+1, bfgsoutX[i]);
		    Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
		  }
		warning("BFGS hit on best individual produced Out of Boundary individual.");
	      }
	    } /* end for loop */
		  
	    /* if we use out of bounds individuals then proceed */
	    /* 0=anything goes, 1: regular; 2: no trespassing! */
	    /* Add new individual to the END because we are doing lexical stuff */
	    if (BoundaryEnforcement==0) {
	      for(i=1;i<=nvars;i++) 
		{
		  population[(pop_size-1)][i]=bfgsoutX[i-1];
		  X[i] = bfgsoutX[i-1];
		}
	      EvaluateLexical(Structure->fn, Structure->rho, 
			      X, nvars, Structure->Lexical, MinMax, LexicalReturn);
	      population[(pop_size-1)][0] = LexicalReturn[0];
	      count = 0;
	      for(i=(nvars+2);i<lexical_end;i++)
		{
		  count++;
		  population[(pop_size-1)][i] = LexicalReturn[count];
		}		      			  

	      /* REDO SORT.  This is inefficient becase we only changed 1 individual*/
	      RlexicalSort(Structure->fnLexicalSort, Structure->rho,
			   population,
			   MinMax, pop_size, nvars, lexical_end, 1);
	    }
	    else if (BoundaryTrigger==0) {
	      for(i=1;i<=nvars;i++) 
		{
		  population[(pop_size-1)][i]=bfgsoutX[i-1];
		  X[i] = bfgsoutX[i-1];
		}
	      EvaluateLexical(Structure->fn, Structure->rho, 
			      X, nvars, Structure->Lexical, MinMax, LexicalReturn);
	      population[(pop_size-1)][0] = LexicalReturn[0];
	      count = 0;
	      for(i=(nvars+2);i<lexical_end;i++)
		{
		  count++;
		  population[(pop_size-1)][i] = LexicalReturn[count];
		}		      			  

	      /* REDO SORT.  This is inefficient becase we only changed 1 individual*/
	      RlexicalSort(Structure->fnLexicalSort, Structure->rho,
			   population,
			   MinMax, pop_size, nvars, lexical_end, 1);
	    }
	  } /*end of LEXICAL bfgs hit */
        else // Transform
          {
            for(j=1; j<=nvars; j++)
              {
                X[j] = bfgsoutX[(j-1)];
              }
            // Call EvaluateTransform now to transform X
            EvaluateTransform(Structure->fn, Structure->rho,
                              X, nvars, Structure->Lexical, MinMax, LexicalReturn);

            /* is the BFGS individual in the bounds? */
            BoundaryTrigger=0; /* outside of bounds ? */
            for (j=1; i<=nvars; i++) {
              if (X[j] < domains[j][1]) {
                BoundaryTrigger=1;
                if (PrintLevel>1)
                  {
                    Rprintf(
                            "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
                    Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n\n",
                            count_gener, j, X[j]);
                    Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
                  }
                warning("BFGS hit on best individual produced Out of Boundary individual.");
              }
              if (X[j] > domains[j][3]) {
                BoundaryTrigger=1;
                if (PrintLevel>1)
                  {
                    Rprintf(
                            "\nNOTE: BFGS hit on best individual produced Out of Boundary individual.\n");
                    Rprintf("NOTE: Generation: %d \t Parameter: %d \t Value: %e\n",
                            count_gener, j, X[j]);
                    Rprintf("NOTE: Fit: %e\n\n", bfgsfit);
                  }
                warning("BFGS hit on best individual produced Out of Boundary individual.");
              }
            } /* end for loop */

            /* if we use out of bounds individuals then proceed */
            /* 0=anything goes, 1: regular; 2: no trespassing! */
            /* Figure out why BoundaryEnforcement==0 and BoundaryTrigerr==0 are separate above */
            if (BoundaryEnforcement == 0 || BoundaryTrigger == 0) {
              if(Structure->Lexical < 2) // Transform but not Lexical
                {
                  for(i=1; i<=nvars; i++)
                    {
                      population[1][i] = X[i];
                    }
                  population[1][0] = LexicalReturn[0]; // from ~45 lines above
                }
              else /* Add new individual to the END because we are doing Transform and Lexical */
                {
                  for(i=1;i<=nvars;i++)
                    {
                      population[(pop_size-1)][i] = X[i];
                    }
                  population[(pop_size-1)][0] = LexicalReturn[0]; // from ~55 lines above
                  count = 0;
                  for(i=(nvars+2);i<lexical_end;i++)
                    {
                      count++;
                      population[(pop_size-1)][i] = LexicalReturn[count];
                    }

                    /* REDO SORT.  This is inefficient becase we only changed 1 individual*/
                    RlexicalSort(Structure->fnLexicalSort, Structure->rho,
                                 population,
                                 MinMax, pop_size, nvars, lexical_end, 1);
                }
            } /* end of boundary enforcment */
          } /*end of TRANSFORM bfgs hit */
      } /* end of UseBFGS */  

      /* check to see if fit is improving */
      if(Structure->Lexical < 2)
	{
	  switch(MinMax)
	    {
	    case 0:
	      if ( (oldfitvalue - SolutionTolerance) > population[1][0]) {
		nochange_gen=0;
		oldfitvalue=population[1][0];
		peak_cnt = count_gener;
	      }
	      else nochange_gen++;
	      break;
	    case 1:
	      if ( (oldfitvalue + SolutionTolerance) < population[1][0]) {
		nochange_gen=0;
		oldfitvalue=population[1][0];
		peak_cnt = count_gener;
	      }
	      else nochange_gen++;	      
	      break;
	    }
	} /*       if(Structure->Lexical < 2) */
      else
	{
	  switch(MinMax)
	    {
	    case 0:
	      LexicalFitsImproving = 0;
	      if ( (oldfitvalue - SolutionTolerance) > population[1][0]) 
		{
		  LexicalFitsImproving = 1;
		} 
	      else
		{
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    if ( (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) {
		      LexicalFitsImproving = 1;
		      break;
		    } /* (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) */
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		} /* else if ( (oldfitvalue - SolutionTolerance) > population[1][0])  */
	      if (LexicalFitsImproving)
		{
		  nochange_gen = 0;
		  peak_cnt = count_gener;
		  oldfitvalue=population[1][0];
		  oldfitvalueVEC[0]=population[1][0];
		  k = 1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    oldfitvalueVEC[k]=population[1][i];
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		}
	      else
		nochange_gen++;
	      break;
	    case 1:
	      LexicalFitsImproving = 0;
	      if ( (oldfitvalue + SolutionTolerance) < population[1][0]) 
		{
		  LexicalFitsImproving = 1;
		} 
	      else
		{
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    if ( (oldfitvalueVEC[k] + SolutionTolerance) < population[1][i] ) {
		      LexicalFitsImproving = 1;
		      break;
		    } /* (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) */
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		} /* else if ( (oldfitvalue - SolutionTolerance) > population[1][0])  */
	      if (LexicalFitsImproving)
		{
		  nochange_gen = 0;
		  peak_cnt = count_gener;
		  oldfitvalue=population[1][0];
		  oldfitvalueVEC[0]=population[1][0];
		  k = 1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    oldfitvalueVEC[k]=population[1][i];
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		}
	      else
		nochange_gen++;
	      break;	      
	    } /* switch(MinMax) */
	} /* else (Structure->Lexical > 2) */

      if(PrintLevel==1)
	{
	  if( nochange_gen==0 )
	    {
	      if(Structure->Lexical > 1)
		{
		  Rprintf("\n%7d \t%e  ", count_gener, population[1][0]);
		  for(j=(nvars+2);j<lexical_end;j++)
		    {
		      Rprintf("%e  ", population[1][j]);
		    }		      
		  Rprintf("\n");	
		}
	      else
		{
		  Rprintf("%7d \t%e\n",
			  count_gener,population[1][0]); 
		  /* fflush(output); */
		}
	    } 
	}

      /* compute and print mean and variance of population */
      if (PrintLevel>1) {
	Rprintf("\nGENERATION: %d\n", count_gener);
	populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);

	if(Structure->Lexical > 1)
	  {
	    Rprintf( "Lexical Fit..... %e  ", population[1][0]);
	    for(j=(nvars+2);j<lexical_end;j++)
	      {
		Rprintf("%e  ", population[1][j]);
	      }		      
	    Rprintf("\n");	    
	  }
	else
	  {
	    Rprintf( "Fitness value... %e\n", population[1][0]);
	    Rprintf( "mean............ %e\n", mean[0]);
	    Rprintf( "variance........ %e\n", var[0]);
	    /*
	      Rprintf( "skewness........ %e\n", skew[i]);
	      Rprintf( "kurtosis........ %e\n", kur[i]);
	    */
	  }

	nnull = pop_size-tobs[0];
	if(nnull > 0)
	  Rprintf( "#null........... %d\n", nnull);
	if(Structure->MemoryUsage==1)
	  Rprintf( "#unique......... %d, #Total UniqueCount: %d\n", 
		  UniqueCount-OldUniqueCount, UniqueCount);
	/* Rprintf( "tobs............ %d\n", tobs[i]); */

	for (i=1; i<=nvars; i++) {
	  Rprintf( "var %d:\n", i);
	  Rprintf( "best............ %e\n", population[1][i]);
	  Rprintf( "mean............ %e\n", mean[i]);
	  Rprintf( "variance........ %e\n", var[i]);
	  /*
	    Rprintf( "skewness........ %e\n", skew[i]);
	    Rprintf( "kurtosis........ %e\n", kur[i]);
	  */
	  nnull = pop_size-tobs[i];
	  if(nnull > 0)
	    Rprintf( "#null........... %d\n", nnull);
	  /* Rprintf( "tobs............ %d\n", tobs[i]); */
	}
      } /* end of printlevel if */
      
      /* if (PrintLevel>0)
	 fflush(output); */

      /* Print the population file */
      if ( PrintLevel == 1 ) {
	if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	  Rprintf("Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	      JaMatrixFree(Memory, MemorySize);

	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 0, nvars+1, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	  free_ivector(parents, 1);

          if(Structure->Lexical > 1 || Structure->Transform == 1)
	    {
	      free(LexicalReturn);
	      free(oldfitvalueVEC);
	    }

	  error("Fatal Error. See output for diagnostic information.");
	}
	print_population((int) pop_size, (int) nvars, (int) count_gener, (int) Structure->Lexical, population, popout);
	fclose(popout);
      } /* end of PrintLevel if */
      if ( PrintLevel>1) {
	if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	  Rprintf("Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);

	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 0, nvars+1, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	  free_ivector(parents, 1);

          if(Structure->Lexical > 1 || Structure->Transform == 1)
	    {
	      free(LexicalReturn);
	      free(oldfitvalueVEC);
	    }

	  error("Fatal Error. See output for diagnostic information.");
	}
	print_population((int) pop_size, (int) nvars, (int) count_gener, (int) Structure->Lexical, population, popout);
	fflush(popout);
	fclose(popout);
      }

      if (nochange_gen > (WaitGenerations)) {
	/* increase the number of WaitGenerations if the gradients are NOT zero! */	  
	if (GradientCheck==0) {
	  if(PrintLevel>0)	  
	    {
	      Rprintf("\n'wait.generations' limit reached.\n");
	      Rprintf("No significant improvement in %d generations.\n", nochange_gen-1);
	      /* fflush(output); */
	    }
	  MaxGenerations = 0;
	  nochange_gen=0;
	}
	else  
	  {
	    for (i=1; i<=nvars; i++)
	      {
		bfgsoutX[i-1]=population[1][i];
	      }
	    if(Structure->UserGradient==0)
	      {	    
		gradient(Structure->fn, Structure->rho,
			 bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, domains);
	      } 
	    else 
	      {
		userGradientfn(Structure->fnGR, Structure->rho, bfgsoutX, grad, nvars);
	      }
	    GradientTrigger = 0;
	    for (i=0; i<nvars; i++) {
	      if (fabs(grad[i]) > SolutionTolerance) {
		GradientTrigger = 1;
		break;
	      }
	    } /* end for loop */
	    if (GradientTrigger==1) {
	      IncreaseGenerations = WaitGenerations;
	      WaitGenerations += IncreaseGenerations;
              if(Structure->BFGSburnin < 0) Structure->BFGSburnin = 0;
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nDoubling 'wait.generations' limit to %d (from %d) ", 
			  WaitGenerations, IncreaseGenerations);
		  Rprintf("because at least one gradient is too large.\n");
		  Rprintf("G[%d]: %e\t Solution Tolerance: %e\n\n", 
			  i+1, grad[i], SolutionTolerance);
		}
	    }
	    else {
	      if(PrintLevel>0)	  
		{
		  Rprintf("\n'wait.generations' limit reached.\n");
		  Rprintf("No significant improvement in %d generations.\n", nochange_gen-1);
		  /* fflush(output); */
		}
	      MaxGenerations = 0;
	      nochange_gen=0;
	    }
	  }/* end else loop */
      } /* end of if (nochange_gen > (WaitGenerations)) { */

      if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 
	{
	  if (HardGenerationLimit==0)
	    {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += IncreaseGenerations;
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because at least one gradient is too large.\n\n");
		}
	    } // if (Structure->HardGenerationLimit==0)
	  else
	    {
	      HardMaximumNumber = 1;
	      warning("Stopped because hard maximum generation limit was hit.\nAt least one gradient is too large.");
	      if(PrintLevel>0)	  
		{
		  Rprintf("\nNOTE: HARD MAXIMUM GENERATION LIMIT HIT\n");
		  Rprintf("        At least one gradient is too large\n");
		}
	    } // else
	} // if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 


      /* increase the number of generations if fitness has been improving */
      if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) ) {
	if (HardGenerationLimit==0)
	  {
	    if (WaitGenerations > MaxGenerations) {
	      IncreaseGenerations = WaitGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because the fitness is still impoving.\n\n");
		}
	    }
	    else {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because the fitness is still improving.\n\n");
		}
	    }
	  } // if (Structure->HardGenerationLimit==0)
	else
	  {
	    if (HardMaximumNumber==0)
	      {
		warning("Stopped because hard maximum generation limit was hit.");
		if(PrintLevel>0)	  
		  {
		    Rprintf("\nNOTE: HARD MAXIMUM GENERATION LIMIT HIT\n");
		  }
	      } /* end of if HardMax */
	  }
      } // if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) )
      
      /* if(PrintLevel>0)	  
	 fflush(output); */
      
    } /* end of do loop */
  /*Increment iteration count and test whether all generations are done*/
  while (++count_gener <= MaxGenerations);

  if(PrintLevel>0)	    
    {
      if(Structure->Lexical > 1)
	{
	  Rprintf("\nSolution Lexical Fitness Value:\n");
	  Rprintf("%e  ", population[1][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[1][j]);
	    }		      
	  Rprintf("\n");
	}
      else
	{
	  Rprintf("\nSolution Fitness Value: %e\n", population[1][0]);
	}
      if (GradientCheck==0 && UseBFGS==0)
	Rprintf("\nParameters at the Solution:\n\n");
      else
	Rprintf("\nParameters at the Solution (parameter, gradient):\n\n");
    }

  /* output data structure */
  Structure->oPeakGeneration=peak_cnt;
  Structure->oGenerations=count_gener-1;

  /* obtain gradients */
  /* print best solution */
  if (GradientCheck==0 && UseBFGS==0)
    {
      for(j = 1; j <= nvars; j++) {
	i = j-1;
	if(PrintLevel>0)	  
	  Rprintf(" X[%2d] :\t%e\n",j,population[1][j]);
	grad[i] = -1.0;
	Results[i] = population[1][j];
	Gradients[i] = grad[i];
      }
    } /* end of if (GradientCheck==0 && UseBFGS==0) */
  else 
    {
      for (i=1; i<=nvars; i++)
	{
	  bfgsoutX[i-1]=population[1][i];
	}
      if(Structure->UserGradient==0)
	{
	  gradient(Structure->fn, Structure->rho,
		   bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, domains);
	} 
      else 
	{
	  userGradientfn(Structure->fnGR, Structure->rho, bfgsoutX, grad, nvars);
	}

      for(j = 1; j <= nvars; j++) {
	i = j-1;
	if(PrintLevel>0)	  
	  Rprintf(" X[%2d] :\t%e\tG[%2d] :\t%e\n",j,population[1][j],j,grad[i]);
	Results[i] = population[1][j];
	Gradients[i] = grad[i];
      }
    } /* end of  else (GradientCheck==0 && UseBFGS==0) */

  Structure->oFitValues[0]=population[1][0];
  if (Structure->Lexical > 1)
    {
      k = 1;
      for (i=(nvars+2);i<lexical_end;i++)  {
	Structure->oFitValues[k]=population[1][i];
	k++;  	  
      }
    } 

  /* free memory */
  /* free populationstats stuff */
  free(mean);
  free(var);
  free(skew);
  free(kur);
  free(tobs);
  
  free(bfgsoutX);
  free(finalhessin);
  free(evalX);
  free(grad);

  if (Structure->MemoryUsage==1)
    JaMatrixFree(Memory, MemorySize);

  JaMatrixFree(population, pop_size+2);
  JaMatrixFree(new_genera, pop_size+2);

  free_matrix(temp, 0, nvars+1, 0);
  free_vector(probab, 1);
  free_vector(t_vec, 1);
  free_vector(cum_probab, 1);
  free_ivector(live, 1);
  free_ivector(parents, 1);

  if(Structure->Lexical > 1 || Structure->Transform == 1)
    {
      free(LexicalReturn);
      free(oldfitvalueVEC);
    }
} /* end optimiztion function */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   sort()                                       */
/*                                                                              */
/*           SYNOPSIS          :   void sort(MinMax, population,pop_size,       */
/*                                 variable)                                    */
/*                                                                              */
/*           DESCRIPTION       :   This function sorts the population, in the   */
/*                                  ascending or the descending order of the    */
/*                                  evaluation function, depending on whether   */
/*                                  it is a maximization or a minimization      */
/*                                  function, respectively.                     */
/*                                                                              */
/*                                  As an alternative, the sortq function below */
/*                                  can be used, That sorting function uses     */
/*                                  the quicksort algorithm.                    */
/*                                                                              */
/*                                                                              */
/*           FUNCTIONS CALLED  :   swap()                                       */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/
void sort(short int MinMax, MATRIX  population, int pop_size,
	  long nvar)
     /*
       short int MinMax;      Tells whether it is a maximizaiton or a minimization function
       int pop_size;          Population size
       MATRIX population;     Array of population
     */
{
  int i,j;


  /*If MinMax is 0 sorts in the descending order, and*/
  /*if it is 1 sorts in the ascending order*/
  /*Sorted in ascending or descending order, based on*/
  /*the evaluated values of each of the agents*/
  switch(MinMax)
    {
    case 0 :
      for(i=1; i<=pop_size; i++)
        for(j=i+1; j<=pop_size; j++)
          if(population[i][nvar] > population[j][nvar])
            swap(&population[i],&population[j]);
      break;

    case 1 :
      for(i=1; i<=pop_size; i++)
        for(j=i+1; j<=pop_size; j++)
          if(population[i][nvar] < population[j][nvar])
            swap(&population[i],&population[j]);
      break;
    }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   swap()                                       */
/*                                                                              */
/*           SYNOPSIS          :   void swap(x,y)                               */
/*                                                                              */
/*           DESCRIPTION       :   This function interchanges the values of     */
/*                                  x and y.                                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   sort()                                       */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



void swap(double **x, double **y)
     /* double **x,**y; */
{
  double *temp;

  temp = *x;
  *x = *y;
  *y = temp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_parent()                                */
/*                                                                              */
/*           SYNOPSIS          :   int find_parent(live,pop_size)               */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the index of the       */
/*                                  agent in the population, which is to be     */
/*                                  chosen for reproduction.                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   irange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

int find_parent(IVECTOR live, int pop_size)
     /*
       int pop_size;    Population size
       IVECTOR live;    Vector containing the number of times each agent
                        is going to reproduce
     */
{
  int i,temp,t1=0,tot=0;

  /*Finding the total number of parents to reproduce*/
  for(i=1; i<=pop_size; i++)
    tot = tot + live[i];
  if(tot==0)
    {
      error("No agents to select");
    }

  /*Choosing one of them randomly*/
  temp = irange_ran(1,tot);

  tot = 0;
  i = 1;
  do{
    if(live[i]!=0)
      t1 = i;
    tot = tot + live[i++];
  }while(tot<temp);

  /*Decrementing the number of times the parent chosen is going to reproduce*/
  // live[t1]--;
  return(t1);
}
/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   assign_probab()                              */
/*                                                                              */
/*           SYNOPSIS          :   void assign_probab(probab,pop_size,Q)        */
/*                                                                              */
/*           DESCRIPTION       :   This function assigns probability of survival*/
/*                                  to each of the agents determined by the     */
/*                                  value provided by the user for the          */
/*                                  probability of the best agnet.              */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void assign_probab(VECTOR probab, int pop_size, double Q)
     /*
       int pop_size;     Population size
       double Q;         The probability of survival of the best agent
       VECTOR probab;    Array to contain the probability of survival
                         of each of the agents
     */
{
  int i;

  /* Q, Q(1-Q)^1, Q(1-Q)^2 ... Q(1-Q)^n */
  for(i=1; i<=pop_size; i++)
    probab[i] = Q * x_pow_y(1-Q,i-1);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   x_pow_y()                                    */
/*                                                                              */
/*           SYNOPSIS          :   double x_pow_y(x,y)                           */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the value of x to the  */
/*                                  power of y.                                 */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   assign_probab()                              */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



double x_pow_y(double x, int y)
     /*
       double x;
       int y;
     */
{
  int i;
  double tot = 1.0;

  for(i=0; i < y; i++)
    tot = tot * x;
  return(tot);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_cum__probab()                           */
/*                                                                              */
/*           SYNOPSIS          :   void find_cum__probab(cum_probab,probab,     */
/*                                                                     pop_size)*/
/*                                                                              */
/*           DESCRIPTION       :   This function finds the cumulative           */
/*                                  probability of each of the agents, from the */
/*                                  individual probability found earlier.       */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/



void find_cum_probab(VECTOR cum_probab, VECTOR probab, int pop_size)
     /*
       int pop_size;     Population size
       VECTOR probab,      Individual probability of survival of each of the agent
       cum_probab;         Cumulative probability of survival of each of the agent
     */
{
  int i;

  cum_probab[1] = probab[1];

  for(i=2; i<=pop_size; i++)
    cum_probab[i] = cum_probab[i-1] + probab[i];

}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_live()                                  */
/*                                                                              */
/*           SYNOPSIS          :   void find_live(cum_probab,live,pop_size,P4+P5*/
/*                                                                              */
/*           DESCRIPTION       :   This function finds the agents from the      */
/*                                  population, who are going to live - those   */
/*                                  who are going to reproduce, which is done   */
/*                                  based on the cumulative probability of      */
/*                                  survival of each of the agents.             */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/********************************************************************************/


void find_live(VECTOR cum_probab, IVECTOR live, int pop_size, int P)
     /*
       VECTOR cum_probab;  Cumulative probability
       IVECTOR live;       Agents that are going to reproduce
       int pop_size,       Population size
       P;                  Total number of parents needed to reproduce
     */
{
  double random;
  int count=0,/*Count of the number of agents chosen to live*/
      i;

  do
    {
      /*Choosing a random cumulative probability*/
      random = frange_ran(0.0,1.0);
      i=0;
      /*Finding the agent with the chosen cumulative probability*/
      do{
        i++;
        }while((random > cum_probab[i]) && (i< pop_size));

      /*Chosing the parent with that probability to reproduce*/
      if(count < P)
        {
          live[i]++;
          count++;
        }
    }while(count < P);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_die()                                   */
/*                                                                              */
/*           SYNOPSIS          :   void find_die(cum_probab,die,pop_size,P4+P5) */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the agents from the      */
/*                                  population, who are going to die.           */
/*                                                                              */
/*           FUNCTIONS CALLED  :   frange_ran()                                 */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/********************************************************************************/



int find_die(VECTOR cum_probab, IVECTOR die, int pop_size)
     /*
       VECTOR cum_probab; Cumulative probability
       IVECTOR die;       Agents that are going to die
       int pop_size;      Population size
     */
{
  double random;
  int i;
  int done = FALSE;

  do
    {
      /*Choosing a random cumulative probability*/
      random = frange_ran(0.0,1.0);
      i=0;
      /*Finding the agent with the chosen cumulative probability*/
      do{
        i++;
        }
      while((random > cum_probab[i]) && (i< pop_size));

      /*Chosing the agent to die*/
      if ((die[pop_size+1-i] == 0) && (i < pop_size))
        done = TRUE;
    }
  while(!done);
  return(pop_size+1-i);
}



void SetRunTimeParameters(struct GND_IOstructure *Structure, 
			  short FirstTime,
			  long *PopSize, long *nvars, long *MaxGenerations, long *WaitGenerations,
			  short *MinMax, short *GradientCheck, short *BoundaryEnforcement, short *UseBFGS,
			  double *SolutionTolerance,
			  long *InstanceNumber, long *P, long *P0, long *P1, long *P2, long *P3, long *P4, long *P5, 
			  long *P6, long *P7, long *P8, short *PrintLevel, 
			  short *HardGenerationLimit)
{
  double tP;
  int i;

  *PopSize=Structure->PopSize;
  *nvars=Structure->nvars;
  *MaxGenerations=Structure->MaxGenerations;
  *WaitGenerations=Structure->WaitGenerations;

  if (FirstTime==1)
    *HardGenerationLimit=Structure->HardGenerationLimit;

  *MinMax=Structure->MinMax;
  *GradientCheck=Structure->GradientCheck;
  *BoundaryEnforcement=Structure->BoundaryEnforcement;
  *UseBFGS=Structure->UseBFGS;
  *InstanceNumber=Structure->InstanceNumber;

  *SolutionTolerance=Structure->SolutionTolerance;
  *PrintLevel=Structure->PrintLevel;

  /* Check to make sure that all operators are positve numbers! */
  if (Structure->P[0] < 0 ) {
    warning("Operator 1 (Cloning) was Assigned an Illegal Value: %d.", Structure->P[0]);
    Structure->P[0]=0.0;
  }
  if (Structure->P[1] < 0 ) {
    warning("Operator 1 (Uniform Mutation) was Assigned an Illegal Value: %d.", Structure->P[1]);
    Structure->P[1]=0.0;
  }
  if (Structure->P[2] < 0 ) {
    warning("Operator 3 (Boundary Mutation) was Assigned an Illegal Value: %d.", Structure->P[2]);
    Structure->P[2]=0;
  }
  if (Structure->P[3] < 0 ) {
    warning("Operator 4 (Non-Uniform Mutation) was Assigned an Illegal Value: %d.", 
	    Structure->P[3]);
    Structure->P[3]=0;
  }
  if (Structure->P[4] < 0 ) {
    warning("Operator 5 (Polytope Crossover) was Assigned an Illegal Value: %d.", 
	    Structure->P[4]);
    Structure->P[4]=0;
  }
  if (Structure->P[5] < 0 ) {
    warning("Operator 6 (Simple Crossover) was Assigned an Illegal Value: %d.", 
	    Structure->P[5]);
    Structure->P[5]=0;
  }
  if (Structure->P[6] < 0 ) {
    warning("Operator 7 (Whole Non-Uniform Mutation) was Assigned an Illegal Value: %d.", 
	    Structure->P[6]);
    Structure->P[6]=0;
  }
  if (Structure->P[7] < 0 ) {
    warning("Operator 8 (Heuristic Crossover) was Assigned an Illegal Value: %d.", 
	    Structure->P[7]);
    Structure->P[7]=0;
  }

  /* if DataType==1 (i.e., integer) we are not giong to use any gradient information etc. */
  if (Structure->DataType==1) {
    *UseBFGS=0;
    *GradientCheck=0;

    if (Structure->P[8] > 0) {
      warning("Operator 9 (Local-Minimum Crossover) was Assigned an Illegal Value: %d\nThis is an illegal value because we are working with integer data.", 
	      Structure->P[8]);
      Structure->P[8]=0;
    } /* end of if */
  }
  else {
    if (Structure->P[8] < 0 ) {
      warning("Operator 9 (Local-Minimum Crossover) was Assigned an Illegal Value: %d.", 
	      Structure->P[8]);
      Structure->P[8]=0;
    }
  } /* end of else */

  /* Let's figure out the number of operators we need.  Move stuff to absolute space */
  tP = 0;
  for (i=0; i<9; i++) {
    tP = tP + Structure->P[i] ;
  }
  
  if (tP > 0) {
    *P0 = Iround(  (Structure->P[0] /  tP) * (double) (*PopSize-2) );
    *P1 = Iround(  (Structure->P[1] /  tP) * (*PopSize-2) );
    *P2 = Iround(  (Structure->P[2] /  tP) * (*PopSize-2) );
    *P3 = Iround(  (Structure->P[3] /  tP) * (*PopSize-2) );
    *P4 = Iround(  (Structure->P[4] /  tP) * (*PopSize-2) );
    *P5 = Iround(  (Structure->P[5] /  tP) * (*PopSize-2) );
    *P6 = Iround(  (Structure->P[6] /  tP) * (*PopSize-2) );
    *P7 = Iround(  (Structure->P[7] /  tP) * (*PopSize-2) );
    *P8 = Iround(  (Structure->P[8] /  tP) * (*PopSize-2) );
  }
  else {
    *P0 = 0;
    *P1 = 0;
    *P2 = 0;
    *P3 = 0;
    *P4 = 0;
    *P5 = 0;
    *P6 = 0;
    *P7 = 0;
    *P8 = 0;
  }

  /* Check to make sure that all operators (i.e., 5, 7) which have to be even numbers are */
  if (fmod((long double) *P5, (long double) 2) > 0.0) {
    if(Structure->PrintLevel>2)
      {
	Rprintf("\nNOTE: Operator 6 (Simple Crossover) may only be started\n");
	Rprintf("NOTE: an even number of times.  I am increasing this operator by one.\n");
      }
    *P5=*P5+1;
  }
  if (fmod((long double) *P7, (long double) 2) > 0.0) {
    if(Structure->PrintLevel>2)
      {
	Rprintf("\nNOTE: Operator 8 (Heuristic Crossover) may only be started\n");
	Rprintf("NOTE: an even number of times.  I am increasing this operator by one.\n");
      }
    *P7=*P7+1;
  }

  /*P is the total number of parents needed for applying all the operators*/
  *P = *P1 + *P2 + *P3 + *P4 + *P5 + *P6 + *P7 + *P8;
  if(*P > *PopSize)
    {
      if(Structure->PrintLevel>0)
	{
	  Rprintf("\nNOTE: The total number of operators greater than population size\n");
	}

      if (fmod((long double) *P+1, (long double) 2) > 0.0) {
	*PopSize = *P+2;
	if(Structure->PrintLevel>0)
	  {
	    Rprintf("NOTE: I'm increasing the population size to %d (operators+2).\n", *PopSize);
	  }
      }
      else {
	*PopSize = *P+1;
	if(Structure->PrintLevel>0)
	  {
	    Rprintf("NOTE: I'm increasing the population size to %d (operators+1).\n", *PopSize);
	  }
      }
    }
  else if ( *P== *PopSize) {
    if(Structure->PrintLevel>0)
      {
	Rprintf("\nNOTE: The total number of operators equal to the population size\n");
      }

    if (fmod( (long double) *P+1, (long double) 2) > 0.0) {
	*PopSize = *P+2;
	if(Structure->PrintLevel>0)
	  {
	    Rprintf("NOTE: I'm increasing the population size to %d (operators+2).\n", *PopSize);
	  }
      }
      else {
	*PopSize = *P+1;
	if(Structure->PrintLevel>0)
	  {
	    Rprintf("NOTE: I'm increasing the population size to %d (operators+1).\n", *PopSize);
	  }
      }
  }

  if (fmod( (long double) *PopSize, (long double) 2) > 0.0) {
    if(Structure->PrintLevel>0)
      {
	Rprintf("NOTE: population size is not an even number\n");
	Rprintf("      increasing population size by 1\n");
      }
    *PopSize=*PopSize+1;
  }

  /* return PopSize and P values to the return data structure */
  *P0 = *PopSize-*P-1;
  Structure->oP[0]=*P0;
  Structure->oP[1]=*P1;
  Structure->oP[2]=*P2;
  Structure->oP[3]=*P3;
  Structure->oP[4]=*P4;
  Structure->oP[5]=*P5;
  Structure->oP[6]=*P6;
  Structure->oP[7]=*P7;
  Structure->oP[8]=*P8;
  Structure->oPopSize=*PopSize;

  if(Structure->PrintLevel>0)
    {
      Rprintf( "\n");
      if (Structure->DataType==1) Rprintf( "Data Type: Integer\n");
      else Rprintf( "Data Type: Floating Point\n");
      
      Rprintf("Operators (code number, name, population) \n");
      Rprintf("\t(1) Cloning........................... \t%d\n", *P0);
      Rprintf("\t(2) Uniform Mutation.................. \t%d\n", *P1);
      Rprintf("\t(3) Boundary Mutation................. \t%d\n", *P2);
      Rprintf("\t(4) Non-Uniform Mutation.............. \t%d\n", *P3);
      Rprintf("\t(5) Polytope Crossover................ \t%d\n", *P4);
      Rprintf("\t(6) Simple Crossover.................. \t%d\n", *P5);
      Rprintf("\t(7) Whole Non-Uniform Mutation........ \t%d\n", *P6);
      Rprintf("\t(8) Heuristic Crossover............... \t%d\n", *P7);
      Rprintf("\t(9) Local-Minimum Crossover........... \t%d\n\n", *P8);
      if (*HardGenerationLimit==0)
	Rprintf("SOFT Maximum Number of Generations: %lu\n", *MaxGenerations);
      else
	Rprintf("HARD Maximum Number of Generations: %lu\n", *MaxGenerations);
      Rprintf("Maximum Nonchanging Generations: %lu\n", *WaitGenerations);
      Rprintf("Population size       : %d\n", *PopSize);
      Rprintf("Convergence Tolerance: %e\n", *SolutionTolerance);
      
      Rprintf( "\n");
      if (*UseBFGS !=0) {
	Rprintf(
		"Using the BFGS Derivative Based Optimizer on the Best Individual Each Generation.\n");
      }
      else {
	Rprintf(
		"Not Using the BFGS Derivative Based Optimizer on the Best Individual Each Generation.\n");
      }
      if (*GradientCheck==0)
	Rprintf("Not Checking Gradients before Stopping.\n");
      else 
	Rprintf("Checking Gradients before Stopping.\n");
      
      if (*BoundaryEnforcement==0) 
	Rprintf("Using Out of Bounds Individuals.\n\n");
      else if (*BoundaryEnforcement==1) 
	Rprintf("Not Using Out of Bounds Individuals But Allowing Trespassing.\n\n");
      else if (*BoundaryEnforcement==2) 
	Rprintf("Not Using Out of Bounds Individuals and Not Allowing Trespassing.\n\n");
    }
    
  /* if(Structure->PrintLevel>0)
     fflush(output); */

} /* End SetOperators */


/********************************************************************************/
/*  JaIntegerOptimization:                                                      */
/*                                                                              */
/*  This function assumes that the X variables are integers.                    */
/*                                                                              */
/*  The cross over operators are different!                                     */
/*                                                                              */
/********************************************************************************/

void JaIntegerOptimization(struct GND_IOstructure *Structure, VECTOR X, 
			     MATRIX domains)
{
  extern struct GND_IOstructure *ExternStructure;
  
  MATRIX new_genera,   /*Temporary storage for the new generation*/
         population,   /*Population of x2 variables*/
         temp;

  VECTOR probab,       /*Probability of agents to die or live*/
         cum_probab,   /*Cumilative probability of agents*/
         t_vec;

  IVECTOR live;
  /* for oper4 */
  IVECTOR parents;


  long count_gener= 1; /*Counter to keep track of the number of generations*/
  unsigned long peak_cnt;

  int                     /*Total number of agents chosen to reproduce*/
    j1,
    j2,
    j3,
    j4,
    j5,
    j6,
    j7,
    j8,
    oper,
    ocnt,
    B,                     /*Parameter for the 3rd operator - nonuniform mutation*/
    STEP,                  /*Parameter for the 5th operator - simple arithmetical crossover*/
    first_live=0,          /*Index of the two parents for crossover parents*/
    second_live=0,
    first_die,             /*Index of the two parents for crossover death*/
    second_die,
    die_now;               /*index of agent to replace in current operation*/

  long i,j, s, k;

  /* for oper 4 */
  int p2use;


  double Q;                   /*Probability of the best agent*/
  FLAG  same;
  double **Jnew;

  double *grad, *evalX, *finalhessin, *bfgsoutX;

  int nochange_gen=0;

  double oldfitvalue=0;

  int IncreaseGenerations;
  short int GradientTrigger=0;
  long InstanceNumber;

  long nvars, MaxGenerations, WaitGenerations, count;
  long pop_size, P, P0, P1, P2, P3, P4, P5, P6, P7, P8;
  short int MinMax, GradientCheck, BoundaryEnforcement, UseBFGS, HardMaximumNumber=0;
  double SolutionTolerance, *Results, *Gradients;
  short PrintLevel, HardGenerationLimit;

  /* Old variables which may change when SetRunTimeParameters is run during a run! */
  long pop_size_old;

  /* Summary Statistics (mean, variance etc) */
  /* double popmean, popvar, popwrk, popstat; */

  /* Population Print population*/
  FILE *popout;
  long *tobs, nnull;
  double *mean, *var, *skew, *kur;

  /* Stuff for the Unique Stuff (how's that for an informative comment! */
  /* A big Matrix which remembers all of our past evaluations. It's
     maximum memory is set in genoud.h */
  extern long Gnvars[MAXINSTANCES];
  double **Memory;
  long MemorySize=0, UniqueCount=0, OldUniqueCount=0;

  /* fine two unique parents count */
  long SameCount, UniquePairs;

  ExternStructure=Structure;

  Results=Structure->oResults;
  Gradients=Structure->oGradients;

  /* Structure Done */
  SetRunTimeParameters(Structure, 1,
		       &pop_size, &nvars, &MaxGenerations, &WaitGenerations,
		       &MinMax, &GradientCheck, &BoundaryEnforcement, &UseBFGS, &SolutionTolerance,
		       &InstanceNumber, &P, &P0, &P1, &P2, &P3, &P4, &P5, &P6, &P7, &P8, 
		       &PrintLevel, &HardGenerationLimit);

  /*Space allocation for all the vectors and matrices involved*/
  long lexical_end = (Structure->Lexical-1)+nvars+2;
  /* population[][0] = fitness value (first)
     population[][1:nvars] = parameter values
     population[][nvars+1] = flag for fitting
     population[][(nvars+2):((Structure->Lexical-1)+nvars+2)] = other fitness for Lexical fitting
  */
  population    = JaMatrixAllocate(pop_size+2, lexical_end);
  new_genera    = JaMatrixAllocate(pop_size+2, lexical_end);

  /* reset population to get rid of odd things being passed to R */
  for(i=1; i<=pop_size; i++)
    {
      for(j=0; j<lexical_end; j++)
	{
	  population[i][j] = 0;
	}
    }

  VECTOR LexicalReturn;
  VECTOR oldfitvalueVEC;
  short int LexicalFitsImproving;
  if(Structure->Lexical > 1)
    {
      LexicalReturn = (double *)  malloc(Structure->Lexical*sizeof(double));  
      oldfitvalueVEC = (double *)  malloc(Structure->Lexical*sizeof(double));  
    }

  temp       = matrix(0,nvars+1,0,nvars);
  probab     = Gvector(1,pop_size);
  t_vec      = Gvector(1,nvars);
  cum_probab = Gvector(1,pop_size);
  live       = ivector(1,pop_size);

  /*for oper4 Find max(2,nvars) parents for crossover operator 4*/
  p2use = nvars > 2 ? nvars : 2;
  parents    = ivector(1,p2use);

  Gnvars[Structure->InstanceNumber]=nvars;

  if (Structure->MemoryUsage==1)
    {
      if (HardGenerationLimit==0)
	MemorySize=3*(MaxGenerations+1)*pop_size+1+pop_size;
      else
	MemorySize=(MaxGenerations+1)*pop_size+1+pop_size;
      
      Memory = JaMatrixAllocate(MemorySize, lexical_end);
    }

  grad = (double *) malloc((nvars)*sizeof(double));
  evalX = (double *) malloc((nvars)*sizeof(double));
  finalhessin = (double *) malloc(((nvars*nvars)+(nvars))*sizeof(double));
  bfgsoutX = (double *) malloc((nvars+1)*sizeof(double));

  /* populationstats variables */
  mean = (double *) malloc((nvars+1)*sizeof(double));
  var = (double *) malloc((nvars+1)*sizeof(double));
  skew = (double *) malloc((nvars+1)*sizeof(double));
  kur = (double *) malloc((nvars+1)*sizeof(double));
  tobs = (long *) malloc((nvars+1)*sizeof(long));

  /* JS: Integer Q was different, why? Q=0.2; */
  Q=0.5;
  B=6;
  STEP=10;

  if(PrintLevel>0)
    {
      switch(MinMax) {
      case 0:
	Rprintf("Minimization Problem.\n");  
	break;
      case 1:
	Rprintf("Maximization Problem.\n");  
	break;
      }
    }

  /*
    if (PrintLevel>2) {
    Rprintf("Parameter B (hardcoded): %d\n", B); 
    Rprintf("Parameter Q (hardcoded): %f\n", Q);
    }
  */

  peak_cnt = 0;

  pop_size_old=0;
  if (Structure->ShareType == 1 || Structure->ShareType == 3) {

    if(PrintLevel>0)
      Rprintf( "Using old population file to initialize new population\n");

    if((popout = fopen(Structure->ProjectPath, "r")) == NULL) {
      Rprintf("         Generating new population\n");
      warning("Unable to open the old project file: %s", Structure->ProjectPath);
    }
    else {
      pop_size_old=ReadPopulation(population, pop_size, nvars, popout, PrintLevel);
      fclose(popout);
      
      for (i=1; i<=pop_size; i++) {
	for (j=1; j<=nvars; j++) {
	  population[i][j] = (int) population[i][j];
	}
      }
      
      if (pop_size_old<2) {
	warning("The old population file appears to be from a different genoud specification.");
	pop_size_old=0;
      }
    }
    if (PrintLevel>1) {
      if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	warning("Unable to open the project file: %s", Structure->ProjectPath);

	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);
	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 0, nvars+1, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	free_ivector(parents, 1);

	if(Structure->Lexical > 1)
	  {
	    free(LexicalReturn);
	    free(oldfitvalueVEC);
	  }

	error("Fatal Error. See output for diagnostic information.");
      }
      fclose(popout);
    }
  } /* end of ShareType 0 */
  else {
    if (PrintLevel>1) {
      if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	warning("Unable to open the project file: %s", Structure->ProjectPath);
	
	/* free populationstats stuff */
	free(mean);
	free(var);
	free(skew);
	free(kur);
	free(tobs);
	
	free(bfgsoutX);
	free(finalhessin);
	free(evalX);
	free(grad);
	
	/* free numeric.c allocations */
	if (Structure->MemoryUsage==1)
	  JaMatrixFree(Memory, MemorySize);
	JaMatrixFree(population, pop_size+2);
	JaMatrixFree(new_genera, pop_size+2);
	
	free_matrix(temp, 0, nvars+1, 0);
	free_vector(probab, 1);
	free_vector(t_vec, 1);
	free_vector(cum_probab, 1);
	free_ivector(live, 1);
	free_ivector(parents, 1);

	if(Structure->Lexical > 1)
	  {
	    free(LexicalReturn);
	    free(oldfitvalueVEC);
	  }
	
	error("Fatal Error. See output for diagnostic information.");
      }
      fclose(popout);
    }
  }

  /* The new initial value matrix: setting a new initial value for every individual */
  if (ExternStructure->nStartingValues > 0) 
    {
      /* Adjust old starting values (from ReadPopulation) so we have enough room for our 
	 starting.values */
      pop_size_old = pop_size_old-ExternStructure->nStartingValues-1;
      if (pop_size_old < 0)
	pop_size_old = 0;
      
      // seed the starting values until we run out of population or starting values!
      j = pop_size_old;

      for(s=0; s<ExternStructure->nStartingValues; s++) {
	j++;
	for(i=1; i<=nvars; i++) {
	  population[j][i] = (int) ExternStructure->StartingValues[s][i-1];
	  population[j][nvars+1] = -1.0;
	}
      } // end of for loop
      pop_size_old = j;

      // randomly add on people if we still have population left over!
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = (int) frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of we have starting values!
  else 
    {
      for(j=pop_size_old+1; j<=pop_size; j++) {
	for(i=1; i<=nvars; i++) { 
	  population[j][i] = (int) frange_ran(domains[i][1], domains[i][3]); 
	  population[j][nvars+1] = -1.0;
	}
      }
    } // end of else


  if (Structure->MemoryUsage==1)
    {
      OldUniqueCount=UniqueCount;

      if (UniqueCount==0)
	UniqueCount = 1;

      UniqueCount = RmemoryMatrixEvaluate(Structure->fnMemoryMatrixEvaluate, Structure->rho,
					  Memory, population,
					  MinMax, pop_size, UniqueCount,
					  nvars, Structure->Lexical, lexical_end);
      
      if ( (UniqueCount+pop_size) >= MemorySize )
	{
	  Structure->MemoryUsage=0;
	  warning("Turned Off MemoryMatrix because memory usage was too great.");
	} /* end of if */
    } // end of Memory based evaluation
  else
    {
      for (i=1; i<=pop_size; i++) 
	{
	  if (population[i][nvars+1]==-1.0 || population[i][nvars+1]==11.0)
	    {
	      for(j=1; j<=nvars; j++)
		X[j] = population[i][j];
	      
	      if (Structure->Lexical < 2)
		{
		  population[i][0] = evaluate(Structure->fn, Structure->rho, X, nvars, MinMax);
		} 
	      else 
		{
		  EvaluateLexical(Structure->fn, Structure->rho, 
				  X, nvars, Structure->Lexical, MinMax, LexicalReturn);
		  
		  population[i][0] = LexicalReturn[0];
		  count = 0;
		  for(j=(nvars+2);j<lexical_end;j++)
		    {
		      count++;
		      population[i][j] = LexicalReturn[count];
		    }		      
		} // else
	    }
	} //end of i loop
    } // end of default evaluation
  
  if(Structure->MemoryUsage!=1)
    {
      /*Sort the initial individuals based on their evaluation function*/
      if (Structure->Lexical < 2)
	{
	  sort(MinMax,population,pop_size,0);
	}
      else
	{
	  /* in eval.cpp because it is like the EvaluateLexical() function */
	  RlexicalSort(Structure->fnLexicalSort, Structure->rho,
		       population,
		       MinMax, pop_size, nvars, lexical_end, 1);
	}
    }

  peak_cnt = count_gener;

  /* since we are on generation 0 */
  oldfitvalue=population[1][0];
  if(Structure->Lexical  > 1)
    {
      oldfitvalueVEC[0]=population[1][0];
      k = 1;
      for (i=(nvars+2);i<lexical_end;i++)  {
	oldfitvalueVEC[k]=population[1][i];
	k++;  
      } /* for (i=(nvars+2);i<lexical_end;i++) */
    }

  /*
  if(PrintLevel>0)
    {
      Rprintf("\nThe best initial individual is:\n");
      print_vector(population[1],1,nvars,output);

      if (Structure->Lexical > 1)
	{
	  Rprintf("\nbest (lexical) fitness:\n");
	  Rprintf("%e  ", population[1][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[1][j]);
	    }		      
	  Rprintf("\n");
	} else {
	Rprintf("\nbest fitness: %e\n", population[1][0]);
      }
      Rprintf("\n");

      if (Structure->Lexical > 1)
	{      
	  Rprintf("The worst (lexical) fitness is:\n");
	  Rprintf("%e  ", population[pop_size][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[pop_size][j]);
	    }		   
	  Rprintf("\n");   	  
	} else {
	Rprintf("The worst fit is: %e\n", 
		population[pop_size][0]);
      }
      Rprintf("\n");
    }
  */

  if(PrintLevel==1)
    {
      if (Structure->Lexical > 1)
	{
	Rprintf("\n\nGeneration#\t    Solution Values (lexical)\n");	  
	Rprintf("\n%7d \t%e  ", 0, population[1][0]);
	for(j=(nvars+2);j<lexical_end;j++)
	  {
	    Rprintf("%e  ", population[1][j]);
	  }		      
	Rprintf("\n");	
	} else {
	Rprintf("\n\nGeneration#\t    Solution Value\n");
	Rprintf("\n%7d \t%e\n", 0, population[1][0]);
      }
    }

  /* compute and print mean and variance of population */
  if (PrintLevel>1) {
    Rprintf("GENERATION: 0 (initializing the population)\n");
    populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);
    
    if(Structure->Lexical > 1)
      {
	Rprintf( "Lexical Fit..... %e  ", population[1][0]);
	for(j=(nvars+2);j<lexical_end;j++)
	  {
	    Rprintf("%e  ", population[1][j]);
	  }		      
	Rprintf("\n");	    
      }
    else
      {
	Rprintf( "Fitness value... %e\n", population[1][0]);
	Rprintf( "mean............ %e\n", mean[0]);
	Rprintf( "variance........ %e\n", var[0]);
	/*
	  Rprintf( "skewness........ %e\n", skew[i]);
	  Rprintf( "kurtosis........ %e\n", kur[i]);
	*/
      }
    nnull = pop_size-tobs[0];
    if(nnull > 0)
      Rprintf( "#null........... %d\n", nnull);    
    if(Structure->MemoryUsage==1)
      Rprintf( "#unique......... %d, #Total UniqueCount: %d\n", 
	      UniqueCount-OldUniqueCount, UniqueCount);
    /* Rprintf( "tobs............ %d\n", tobs[i]); */
    
    for (i=1; i<=nvars; i++) {
      Rprintf( "var %d:\n", i);
      Rprintf( "best............ %e\n", population[1][i]);
      Rprintf( "mean............ %e\n", mean[i]);
      Rprintf( "variance........ %e\n", var[i]);
      /*
	Rprintf( "skewness........ %e\n", skew[i]);
	Rprintf( "kurtosis........ %e\n", kur[i]);
      */
      nnull = pop_size-tobs[i];
      if(nnull > 0)
	Rprintf( "#null........... %d\n", nnull);
      /* Rprintf( "tobs............ %d\n", tobs[i]); */
    }
  } /* end of printlevel if */
  
  /* if(PrintLevel>0)
     fflush(output); */

  /* Print the population file */
  if ( PrintLevel == 1 ) {
    if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
      Rprintf("Unable to open the project file: %s", 
	      Structure->ProjectPath);
      
      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);
      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 0, nvars+1, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      free_ivector(parents, 1);

      if(Structure->Lexical > 1)
	{
	  free(LexicalReturn);
	  free(oldfitvalueVEC);
	}
      
      error("Fatal Error. See output for diagnostic information.");
    }
    print_population((int) pop_size, (int) nvars, 0, (int) Structure->Lexical, population, popout);
    fclose(popout);
  } /* end of PrintLevel if */
  if ( PrintLevel>1) {
    if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
      Rprintf("Unable to open the project file: %s", 
	      Structure->ProjectPath);
      
      /* free populationstats stuff */
      free(mean);
      free(var);
      free(skew);
      free(kur);
      free(tobs);
      
      free(bfgsoutX);
      free(finalhessin);
      free(evalX);
      free(grad);
      
      /* free numeric.c allocations */
      if (Structure->MemoryUsage==1)
	JaMatrixFree(Memory, MemorySize);
      JaMatrixFree(population, pop_size+2);
      JaMatrixFree(new_genera, pop_size+2);
      
      free_matrix(temp, 0, nvars+1, 0);
      free_vector(probab, 1);
      free_vector(t_vec, 1);
      free_vector(cum_probab, 1);
      free_ivector(live, 1);
      free_ivector(parents, 1);

      if(Structure->Lexical > 1)
	{
	  free(LexicalReturn);
	  free(oldfitvalueVEC);
	}
      
      error("Fatal Error. See output for diagnostic information.");
    }
    print_population((int) pop_size, (int) nvars, 0, (int) Structure->Lexical, population, popout);
    fflush(popout);
    fclose(popout);
  }

  /* Interrupt setup.  Let's print a nice message to recover the best
     solution so far if at least generation 0 has been run */
  if (PrintLevel > 0 & (strcmp(Structure->ProjectPath, "/dev/null")!=0))
    setVar(install("interrupted"), ScalarLogical(1), Structure->rho);
  
  /*Assigning probability of survival for each of the agent, with the*/
  /*probability provided by the user for the best agent*/
  assign_probab(probab,pop_size,Q); 

  /*Finding the cumulative probability of the agents*/
  find_cum_probab(cum_probab,probab,pop_size);

  /*Reproducing and evaluating for the total number of generations times*/
  do
    {

      /*Initializing the live vector*/
      for(j=1; j<=pop_size; j++)
        {
          live[j] = 0;
          for(i=0; i<lexical_end; i++)
            new_genera[j][i] = population[j][i];
	  new_genera[j][nvars+1]=0;
        }

      /*Finding the agents that will die and the agents that will reproduce*/
      find_live(cum_probab,live,pop_size,P);
      /* set die_now counter to start replacements with the worst agent.
         use of die_now is okay if the entire population (except the previous
         best) is to be replaced in each generation */
      die_now = pop_size;

      j1=j2=j3=j4=j5=j6=j7=j8=0;

      /* This was causing a difference to appear between MemoryMatrix and !MemoryMatrix runs see oper5 and oper7
      UniquePairs= UniqueCount-OldUniqueCount;
      UniquePairs= (int) (0.5*(UniquePairs*UniquePairs-UniquePairs));
      if ( MAX_OPER_UNIQUE_TRY < UniquePairs)
	UniquePairs = MAX_OPER_UNIQUE_TRY;
      */
      UniquePairs = MAX_OPER_UNIQUE_TRY;

      /* main operator loop */
      while(j1+j2+j3+j4+j4+j5+j5+j6+j7+j7 < P)
        {
          oper = irange_ran(1,7);
          switch (oper)
            {
	    case 1:
	      /* JS Description: Uniform Mutation */
	      /*Applying the first operator, uniform mutation*/
	      if (j1 < P1)
		{
		  /*Find one parent for mutation operator 1*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    error("No agents to be replaced\n");
		  }
		  
		  new_genera[die_now][nvars+1] = 1.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
		    JaIntegerOper1(t_vec,domains,nvars);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j1++;
		}
	      break;
	    case 2:
	      /* JS Description: Boundary Mutation */
	      /*Applying the second operator, boundary mutation*/
	      if (j2 < P2)
		{
		  /*Find one parent for mutation operator 2*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    error( "No agents to be replaced\n");
		  }
		  
		  new_genera[die_now][nvars+1] = 2.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  JaIntegerOper2(t_vec,domains,nvars);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j2++;
		}
	      break;
	    case 3:
	      /* JS Description: Non-uniform Mutation */
	      /*Applying the third operator, non-uniform mutation*/
	      if (j3 < P3)
		{
		  /*Find one parent for mutation operator 3*/
		  first_live  = find_parent(live,pop_size);
		  live[first_live]--;
		  /* check that agent to replace is in range */
		  if (die_now < 2) {
		    error( "No agents to be replaced\n");
		  }
		  
		  new_genera[die_now][nvars+1] = 3.0;
		  for(i=1; i<=nvars; i++)
		    t_vec[i] = population[first_live][i];
		  for (ocnt = irange_ran(1,nvars); ocnt>0; ocnt--)
		    JaIntegerOper3(t_vec,domains,nvars,MaxGenerations,count_gener,B);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i] = t_vec[i];
		  die_now--;
		  j3++;
		}
	      break;
	      
	    case 4:
	      /*Applying the fourth operator, GENOUD Polytope Crossover */
	      if (j4 < (int) P4)
		{
		  /*Find max(2,nvars) parents for crossover operator 4*/
		  for (i=1; i<p2use; i++) {
		    parents[i] = find_parent(live,pop_size);
		    live[parents[i]]++;  /* no decr. first p2use-1 parents */
		  }
		  parents[p2use] = find_parent(live,pop_size);
		  /* check that agents to replace are in range */
		  if (die_now < 2) {
		    error("No agents to be replaced\n");
		  }
		  new_genera[die_now][nvars+1]  = 4.0;
		  for(j=1; j<=p2use; j++)
		    for(i=1; i<=nvars; i++)
		      temp[j][i] = population[parents[j]][i];
		  JaIntegeroper4(temp,p2use,nvars,domains);
		  for(i=1; i<=nvars; i++)
		    new_genera[die_now][i]  = temp[1][i];
		  die_now--;
		  j4++;
		}		
	      break;
	    case 5:
	      /* JS Description: Simple Crossover
		 Applying the fifth operator, simple arithmetical crossover*/
	      if (j5 < (int) P5/2)
		{
		  /*Find two distinct parents for crossover operator 5*/
		  same = TRUE;
		  SameCount=0;
		  while (same==TRUE) {
		    SameCount++;
		    
		    first_live  = find_parent(live,pop_size);
		    second_live = find_parent(live,pop_size);

		    if (SameCount >= (UniquePairs) ) 
		      break;

		    for(i=1; i<=nvars; i++)
		      {
			if ((int) population[first_live][i] != (int) population[second_live][i])
			  {
			    same = FALSE;
			    break;
			  }
		      }
		  } /* end of while same==TRUE loop */
		  /* check that agents to replace are in range */
		  if (die_now < 3) {
		    error("Not enough agents to be replaced\n");
		  }
		  live[first_live]--;
		  live[second_live]--;
		  first_die   = die_now-- ;
		  second_die  = die_now-- ;
		  new_genera[first_die][nvars+1]  = 5.0;
		  new_genera[second_die][nvars+1] = 5.0;
		  if (!same)
		    {
		      for(i=1; i<=nvars; i++)
			{
			  temp[1][i] = population[first_live][i];
			  temp[2][i] = population[second_live][i];
			}
		      JaIntegerOper5(temp[1],temp[2],STEP,domains,nvars);
		      for(i=1; i<=nvars; i++)
			{
			  new_genera[first_die][i]  = temp[1][i];
			  new_genera[second_die][i] = temp[2][i];
			}
		    }
		  else {
		    /* copy agent chosen twice into two new indivs */
		    for(i=1; i<=nvars; i++) {
		      new_genera[first_die][i]  = 
			      population[first_live][i];
		      new_genera[second_die][i] = 
			population[second_live][i];
		    }
		  }
		  j5++;
		}
	      break;
	    case 6:
	      /* JS Description: Whole Non-uniform Mutation */
	      /*Applying the sixth operator, whole non-uniform mutation*/
                    if (j6 < P6)
                      {
                        /*Find one parent for mutation operator 6*/
                        first_live  = find_parent(live,pop_size);
			live[first_live]--;
			/* check that agent to replace is in range */
			if (die_now < 2) {
			  error( "No agents to be replaced\n");
			}
			
                        new_genera[die_now][nvars+1] = 6.0;
                        for(i=1; i<=nvars; i++)
                          t_vec[i] = population[first_live][i];
                        JaIntegerOper6(t_vec,domains,nvars,MaxGenerations,count_gener,B);
                        for(i=1; i<=nvars; i++)
                          new_genera[die_now][i] = t_vec[i];
			die_now--;
                        j6++;
                      }
                    break;
	    case 7:
	      /* JS Description: Heuristic Crossover */
	      /*Applying the seventh operator*/
	      if (j7 < (int) P7/2)
		{
		  /*Find two distinct parents for operator 7*/
		  same = TRUE;
		  SameCount=0;
		  while (same==TRUE) {
		    SameCount++;
		    first_live  = find_parent(live,pop_size);
		    second_live = find_parent(live,pop_size);

		    if (SameCount >= (UniquePairs) ) 
		      break;

		    for(i=1; i<=nvars; i++)
		      {
			if ((int) population[first_live][i] != (int) population[second_live][i])
			  {
			    same = FALSE;
			    break;
			  }
		      }
		  } /* end of while same==TRUE loop */
		  /* check that agents to replace are in range */
		  if (die_now < 3) {
		    error("Not enough agents to be replaced\n");
		  }
		  live[first_live]--;
		  live[second_live]--;
		  first_die   = die_now-- ;
		  second_die  = die_now-- ;
		  new_genera[first_die][nvars+1]  = 7.0;
		  new_genera[second_die][nvars+1] = 7.0;
		  if (!same) {
		    if (first_live < second_live)
		      /* first agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[first_live][i];
			temp[1][i] = population[second_live][i];
			    }
		    else
		      /* second agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[second_live][i];
			temp[1][i] = population[first_live][i];
		      }
		    JaIntegerOper7(temp[1],temp[2],domains,nvars);
		    for(i=1; i<=nvars; i++)
		      new_genera[first_die][i]  = temp[1][i];
		    if (first_live < second_live)
		      /* first agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[first_live][i];
			temp[1][i] = population[second_live][i];
			    }
		    else
		      /* second agent is better agent */
		      for(i=1; i<=nvars; i++) {
			temp[2][i] = population[second_live][i];
			temp[1][i] = population[first_live][i];
		      }
		    JaIntegerOper7(temp[1],temp[2],domains,nvars);
		    for(i=1; i<=nvars; i++)
		      new_genera[second_die][i]  = temp[1][i];
		  }
		  else {
		    /* copy agent chosen twice into two new indivs */
		    for(i=1; i<=nvars; i++) {
			    new_genera[first_die][i]  = 
			      population[first_live][i];
			    new_genera[second_die][i] = 
			      population[second_live][i];
		    }
		  }
		  j7++;
		}
	    }
	}
      
      /*Replace the population with the new generation */
      Jnew = new_genera;
      new_genera = population;
      population = Jnew;

      if (Structure->MemoryUsage==1)
	{
	  OldUniqueCount=UniqueCount;

	  UniqueCount = RmemoryMatrixEvaluate(Structure->fnMemoryMatrixEvaluate, Structure->rho,
					      Memory, population,
					      MinMax, pop_size, UniqueCount,
					      nvars, Structure->Lexical, lexical_end);	  	  

	  if ( (UniqueCount+pop_size) >= MemorySize )
	    {
	      Structure->MemoryUsage=0;
	      warning("Turned Off MemoryMatrix because memory usage was too great.");
	    } /* end of if */
	} // end of MemoryUsage==1
        else
	  {
	    for (i=1; i<=pop_size; i++) 
	      {
		if (population[i][nvars+1]!=0)
		  {
		    for(j=1; j<=nvars; j++)
		      X[j] = population[i][j];
		    
		    if (Structure->Lexical < 2)
		      {
			population[i][0] = evaluate(Structure->fn, Structure->rho, X, nvars, MinMax);
		      } 
		    else 
		      {
			EvaluateLexical(Structure->fn, Structure->rho, 
					X, nvars, Structure->Lexical, MinMax, LexicalReturn);
			population[i][0] = LexicalReturn[0];
			count = 0;
			for(j=(nvars+2);j<lexical_end;j++)
			  {
			    count++;
			    population[i][j] = LexicalReturn[count];
			  }		      			  
		      }
		  }
	      } //end of i loop
	} //end of default evaluation scheme

      if(Structure->MemoryUsage!=1)
	{
	  /*Sort the new population based on their evaluation function*/
	  if (Structure->Lexical < 2)
	    {
	      sort(MinMax,population,pop_size,0);
	    }
	  else
	    {
	      /* in eval.cpp because it is like the EvaluateLexical() function */
	      RlexicalSort(Structure->fnLexicalSort, Structure->rho,
			   population,
			   MinMax, pop_size, nvars, lexical_end, 1);
	    }
	}

      /* check to see if fit is improving */
      if(Structure->Lexical < 2)
	{
	  switch(MinMax)
	    {
	    case 0:
	      if ( (oldfitvalue - SolutionTolerance) > population[1][0]) {
		nochange_gen=0;
		oldfitvalue=population[1][0];
		peak_cnt = count_gener;
	      }
	      else nochange_gen++;
	      break;
	    case 1:
	      if ( (oldfitvalue + SolutionTolerance) < population[1][0]) {
		nochange_gen=0;
		oldfitvalue=population[1][0];
		peak_cnt = count_gener;
	      }
	      else nochange_gen++;	      
	      break;
	    }
	} /*       if(Structure->Lexical < 2) */
      else
	{
	  switch(MinMax)
	    {
	    case 0:
	      LexicalFitsImproving = 0;
	      if ( (oldfitvalue - SolutionTolerance) > population[1][0]) 
		{
		  LexicalFitsImproving = 1;
		} 
	      else
		{
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    if ( (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) {
		      LexicalFitsImproving = 1;
		      break;
		    } /* (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) */
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		} /* else if ( (oldfitvalue - SolutionTolerance) > population[1][0])  */
	      if (LexicalFitsImproving)
		{
		  nochange_gen = 0;
		  peak_cnt = count_gener;
		  oldfitvalue=population[1][0];
		  oldfitvalueVEC[0]=population[1][0];
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    oldfitvalueVEC[k]=population[1][i];
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		}
	      else
		nochange_gen++;
	      break;
	    case 1:
	      LexicalFitsImproving = 0;
	      if ( (oldfitvalue + SolutionTolerance) < population[1][0]) 
		{
		  LexicalFitsImproving = 1;
		} 
	      else
		{
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    if ( (oldfitvalueVEC[k] + SolutionTolerance) < population[1][i] ) {
		      LexicalFitsImproving = 1;
		      break;
		    } /* (oldfitvalueVEC[k] - SolutionTolerance) > population[1][i] ) */
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		} /* else if ( (oldfitvalue - SolutionTolerance) > population[1][0])  */
	      if (LexicalFitsImproving)
		{
		  nochange_gen = 0;
		  peak_cnt = count_gener;
		  oldfitvalue=population[1][0];
		  oldfitvalueVEC[0]=population[1][0];
		  k=1;
		  for (i=(nvars+2);i<lexical_end;i++)  {
		    oldfitvalueVEC[k]=population[1][i];
		    k++;  
		  } /* for (i=(nvars+2);i<lexical_end;i++) */
		}
	      else
		nochange_gen++;
	      break;	      
	    } /* switch(MinMax) */
	} /* else (Structure->Lexical > 2) */


      if(PrintLevel==1)
	{
	  if( nochange_gen==0 )
	    {
	      if(Structure->Lexical > 1)
		{
		  Rprintf("\n%7d \t%e  ", count_gener, population[1][0]);
		  for(j=(nvars+2);j<lexical_end;j++)
		    {
		      Rprintf("%e  ", population[1][j]);
		    }		      
		  Rprintf("\n");	
		}
	      else
		{
		  Rprintf("%7d \t%e\n",
			  count_gener,population[1][0]); 
		  /* fflush(output); */
		}
	    } 
	}

      /* compute and print mean and variance of population */
      if (PrintLevel>1) {
	Rprintf("\nGENERATION: %d\n", count_gener);
	populationstats(population, pop_size, nvars, mean, var, skew, kur, tobs);

	if(Structure->Lexical > 1)
	  {
	    Rprintf( "Lexical Fit..... %e  ", population[1][0]);
	    for(j=(nvars+2);j<lexical_end;j++)
	      {
		Rprintf("%e  ", population[1][j]);
	      }		      
	    Rprintf("\n");	    
	  }
	else
	  {
	    Rprintf( "Fitness value... %e\n", population[1][0]);
	    Rprintf( "mean............ %e\n", mean[0]);
	    Rprintf( "variance........ %e\n", var[0]);
	    /*
	      Rprintf( "skewness........ %e\n", skew[i]);
	      Rprintf( "kurtosis........ %e\n", kur[i]);
	    */
	  }

	nnull = pop_size-tobs[0];
	if(nnull > 0)
	  Rprintf( "#null........... %d\n", nnull);
	if(Structure->MemoryUsage==1)
	  Rprintf( "#unique......... %d, #Total UniqueCount: %d\n", 
		  UniqueCount-OldUniqueCount, UniqueCount);
	/* Rprintf( "tobs............ %d\n", tobs[i]); */

	for (i=1; i<=nvars; i++) {
	  Rprintf( "var %d:\n", i);
	  Rprintf( "best............ %e\n", population[1][i]);
	  Rprintf( "mean............ %e\n", mean[i]);
	  Rprintf( "variance........ %e\n", var[i]);
	  /*
	    Rprintf( "skewness........ %e\n", skew[i]);
	    Rprintf( "kurtosis........ %e\n", kur[i]);
	  */
	  nnull = pop_size-tobs[i];
	  if(nnull > 0)
	    Rprintf( "#null........... %d\n", nnull);
	  /* Rprintf( "tobs............ %d\n", tobs[i]); */
	}
      } /* end of printlevel if */
      
      /* if (PrintLevel>0)
	 fflush(output); */

      /* Print the population file */
      if ( PrintLevel == 1 ) {
	if((popout = fopen(Structure->ProjectPath, "w")) == NULL) {
	  Rprintf("Unable to open the project file: %s", 
		  Structure->ProjectPath);

	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);
	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 0, nvars+1, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	  free_ivector(parents, 1);

	  if(Structure->Lexical > 1)
	    {
	      free(LexicalReturn);
	      free(oldfitvalueVEC);
	    }

	  error("Fatal Error. See output for diagnostic information.");
	}
	print_population((int) pop_size, (int) nvars, (int) count_gener, (int) Structure->Lexical, population, popout);
	fclose(popout);
      } /* end of PrintLevel if */
      if ( PrintLevel>1) {
	if((popout = fopen(Structure->ProjectPath, "a")) == NULL) {
	  Rprintf("Unable to open the project file: %s", 
		  Structure->ProjectPath);
	  
	  /* free populationstats stuff */
	  free(mean);
	  free(var);
	  free(skew);
	  free(kur);
	  free(tobs);
	  
	  free(bfgsoutX);
	  free(finalhessin);
	  free(evalX);
	  free(grad);
	  
	  /* free numeric.c allocations */
	  if (Structure->MemoryUsage==1)
	    JaMatrixFree(Memory, MemorySize);
	  JaMatrixFree(population, pop_size+2);
	  JaMatrixFree(new_genera, pop_size+2);
	  
	  free_matrix(temp, 0, nvars+1, 0);
	  free_vector(probab, 1);
	  free_vector(t_vec, 1);
	  free_vector(cum_probab, 1);
	  free_ivector(live, 1);
	  free_ivector(parents, 1);

	  if(Structure->Lexical > 1)
	    {
	      free(LexicalReturn);
	      free(oldfitvalueVEC);
	    }
	  
	  error("Fatal Error. See output for diagnostic information.");
	}
	print_population((int) pop_size, (int) nvars, (int) count_gener, (int) Structure->Lexical, population, popout);
	fflush(popout);
	fclose(popout);
      }
      
      if (nochange_gen > (WaitGenerations)) {
	/* increase the number of WaitGenerations if the gradients are NOT zero! */	  
	if (GradientCheck==0) {
	  if(PrintLevel>0)	  
	    {
	      Rprintf("\n'wait.generations' limit reached.\n");
	      Rprintf("No significant improvement in %d generations.\n", nochange_gen-1);
	      /* fflush(output); */
	    }
	  MaxGenerations = 0;
	  nochange_gen=0;
	}
	else  
	  {
	    for (i=1; i<=nvars; i++)
	      {
		bfgsoutX[i-1]=population[1][i];
	      }
	    if(Structure->UserGradient==0)
	      {	    
		gradient(Structure->fn, Structure->rho,
			 bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, domains);
	      } 
	    else 
	      {
		userGradientfn(Structure->fnGR, Structure->rho, bfgsoutX, grad, nvars);
	      }
	    GradientTrigger = 0;
	    for (i=0; i<nvars; i++) {
	      if (fabs(grad[i]) > SolutionTolerance) {
		GradientTrigger = 1;
		break;
	      }
	    } /* end for loop */
	    if (GradientTrigger==1) {
	      IncreaseGenerations = WaitGenerations;
	      WaitGenerations += IncreaseGenerations;
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nDoubling 'wait.generations' limit to %d (from %d) ", 
			  WaitGenerations, IncreaseGenerations);
		  Rprintf("because at least one gradient is too large.\n");
		  Rprintf("G[%d]: %e\t Solution Tolerance: %e\n\n", 
			  i+1, grad[i], SolutionTolerance);
		}
	    }
	    else {
	      if(PrintLevel>0)	  
		{
		  Rprintf("\n'wait.generations' limit reached.\n");
		  Rprintf("No significant improvement in %d generations.\n", nochange_gen-1);
		  /* fflush(output); */
		}
	      MaxGenerations = 0;
	      nochange_gen=0;
	    }
	  }/* end else loop */
      } /* end of if (nochange_gen > (WaitGenerations)) { */
      
      if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 
	{
	  if (HardGenerationLimit==0)
	    {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += IncreaseGenerations;
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because at least one gradient is too large.\n\n");
		}
	    } // if (Structure->HardGenerationLimit==0)
	  else
	    {
	      HardMaximumNumber = 1;
	      warning("Stopped because hard maximum generation limit was hit.\nAt least one gradient is too large.");
	      if(PrintLevel>0)	  
		{
		  Rprintf("\nNOTE: HARD MAXIMUM GENERATION LIMIT HIT\n");
		  Rprintf("        At least one gradient is too large\n");
		}
	    } // else
	} // if ( (count_gener == MaxGenerations) && (GradientTrigger==1) ) 


      /* increase the number of generations if fitness has been improving */
      if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) ) {
	if (HardGenerationLimit==0)
	  {
	    if (WaitGenerations > MaxGenerations) {
	      IncreaseGenerations = WaitGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because the fitness is still impoving.\n\n");
		}
	    }
	    else {
	      IncreaseGenerations = MaxGenerations;
	      MaxGenerations += (int) (IncreaseGenerations);
	      if(PrintLevel>0)	  
		{
		  Rprintf(
			  "\nIncreasing 'max.generations' limit by %d generations to %d ", 
			  IncreaseGenerations, MaxGenerations);
		  Rprintf("because the fitness is still improving.\n\n");
		}
	    }
	  } // if (Structure->HardGenerationLimit==0)
	else
	  {
	    if (HardMaximumNumber==0)
	      {
		warning("Stopped because hard maximum generation limit was hit.");
		if(PrintLevel>0)	  
		  {
		    Rprintf("\nNOTE: HARD MAXIMUM GENERATION LIMIT HIT\n");
		  }
	      } /* end of if HardMax */		
	  }
      } // if ( (count_gener == MaxGenerations) &&  (nochange_gen < WaitGenerations) )
      
      /* if(PrintLevel>0)	  
	 fflush(output); */
      
    } /* end of do loop */
  /*Increment iteration count and test whether all generations are done*/
  while (++count_gener <= MaxGenerations);
  
  if(PrintLevel>0)	    
    {
      if(Structure->Lexical > 1)
	{
	  Rprintf("\nSolution Lexical Fitness Value:\n");
	  Rprintf("%e  ", population[1][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      Rprintf("%e  ", population[1][j]);
	    }		      
	  Rprintf("\n");
	}
      else
	{
	  Rprintf("\nSolution Fitness Value: %e\n", population[1][0]);
	}
      if (GradientCheck==0 && UseBFGS==0)
	Rprintf("\nParameters at the Solution:\n\n");
      else
	Rprintf("\nParameters at the Solution (parameter, gradient):\n\n");
    }

  /* output data structure */
  Structure->oPeakGeneration=peak_cnt;
  Structure->oGenerations=count_gener-1;

  /* obtain gradients */
  /* print best solution */
  if (GradientCheck==0 && UseBFGS==0)
    {
      for(j = 1; j <= nvars; j++) {
	i = j-1;
	if(PrintLevel>0)	  
	  Rprintf(" X[%2d] :\t%e\n",j,population[1][j]);
	grad[i] = -1.0;
	Results[i] = population[1][j];
	Gradients[i] = grad[i];
      }
    } /* end of if (GradientCheck==0 && UseBFGS==0) */
  else 
    {
      for (i=1; i<=nvars; i++)
	{
	  bfgsoutX[i-1]=population[1][i];
	}
      if(Structure->UserGradient==0)
	{
	  gradient(Structure->fn, Structure->rho,
		   bfgsoutX, grad, nvars, MinMax, BoundaryEnforcement, domains);
	} 
      else 
	{
	  userGradientfn(Structure->fnGR, Structure->rho, bfgsoutX, grad, nvars);
	}

      for(j = 1; j <= nvars; j++) {
	i = j-1;
	if(PrintLevel>0)	  
	  Rprintf(" X[%2d] :\t%e\tG[%2d] :\t%e\n",j,population[1][j],j,grad[i]);
	Results[i] = population[1][j];
	Gradients[i] = grad[i];
      }
    } /* end of  else (GradientCheck==0 && UseBFGS==0) */

  Structure->oFitValues[0]=population[1][0];
  if (Structure->Lexical > 1)
    {
      k = 1;
      for (i=(nvars+2);i<lexical_end;i++)  {
	Structure->oFitValues[k]=population[1][i];
	k++;  	  
      }
    } 

  /* free memory */
  /* free populationstats stuff */
  free(mean);
  free(var);
  free(skew);
  free(kur);
  free(tobs);
  
  free(bfgsoutX);
  free(finalhessin);
  free(evalX);
  free(grad);

  if (Structure->MemoryUsage==1)
    JaMatrixFree(Memory, MemorySize);

  JaMatrixFree(population, pop_size+2);
  JaMatrixFree(new_genera, pop_size+2);

  free_matrix(temp, 0, nvars+1, 0);
  free_vector(probab, 1);
  free_vector(t_vec, 1);
  free_vector(cum_probab, 1);
  free_ivector(live, 1);
  free_ivector(parents, 1);

  if(Structure->Lexical > 1)
    {
      free(LexicalReturn);
      free(oldfitvalueVEC);
    }
} /* end JaIntegerOptimization */



/********************************************************************************/
/*  JaIntegerSort():                                                            */
/*                                                                              */
/*  This function sorts a double** on an integer basis.                         */
/*  The function also assumes that the double** is indexed from 1 in its rows   */
/*  and from zero in its columns.                                               */
/*                                                                              */
/********************************************************************************/

void JaIntegerSort(double **InMatrix, long n, long k)
{
  /* extern int JaIntegerCMP(); */
  long i, j;
  double **Tmp;
  extern long Gnvars[MAXINSTANCES];
  long nvars;

  Tmp = JaMatrixAllocate(n, k);

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      Tmp[i-1][j] = InMatrix[i][j];
    }
  }

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      InMatrix[i][j] = Tmp[i-1][j];
    }
  }

  JaMatrixFree(Tmp, n);
} /* end of JaIntegerSort */


/********************************************************************************/
/*  JaDoubleSort():                                                             */
/*                                                                              */
/*  This function sorts a double** on an double  basis.                         */
/*  The function also assumes that the double** is indexed from 1 in its rows   */
/*  and from zero in its columns.                                               */
/*                                                                              */
/********************************************************************************/

void JaDoubleSort(double **InMatrix, long n, long k)
{
  /* extern int JaDoubleCMP(); */
  long i, j;
  double **Tmp;
  extern long Gnvars[MAXINSTANCES];
  long nvars;

  Tmp = JaMatrixAllocate(n, k);

  nvars=Gnvars[ExternStructure->InstanceNumber];

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      Tmp[i-1][j] = InMatrix[i][j];
    }
  }

  for (i=1; i<=n; i++) {
    for (j=0; j<k; j++) {
      InMatrix[i][j] = Tmp[i-1][j];
    }
  }

  JaMatrixFree(Tmp, n);
} /* end of JaDoubleSort */



