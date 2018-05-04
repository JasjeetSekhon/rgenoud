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

*/


#include "genoud.h"


/********************************************************************************/
/*  ReadPopulation():                                                           */
/*                                                                              */
/*  This function reads in an old population file and initializes the           */
/*  newpopulation with it.                                                      */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

long ReadPopulation(double **Data, long NewPopSize, long NewVars, FILE *fp, short PrintLevel)
{
  char ctmp[MAXPATH];
  int generation, PopSize, nvars, UsePopSize, FitVals;
  int i, j, ltmp, fint;
  double **OldData;
  short trip=0;

  /* This reads the "Generations:" name */

  while(!(feof(fp))) {

    /* pos = ftell(fp); */

    fint = fscanf(fp, "%s", ctmp);
    fint = fscanf(fp, " %d", &generation);

    if(PrintLevel>0)
      Rprintf( "Generation: %d\n", generation);
    /* This reads the "Population" name */
    fint = fscanf(fp, "%s", ctmp);
    /* This reads the "Size:" name */
    fint = fscanf(fp, "%s", ctmp);
    fint = fscanf(fp, " %d", &PopSize);

    if( ((PrintLevel>0) & (trip==0)) )
      Rprintf( "Population Size: %d\n", PopSize); 

    fint = fscanf(fp, "%s", ctmp);     /* reads "Fit" */
    fint = fscanf(fp, "%s", ctmp);     /* reads "Values:" */
    fint = fscanf(fp, "%d", &FitVals); /* reads number of fit values */

    if(FitVals > 1)
      warning("Reading an existing population file is not supported for Fit Values != 1");

    /* This reads the "Variables:" name */
    fint = fscanf(fp, "%s", ctmp);
    fint = fscanf(fp, " %d", &nvars);

    if( ((PrintLevel>0) & (trip==0)) )
      Rprintf( "Number of Variables: %d\n", nvars); 

    if (trip==0) {
      if (nvars!=NewVars) return(0);

      OldData = JaMatrixAllocate(PopSize+2, nvars+2);
      trip++;
    }
    
    /* loop over the main data part */
    for (i=1; i<=PopSize; i++) {
      fint = fscanf(fp,"%d",&ltmp);
      for (j=0; j<=nvars; j++) {
	fint = fscanf(fp,"%lf", &OldData[i][j]);
      }
    }

  }

  /* Map OldData to Data */
  if (NewPopSize < PopSize) UsePopSize = NewPopSize;
  else UsePopSize=PopSize;
  for (i=1; i<=UsePopSize; i++) {
    Data[i][nvars+1] = 0.0;
    for (j=0; j<=nvars; j++) {
      Data[i][j] = OldData[i][j];
    }
  }

  /* let's print the population file */
  if(PrintLevel>1)
    {
      Rprintf( "\nRead in Population. Used Population Size: %d\n", UsePopSize); 
      for (i=1; i<=UsePopSize; i++) {
	Rprintf( "%d \t", i); 
	for (j=0; j<=nvars; j++) {
	  Rprintf( "%e \t", Data[i][j]);
	}
	Rprintf( "\n");
      }
      Rprintf( "\n");
      /* fflush(output); */
    }

  JaMatrixFree(OldData, PopSize);
  return(PopSize);
} /* end Read Population */



/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_domains()                              */
/*                                                                              */
/*           SYNOPSIS          :   void print_domains(equal,t_equ)              */
/*                                                                              */
/*           DESCRIPTION       :   This function prints the matrix passed, on to*/
/*                                  the standard output, in the format of       */
/*                                  domains.                                    */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void print_domains(MATRIX equal, int t_equ, short DataType)
     /*
       MATRIX equal;   the domains matrix, with the upper and lower limits
       int t_equ;      *the total number of domains
     */
{
  int i,j;

  Rprintf("Domains:\n");
  //Integer
  if (DataType==1) 
  {
      for(i=1; i<=t_equ; i++)
      {
	  for(j=1; j<=3; j++)
	  {
	      if(j == 2)
		  Rprintf("  <=  X%-2d  <=   ",(int)equal[i][j]);
	      else
		  Rprintf(" %d ",(int) equal[i][j]);
	  }
	  Rprintf("\n");
      }
  } else {
      for(i=1; i<=t_equ; i++)
      {
	  for(j=1; j<=3; j++)
	  {
	      if(j == 2)
		  Rprintf("  <=  X%-2d  <=   ",(int)equal[i][j]);
	      else
		  Rprintf(" %e ",equal[i][j]);
	  }
	  Rprintf("\n");
      }
  }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   print_population()                           */
/*                                                                              */
/********************************************************************************/

void print_population(int popsize, int nvars, int generation, int lexical, double **foo, FILE *out)
{
  int i,j;

  if (lexical < 2)
    {
      fprintf(out,"Generation: %d \t Population Size: %d \t Fit Values: 1 \t Variables: %d\n\n", 
	      generation, popsize, nvars);
      for(i = 1; i <= popsize; i++)
	{
	  fprintf(out,"%d \t %e \t",i, foo[i][0]);
	  for (j = 1; j <= nvars; j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }
	  fprintf(out,"\n");
	}
      fprintf(out,"\n\n");
    }
  else 
    {
      long lexical_end = lexical-1+nvars+2;

      fprintf(out,"Generation: %d \t Population Size: %d \t Fit Values: %d \t Variables: %d\n\n", 
	      generation, popsize, lexical, nvars);
      for(i = 1; i <= popsize; i++)
	{
	  fprintf(out,"%d \t ", i);

	  /* print lexical fit values */
	  fprintf(out,"%e \t ",foo[i][0]);
	  for(j=(nvars+2);j<lexical_end;j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }		      
	  
	  /* print variables */
	  for (j = 1; j <= nvars; j++)
	    {
	      fprintf(out,"%e \t ",foo[i][j]);
	    }
	  fprintf(out,"\n");
	}
      fprintf(out,"\n\n");
    }
} /* end */


