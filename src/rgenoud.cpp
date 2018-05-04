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

  August 3, 2009
*/

#include "genoud.h"

void genoud(struct GND_IOstructure *Structure);

extern "C" 
{

  // mkanswer
  SEXP mkans(double *oFitValues, double *oResults, double *oGradients, long *oP, long oGenerations,
	     long oPeakGeneration, long oPopSize, long nvars, long lexical)
  { 
    SEXP ans;
    long length, i, indx, operators;
    
    operators=9;
    length= lexical + (nvars*2) + 3 + operators;

    PROTECT(ans=allocVector(REALSXP,length));
    REAL(ans)[0] = (double) oGenerations;
    REAL(ans)[1] = (double) oPeakGeneration;
    REAL(ans)[2] = (double) oPopSize;
    indx = 2;
    // include fit values
    for (i=0; i<lexical; i++) {
      indx++;
      REAL(ans)[indx] = oFitValues[i];
    }    
    // include results
    for (i=0; i<nvars; i++) {
      indx++;
      REAL(ans)[indx] = oResults[i];
    }
    // include gradients
    for (i=0; i<nvars; i++) {
      indx++;
      REAL(ans)[indx] = oGradients[i];
    }
    // include the actual operator count
    for (i=0; i<operators; i++) {
      indx++;
      REAL(ans)[indx] = oP[i];
    }
    UNPROTECT(1);

    return(ans);
  } // end of mkans


  double genoud_optim(SEXP fn_optim, SEXP rho, double *X, long parameters)
  {
    SEXP ans, R_fcall, x;
    double fit;
    long i;

    PROTECT(x = allocVector(REALSXP, parameters));

    for (i=0; i<parameters; i++)
      {
	REAL(x)[i] = X[i];
      }

    PROTECT(R_fcall = lang2(fn_optim, R_NilValue));
    SETCADR(R_fcall, x);

    ans = eval(R_fcall, rho);
    fit = REAL(ans)[0];

    for(i=0; i<parameters; i++)
      {
	X[i] = REAL(ans)[i+1];
      }

    UNPROTECT(2);
    return(fit);
  } // end of genoud_optim()



  SEXP rgenoud(SEXP fn, SEXP rho,
	       SEXP nvars, SEXP pop_size, SEXP max_generations, SEXP wait_generations,
	       SEXP n_starting_values, SEXP starting_values,
	       SEXP P, SEXP Domains, 
	       SEXP max, SEXP gradient_check, SEXP boundary_enforcement,
	       SEXP solution_tolerance, SEXP BFGS, SEXP data_type_int,
	       SEXP provide_seeds, SEXP unif_seed, SEXP int_seed,
	       SEXP print_level, SEXP share_type, SEXP instance_number,
	       SEXP MemoryMatrix, SEXP Debug,
	       SEXP output_path, SEXP output_type, SEXP project_path,
	       SEXP hard_generation_limit,
	       SEXP fn_optim, 
	       SEXP lexical, SEXP fnLexicalSort, SEXP fnMemoryMatrixEvaluate,
	       SEXP RuserGradient, SEXP fnGR,
               SEXP RP9mix, SEXP BFGSburnin, SEXP transform)
  {

    SEXP ret;
    long parameters, i, j;

    double *FitValues, *Results, *Gradients;

    if(!isEnvironment(rho)) 
      error ("`rho' should be an environment");

    parameters = asInteger(nvars);

    // setup GENOUD
    struct GND_IOstructure *MainStructure;
    MainStructure = (struct GND_IOstructure *) malloc(sizeof(struct GND_IOstructure));

    double **domains;
    domains = (double **) malloc(parameters*sizeof(double));
    for (i=0; i<parameters; i++) {
      domains[i] = (double *) malloc(2*sizeof(double));
    }

    for (j=0; j<2; j++) {
      for (i=0; i<parameters; i++) {
	domains[i][j] = REAL(Domains)[i + j*parameters];
      }
    }

    // starting values
    double **StartingValues;
    int nStartingValues;
    nStartingValues = asInteger(n_starting_values);
    if (nStartingValues > 0) {
      /* need to free a matrix of StaringValues below */
      StartingValues = (double **) malloc(nStartingValues*sizeof(double));
      for (i=0; i<nStartingValues; i++) {
	StartingValues[i] = (double *) malloc(parameters*sizeof(double));
        for(j=0; j<parameters; j++) 
	  StartingValues[i][j] = REAL(starting_values)[(i * parameters + j)];
      }
    }

    MainStructure->fn=fn;
    MainStructure->rho=rho;
    MainStructure->fnLexicalSort=fnLexicalSort;
    MainStructure->fnMemoryMatrixEvaluate=fnMemoryMatrixEvaluate;
    MainStructure->fnGR=fnGR;
    MainStructure->fn_optim=fn_optim;
    MainStructure->Lexical=asInteger(lexical);
    MainStructure->UserGradient=asInteger(RuserGradient);
    MainStructure->nvars=parameters;
    MainStructure->PopSize=asInteger(pop_size);
    MainStructure->MaxGenerations=asInteger(max_generations);
    MainStructure->WaitGenerations=asInteger(wait_generations);
    MainStructure->HardGenerationLimit=asInteger(hard_generation_limit);
    MainStructure->nStartingValues=nStartingValues;
    MainStructure->StartingValues=StartingValues;
    MainStructure->P[0]=REAL(P)[0];
    MainStructure->P[1]=REAL(P)[1];
    MainStructure->P[2]=REAL(P)[2];
    MainStructure->P[3]=REAL(P)[3];
    MainStructure->P[4]=REAL(P)[4];
    MainStructure->P[5]=REAL(P)[5];
    MainStructure->P[6]=REAL(P)[6];
    MainStructure->P[7]=REAL(P)[7];
    MainStructure->P[8]=REAL(P)[8];
    MainStructure->Domains=domains;
    MainStructure->MinMax=asInteger(max);
    MainStructure->GradientCheck=asInteger(gradient_check);
    MainStructure->BoundaryEnforcement=asInteger(boundary_enforcement);
    MainStructure->SolutionTolerance=asReal(solution_tolerance);
    MainStructure->UseBFGS=asInteger(BFGS);

    MainStructure->MemoryUsage=asInteger(MemoryMatrix);
    MainStructure->Debug=asInteger(Debug);

    MainStructure->InstanceNumber=asInteger(instance_number);

    MainStructure->ProvideSeeds=asInteger(provide_seeds);
    MainStructure->UnifSeed=asInteger(unif_seed);
    MainStructure->IntSeed=asInteger(int_seed);
    MainStructure->PrintLevel=asInteger(print_level);
    MainStructure->DataType=asInteger(data_type_int);

    /* 
       Share Type:
       (0) no reading of the existing project file and no looking at the public population file
       (1) reading of any existing project file, but no examining of public population file
       (2) NO reading of any existing project file but examination of public population file
       (3) BOTH reading of any existing project file AND examination of public population file
    */
    MainStructure->ShareType=asInteger(share_type);

    //Paths
    char OutputPath[1000], ProjectPath[1000];
    strcpy(OutputPath,STRING_VALUE(output_path));
    strcpy(ProjectPath,STRING_VALUE(project_path));
    MainStructure->OutputPath=OutputPath;
    MainStructure->ProjectPath=ProjectPath;
    MainStructure->OutputType=asInteger(output_type);

    /* output data structures */
    FitValues = (double *) malloc(MainStructure->Lexical*sizeof(double));  
    Results = (double *)  malloc(parameters*sizeof(double));
    Gradients = (double *)  malloc(parameters*sizeof(double));
  
    MainStructure->oFitValues=FitValues;
    MainStructure->oResults=Results;
    MainStructure->oGradients=Gradients;

    /* from setupGenoud */
    /* output data structures */
    MainStructure->oGenerations=0;
    MainStructure->oPeakGeneration=0;
    MainStructure->oPopSize=0;
    MainStructure->ThreadNumber=0;

    /* Operator Options */
    MainStructure->P9mix=asReal(RP9mix);
    MainStructure->BFGSburnin=asInteger(BFGSburnin);

    /* Transform Related Variables */
    /* whichFUN == 3 implies EvaluateTransform should be called */
    /* whichFUN == 2 implies EvaluateLexical should be called */
    /* whichFUN == 1 implies evaluate should be called */
    MainStructure->Transform=asInteger(transform);
    if(MainStructure->Transform == 1)
        MainStructure->whichFUN = 3;
    else if(MainStructure->Lexical > 1)
        MainStructure->whichFUN = 2;
    else
        MainStructure->whichFUN = 1;

    genoud(MainStructure);

    ret = mkans(MainStructure->oFitValues,
		MainStructure->oResults, MainStructure->oGradients, MainStructure->oP,
		MainStructure->oGenerations, MainStructure->oPeakGeneration,
		MainStructure->oPopSize, MainStructure->nvars, MainStructure->Lexical);

    // Free memory
    free(MainStructure);
    for (i=0; i<parameters; i++) 
      free(domains[i]);
    free(domains);
    free(Results);
    free(Gradients);
    free(FitValues);

    if (nStartingValues > 0) {
      for (i=0; i<nStartingValues; i++) 
	free(StartingValues[i]);
      free(StartingValues);
    }

    return(ret);
  } // end of rgenoud()
  
} // end of extern "C"




