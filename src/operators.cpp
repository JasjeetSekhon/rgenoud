/*

  RGENOUD

  Walter R. Mebane, Jr.
  University of Michigan
  http://www-personal.umich.edu/~wmebane/
  <wmebane@umich.edu>

  Jasjeet Singh Sekhon 
  UC Berkeley
  http://sekhon.polisci.berkeley.edu
  <sekhon@berkeley.edu>

  May 6, 2013

*/

#include "genoud.h"

extern "C" 
{
  double genoud_optim(SEXP fn_optim, SEXP rho, double *X, long parameters);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper1()                                      */
/*                                 Uniform Mutation                             */
/*                                                                              */
/*           SYNOPSIS          :   void oper1(parent,domains,nvars)             */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, uniform mutation.             */
/*                                                                              */
/********************************************************************************/

void oper1(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
{
  int comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  double tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;

      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Find a random value between the lower and the upper limits, to substitute*/
      /*for the old value*/
      tmp = frange_ran(llim,ulim);

      if ( parent[comp] != tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper2()                                      */
/*                                 Boundary Mutatation                          */
/*                                 No Uniqueness checking here                  */
/*                                 Don't use this oper often!                   */
/*                                                                              */
/********************************************************************************/


void oper2(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
{
  int comp;
  double llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  double tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;
      
      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_range(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Replace either the lower limit or the upper limit at random,*/
      /*for the old value*/
      tmp = (flip() == TAIL) ? llim : ulim;

      if ( tmp != parent[comp])
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */
  
  parent[comp] = tmp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper3()                                      */
/*                                                                              */
/*           SYNOPSIS          :   void oper3(parent,fin_mat,r,c,T,t,B)         */
/*                                                                              */
/********************************************************************************/

void oper3(VECTOR parent, double **domains, int nvars, int T, int t, int B)
  /* VECTOR parent; */
  /* unsigned long T;   Total number of generations*/
  /* unsigned long t;   Current generation number*/
  /* int B; */
{
  int comp;
  double llim,ulim;

  FLAG SAME;
  double tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    comp = irange_ran(1,nvars);
    find_range(&llim,&ulim,comp,domains,nvars,parent);
    
    /*From the current value of the component to be mutated, chooose at random*/
    /*whether to mutate with a lesser value or a greater value*/
    /*Then find a value lesser or greater than the original value from the*/
    /*function get_f()*/
    tmp = (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
      parent[comp]+get_F(T,t,ulim-parent[comp],B);

    if ( parent[comp] != tmp)
      SAME=FALSE;
    else if (count >= MAX_OPER_UNIQUE_TRY)
      SAME=FALSE;
  } /* end of while */

  parent[comp] = tmp;
}



/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper4()                                      */
/*                                 Polytope Crossover                           */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void oper4(MATRIX p, int p2use, int nvars)
  /* int p The parents chosen for crossover */
  /* p2use;     number of parents (rows) in p */
  /* int nvars Length of the parameter vector (cols in p) */
{
  double *A, sum;
  int    i,k;
  
  A = (double *) malloc((p2use+1)*sizeof(double));
  
  sum=0.0;
  for (k=1; k<=p2use; k++) {
    do
      A[k] = frange_ran(0.0,1.0);
    while (A[k]==0.0);                   /* insure A[k] is above 0.0 */
    sum += A[k];
  }
  sum = 1.0/sum;
  for (k=1; k<=p2use; k++) {    /* rescale A[k] to sum to 1.0 */
    A[k] *= sum;
  }
  
  for(i=1; i<=nvars; i++) {
    sum = p[1][i] * A[1];
    for (k=2; k<=p2use; k++)
      sum += p[k][i] * A[k];
    p[1][i] = sum;
  }
  
  free(A);
} /* end of oper4 */

#ifdef NEVERDEFINED
/* This is the arithmetic crossover operator */
void oper4(VECTOR p1, VECTOR p2, int nvars)
     /* VECTOR p1,p2;  The two parents chosen for crossover*/
     /* int nvars;   Length of the vector*/
{
  double **child;
  long    i;
  double  A;

  FLAG SAME;
  long count, tcount;

  child = JaMatrixAllocate(3, nvars+1); 

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    do
      A = frange_ran(0.0,1.0);
    while (A==0);                   /* insure A is above 0 */
    
    for(i=1; i<=nvars; i++)
      {
	child[1][i] = p1[i] * A + p2[i] * (1.0-A);
	child[2][i] = p2[i] * A + p1[i] * (1.0-A);
      }

    if (count >= MAX_OPER_UNIQUE_TRY)
      {
	SAME=FALSE;
	break;
      }
    
    /* Are the two new individuals unique? */    
    tcount=0;
    for (i=1; i<=nvars; i++) {
      if ( (int) child[1][i] != (int) p1[i] )
	tcount++;
      
      if ( (int) child[2][i] != (int) p2[i] )
	tcount++;
    } /* end of i loop */
    
    if (tcount==(nvars*2)) SAME=FALSE;
  } /* end of while */

  for(i=1; i<=nvars; i++)
    {
      p1[i] = child[1][i];
      p2[i] = child[2][i];
    }
  
  JaMatrixFree(child, 3);  
} /* end of oper4() */
#endif


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper5()                                      */
/*                                 Simple Crossover                             */
/*                                                                              */
/********************************************************************************/

void oper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars)
     /* VECTOR p1,p2;   *The two parents for crossing over*/
     /* int    STEP;    *Parameter for the crossover*/
{
  MATRIX child;
  FLAG BFLAG1 = FALSE,/*Check to see if the newly created vectors satisfies the*/
       BFLAG2 = FALSE;/*set of constraints*/
  int i,n=1,cut;

  /* unique check variables */
  FLAG SAME;
  int count, tcount, ccount;

  child = matrix(1,2,1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;

      /*Get a random spot on the vector for crossover*/
      cut = irange_ran(1,nvars);
      /*Copy the parent vectors on to the child vectors*/
      for(i=1; i<=cut; i++)
	{
	  child[1][i] = p1[i];
	  child[2][i] = p2[i];
	}
      do
	{
	  /*Cross the two vectors*/
	  ccount = 0;
	  for(i=cut + 1; i<=nvars; i++)
	    {
	      child[1][i] = p1[i] * (double)n/(double)STEP + p2[i] * (1.0-(double)n/(double)STEP);
	      child[2][i] = p2[i] * (double)n/(double)STEP + p1[i] * (1.0-(double)n/(double)STEP);
	      ccount++;
	    }
	  
	  /*Check to see if they satisfy the constraints*/
	  BFLAG1 = InBounds(child[1],domains,nvars);
	  BFLAG2 = InBounds(child[2],domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then generate another*/
	  /*set of crossed over values*/
	}while((n<=STEP) && ((BFLAG1 == FALSE) || (BFLAG2 == FALSE)));
      
      /* Are the two new individuals unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}

      tcount=0;
      for (i=cut+1; i<=nvars; i++) {
	if ( child[1][i] != p1[i] )
	  tcount++;
	
	if ( child[2][i] != p2[i] )
	  tcount++;
      } /* end of i loop */

      if (tcount>=(ccount*2)) SAME=FALSE;

    } /* end of while (SAME==TRUE) */

  if (BFLAG1==TRUE && BFLAG2==TRUE)
    {
      for(i=1; i<=nvars; i++)
	{
	  p1[i] = child[1][i];
	  p2[i] = child[2][i];
	}
    }
  
  free_matrix(child,1,2,1);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper6()                                      */
/*                                 Whole Non-Uniform Mutation                   */
/*                                                                              */
/********************************************************************************/


void oper6(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	MATRIX fin_mat;
	unsigned long T;    Total number of generations
	unsigned long t;    Current generation number
	int B; */
{
  int  i;
  double llim,ulim;

  /* unique check variables */
  FLAG SAME;
  long count;
  double tmp=0;

  count=0;
  SAME=TRUE;

  while (SAME==TRUE)
    {
      for (i=1; i<=nvars; i++)
	{
	  count++;
	  find_range(&llim,&ulim,i,domains,nvars,parent);
	  
	  /*From the current value of the component to be mutated, chooose at random*/
	  /*whether to mutate with a lesser value or a greater value*/
	  /*Then find a value lesser or greater than the original value from the*/
	  /*function get_f()*/
	  tmp = (flip() == TAIL) ? parent[i]-get_F(T,t,parent[i]-llim,B) :
	    parent[i]+get_F(T,t,ulim-parent[i],B);

	  if ( parent[i] != tmp)
	    SAME=FALSE;
	  else if (count >= MAX_OPER_UNIQUE_TRY)
	    SAME=FALSE;

	  parent[i] = tmp;
	}//end of for loop
    } // end of while loop 
}//oper6


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper7()                                      */
/*                                 Heuristic Crossover                          */
/*                                                                              */
/********************************************************************************/

void oper7(VECTOR p1, VECTOR p2, double **domains, int nvars)
{
  VECTOR child;
  FLAG BFLAG = FALSE;/*Check to see if the newly created vector satisfies the*/
                      /*set of constraints*/
  int i,n=2,tries=MAX_OPER_UNIQUE_TRY;
  double A;

  /* unique check variables */
  FLAG SAME;
  long count;

  child = Gvector(1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      
      do
	{
	  A = frange_ran(0.0,1.0);
	  for(i=1; i<=nvars; i++)
	    child[i] = (  A * (p2[i] - p1[i]) + p2[i] );
	  
	  /*Check to see if it satisfies the constraints*/
	  BFLAG = InBounds(child,domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then try again */
	}
      while((n<=tries) && (BFLAG == FALSE));

      /* Is the new individual unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}
      
      for (i=1; i<=nvars; i++) {
	if ( child[i] != p1[i] ) {
	  SAME=FALSE;
	  break;
	}
      }

    } /* end of while SAME loop */

  if (BFLAG==TRUE)
    {
      for(i=1; i<=nvars; i++)
	p1[i] = child[i];
    }


  free_vector(child,1);
} /* end of oper7() */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_range()                                 */
/*                                                                              */
/*           SYNOPSIS          :   void find_range(llim,ulim,domains,nvars,     */
/*                                                 parent                       */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the upper and lower      */
/*                                  limits, within which the mutation should    */
/*                                  occur.                                      */
/*                                                                              */
/********************************************************************************/

void find_range(double *llim, double *ulim, int comp, double **domains, int nvars, VECTOR parent)
     /* double *llim,*ulim; Upper and lower limits*/
     /* int comp;           Component of the vector to be mutated*/
     /* VECTOR parent;      The vector with the values of the variables*/
{
  double A, B;

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *llim = (A*domains[comp][1]) + B* parent[comp];

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *ulim = B * parent[comp] + (A*domains[comp][3]);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_rangeInt()                                 */
/*                                                                              */
/*           SYNOPSIS          :   void find_range(llim,ulim,domains,nvars,     */
/*                                                 parent                       */
/*                                                                              */
/*           DESCRIPTION       :   This function finds the upper and lower      */
/*                                  limits, within which the mutation should    */
/*                                  occur.                                      */
/*                                                                              */
/********************************************************************************/

void find_rangeInt(int *llim, int *ulim, int comp, double **domains, int nvars, VECTOR parent)
     /* int *llim,*ulim; Upper and lower limits*/
     /* int comp;           Component of the vector to be mutated*/
     /* VECTOR parent;      The vector with the values of the variables*/
{
  double A, B;

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *llim = (int) ( (A*domains[comp][1]) + B* parent[comp]);
  if( (int) domains[comp][1] > *llim)
    *llim = (int) domains[comp][1];

  A = frange_ran(0.0,1.0);
  B = 1.0 - A;

  *ulim = (int) ( B * parent[comp] + (A*domains[comp][3]) );
  if( (int) domains[comp][3] < *ulim)
    *ulim = (int) domains[comp][3];

}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   get_F()                                      */
/*                                                                              */
/*           SYNOPSIS          :   double get_F(T,t,y,B)                         */
/*                                                                              */
/*           DESCRIPTION       :   This function returns the double value which  */
/*                                  is the evaluated value of the function,     */
/*                                  needed for the operators 3 and 6            */
/*                                                                              */
/********************************************************************************/

double get_F(int T, int t, double y, int B)
     /*
       unsigned long t,T;
       int           B;
       double         y;
     */
{
  double factor;

  factor =  (double) pow(1.0 - (double)t/(double)T,(double)B);
  factor = factor * frange_ran(0.0,1.0);
  if (factor < 0.00001)
    factor = 0.00001;
  return(y * factor);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   oper8()                                      */
/*                                 Local-Minimum Crossover: bfgs                */
/*                                                                              */
/*           SYNOPSIS          :   void oper8(parent,fin_mat,rc)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, BFGS.                         */
/*                                                                              */
/*                                                                              */
/*           CALLING FUNCITONS :   optimization()                               */
/*                                                                              */
/********************************************************************************/

void oper8(SEXP fn_optim, SEXP rho,
	   VECTOR parent, MATRIX domains, 
	   double SolutionTolerance, long nvars, 
	   short BoundaryEnforcement, 
	   short PrintLevel,
	   double mix)
{

  double *parm, *work;
  long i, j, btest;
  double bfgsfit;
  double A, B;

  parm  = (double *) malloc((nvars+1)*sizeof(double)); 
  work  = (double *) malloc((nvars+1)*sizeof(double));

  if( mix < 0)
    {
      A = frange_ran(0.0,1.0);
    } 
  else 
    {
      A = mix;
    }
  B = 1.0 - A;

  for (i=0; i<nvars; i++) {
    parm[i] = parent[i+1];
  }

  bfgsfit = genoud_optim(fn_optim, rho, parm, nvars);

  if (BoundaryEnforcement==0) {
    for(i=1; i<=nvars; i++) {
      parent[i] = A * parm[i-1] + B * parent[i];
    }
  }
  else {
    for (j=0; j<20; j++) 
      {
	btest = 0;
	for (i=1; i<=nvars; i++) 
	  {
	    work[i] = A * parm[i-1] + B * parent[i]; 

	    btest = (domains[i][1] > work[i]) || (work[i] > domains[i][3]) ;
	    /* shrink point until all parameters are in bounds */
	    if (btest) 
	      {
		if(PrintLevel > 1)
		  {
		    Rprintf("NOTE: killing out-of-bounds individual created by bfgs oper(9). fit:%10.8lf\n",bfgsfit);
		    Rprintf("NOTE: oper(9) Parameter: %d \t Value: %e\n\n", i, work[i]);
		  }
		warning("killed out-of-bounds individual created by bfgs oper(9)");
	      }
	  }
	if (btest==0) break;
	A *= 0.5 ;
	B = 1.0 - A;
      }
    if (j<20) 
      { 
	/* leave parent unchanged if not in boundary after 20 halvings */
	for (i=1; i<=nvars; i++) 
	  {
	    parent[i] = work[i];
	  }
      }
  }

  free(work);
  free(parm);

  return ;
}

/********************************************************************************/
/* Integer Operators!                                                           */
/*                                                                              */
/********************************************************************************/

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper1()                             */
/*                                 Uniform Mutation                             */
/*                                                                              */
/*           SYNOPSIS          :   void oper1(parent,fin_mat,rc)                */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, uniform mutation.             */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper1(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
     /* INDEX rc;       Row and column of the final matrix*/
{
  long comp;
  int llim,ulim;/*Lower and Upper limits of the value to be mutated*/
  FLAG SAME;
  int tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;
      comp = irange_ran(1,nvars);
        
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_rangeInt(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Find a random value between the lower and the upper limits, to substitute*/
      /*for the old value*/
      tmp =  irange_ran(llim,ulim);

      if ( (int) parent[comp] != (int) tmp)
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = (int) tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper2()                             */
/*                                 Boundary Mutatation                          */
/*                                 No Uniqueness checking here                  */
/*                                 Don't use this oper often!                   */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper2(VECTOR parent, double **domains, int nvars)
     /* VECTOR parent;  The parent vector*/
     /* MATRIX fin_mat; The final matrix*/
     /* INDEX rc;       Row and column of the final matrix*/
{
  int comp;
  int llim,ulim;/*Lower and Upper limits of the value to be mutated*/

  FLAG SAME;
  int tmp;
  long count;

  count=0;
  SAME=TRUE;
  while (SAME==TRUE) 
    {
      count++;

      comp = irange_ran(1,nvars);
      
      /*Finding the lower and upper limits between which the values are to be mutated*/
      find_rangeInt(&llim,&ulim,comp,domains,nvars,parent);
      
      /*Replace either the lower limit or the upper limit at random,*/
      /*for the old value*/
      tmp =  (int) (flip() == TAIL) ? llim : ulim;

      if ( (int) tmp != (int) parent[comp])
	SAME=FALSE;
      else if (count >= MAX_OPER_UNIQUE_TRY)
	SAME=FALSE;
    } /* end of while */

  parent[comp] = (int) tmp;
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegeroper3()                             */
/*                                 Non-Uniform Mutation                         */
/*                                                                              */
/*           SYNOPSIS          :   VECTOR oper3(parent,fin_mat,r,c,T,t,B)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns a new vector generated */
/*                                  from the parent vector, after applying      */
/*                                  the operator, non-uniform mutation.         */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper3(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	MATRIX fin_mat;
	INDEX rc; */
     /* unsigned long T;   Total number of generations*/
     /* unsigned long t;   Current generation number*/
     /* int B; */
{
  int comp;
  int llim,ulim;

  FLAG SAME;
  int tmp;
  long count;

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;
    comp = irange_ran(1,nvars);

    find_rangeInt(&llim,&ulim,comp,domains,nvars,parent);

    /*From the current value of the component to be mutated, chooose at random*/
    /*whether to mutate with a lesser value or a greater value*/
    /*Then find a value lesser or greater than the original value from the*/
    /*function get_f()*/
    tmp = 
      (int) (flip() == TAIL) ? parent[comp]-get_F(T,t,parent[comp]-llim,B) :
      parent[comp]+get_F(T,t,ulim-parent[comp],B);

    if ( (int) parent[comp] != (int) tmp)
      SAME=FALSE;
    else if (count >= MAX_OPER_UNIQUE_TRY)
      SAME=FALSE;
  } /* end of while */
  parent[comp] = (int) tmp;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegeroper4()                             */
/*                                 Polytope Crossover                           */
/*                                                                              */
/*                                                                              */
/********************************************************************************/

void JaIntegeroper4(MATRIX p, int p2use, int nvars, double **domains)
  /* int p The parents chosen for crossover */
  /* p2use;     number of parents (rows) in p */
  /* int nvars Length of the parameter vector (cols in p) */
{
  double *A, sum;
  int    i,k;

  A = (double *) malloc((p2use+1)*sizeof(double));

  sum=0.0;
  for (k=1; k<=p2use; k++) {
    do
      A[k] = frange_ran(0.0,1.0);
    while (A[k]==0.0);                   /* insure A[k] is above 0.0 */
    sum += A[k];
  }
  sum = 1.0/sum;
  for (k=1; k<=p2use; k++) {    /* rescale A[k] to sum to 1.0 */
    A[k] *= sum;
  }

  for(i=1; i<=nvars; i++) {
    sum = p[1][i] * A[1];
    for (k=2; k<=p2use; k++)
      sum += p[k][i] * A[k];

    p[1][i] = (int) sum;

    if( (int) domains[i][1] > (int) p[1][i])
      p[1][i] = (int) domains[i][1];
    if( (int) domains[i][3] < (int) p[1][i])
      p[1][i] = (int) domains[i][3];
  }

  free(A);
}


#ifdef NEVERDEFINED
void JaIntegerOper4(VECTOR p1, VECTOR p2, int nvars)
     /* VECTOR p1,p2;  The two parents chosen for crossover*/
     /* int nvars;   Length of the vector*/
{
  double **child;
  long    i;
  double  A;

  FLAG SAME;
  long count, tcount;

  child = JaMatrixAllocate(3, nvars+1); 

  count=0;

  SAME=TRUE;
  while (SAME==TRUE) {
    count++;

    do
      A = frange_ran(0.0,1.0);
    while (A==0);                   /* insure A is above 0 */
    
    for(i=1; i<=nvars; i++)
      {
	child[1][i] = ( p1[i] * A + p2[i] * (1.0-A) );
	child[2][i] = ( p2[i] * A + p1[i] * (1.0-A) );
      }

    if (count >= MAX_OPER_UNIQUE_TRY)
      {
	SAME=FALSE;
	break;
      }

    /* Are the two new individuals unique? */    
    tcount=0;
    for (i=1; i<=nvars; i++) {
      if ( (int) child[1][i] != (int) p1[i] )
	tcount++;

      if ( (int) child[2][i] != (int) p2[i] )
	tcount++;
    } /* end of i loop */
	 
    if (tcount==(nvars*2)) SAME=FALSE;
  } /* end of SAME while */

    for(i=1; i<=nvars; i++)
      {
	p1[i] = (int) child[1][i];
	p2[i] = (int) child[2][i];
      }


  JaMatrixFree(child, 3);  
} /* end of JaIntegerOper4 */
#endif


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper5()                             */
/*                                 Simple Crossover                             */
/*                                                                              */
/*           SYNOPSIS          :   void oper5(p1,p2,STEP,rc,fin_mat,X,x2)       */
/*                                                                              */
/*           DESCRIPTION       :   This function returns two new vectors        */
/*                                  generated after simple arithmetical         */
/*                                  crossover, from the two parent vectors.     */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars)
     /* VECTOR p1,p2;   *The two parents for crossing over*/
     /* INDEX rc;       *Row and column of the final matrix*/
     /* MATRIX fin_mat; *The final matrix*/
     /* int    STEP;    *Parameter for the crossover*/
{
  MATRIX child;
  FLAG BFLAG1 = FALSE,/*Check to see if the newly created vectors satisfies the*/
       BFLAG2 = FALSE;/*set of constraints*/
  int i,n=1,cut;

  /* unique check variables */
  FLAG SAME;
  int count, tcount, ccount;


  child = matrix(1,2,1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      /*Get a random spot on the vector for crossover*/
      cut = irange_ran(1,nvars);
      /*Copy the parent vectors on to the child vectors*/
      for(i=1; i<=cut; i++)
	{
	  child[1][i] = p1[i];
	  child[2][i] = p2[i];
	}
      do
	{
	  /*Cross the two vectors*/
	  ccount = 0;
	  for(i=cut + 1; i<=nvars; i++)
	    {
	      child[1][i] = p1[i] * (double)n/(double)STEP + p2[i] * (1.0-(double)n/(double)STEP);
	      child[2][i] = p2[i] * (double)n/(double)STEP + p1[i] * (1.0-(double)n/(double)STEP);
	      ccount++;
	    }
	  
	  /*Check to see if they satisfy the constraints*/
	  BFLAG1 = InBounds(child[1],domains,nvars);
	  BFLAG2 = InBounds(child[2],domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then generate another*/
	  /*set of crossed over values*/
	}while((n<=STEP) && ((BFLAG1 == FALSE) || (BFLAG2 == FALSE)));

      /* Are the two new individuals unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}

      tcount=0;
      for (i=cut+1; i<=nvars; i++) {
	if ( (int) child[1][i] != (int) p1[i] )
	  tcount++;
	
	if ( (int) child[2][i] != (int) p2[i] )
	  tcount++;
      } /* end of i loop */

      if (tcount>=(ccount*2)) SAME=FALSE;
      
    } /* end of while (SAME==TRUE); */

  if (BFLAG1==TRUE && BFLAG2==TRUE)
    {
      for(i=1; i<=nvars; i++)
	{
	  p1[i] = (int) child[1][i];
	  p2[i] = (int) child[2][i];
	}
    }

  free_matrix(child,1,2,1);
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper6()                             */
/*                                 Whole Non-Uniform Mutation                   */
/*                                                                              */
/********************************************************************************/


void JaIntegerOper6(VECTOR parent, double **domains, int nvars, int T, int t, int B)
     /* VECTOR parent;
	unsigned long T;    Total number of generations
	unsigned long t;    Current generation number
	int B; */
{
  int  i;
  int llim,ulim;

  /* unique check variables */
  FLAG SAME;
  long count;
  int tmp;

  count=0;
  SAME=TRUE;

  while (SAME==TRUE)
    {
      for (i=1; i<=nvars; i++)
	{
	  count++;
	  find_rangeInt(&llim,&ulim,i,domains,nvars,parent);
	  
	  /*From the current value of the component to be mutated, chooose at random*/
	  /*whether to mutate with a lesser value or a greater value*/
	  /*Then find a value lesser or greater than the original value from the*/
	  /*function get_f()*/
	  tmp = (int) (flip() == TAIL) ? parent[i]-get_F(T,t,parent[i]-llim,B) :
	    parent[i]+get_F(T,t,ulim-parent[i],B);

	  if ( (int) parent[i] != (int) tmp)
	    SAME=FALSE;
	  else if (count >= MAX_OPER_UNIQUE_TRY)
	    SAME=FALSE;
	  
	  parent[i] = (int) tmp;
	}//end of for loop
    } // end of while loop 
}//oper6

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   JaIntegerOper7()                             */
/*                                 Heuristic Crossover                          */
/*                                                                              */
/********************************************************************************/

void JaIntegerOper7(VECTOR p1, VECTOR p2, double **domains, int nvars)
{
  VECTOR child;
  FLAG BFLAG = FALSE;/*Check to see if the newly created vector satisfies the*/
                      /*set of constraints*/
  int i,n=2,tries=MAX_OPER_UNIQUE_TRY;
  double A;

  /* unique check variables */
  FLAG SAME;
  long count;

  child = Gvector(1,nvars);

  count=0;
  SAME=TRUE;
  while (SAME==TRUE)
    {
      count++;
      
      do
	{
	  A = frange_ran(0.0,1.0);
	  for(i=1; i<=nvars; i++)
	    child[i] = (int) (  A * (p2[i] - p1[i]) + p2[i] );
	  
	  /*Check to see if it satisfies the constraints*/
	  BFLAG = InBounds(child,domains,nvars);
	  n++;
	  /*If the constraints not satisfied, then try again */
	}
      while((n<=tries) && (BFLAG == FALSE));

      /* Is the new individual unique? */
      if (count >= MAX_OPER_UNIQUE_TRY)
	{
	  SAME=FALSE;
	  break;
	}
      
      for (i=1; i<=nvars; i++) {
	if ( (int) child[i] != (int) p1[i] ) {
	  SAME=FALSE;
	  break;
	}
      }

    } /* end of while SAME loop */

  if (BFLAG==TRUE) 
    {
      for(i=1; i<=nvars; i++)
	p1[i] = (int) child[i];
    }


  free_vector(child,1);
} /* end of JaIntegerOper7() */


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   InBounds(child, domains, nvars)              */
/*                                                                              */
/*           DESCRIPTION       :   This function returns TRUE or FALSE depending*/
/*                                  on whether the vector passed satisfies all  */
/*                                  the constraints or not.                     */
/*                                                                              */
/********************************************************************************/
FLAG InBounds(VECTOR child, double **domains, int nvars)
     /* VECTOR child;   The vector to be checked for violation of constriants*/
{
  int i;
  

  for(i=1; i<=nvars; i++)
    {
      if( (child[i] < domains[i][1]) || (child[i] > domains[i][3]) )
	return(FALSE);
    }

  return(TRUE);
}


