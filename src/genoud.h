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

#include<R.h>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<time.h>
#include<string.h>
#include<stdarg.h>
#include<Rdefines.h>

extern "C"
{
  /* function definitions */
  double evaluate(SEXP fn, SEXP rho, double *X, long nvars, short int MinMax);
  void EvaluateLexical(SEXP fn, SEXP rho,
		       double *X, long nvars, long lexical, short int MinMax, double *ret);
  void EvaluateTransform(SEXP fn, SEXP rho,
		       double *X, long nvars, long lexical, short int MinMax, double *ret);  
} /*end of extern C */

#define M(ROW,COL,NCOLS) (((ROW)*(NCOLS))+(COL))
#define EVALUATE -645271937
#define DOUBLEMAX DOUBLE_XMAX
#define MAXPATH 1000
#define MAXTHREADS 20
#define MAXINSTANCES 20
#define ERROR_CODE -99999
#define MAX_OPER_UNIQUE_TRY 1000

#define TRUE 1
#define FALSE 0

#define DOS_SYS  FALSE           /* set to true for dos, false for unix */
#define UNIX_SYS TRUE            /* set to true for unix, false for dos */

#define flip()  ((int) ((newrand()*(long)2)/(long) 65535))
#define MIN -32768
#define MAX 32768
#define HEAD 1
#define TAIL 0
#define TRIES 1000
#define MULT 25173
#define INCR 13849
#define MOD ((long int) 65536)
#define SHUFFLE 256   /* size of random number shuffle array */
#define EST_OFFSET 0  /* offset in hours from Eastern Standard Time zone)  */

#define NOTDONE_ADD 0.25
/* number of generations with no changes to be treated as convergence at
   generation limit */
#define NOTDONE_LIM 50

typedef double **MATRIX;
typedef double *VECTOR;
typedef int **IMATRIX;
typedef int *IVECTOR;
typedef int FLAG;
typedef int TOSS;
typedef struct {int r; int c;}INDEX;



struct GND_IOstructure
{
  /* --- Basic Input Parameters ---- */
  SEXP          fn;
  SEXP          rho;
  SEXP          fnLexicalSort;
  SEXP          fnMemoryMatrixEvaluate;
  SEXP          fnGR;
  SEXP          fn_optim;
  long		nvars;
  long		PopSize;
  long		MaxGenerations;
  long		WaitGenerations; 
  double	P[9];		  /* Operators */
  double	**Domains;
  short		MinMax;
  short		GradientCheck;
  short		BoundaryEnforcement;  /* 0=anything goes, 1: regular; 2: no trespassing! */
  double	SolutionTolerance;
  long		ThreadNumber;	/* indexed from zero */
  long          InstanceNumber; /* indexed from zero, the number of parallel runs */
  short		UseBFGS;        /* Use BFGS on the Best Individual 1= True, 0=False */
  short         DataType;       /* 1== integer, everything else equals float */
  short         MemoryUsage;    /* 0=minimize, 1=normal */
  short         Debug;          /* T, F */
  short         HardGenerationLimit; // T, F

  /* Starting Values (if we want to provide any) */
  double        **StartingValues; /* a matrix of starting values (each set consists of a row) */
  long          nStartingValues;  /* number of starting values */

  /* Random Number Stuff */
  short         ProvideSeeds; /* 0: no, 1: yes */
  long          UnifSeed;
  long          IntSeed;

  /* --- Ouput Diagnostics --- */
  double	*oFitValues; 
  double	*oResults; 
  double	*oGradients;
  long		oP[9];				/* operators used */
  long		oGenerations;
  long		oPeakGeneration;
  long		oPopSize;

  /* Output Files */
  char*		OutputPath;
  char*         ProjectPath;

  /* I/O types */
  short		OutputType;
  short         PrintLevel;

  /* Parallel Processing Stuff */
  short         ShareType;

  /* lexical sorting */
  long          Lexical;

  short int     UserGradient;

  /* Operator Options */
  double        P9mix;
  int           BFGSburnin;

  /* Transform Related Variables */
  short int     whichFUN;  /* 1, 2, or 3 corresponding to which evaluate function to call   */
  short int     Transform; /* 0 or 1 indicating whether transformed parameters are returned */
};

/* bfgs.c */
void dfgsmin(SEXP fn, SEXP rho,
	     double *p, int n, double gtol, int *iter, double *fret, double *hessian,
	     short int MinMax, short int BoundaryEnforcement, long InstanceNumber,
	     double **Domains, short PrintLevel, FILE *output);

/* change_order.c file */
void get_var_order(IVECTOR tot, IVECTOR cart, IMATRIX var_order);
void find_x1_x2(int tot, IMATRIX var_order, IVECTOR x1, IVECTOR x2);
void find_ac1_ac2(int t1, int t2, int t3, IVECTOR x1, IVECTOR x2, MATRIX mat, MATRIX ac1, MATRIX ac2);
void find_lu1_lu2(IVECTOR tot, IVECTOR x1, IVECTOR x2, VECTOR dom, VECTOR dom1, VECTOR dom2);
void find_limits(int tot, MATRIX domains, VECTOR llim, VECTOR ulim);
void find_new_in_eq(VECTOR a1b, MATRIX a1a2, VECTOR ll, VECTOR ul, INDEX rc, MATRIX newin);
void find_org_in_eq(VECTOR a1_b, MATRIX a1_a2, VECTOR vec_d, MATRIX c1, MATRIX c2, int c1row,
		    INDEX a1a2, MATRIX org_ineq);
void initialize(MATRIX mat, INDEX rc);
void find_final_mat1(VECTOR l2, VECTOR u2, MATRIX finmat, int row, int col);
void find_final_mat2(MATRIX newin, int r, int c, int finr, MATRIX finmat);
void find_final_mat3(MATRIX orgin, int r, int c, int finr, MATRIX finmat);

/* evaluate.c */
void optimization(struct GND_IOstructure *Structure, VECTOR X, 
		    MATRIX domains);
void sort(short int MinMax, MATRIX  population, int pop_size,
	  long nvar);
void swap(double **x, double **y);
int find_parent(IVECTOR live, int pop_size);
void assign_probab(VECTOR probab, int pop_size, double Q);
double x_pow_y(double x, int y);
void find_cum_probab(VECTOR cum_probab, VECTOR probab, int pop_size);
void find_live(VECTOR cum_probab, IVECTOR live, int pop_size, int P);
int find_die(VECTOR cum_probab, IVECTOR die, int pop_size);
void SetRunTimeParameters(struct GND_IOstructure *Structure, 
			  short FirstTime,
			  long *PopSize, long *nvars, long *MaxGenerations, long *WaitGenerations,
			  short *MinMax, short *GradientCheck, short *BoundaryEnforcement, short *UseBFGS,
			  double *SolutionTolerance,
			  long *InstanceNumber, long *P, long *P0, long *P1, long *P2, long *P3, long *P4, long *P5, 
			  long *P6, long *P7, long *P8, short *PrintLevel, 
			  short *HardGenerationLimit);
void JaIntegerOptimization(struct GND_IOstructure *Structure, VECTOR X, 
			     MATRIX domains);
void JaIntegerSort(double **InMatrix, long n, long k);
int JaIntegerCMP(double **a, double **b) ;
void JaDoubleSort(double **InMatrix, long n, long k);
int JaDoubleCMP(double **a, double **b) ;

/* frange_ran.c */
double newunif(void);
double frange_ran(double llim, double ulim);
unsigned int randint(void);
unsigned int newrand(void);

/* math.c */
/* Not needed here.  In here for completeness! */
void add(double *in1, double *in2, double *out, int row, int col);
void copy(double *in, double *target, int row, int col);
void multi(double *in1, double *in2, double *out,
	   int row1, int col1, int row2, int col2, int outrowcol[2]);
void scalarmulti(double scalar, double *in1, double *out, int row, int col) ;
void scalarmultioffdiag(double scalar, double *in1, double *out, int row, int col) ;
void subtract(double *in1, double *in2, double *out, int row, int col);
double trace(double *a, int n);
void transpose(double *orig_matrix, double *t_matrix, int orig_rows, int orig_columns);
void copy_matrix(MATRIX mat1, MATRIX mat2, int lr, int ur, int lc, int uc);
int Iround(double in);
void samplestats(double **obsdata, int numobsv, int novarsv, int weightflag, 
		 double *weightdata);
void populationstats(double **popdata, int numobsv, int novarsv, 
		     double *mean, double *var, double *skew, double *kur,
		     long *tobs);

/* multiply.c */
void mmprod(int m, int nm, int n, MATRIX mul_cm, MATRIX mul_am, MATRIX mul_bm);
void mvprod(int m, int nm, VECTOR cm, MATRIX am, VECTOR bm);

/* numerics.c */
double **JaMatrixAllocate(long n, long k);
void JaMatrixFree(double **M, long k);
short **JaShortMatrixAllocate(long nobs, long nvars);
void JaShortMatrixFree(double **M, long nobs);
MATRIX matrix(int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
double *Gvector(int nl, int nh);
int **imatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
void free_vector(double *v, int nl);
void free_ivector(int *v, int nl);
void free_matrix(double **m, int nrl, int nrh, int ncl);
void free_imatrix(int **m, int nrl, int nrh, int ncl);


/* operators.c file */
void oper1(VECTOR parent, double **domains, int nvars);
void oper2(VECTOR parent, double **domains, int nvars);
void oper3(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void oper4(MATRIX p, int p2use, int nvars);
void oper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars);
void oper6(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void oper7(VECTOR p1, VECTOR p2, double **domains, int nvars);
void oper8(SEXP fn, SEXP rho,
	   VECTOR parent, MATRIX domains, 
	   double SolutionTolerance, long nvars, 
	   short BoundaryEnforcement, 
	   short PrintLevel, double mix);
void find_range(double *llim, double *ulim, int comp, double **domains, int nvars, VECTOR parent);
void find_rangeInt(int *llim, int *ulim, int comp, double **domains, int nvars, VECTOR parent);
int irange_ran(int llim, int ulim);
double get_F(int T, int t, double y, int B);
void JaIntegerOper1(VECTOR parent, double **domains, int nvars);
void JaIntegerOper2(VECTOR parent, double **domains, int nvars);
void JaIntegerOper3(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void JaIntegeroper4(MATRIX p, int p2use, int nvars, double **domains);
void JaIntegerOper5(VECTOR p1, VECTOR p2, int STEP, double **domains, int nvars);
void JaIntegerOper6(VECTOR parent, double **domains, int nvars, int T, int t, int B);
void JaIntegerOper7(VECTOR p1, VECTOR p2, double **domains, int nvars);
FLAG InBounds(VECTOR child, double **domains, int nvars);

/*print_format.c */
long ReadPopulation(double **Data, long NewPopSize, long NewVars, FILE *fp, short PrintLevel);
void print_domains(MATRIX equal, int t_equ, short DataType);
void print_population(int popsize, int nvars, int generation, int
		      lexical, double **foo, FILE *out);


