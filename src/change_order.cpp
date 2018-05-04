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

  $Header: /home/jsekhon/xchg/genoud/rgenoud.distribution/sources/RCS/change_order.cpp,v 2.15 2005/10/29 06:14:44 jsekhon Exp jsekhon $

*/


#include "genoud.h"


void get_var_order(IVECTOR tot, IVECTOR cart, IMATRIX var_order)
     /* IVECTOR tot,       array with total number of variables and equalities
        cart;              array with the subscript of the variables to be eliminated
	IMATRIX var_order; matrix with the variables in one column and a tag in the 
	                   other to identify as to whether it is to be eliminated 
			   or not
     */
{
  int i;

  for(i=1; i<=tot[0]; i++)
    {
      var_order[i][1] = i;
      var_order[i][2] = 0;
    }

  for(i=1; i<=tot[1]; i++)
    var_order[cart[i]][2] = 1;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_x1_x2()                                 */
/*                                                                              */
/*           SYNOPSIS          :   void find_x1_x2(tot,var_order,eq_co,x1,x2)   */
/*                                                                              */
/*           DESCRIPTION       :   This function splits the original vector of  */
/*                                  variables into x1 and x2, where x1 consists */
/*                                  of the 'p' variables to be eliminated and   */
/*                                  x2 consists of the remaining variables      */
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




void find_x1_x2(int tot, IMATRIX var_order, IVECTOR x1, IVECTOR x2)
     /*
       int tot;            total number of variables
       IMATRIX var_order;  array of variables with tag identifying them 
                           to be eliminated
       IVECTOR x1,         array of variables to be eliminated
       x2;                 array of remaining variables
     */
{
  int i,j=1,k=1;

  for(i=1; i<=tot; i++)
    {
      if(var_order[i][2] == 1)
        x1[j++] = var_order[i][1];
      else
        x2[k++] = var_order[i][1];
    }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_ac1_ac2()                               */
/*                                                                              */
/*           SYNOPSIS          :   void find_ac1_ac2(t1,t2,t3,x1,x2,mat,ac1,ac2)*/
/*                                                                              */
/*           DESCRIPTION       :   This function splits the original equality   */
/*                                  or the inequality matrix into two matrices; */
/*                                  the p-equalities or the inequalities and the*/
/*                                  remaining part of the matrix                */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void find_ac1_ac2(int t1, int t2, int t3, IVECTOR x1, IVECTOR x2, MATRIX mat, MATRIX ac1, MATRIX ac2)
     /* 
	int t1,t2,t3;
	IVECTOR x1,x2;   the variables corresponding to the split matrices
	MATRIX  mat,     the original matrix to be split
        ac1,ac2;         the split matrices
     */
{
  int i,j;

  for(i=1; i<=t1; i++)
    for(j=1; j<=t2; j++)
      ac1[j][i] = mat[j][x1[i]];

  for(i=1; i<=t3; i++)
    for(j=1; j<=t2; j++)
      ac2[j][i] = mat[j][x2[i]];
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_lu1_lu2()                               */
/*                                                                              */
/*           SYNOPSIS          :   void find_lu1_lu2_(tot,x1,x2,dom,dom1,dom2)  */
/*                                                                              */
/*           DESCRIPTION       :   This function splits the lower or the upper  */
/*                                  bounds of the total domain constraints into */
/*                                  two groups, one with the domain constraints */
/*                                  for the variables to be eliminated and the  */
/*                                  other one for the remaining variables       */
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




void find_lu1_lu2(IVECTOR tot, IVECTOR x1, IVECTOR x2, VECTOR dom, VECTOR dom1, VECTOR dom2)
     /*
       IVECTOR tot,x1,x2;
       VECTOR  dom,        the original array of the lower or the upper bounds
       dom1,dom2;          the original bounds split in to two groups
     */
{
  int i;

  for(i=1; i<=tot[1]; i++)
      dom1[i] = dom[x1[i]];

  for(i=1; i<=tot[0]-tot[1]; i++)
      dom2[i] = dom[x2[i]];
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_limits()                                */
/*                                                                              */
/*           SYNOPSIS          :   void find_limits(tot,domains,llim,ulim)      */
/*                                                                              */
/*           DESCRIPTION       :   This function forms seperate arrays for the  */
/*                                  lower and upper limits of the domains from  */
/*                                  the original domain constraints read already*/
/*                                  from the input file                         */
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




void find_limits(int tot, MATRIX domains, VECTOR llim, VECTOR ulim)
     /*
       int tot;
       MATRIX domains;     matrix containing the domain variables and the limits
       VECTOR llim,ulim;     vectors of lower and upper limits
     */
{
  int i;

  for(i=1; i<=tot; i++)
    {
      llim[i] = domains[i][1];
      ulim[i] = domains[i][3];
    }
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_new_in_eq()                             */
/*                                                                              */
/*           SYNOPSIS          :   void find_new_in_eq(a1b.a1a2,ll,ul,rc,newin) */
/*                                                                              */
/*           DESCRIPTION       :   This function converts the original          */
/*                                  equality constraints into domain            */
/*                                  constraints eliminating the p-equalities    */
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




void find_new_in_eq(VECTOR a1b, MATRIX a1a2, VECTOR ll, VECTOR ul, INDEX rc, MATRIX newin)
     /*
       VECTOR a1b,    the product of a-inverse and b
       ll,ul;         upper and lower limits of the domain constraints
       MATRIX a1a2,   the products of a-inverse and matrix a2
       newin;         the final matrix with domain constraints
       INDEX rc;
     */
{
  int i,j;

  for(i=1; i<=rc.r; i++)
    for(j=1; j<=rc.c; j++)
      if(j==1)
        newin[i][j] = ll[i] - a1b[i];/*eliminating the constrants from the*/
      else if(j==rc.c)                      /*equations in the domain constraints*/
        newin[i][j] = ul[i] - a1b[i];
      else
        newin[i][j] = 0 - a1a2[i][j-1];
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_org_in_eq()                             */
/*                                                                              */
/*           SYNOPSIS          :   void find_org_in_eq(a1_b,a1_a2,vec_d,c1,c2,  */
/*                                                          c1row,a1a2,org_ineq)*/
/*                                                                              */
/*           DESCRIPTION       :   This function converts the original          */
/*                                  inequality constraints into domain          */
/*                                  constraints, with the variables other than  */
/*                                  the p-variables eliminated                  */
/*                                                                              */
/*           FUNCTIONS CALLED  :   matrix().                                    */
/*                                 mmprod(),                                    */
/*                                 mvprod(),                                    */
/*                                 Gvector(),                                    */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void find_org_in_eq(VECTOR a1_b, MATRIX a1_a2, VECTOR vec_d, MATRIX c1, MATRIX c2, int c1row,
		    INDEX a1a2, MATRIX org_ineq)
     /*
       VECTOR vec_d,       the RHS constant vector
       a1_b;               product of a1-inverse and b
       MATRIX a1_a2,       product of a1-inverse and a2
       org_ineq,           the converted inequalities into domains
       c1,c2;              p_inequalities and the remaining inequalities
       int c1row;
       INDEX a1a2;         rows and columns of the matrices
     */
{
  int i,j;
  VECTOR temp;
  MATRIX mat;

  temp = Gvector(1,c1row);
  mat = matrix(1,c1row,1,a1a2.c-1);

  mvprod(c1row,a1a2.r,temp,c1,a1_b);/*matrix, vector product C1.invA1.b*/
  mmprod(c1row,a1a2.r,a1a2.c-1,mat,c1,a1_a2);/*matrix, matrix product*/

  for(i=1; i<=c1row; i++)
    for(j=1; j<=a1a2.c; j++)
      {
        if (j==a1a2.c)
          org_ineq[i][j] = vec_d[i] - temp[i];
        else
          org_ineq[i][j] = c2[i][j] - mat[i][j];
      }
  free_vector(temp,1);
  free_matrix(mat,1,c1row,1);
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   initialize()                                 */
/*                                                                              */
/*           SYNOPSIS          :   void initialize(mat,rc)                      */
/*                                                                              */
/*           DESCRIPTION       :   This function initializes all the components */
/*                                  of the given matrix to zero                 */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void initialize(MATRIX mat, INDEX rc)
     /*
       MATRIX mat;
       INDEX rc;
     */
{
  int i,j;

  for(i=1; i<=rc.r; i++)
    for(j=1; j<=rc.c; j++)
      mat[i][j] = 0.0;
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_final_mat1()                            */
/*                                                                              */
/*           SYNOPSIS          :   void find_final_mat1(x2,l2,u2,finmat,row,col)*/
/*                                                                              */
/*           DESCRIPTION       :   This function copies the remaining original  */
/*                                  domain constraints on to the final matrix   */
/*                                  to be printed                               */
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




void find_final_mat1(VECTOR l2, VECTOR u2, MATRIX finmat, int row, int col)
     /*
       MATRIX finmat;
       VECTOR l2,u2;
       int row,col;
     */
{
  int i,j=2;

  for(i=1; i<=row; i++)
    {
      finmat[i][1] = l2[i];
      finmat[i][col] = u2[i];
      finmat[i][j++] = 1.0;
    }
}

/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_final_mat2()                            */
/*                                                                              */
/*           SYNOPSIS          :   void find_final_mat2(newin,r,c,finr,finmat) */
/*                                                                              */
/*           DESCRIPTION       :   This function appends the new inequalities   */
/*                                  got from the original equalities, on to the */
/*                                  final matrix to be printed                  */
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




void find_final_mat2(MATRIX newin, int r, int c, int finr, MATRIX finmat)
     /*
       MATRIX newin,finmat;
       int r,c,finr;
     */
{
  int i,j;

  for(i=1; i<=r; i++)
    {
      for(j=1; j<=c; j++)
        finmat[finr][j] = newin[i][j];
      finr++;
    }
}


/********************************************************************************/
/*                                                                              */
/*           FUNCTION NAME     :   find_final_mat3()                            */
/*                                                                              */
/*           SYNOPSIS          :   void find_final_mat3(orgin,r,c,finr,finmat)  */
/*                                                                              */
/*           DESCRIPTION       :   This function appends the inequalities, which*/
/*                                  were converted into domain constraints, on  */
/*                                  to the final matrix to be printed           */
/*                                                                              */
/*           FUNCTIONS CALLED  :   None                                         */
/*                                                                              */
/*           CALLING FUNCITONS :   main()                                       */
/*                                                                              */
/*                                                                              */
/*                                                                              */
/*           REV            DATE            BY           DESCRIPTION            */
/*           ---            ----            --           -----------            */
/*                                                                              */
/*                                                                              */
/********************************************************************************/




void find_final_mat3(MATRIX orgin, int r, int c, int finr, MATRIX finmat)
     /*
       MATRIX orgin,finmat;
       int r,c,finr;
     */
{
  int i,j;

  for(i=1; i<=r; i++)
    {
      finmat[finr][1] = MIN;
      for(j=1; j<=c; j++)
        finmat[finr][j+1] = orgin[i][j];
      finr++;
    }
}
