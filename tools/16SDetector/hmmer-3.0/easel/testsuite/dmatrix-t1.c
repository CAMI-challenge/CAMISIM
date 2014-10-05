/* dmatrix-t1
 * Test of LUP decomposition.
 * 
 * 
 * CVS $Id$
 * SRE, Mon Jul 12 13:15:16 2004 [St. Louis]
 */

#include <stdio.h>
#include <stdlib.h>
#include <easel/esl_core.h>
#include <easel/esl_dmatrix.h>

#define SIZE     4
#define ALPHABET "1234"

/* This example is from Cormen, Leiserson, Rivest;
 * p.753 of second edition gives L,U at each step.
 * 
 * Final result should be
 * 
 *     0 0 1 0    2   0   2 0.6        1   0   0   0     5   5   4    2
 *     1 0 0 0    3   3   4  -2      0.4   1   0   0     0  -2 0.4 -0.2
 *     0 0 0 1    5   5   4   2     -0.2 0.5   1   0     0   0   4 -0.5
 *     0 1 0 0   -1  -2 3.4  -1      0.6   0 0.4   1     0   0   0   -3
 *        P            A           =        L                  U
 */ 
double testdata[SIZE][SIZE] = { 
  { 2., 0., 2., 0.6, },
  { 3., 3., 4., -2., },
  { 5., 5., 4., 2., },
  { -1.,-2.,3.4,-1., },
};


int
main(int argc, char **argv)
{
  int               verbose;
  ESL_DMATRIX      *A;
  ESL_DMATRIX      *LU;
  ESL_DMATRIX      *L;
  ESL_DMATRIX      *U;
  ESL_DMATRIX      *PA;
  ESL_DPERMUTATION *P;
  int               i,j;
  int               status;

  verbose = (argc > 1) ? TRUE:FALSE;

  A  = esl_dmx_Alloc(SIZE, SIZE);
  LU = esl_dmx_Alloc(SIZE, SIZE);
  L  = esl_dmx_Alloc(SIZE, SIZE);
  U  = esl_dmx_Alloc(SIZE, SIZE);
  P  = esl_permutation_Alloc(SIZE);
  PA = esl_dmx_Alloc(SIZE, SIZE);

  /* Set up matrix A.
   */
  for (i = 0; i < SIZE; i++)
    for (j = 0; j < SIZE; j++)
      A->mx[i][j] = testdata[i][j];

  /* Make a copy in LU, then compute PA = LU
   */
  esl_dmx_Copy(A, LU);
  esl_dmx_LUP_decompose(LU, P);

  /* Permute A to get PA...
   */
  esl_permute_PA(P, A, PA);
  
  /* Calculate the product LU by multiplying L,U together again...
   */
  esl_dmx_LU_separate(LU, L, U);
  esl_dmx_Multiply(L, U, LU);

  /* Print stuff for inspection, if we're supposed to.
   */
  if (verbose)
    {
      printf("\nA:\n");
      esl_dmx_fprintf_alphalabeled(stdout, A, ALPHABET);

      printf("\nL:\n");
      esl_dmx_fprintf_alphalabeled(stdout, L, ALPHABET);

      printf("\nU:\n");
      esl_dmx_fprintf_alphalabeled(stdout, U, ALPHABET);

      printf("\nP:\n");
      esl_permutation_fprintf_numlabeled(stdout, P);

      printf("\nPA:\n");
      esl_dmx_fprintf_alphalabeled(stdout, PA, ALPHABET);

      printf("\nLU:\n");
      esl_dmx_fprintf_alphalabeled(stdout, LU, ALPHABET);
    }

  /* The test: PA, LU oughta be equal.
   */
  if (! esl_dmx_MatricesEqual(PA, LU, 0.001)) 
    {
      fprintf(stderr, "fail: PA,LU not equal\n");
      exit(ESL_ETESTFAIL);
    }
  else
    status = ESL_OK;
  
  esl_dmx_Free(A);
  esl_dmx_Free(LU);
  esl_dmx_Free(L);
  esl_dmx_Free(U);
  esl_dmx_Free(PA);
  esl_permutation_Free(P);
  return ESL_OK;
}
