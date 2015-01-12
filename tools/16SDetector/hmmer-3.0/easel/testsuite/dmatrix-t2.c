/* dmatrix-t2
 * Test of matrix inversion.
 * 
 * 
 * CVS $Id$
 * SRE, Mon Jul 12 14:37:46 2004
 */

#include <stdio.h>
#include <stdlib.h>
#include <easel/esl_core.h>
#include <easel/esl_dmatrix.h>

#define SIZE      4
#define ALPHABET  "1234"

int
main(int argc, char **argv)
{
  int               verbose;
  int               status;
  ESL_DMATRIX      *A;
  ESL_DMATRIX      *Ai;

  verbose = (argc > 1) ? TRUE:FALSE;

  A  = esl_dmx_Alloc(SIZE, SIZE);
  Ai = esl_dmx_Alloc(SIZE, SIZE);

  /* Make A the identity matrix.
   */
  esl_dmx_SetIdentity(A);

  /* Invert it.
   * This is a trivial test; the code has to be horribly
   * broken to fail.
   *  I^-1 = I
   */
  esl_dmx_Invert(A, Ai);

  /* Print stuff for inspection, if we're supposed to.
   */
  if (verbose)
    {
      printf("\nA:\n");
      esl_dmx_fprintf_alphalabeled(stdout, A, ALPHABET);

      printf("\nA^{-1}:\n");
      esl_dmx_fprintf_alphalabeled(stdout, Ai, ALPHABET);
    }

  /* The test: A, Ai oughta both be I.
   */
  if (! esl_dmx_MatricesEqual(A, Ai, 0.001)) 
    {
      printf("FAIL\n");
      status = ESL_ETESTFAIL;
    }
  else
    {
      printf("ok\n");
      status = ESL_OK;
    }

  esl_dmx_Free(A);
  esl_dmx_Free(Ai);
  return status;
}
