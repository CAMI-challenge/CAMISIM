/* Routines for manipulating sequence alignment score matrices,
 * such as the BLOSUM and PAM matrices.
 * 
 * Contents:
 *   1. The ESL_SCOREMATRIX object.
 *   2. Reading/writing score matrices.
 *   3. Interpreting score matrices probabilistically.
 *   4. Utility programs.
 *   5. Unit tests.
 *   6. Test driver.
 *   7. Example program.
 *   8. License and copyright.
 * 
 * SRE, Mon Apr  2 08:25:05 2007 [Janelia]
 * SVN $Id: esl_scorematrix.c 337 2009-05-12 02:13:02Z eddys $
 */

#include "esl_config.h"

#include <string.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_fileparser.h"
#include "esl_rootfinder.h"
#include "esl_ratematrix.h"
#include "esl_scorematrix.h"
#include "esl_vectorops.h"

/*****************************************************************
 * 1. The ESL_SCOREMATRIX object
 *****************************************************************/

/* Function:  esl_scorematrix_Create()
 * Synopsis:  Create an <ESL_SCOREMATRIX>.
 * Incept:    SRE, Mon Apr  2 08:38:10 2007 [Janelia]
 *
 * Purpose:   Allocates a score matrix for alphabet <abc>, initializes
 *            all scores to zero.
 *
 * Args:      abc   - pointer to digital alphabet 
 *
 * Returns:   a pointer to the new object.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SCOREMATRIX *
esl_scorematrix_Create(const ESL_ALPHABET *abc)
{
  int status;
  int i;
  ESL_SCOREMATRIX *S = NULL;

  ESL_ALLOC(S, sizeof(ESL_SCOREMATRIX));
  S->s          = NULL;
  S->K          = abc->K;
  S->Kp         = abc->Kp;
  S->isval      = NULL;
  S->abc_r      = abc;
  S->nc         = 0;
  S->outorder   = NULL;
  S->name       = NULL;
  S->path       = NULL;

  ESL_ALLOC(S->s, sizeof(int *) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->s[i] = NULL;
  ESL_ALLOC(S->isval, sizeof(char) * abc->Kp);
  for (i = 0; i < abc->Kp; i++) S->isval[i] = FALSE;
  ESL_ALLOC(S->outorder, sizeof(char) * abc->Kp);
  S->outorder[0] = '\0';		/* init to empty string. */

  ESL_ALLOC(S->s[0], sizeof(int) * abc->Kp * abc->Kp);
  for (i = 1; i < abc->Kp; i++) S->s[i] = S->s[0] + abc->Kp * i;

  for (i = 0; i < abc->Kp*abc->Kp; i++) S->s[0][i] = 0;
  return S;

 ERROR:
  esl_scorematrix_Destroy(S);
  return NULL;
}

/* Function:  esl_scorematrix_SetIdentity()
 * Synopsis:  Set matrix to +1 match, 0 mismatch.
 * Incept:    SRE, Mon Apr 16 20:17:00 2007 [Janelia]
 *
 * Purpose:   Sets score matrix <S> to be +1 for a match, 
 *            0 for a mismatch. <S> may be for any alphabet.
 *            
 *            Rarely useful in real use, but may be useful to create
 *            simple examples (including debugging).
 *
 * Returns:   <eslOK> on success, and the scores in <S> are set.
 */
int
esl_scorematrix_SetIdentity(ESL_SCOREMATRIX *S)
{
  int a;
  int x;

  for (a = 0; a < S->abc_r->Kp*S->abc_r->Kp; a++) S->s[0][a] = 0;
  for (a = 0; a < S->K; a++)                      S->s[a][a] = 1;

  for (x = 0;           x < S->K;  x++)      S->isval[x] = TRUE;
  for (x = S->abc_r->K; x < S->Kp; x++)      S->isval[x] = FALSE;
  
  strncpy(S->outorder, S->abc_r->sym, S->K);  
  S->outorder[S->K] = '\0';
  S->nc             = S->K;
  return eslOK;
}

/* Function:  esl_scorematrix_SetBLOSUM62
 * Synopsis:  Set matrix to BLOSUM62 scores.
 * Incept:    SRE, Tue Apr  3 13:22:03 2007 [Janelia]
 *
 * Purpose:   Set the 20x20 canonical residue scores in an 
 *            allocated amino acid score matrix <S> to BLOSUM62
 *            scores \citep{Henikoff92}.
 *
 * Returns:   <eslOK> on success, and the scores in <S> are set.
 * 
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S)
{
  int x,y;
  static int blosum62[29][29] = {
    /*  A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    -    B    J    Z    O    U    X    *    ~  */
    {   4,   0,  -2,  -1,  -2,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -1,   1,   0,   0,  -3,  -2,   0,  -2,   0,  -1,   0,   0,   0,  -4,   0,  }, /* A */
    {   0,   9,  -3,  -4,  -2,  -3,  -3,  -1,  -3,  -1,  -1,  -3,  -3,  -3,  -3,  -1,  -1,  -1,  -2,  -2,   0,  -3,   0,  -3,   0,   0,  -2,  -4,   0,  }, /* C */
    {  -2,  -3,   6,   2,  -3,  -1,  -1,  -3,  -1,  -4,  -3,   1,  -1,   0,  -2,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,  -4,   0,  }, /* D */
    {  -1,  -4,   2,   5,  -3,  -2,   0,  -3,   1,  -3,  -2,   0,  -1,   2,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,  -4,   0,  }, /* E */
    {  -2,  -2,  -3,  -3,   6,  -3,  -1,   0,  -3,   0,   0,  -3,  -4,  -3,  -3,  -2,  -2,  -1,   1,   3,   0,  -3,   0,  -3,   0,   0,  -1,  -4,   0,  }, /* F */
    {   0,  -3,  -1,  -2,  -3,   6,  -2,  -4,  -2,  -4,  -3,   0,  -2,  -2,  -2,   0,  -2,  -3,  -2,  -3,   0,  -1,   0,  -2,   0,   0,  -1,  -4,   0,  }, /* G */
    {  -2,  -3,  -1,   0,  -1,  -2,   8,  -3,  -1,  -3,  -2,   1,  -2,   0,   0,  -1,  -2,  -3,  -2,   2,   0,   0,   0,   0,   0,   0,  -1,  -4,   0,  }, /* H */
    {  -1,  -1,  -3,  -3,   0,  -4,  -3,   4,  -3,   2,   1,  -3,  -3,  -3,  -3,  -2,  -1,   3,  -3,  -1,   0,  -3,   0,  -3,   0,   0,  -1,  -4,   0,  }, /* I */
    {  -1,  -3,  -1,   1,  -3,  -2,  -1,  -3,   5,  -2,  -1,   0,  -1,   1,   2,   0,  -1,  -2,  -3,  -2,   0,   0,   0,   1,   0,   0,  -1,  -4,   0,  }, /* K */
    {  -1,  -1,  -4,  -3,   0,  -4,  -3,   2,  -2,   4,   2,  -3,  -3,  -2,  -2,  -2,  -1,   1,  -2,  -1,   0,  -4,   0,  -3,   0,   0,  -1,  -4,   0,  }, /* L */
    {  -1,  -1,  -3,  -2,   0,  -3,  -2,   1,  -1,   2,   5,  -2,  -2,   0,  -1,  -1,  -1,   1,  -1,  -1,   0,  -3,   0,  -1,   0,   0,  -1,  -4,   0,  }, /* M */
    {  -2,  -3,   1,   0,  -3,   0,   1,  -3,   0,  -3,  -2,   6,  -2,   0,   0,   1,   0,  -3,  -4,  -2,   0,   3,   0,   0,   0,   0,  -1,  -4,   0,  }, /* N */
    {  -1,  -3,  -1,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -2,  -2,   7,  -1,  -2,  -1,  -1,  -2,  -4,  -3,   0,  -2,   0,  -1,   0,   0,  -2,  -4,   0,  }, /* P */
    {  -1,  -3,   0,   2,  -3,  -2,   0,  -3,   1,  -2,   0,   0,  -1,   5,   1,   0,  -1,  -2,  -2,  -1,   0,   0,   0,   3,   0,   0,  -1,  -4,   0,  }, /* Q */
    {  -1,  -3,  -2,   0,  -3,  -2,   0,  -3,   2,  -2,  -1,   0,  -2,   1,   5,  -1,  -1,  -3,  -3,  -2,   0,  -1,   0,   0,   0,   0,  -1,  -4,   0,  }, /* R */
    {   1,  -1,   0,   0,  -2,   0,  -1,  -2,   0,  -2,  -1,   1,  -1,   0,  -1,   4,   1,  -2,  -3,  -2,   0,   0,   0,   0,   0,   0,   0,  -4,   0,  }, /* S */
    {   0,  -1,  -1,  -1,  -2,  -2,  -2,  -1,  -1,  -1,  -1,   0,  -1,  -1,  -1,   1,   5,   0,  -2,  -2,   0,  -1,   0,  -1,   0,   0,   0,  -4,   0,  }, /* T */
    {   0,  -1,  -3,  -2,  -1,  -3,  -3,   3,  -2,   1,   1,  -3,  -2,  -2,  -3,  -2,   0,   4,  -3,  -1,   0,  -3,   0,  -2,   0,   0,  -1,  -4,   0,  }, /* V */
    {  -3,  -2,  -4,  -3,   1,  -2,  -2,  -3,  -3,  -2,  -1,  -4,  -4,  -2,  -3,  -3,  -2,  -3,  11,   2,   0,  -4,   0,  -3,   0,   0,  -2,  -4,   0,  }, /* W */
    {  -2,  -2,  -3,  -2,   3,  -3,   2,  -1,  -2,  -1,  -1,  -2,  -3,  -1,  -2,  -2,  -2,  -1,   2,   7,   0,  -3,   0,  -2,   0,   0,  -1,  -4,   0,  }, /* Y */
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  }, /* - */
    {  -2,  -3,   4,   1,  -3,  -1,   0,  -3,   0,  -4,  -3,   3,  -2,   0,  -1,   0,  -1,  -3,  -4,  -3,   0,   4,   0,   1,   0,   0,  -1,  -4,   0,  }, /* B */
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  }, /* J */
    {  -1,  -3,   1,   4,  -3,  -2,   0,  -3,   1,  -3,  -1,   0,  -1,   3,   0,   0,  -1,  -2,  -3,  -2,   0,   1,   0,   4,   0,   0,  -1,  -4,   0,  }, /* Z */
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  }, /* O */
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  }, /* U */
    {   0,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,  -1,  -1,   0,   0,  -1,  -2,  -1,   0,  -1,   0,  -1,   0,   0,  -1,  -4,   0,  }, /* X */
    {  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,   0,  -4,   0,  -4,   0,   0,  -4,   1,   0,  }, /* * */
    {   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  }, /* ~ */
  };
  /* The BLOSUM62 background frequencies are the actual frequencies used to create
   * the matrix in 1992. */
  /*                           A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y */
  /* double blosum62f[20] = { 0.074, 0.025, 0.054, 0.054, 0.047, 0.074, 0.026, 0.068, 0.058, 0.099, 0.025, 0.045, 0.039, 0.034, 0.052, 0.057, 0.051, 0.073, 0.013, 0.032 };
   */

  for (x = 0;           x < S->K;  x++)      S->isval[x] = TRUE;
  for (x = S->abc_r->K; x < S->Kp; x++)      S->isval[x] = FALSE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'B'); S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'Z'); S->isval[x] = TRUE;
  x = esl_abc_DigitizeSymbol(S->abc_r, 'X'); S->isval[x] = TRUE;
    
  for (x = 0; x < S->Kp; x++)
    for (y = 0; y < S->Kp; y++)
      S->s[x][y] = blosum62[x][y];

  /* Bookkeeping necessary to be able to reproduce BLOSUM62 output format exactly, if we need to Write() */
  strcpy(S->outorder, "ARNDCQEGHILKMFPSTWYVBZX*");
  S->nc         = strlen(S->outorder);

  if (esl_strdup("BLOSUM62", -1, &(S->name)) != eslOK) return eslEMEM;
  return eslOK;
}


/* Function:  esl_scorematrix_SetWAG()
 * Synopsis:  Set matrix using the WAG evolutionary model.           
 * Incept:    SRE, Thu Apr 12 13:23:28 2007 [Janelia]
 *
 * Purpose:   Parameterize an amino acid score matrix <S> using the WAG
 *            rate matrix \citep{WhelanGoldman01} as the underlying
 *            evolutionary model, at a distance of <t>
 *            substitutions/site, with scale factor <lambda>.
 *
 * Args:      S      - score matrix to set parameters in. Must be created for
 *                     an amino acid alphabet.
 *            lambda - scale factor for scores     
 *            t      - distance to exponentiate WAG to, in substitutions/site         
 *                 
 * Returns:   <eslOK> on success, and the 20x20 residue scores in <S> are set.
 *
 * Throws:    <eslEINVAL> if <S> isn't an allocated amino acid score matrix.
 *            <eslEMEM> on allocation failure.
 */
int
esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, double lambda, double t)
{
  int status;
  int i,j;
  ESL_DMATRIX *Q = NULL;
  ESL_DMATRIX *P = NULL;
  static double wagpi[20];

  if (S->K != 20) ESL_EXCEPTION(eslEINVAL, "Must be using an amino acid alphabet (K=20) to make WAG-based matrices");

  if (( Q = esl_dmatrix_Create(20, 20)) == NULL)  goto ERROR;
  if (( P = esl_dmatrix_Create(20, 20)) == NULL)  goto ERROR;
  if ( esl_composition_WAG(wagpi)       != eslOK) goto ERROR;
  if ( esl_rmx_SetWAG(Q, wagpi)         != eslOK) goto ERROR;
  if ( esl_dmx_Exp(Q, t, P)             != eslOK) goto ERROR;

  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      P->mx[i][j] *= wagpi[i];	/* P_ij = P(j|i) pi_i */
  
  esl_scorematrix_SetFromProbs(S, lambda, P, wagpi, wagpi);

  if ((status = esl_strdup("WAG", -1, &(S->name))) != eslOK) goto ERROR;

  esl_dmatrix_Destroy(Q);
  esl_dmatrix_Destroy(P);
  return eslOK;

 ERROR:
  if (Q != NULL) esl_dmatrix_Destroy(Q);
  if (Q != NULL) esl_dmatrix_Destroy(P);
  return status;
}


/* Function:  esl_scorematrix_SetFromProbs()
 * Synopsis:  Set matrix from target and background probabilities.
 * Incept:    SRE, Wed Apr 11 17:37:45 2007 [Janelia]
 *
 * Purpose:   Sets the scores in a new score matrix <S> from target joint
 *            probabilities in <P>, query background probabilities <fi>, and 
 *            target background probabilities <fj>, with scale factor <lambda>:
 *                 $s_{ij} = \frac{1}{\lambda} \frac{p_{ij}}{f_i f_j}$.
 *                 
 *            Size of everything must match the canonical alphabet
 *            size in <S>. That is, <S->abc->K> is the canonical
 *            alphabet size of <S>; <P> must contain $K times K$
 *            probabilities $P_{ij}$, and <fi>,<fj> must be vectors of
 *            K probabilities. All probabilities must be nonzero.
 *            
 * Args:      S      - score matrix to set scores in
 *            lambda - scale factor     
 *            P      - matrix of joint probabilities P_ij (KxK)
 *            fi     - query background probabilities (0..K-1)
 *            fj     - target background probabilities 
 *
 * Returns:   <eslOK> on success, and <S> contains the calculated score matrix.
 */
int
esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, double lambda, const ESL_DMATRIX *P, const double *fi, const double *fj)
{
  int    i,j;
  double sc;
  
  for (i = 0; i < S->abc_r->K; i++)
    for (j = 0; j < S->abc_r->K; j++)
      {
	sc = log(P->mx[i][j] / (fi[i] * fj[j])) / lambda;
	S->s[i][j] = (int) (sc + (sc>0 ? 0.5 : -0.5)); /* that's rounding to the nearest integer */
      }

  for (i = 0; i < S->abc_r->K; i++)
    S->isval[i] = TRUE;
  S->nc = S->abc_r->K;

  strncpy(S->outorder, S->abc_r->sym, S->abc_r->K);
  S->outorder[S->nc] = '\0';
  return eslOK;
}


/* Function:  esl_scorematrix_Copy()
 * Synopsis:  Copy <src> matrix to <dest>.
 * Incept:    SRE, Tue May 15 10:24:20 2007 [Janelia]
 *
 * Purpose:   Copy <src> score matrix into <dest>. Caller
 *            has allocated <dest> for the same alphabet as
 *            <src>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINCOMPAT> if <dest> isn't allocated for
 *            the same alphabet as <src>.
 *            <eslEMEM> on allocation error.
 */
int
esl_scorematrix_Copy(const ESL_SCOREMATRIX *src, ESL_SCOREMATRIX *dest)
{
  int i,j;
  int status;

  if (src->abc_r->type != dest->abc_r->type || src->K != dest->K || src->Kp != dest->Kp)
    ESL_EXCEPTION(eslEINCOMPAT, "source and dest score matrix types don't match");

  for (i = 0; i < src->K; i++)
    for (j = 0; j < src->K; j++)
      dest->s[i][j] = src->s[i][j];
  for (i = 0; i < src->Kp; i++)
    dest->isval[i] = src->isval[i];
  dest->nc = src->nc;
  for (i = 0; i < src->nc; i++)
    dest->outorder[i] = src->outorder[i];
  dest->outorder[dest->nc] = '\0';

  if ((status = esl_strdup(src->name, -1, &(dest->name))) != eslOK) return status;
  if ((status = esl_strdup(src->path, -1, &(dest->path))) != eslOK) return status;
  return eslOK;
}

/* Function:  esl_scorematrix_Clone()
 * Synopsis:  Allocate a duplicate of a matrix. 
 * Incept:    SRE, Tue May 15 10:24:20 2007 [Janelia]
 *
 * Purpose:   Allocates a new matrix and makes it a duplicate
 *            of <S>. Return a pointer to the new matrix.
 *
 * Throws:    <NULL> on allocation failure.
 */
ESL_SCOREMATRIX *
esl_scorematrix_Clone(const ESL_SCOREMATRIX *S)
{
  ESL_SCOREMATRIX *dup = NULL;

  if ((dup = esl_scorematrix_Create(S->abc_r)) == NULL)  return NULL;
  if (esl_scorematrix_Copy(S, dup)             != eslOK) { esl_scorematrix_Destroy(dup); return NULL; }
  return dup;
}


/* Function:  esl_scorematrix_Compare()
 * Synopsis:  Compare two matrices for equality.
 * Incept:    SRE, Tue Apr  3 14:17:12 2007 [Janelia]
 *
 * Purpose:   Compares two score matrices. Returns <eslOK> if they 
 *            are identical, <eslFAIL> if they differ. Every aspect
 *            of the two matrices is compared.
 *            
 *            The annotation (name, filename path) are not
 *            compared; we may want to compare an internally
 *            generated scorematrix to one read from a file.
 */
int
esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2)
{
  int a,b;

  if (strcmp(S1->outorder, S2->outorder) != 0) return eslFAIL;
  if (S1->nc         != S2->nc)                return eslFAIL;
  
  for (a = 0; a < S1->nc; a++)
    if (S1->isval[a] != S2->isval[a])          return eslFAIL;
  
  for (a = 0; a < S1->Kp; a++)
    for (b = 0; b < S1->Kp; b++)
      if (S1->s[a][b] != S2->s[a][b]) return eslFAIL;

  return eslOK;
}

/* Function:  esl_scorematrix_CompareCanon()
 * Synopsis:  Compares scores of canonical residues for equality.
 * Incept:    SRE, Tue May 15 11:00:38 2007 [Janelia]
 *
 * Purpose:   Compares the scores of canonical residues in 
 *            two score matrices <S1> and <S2> for equality.
 *            Returns <eslOK> if they are identical, <eslFAIL> 
 *            if they differ. Peripheral aspects of the scoring matrices
 *            having to do with noncanonical residues, output
 *            order, and suchlike are ignored.
 */
int
esl_scorematrix_CompareCanon(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2)
{
  int a,b;

  for (a = 0; a < S1->K; a++)
    for (b = 0; b < S1->K; b++)
      if (S1->s[a][b] != S2->s[a][b]) return eslFAIL;
  return eslOK;
}



/* Function:  esl_scorematrix_Max()
 * Synopsis:  Returns maximum value in score matrix.
 * Incept:    SRE, Thu Apr 12 18:04:35 2007 [Janelia]
 *
 * Purpose:   Returns the maximum value in score matrix <S>.
 */
int
esl_scorematrix_Max(const ESL_SCOREMATRIX *S)
{
  int i,j;
  int max = S->s[0][0];

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      if (S->s[i][j] > max) max = S->s[i][j];
  return max;
}

/* Function:  esl_scorematrix_Min()
 * Synopsis:  Returns minimum value in score matrix.
 * Incept:    SRE, Thu Apr 12 18:06:50 2007 [Janelia]
 *
 * Purpose:   Returns the minimum value in score matrix <S>.
 */
int
esl_scorematrix_Min(const ESL_SCOREMATRIX *S)
{
  int i,j;
  int min = S->s[0][0];

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      if (S->s[i][j] < min) min = S->s[i][j];
  return min;
}


/* Function:  esl_scorematrix_IsSymmetric()
 * Synopsis:  Returns <TRUE> for symmetric matrix.
 * Incept:    SRE, Sat May 12 18:17:17 2007 [Janelia]
 *
 * Purpose:   Returns <TRUE> if matrix <S> is symmetric,
 *            or <FALSE> if it's not.
 */
int
esl_scorematrix_IsSymmetric(const ESL_SCOREMATRIX *S)
{
  int i,j;

  for (i = 0; i < S->K; i++)
    for (j = i; j < S->K; j++)
      if (S->s[i][j] != S->s[j][i]) return FALSE;
  return TRUE;
}



/* Function:  esl_scorematrix_Destroy()
 * Synopsis:  Frees a matrix.
 * Incept:    SRE, Mon Apr  2 08:46:44 2007 [Janelia]
 *
 * Purpose:   Frees a score matrix.
 */
void
esl_scorematrix_Destroy(ESL_SCOREMATRIX *S)
{
  if (S == NULL) return;
  if (S->s != NULL) {
    if (S->s[0] != NULL) free(S->s[0]);
    free(S->s);
  }
  if (S->isval    != NULL) free(S->isval);
  if (S->outorder != NULL) free(S->outorder);
  if (S->name     != NULL) free(S->name);
  if (S->path     != NULL) free(S->path);
  free(S);
  return;
}





/*****************************************************************
 * 2. Reading/writing score matrices.
 *****************************************************************/

/* Function:  esl_sco_Read()
 * Synopsis:  Read a standard matrix input file.
 * Incept:    SRE, Mon Apr  2 08:26:40 2007 [Janelia]
 *
 * Purpose:   Given a pointer <efp> to an open file parser for a file
 *            containing a score matrix (such as a PAM or BLOSUM
 *            matrix), parse the file and create a new score matrix
 *            object. The scores are expected to be for the alphabet
 *            <abc>. 
 *            
 *            The score matrix file is in the format that BLAST or
 *            FASTA use. The first line is a header contains N
 *            single-letter codes for the residues. Each of N
 *            subsequent rows optionally contains a residue row label
 *            (in the same order as the columns), followed by N
 *            residue scores.  (Older matrix files do not contain the
 *            leading row label; newer ones do.) The residues may
 *            appear in any order. They must minimally include the
 *            canonical K residues (K=4 for DNA, K=20 for protein),
 *            and may also contain none, some, or all degeneracy
 *            codes. Any other residue code that is not in the Easel
 *            digital alphabet (including, in particular, the '*' code
 *            for a stop codon) is ignored by the parser.
 *
 * Returns:   <eslOK> on success, and <ret_S> points to a newly allocated 
 *            score matrix. 
 *
 *            Returns <eslEFORMAT> on parsing error; in which case, <ret_S> is
 *            returned <NULL>, and <efp->errbuf> contains an informative
 *            error message.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_sco_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S)
{
  int status;
  ESL_SCOREMATRIX *S     = NULL;
  int             *map   = NULL; /* maps col/row index to digital alphabet x */
  char            *tok;
  int              toklen;
  int              c, x;
  int              row,col;

  /* Allocate the matrix
   */
  if ((S = esl_scorematrix_Create(abc)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Make sure we've got the comment character set properly in the fileparser.
   * Score matrices use #.
   */
  esl_fileparser_SetCommentChar(efp, '#');

  /* Look for the first non-blank, non-comment line in the file.  That line
   * gives us the single-letter codes in the order that the file's using.
   */
  if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "file appears to be empty");

  /* Read the characters: count them and store them in order in label[0..nc-1].
   * nc cannot exceed Kp+1 in our expected alphabet (+1, for the stop character *)
   */
  S->nc = 0;
  while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) == eslOK)
    {
      if (S->nc >= abc->Kp) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Header contains more residues than expected for alphabet");
      if (toklen != 1)      ESL_XFAIL(eslEFORMAT, efp->errbuf, "Header can only contain single-char labels; %s is invalid", tok);
      S->outorder[S->nc++] = *tok;
    }
  if (status != eslEOL) ESL_XFAIL(status, efp->errbuf, "Unexpected failure of esl_fileparser_GetTokenOnLine()");
  S->outorder[S->nc] = '\0';	/* NUL terminate */
  
  /* Verify that these labels for the score matrix seem plausible, given our alphabet.
   * This sets S->isval array: which residues we have scores for.
   * It also sets the map[] array, which maps coord in label[] to x in alphabet.
   */
  ESL_ALLOC(map, sizeof(int) * S->nc);
  for (c = 0; c < S->nc; c++)
    {
      if (esl_abc_CIsValid(abc, S->outorder[c])) 
	{  
	  x = esl_abc_DigitizeSymbol(abc, S->outorder[c]);
	  map[c] = x;
	  S->isval[x] = TRUE;
	}
      else
	ESL_XFAIL(eslEFORMAT, efp->errbuf, "Don't know how to deal with residue %c in matrix file", S->outorder[c]);
    }
  for (x = 0; x < abc->K; x++)
    if (! S->isval[x]) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Expected to see a column for residue %c", abc->sym[x]);


  /* Read nc rows, one at a time;
   * on each row, read nc+1 or nc tokens, of which nc are scores (may lead with a label or not)
   */
  for (row = 0; row < S->nc; row++)
    {
      if ((status = esl_fileparser_NextLine(efp)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of lines in file");
      for (col = 0; col < S->nc; col++)
	{
	  if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Unexpectedly ran out of fields on line");
	  if (col == 0 && *tok == S->outorder[row]) { col--; continue; } /* skip leading label */

	  S->s[map[row]][map[col]] = atoi(tok);
	}
      if ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslEOL)  ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many fields on line");
    }
  if ((status = esl_fileparser_NextLine(efp)) != eslEOF) ESL_XFAIL(eslEFORMAT, efp->errbuf, "Too many lines in file");
  

  /* Annotate the score matrix */
  if ((status = esl_strdup  (efp->filename, -1,    &(S->path))) != eslOK) goto ERROR;
  if ((status = esl_FileTail(efp->filename, FALSE, &(S->name))) != eslOK) goto ERROR;

  free(map);
  *ret_S = S;
  return eslOK;

 ERROR:
  esl_scorematrix_Destroy(S);
  if (map != NULL) free(map);
  *ret_S = NULL;
  return status;
}

/* Function:  esl_sco_Write()
 * Synopsis:  Write a BLAST-compatible score matrix file.
 * Incept:    SRE, Tue Apr  3 13:55:10 2007 [Janelia]
 *
 * Purpose:   Writes a score matrix <S> to an open stream <fp>, in 
 *            format compatible with BLAST, FASTA, and other common
 *            sequence alignment software.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_sco_Write(FILE *fp, const ESL_SCOREMATRIX *S)
{
  int a,b;			
  int x,y;
  int nc = S->nc;
  
  /* The header line, with column labels for residues */
  fprintf(fp, "  ");
  for (a = 0; a < nc; a++) fprintf(fp, "  %c ", S->outorder[a]);
  fprintf(fp, "\n");
  
  /* The data */
  for (a = 0; a < nc; a++)
    {
      fprintf(fp, "%c ", S->outorder[a]);
      for (b = 0; b < nc; b++)
	{
	  x = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[a]);
	  y = esl_abc_DigitizeSymbol(S->abc_r, S->outorder[b]);
	  fprintf(fp, "%3d ", S->s[x][y]);
	}
      fprintf(fp, "\n");
    }
  return eslOK;
}

/*****************************************************************
 * 3. Interpreting score matrices probabilistically.
 *****************************************************************/ 

static int set_degenerate_probs(const ESL_ALPHABET *abc, ESL_DMATRIX *P, double *fi, double *fj);


struct lambda_params {
  const double *fi;
  const double *fj;
  const ESL_SCOREMATRIX *S;
};

static int
lambda_fdf(double lambda, void *params, double *ret_fx, double *ret_dfx)
{
  struct lambda_params *p = (struct lambda_params *) params;
  int    i,j;
  double tmp;
  
  *ret_fx  = 0.;
  *ret_dfx = 0.;
  for (i = 0; i < p->S->K; i++)
    for (j = 0; j < p->S->K; j++)
      {
	tmp      = p->fi[i] * p->fj[j] * exp(lambda * (double) p->S->s[i][j]);
	*ret_fx  += tmp;
	*ret_dfx += tmp * (double) p->S->s[i][j];
      }
  *ret_fx -= 1.0;
  return eslOK;
}

/* Function:  esl_sco_ProbifyGivenBG()
 * Synopsis:  Obtain $P_{ij}$ for matrix with known $\lambda$ and background. 
 * Incept:    SRE, Thu Apr 12 17:46:20 2007 [Janelia]
 *
 * Purpose:   Given a score matrix <S> and known query and target
 *            background frequencies <fi> and <fj>, calculate scale
 *            <lambda> and implicit target probabilities \citep{Altschul01}. 
 *            Optionally returns either (or both) in <opt_lambda> and <opt_P>.
 *
 *            The implicit target probabilities are returned in a
 *            newly allocated $K \times K$ <ESL_DMATRIX>, over only
 *            the canonical (typically 4 or 20) residues in the
 *            residue alphabet.
 *            
 * Args:      S          - score matrix
 *            fi         - background frequencies for sequence i
 *            fj         - background frequencies for sequence j
 *            opt_lambda - optRETURN: calculated $\lambda$ parameter
 *            opt_P      - optRETURN: implicit target probabilities $p_{ij}$; a KxK DMATRIX.                  
 *
 * Returns:   <eslOK> on success, <*ret_lambda> contains the
 *            calculated $\lambda$ parameter, and <*ret_P> points to
 *            the target probability matrix (which is allocated here,
 *            and must be free'd by caller with <esl_dmatrix_Destroy(*ret_P)>.
 *            
 * Throws:    <eslEMEM> on allocation error; 
 *            <eslEINVAL> if matrix is invalid and has no solution for $\lambda$;
 *            <eslENOHALT> if the solver fails to find $\lambda$.
 *            In these cases, <*ret_lambda> is 0.0, and <*ret_P> is <NULL>. 
 */
int
esl_sco_ProbifyGivenBG(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, 
		       double *opt_lambda, ESL_DMATRIX **opt_P)
{
  int    status;
  ESL_ROOTFINDER *R = NULL;
  ESL_DMATRIX    *P = NULL;
  struct lambda_params p;
  double lambda_guess;
  double lambda;
  int    i,j;
  double fx, dfx;

  /* First, solve for lambda by rootfinding.
   */
  /* Set up the data passed to the lambda_fdf function. */
  p.fi = fi;
  p.fj = fj;
  p.S  = S;

  /* Bracket the root.
   * It's important that we come at the root from the far side, where
   * f(lambda) is positive; else we may identify the root we don't want
   * at lambda=0.
   */
  lambda_guess = 1. / (double) esl_scorematrix_Max(S);
  for (; lambda_guess < 50.; lambda_guess *= 2.0) {
    lambda_fdf(lambda_guess, &p, &fx, &dfx);
    if (fx > 0) break;
  }
  if (fx <= 0) ESL_EXCEPTION(eslEINVAL, "Failed to bracket root for solving lambda");

  /* Create a solver and find lambda by Newton/Raphson */
  if ((    R   = esl_rootfinder_CreateFDF(lambda_fdf, &p) )         == NULL) { status = eslEMEM; goto ERROR; }
  if (( status = esl_root_NewtonRaphson(R, lambda_guess, &lambda))  != eslOK) goto ERROR;
  
  /* Now, given solution for lambda, calculate P
   */
  if (opt_P != NULL) 
    {
      if ((P = esl_dmatrix_Create(S->Kp, S->Kp)) == NULL) { status = eslEMEM; goto ERROR; }
      for (i = 0; i < S->K; i++)
	for (j = 0; j < S->K; j++)
	  P->mx[i][j] = fi[i] * fj[j] * exp(lambda * (double) S->s[i][j]);
      set_degenerate_probs(S->abc_r, P, NULL, NULL);
    }

  esl_rootfinder_Destroy(R);
  if (opt_lambda != NULL) *opt_lambda = lambda;
  if (opt_P      != NULL) *opt_P      = P;  
  return eslOK;

 ERROR:
  if (R != NULL) esl_rootfinder_Destroy(R);
  if (opt_lambda != NULL) *opt_lambda = 0.;
  if (opt_P      != NULL) *opt_P      = NULL;
  return status;


}



/* This section is an implementation of one of the ideas in
 * Yu and Altschul, PNAS 100:15688, 2003 [YuAltschul03]:
 * Given a valid score matrix, calculate its probabilistic
 * basis (P_ij, f_i, f_j, and lambda), on the assumption that
 * the background probabilities are the marginals of P_ij.
 */
struct yualtschul_params {
  ESL_DMATRIX *S;   /* pointer to the KxK score matrix w/ values cast to doubles */		
  ESL_DMATRIX *M;   /* not a param per se: alloc'ed storage for M matrix provided to the objective function */
  ESL_DMATRIX *Y;   /* likewise, alloc'ed storage for Y (M^-1) matrix provided to obj function */
};

/* yualtschul_func()
 *
 * This is the objective function we try to find a root of. 
 * Its prototype is dictated by the esl_rootfinder API.
 */
static int
yualtschul_func(double lambda, void *params, double *ret_fx)
{
  int status;
  struct yualtschul_params *p = (struct yualtschul_params *) params;
  ESL_DMATRIX  *S = p->S;
  ESL_DMATRIX  *M = p->M;
  ESL_DMATRIX  *Y = p->Y;
  int i,j;

  /* the M matrix has entries M_ij = e^{lambda * s_ij} */
  for (i = 0; i < S->n; i++)
    for (j = 0; j < S->n; j++)
      M->mx[i][j] = exp(lambda * S->mx[i][j]);

  /* the Y matrix is the inverse of M */
  if ((status = esl_dmx_Invert(M, Y)) != eslOK) return status;

  /* We're trying to find the root of \sum_ij Y_ij - 1 = 0 */
  *ret_fx = esl_dmx_Sum(Y) - 1.;
  return eslOK;
}

/* yualtschul_engine()
 *
 * This function backcalculates the probabilistic basis for a score
 * matrix S, when S is a double-precision matrix. Providing this
 * as a separate "engine" and writing esl_sco_Probify()
 * as a wrapper around it allows us to separately test inaccuracy
 * due to numerical performance of our linear algebra, versus 
 * inaccuracy due to integer roundoff in integer scoring matrices.
 * 
 * It is not uncommon for this to fail when S is derived from
 * integer scores. Because the scores may have been provided by the
 * user, and this may be our first chance to detect the "user error"
 * of an invalid matrix, this engine returns <eslENORESULT> as a normal error
 * if it can't reach a valid solution.
 */
static int 
yualtschul_engine(ESL_DMATRIX *S, ESL_DMATRIX *P, double *fi, double *fj, double *ret_lambda)
{
  int status;
  ESL_ROOTFINDER *R = NULL;
  struct yualtschul_params p;
  double lambda;
  double xl, xr;
  double fx;
  int    i,j;

  /* Set up a bisection method to find lambda */
  p.S = S;
  p.M = p.Y = NULL;
  if ((p.M = esl_dmatrix_Create(S->n, S->n))           == NULL) { status = eslEMEM; goto ERROR; }
  if ((p.Y = esl_dmatrix_Create(S->n, S->n))           == NULL) { status = eslEMEM; goto ERROR; }
  if ((R = esl_rootfinder_Create(yualtschul_func, &p)) == NULL) { status = eslEMEM; goto ERROR; }

  /* Need a reasonable initial guess for lambda; if we use extreme
   * lambda guesses, we'll introduce numeric instability in the
   * objective function, and may even blow up the values of e^{\lambda
   * s_ij} in the M matrix. Appears to be safe to start with lambda on
   * the order of 2/max(s_ij).
   */
  xr = 1. / esl_dmx_Max(S);
  
  /* Identify suitable brackets on lambda. */
  for (xl = xr; xl > 1e-10; xl /= 1.6) {
    if ((status = yualtschul_func(xl, &p, &fx))  != eslOK) goto ERROR;
    if (fx > 0.) break;
  }
  if (fx <= 0.) { status = eslENORESULT; goto ERROR; }

  for (; xr < 100.; xr *= 1.6) {
    if ((status = yualtschul_func(xr, &p, &fx))  != eslOK) goto ERROR;
    if (fx < 0.) break;
  }
  if (fx >= 0.) { status = eslENORESULT; goto ERROR; }

  /* Find lambda by bisection */
  if (esl_root_Bisection(R, xl, xr, &lambda) != eslOK)     goto ERROR;

  /* Find fi, fj from Y: fi are column sums, fj are row sums */
  for (i = 0; i < S->n; i++) {
    fi[i] = 0.;
    for (j = 0; j < S->n; j++) fi[i] += p.Y->mx[j][i];
  }
  for (j = 0; j < S->n; j++) {
    fj[j] = 0.;
    for (i = 0; i < S->n; i++) fj[j] += p.Y->mx[j][i];
  }

  /* Find p_ij */
  for (i = 0; i < S->n; i++) 
    for (j = 0; j < S->n; j++)
      P->mx[i][j] = fi[i] * fj[j] * p.M->mx[i][j];

  *ret_lambda = lambda;
  esl_dmatrix_Destroy(p.M);
  esl_dmatrix_Destroy(p.Y);
  esl_rootfinder_Destroy(R);
  return eslOK;

 ERROR:
  if (p.M != NULL) esl_dmatrix_Destroy(p.M);
  if (p.Y != NULL) esl_dmatrix_Destroy(p.Y);
  if (R   != NULL) esl_rootfinder_Destroy(R);
  return status;
}

/* Function:  esl_sco_Probify()
 * Synopsis:  Calculate the probabilistic basis of a score matrix.
 * Incept:    SRE, Wed Apr 11 07:56:44 2007 [Janelia]
 *
 * Purpose:   Reverse engineering of a score matrix: given a "valid"
 *            substitution matrix <S>, obtain implied joint
 *            probabilities $p_{ij}$, query composition $f_i$, target
 *            composition $f_j$, and scale $\lambda$, by assuming that
 *            $f_i$ and $f_j$ are the appropriate marginals of $p_{ij}$.
 *            Optionally return any or all of these solutions in
 *            <*opt_P>, <*opt_fi>, <*opt_fj>, and <*opt_lambda>.
 *            
 *            The calculation is run only on canonical residue scores
 *            $0..K-1$ in S, to calculate joint probabilities for all
 *            canonical residues. Joint and background probabilities 
 *            involving degenerate residues are then calculated by
 *            appropriate marginalizations.
 *
 *            This implements an algorithm described in
 *            \citep{YuAltschul03}.
 *            
 *            This algorithm works fine in principle, but when it is
 *            applied to rounded integer scores with small dynamic
 *            range (the typical situation for score matrices) it may
 *            fail due to integer roundoff error. It works best for
 *            score matrices built using small values of $\lambda$. Yu
 *            and Altschul use $\lambda = 0.00635$ for BLOSUM62, which
 *            amounts to scaling default BLOSUM62 up 50-fold. It
 *            happens that default BLOSUM62 (which was created with
 *            lambda = 0.3466, half-bits) can be successfully reverse
 *            engineered (albeit with some loss of accuracy;
 *            calculated lambda is 0.3240) but other common matrices
 *            may fail. This failure results in a normal returned
 *            error of <eslENORESULT>. 
 *            
 * Args:      S          - score matrix 
 *            opt_P      - optRETURN: Kp X Kp matrix of implied target probs $p_{ij}$
 *            opt_fi     - optRETURN: vector of Kp $f_i$ background probs, 0..Kp-1
 *            opt_fj     - optRETURN: vector of Kp $f_j$ background probs, 0..Kp-1
 *            opt_lambda - optRETURN: calculated $\lambda$ parameter
 *
 * Returns:   <eslOK> on success, and <opt_P>, <opt_fi>, <opt_fj>, and <opt_lambda>
 *            point to the results (for any of these that were passed non-<NULL>).
 *
 *            <opt_P>, <opt_fi>, and <opt_fj>, if requested, are new
 *            allocations, and must be freed by the caller.
 *            
 *            Returns <eslENORESULT> if the algorithm fails to determine a valid solution.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J1/35.
 */
int
esl_sco_Probify(const ESL_SCOREMATRIX *S, ESL_DMATRIX **opt_P, double **opt_fi, double **opt_fj, double *opt_lambda)
{
  int status;
  ESL_DMATRIX  *Sd  = NULL;
  ESL_DMATRIX  *P   = NULL;
  double       *fi  = NULL;
  double       *fj  = NULL;
  double        lambda;
  int i,j;

  if (( Sd = esl_dmatrix_Create(S->K,  S->K))  == NULL) {status = eslEMEM; goto ERROR; }
  if (( P  = esl_dmatrix_Create(S->Kp, S->Kp)) == NULL) {status = eslEMEM; goto ERROR; }
  ESL_ALLOC(fi, sizeof(double) * S->Kp);
  ESL_ALLOC(fj, sizeof(double) * S->Kp);

  /* Construct a double-precision dmatrix from S.
   * I've tried integrating over the rounding uncertainty by
   * averaging over trials with values jittered by +/- 0.5,
   * but it doesn't appear to help much, if at all.
   */
  for (i = 0; i < S->K; i++) 
    for (j = 0; j < S->K; j++)
      Sd->mx[i][j] = (double) S->s[i][j];

  /* Reverse engineer the doubles */
  if ((status = yualtschul_engine(Sd, P, fi, fj, &lambda)) != eslOK) goto ERROR;

  /* Set the degenerate probabilities by appropriate sums */
  set_degenerate_probs(S->abc_r, P, fi, fj);
      
  /* Done. */
  esl_dmatrix_Destroy(Sd);
  if (opt_P      != NULL) *opt_P      = P;       else esl_dmatrix_Destroy(P);
  if (opt_fi     != NULL) *opt_fi     = fi;      else free(fi);
  if (opt_fj     != NULL) *opt_fj     = fj;      else free(fj);
  if (opt_lambda != NULL) *opt_lambda = lambda;
  return eslOK;

 ERROR:
  if (Sd  != NULL) esl_dmatrix_Destroy(Sd);
  if (P   != NULL) esl_dmatrix_Destroy(P);
  if (fi  != NULL) free(fi);
  if (fj  != NULL) free(fj);
  if (opt_P      != NULL) *opt_P      = NULL;
  if (opt_fi     != NULL) *opt_fi     = NULL;
  if (opt_fj     != NULL) *opt_fj     = NULL;
  if (opt_lambda != NULL) *opt_lambda = 0.;
  return status;
}



/* Function:  esl_sco_RelEntropy()
 * Synopsis:  Calculates relative entropy of a matrix.
 * Incept:    SRE, Sat May 12 18:14:02 2007 [Janelia]
 *
 * Purpose:   Calculates the relative entropy of score matrix <S> in
 *            bits, given its background distributions <fi> and <fj> and
 *            its scale <lambda>.
 *
 * Args:      S          - score matrix
 *            fi         - background freqs for sequence i
 *            fj         - background freqs for sequence j
 *            lambda     - scale factor $\lambda$ for <S>
 *            ret_D      - RETURN: relative entropy.
 * 
 * Returns:   <eslOK> on success, and <ret_D> contains the relative
 *            entropy.
 *
 * Throws:    <eslEMEM> on allocation error. 
 *            <eslEINVAL> if the implied $p_{ij}$'s don't sum to one,
 *            probably indicating that <lambda> was not the correct
 *            <lambda> for <S>, <fi>, and <fj>.
 *            In either exception, <ret_D> is returned as 0.0.
 */
int
esl_sco_RelEntropy(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, double lambda, double *ret_D)
{
  int    status;
  double pij;
  double sum = 0.;
  int    i,j;
  double D = 0;

  for (i = 0; i < S->K; i++)
    for (j = 0; j < S->K; j++)
      {
	pij  = fi[i] * fj[j] * exp(lambda * (double) S->s[i][j]);
	sum += pij;
	if (pij > 0.) D += pij * log(pij / (fi[i] * fj[j]));
	
      }
  if (esl_DCompare(sum, 1.0, 1e-3) != eslOK) 
    ESL_XEXCEPTION(eslEINVAL, "pij's don't sum to one: bad lambda?");

  D /= eslCONST_LOG2;
  *ret_D = D;
  return eslOK;

 ERROR:
  *ret_D = 0.;
  return status;
}


/* Input: P->mx[i][j] are joint probabilities p_ij for the canonical alphabet 0..abc->K-1,
 *        but P matrix is allocated for Kp X Kp
 * 
 * Fill in [i][j'=K..Kp-1], [i'=K..Kp-1][j], and [i'=K..Kp-1][j'=K..Kp-1] for degeneracies i',j'
 * Any p_ij involving a gap (K), nonresidue (Kp-2), or missing data (Kp-1) character is set to 0.0 by convention.
 *
 * Don't assume symmetry. 
 * 
 * If <fi> or <fj> background probability vectors are non-<NULL>, set them too.
 * (Corresponding to the assumption of background = marginal probs, rather than
 *  background being given.)
 */
static int
set_degenerate_probs(const ESL_ALPHABET *abc, ESL_DMATRIX *P, double *fi, double *fj)
{
  int i,j,ip,jp;

  for (i = 0; i < abc->K; i++)
    {
      P->mx[i][abc->K] = 0.0;
      for (jp = abc->K+1; jp < abc->Kp; jp++)
	{
	  P->mx[i][jp] = 0.0;
	  for (j = 0; j < abc->K; j++)
	    if (abc->degen[jp][j]) P->mx[i][jp] += P->mx[i][j];
	}
      P->mx[i][abc->Kp-2] = 0.0;
      P->mx[i][abc->Kp-1] = 0.0;
    }

  esl_vec_DSet(P->mx[abc->K],    abc->Kp, 0.0); /* gap row */

  for (ip = abc->K+1; ip < abc->Kp-2; ip++)
    {
      for (j = 0; j < abc->K; j++)      
	{
	  P->mx[ip][j] = 0.0;
	  for (i = 0; i < abc->K; i++)
	    if (abc->degen[ip][i]) P->mx[ip][j] += P->mx[i][j];
	}
      P->mx[ip][abc->K] = 0.0;

      for (jp = abc->K+1; jp < abc->Kp; jp++)      
	{
	  P->mx[ip][jp] = 0.0;
	  for (j = 0; j < abc->K; j++)
	    if (abc->degen[jp][j]) P->mx[ip][jp] += P->mx[ip][j];
	}
      P->mx[ip][abc->Kp-2] = 0.0;      
      P->mx[ip][abc->Kp-1] = 0.0;      
    }

  esl_vec_DSet(P->mx[abc->Kp-2], abc->Kp, 0.0); /* nonresidue data ~ row   */
  esl_vec_DSet(P->mx[abc->Kp-1], abc->Kp, 0.0); /* missing data ~ row   */

  if (fi != NULL) { /* fi[i'] = p(i',X) */
    fi[abc->K] = 0.0;
    for (ip = abc->K+1; ip < abc->Kp-2; ip++) fi[ip] = P->mx[ip][abc->Kp-3];
    fi[abc->Kp-2] = 0.0;
    fi[abc->Kp-1] = 0.0;
  }

  if (fj != NULL) { /* fj[j'] = p(X,j')*/
    fj[abc->K] = 0.0;
    for (jp = abc->K+1; jp < abc->Kp-2; jp++) fj[jp] = P->mx[abc->Kp-3][jp];
    fj[abc->Kp-2] = 0.0;
    fj[abc->Kp-1] = 0.0;
  }

  return eslOK;
}


/*****************************************************************
 * 4. Utilities
 *****************************************************************/ 

/* Reformat a score matrix file, canonical residues only, into
 * Easel internal digital alphabet order, suitable for making 
 * a static data structure.
 */
#ifdef eslSCOREMATRIX_UTILITY1
/* 
    gcc -g -Wall -o utility -I. -L. -DeslSCOREMATRIX_UTILITY1 esl_scorematrix.c -leasel -lm
    ./utility BLOSUM62
*/
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_scorematrix.h"
#include "esl_fileparser.h"

int
main(int argc, char **argv)
{
  char *infile = argv[1];
  ESL_ALPHABET    *abc;
  ESL_FILEPARSER  *efp;
  ESL_SCOREMATRIX *S;
  int x,y;

  abc = esl_alphabet_Create(eslAMINO);

  if (esl_fileparser_Open(infile, NULL, &efp) != eslOK) esl_fatal("Failed to open %s\n", infile);
  if (esl_scorematrix_Read(efp, abc, &S)      != eslOK) esl_fatal("parse failed: %s", efp->errbuf);

  for (x = 0; x < abc->Kp; x++) {
    printf("{ ");
    for (y = 0; y < abc->Kp; y++)
      printf("%3d, ", S->s[x][y]);
    printf(" },\n");
  }
  
  esl_scorematrix_Destroy(S);
  esl_fileparser_Close(efp);
  esl_alphabet_Destroy(abc);
  return eslOK;
}
#endif /*eslSCOREMATRIX_UTILITY1*/




/* Utility 2: joint or conditional probabilities from BLOSUM62 (depending on how compiled)
 */
#ifdef eslSCOREMATRIX_UTILITY2
/* 
    gcc -g -Wall -o utility2 -I. -L. -DeslSCOREMATRIX_UTILITY2 esl_scorematrix.c -leasel -lm
    ./utility2
*/
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_scorematrix.h"

int
main(int argc, char **argv)
{
  ESL_ALPHABET    *abc      = esl_alphabet_Create(eslAMINO);
  ESL_SCOREMATRIX *S        = esl_scorematrix_Create(abc);
  ESL_DMATRIX     *Q        = NULL;
  double          *fa       = NULL;
  double          *fb       = NULL;
  double           slambda;
  int              a,b;

  esl_scorematrix_SetBLOSUM62(S);
  
  esl_sco_Probify(S, &Q, &fa, &fb, &slambda);
#if 0
  for (a = 0; a < abc->K; a++)
    for (b = 0; b < abc->K; b++)
      Q->mx[a][b] /= fa[a];	/* Q->mx[a][b] is now P(b | a) */
#endif
  
  esl_dmatrix_Dump(stdout, Q, abc->sym, abc->sym);
  
  esl_dmatrix_Destroy(Q);
  esl_scorematrix_Destroy(S);
  esl_alphabet_Destroy(abc);
  return eslOK;
}
#endif /*eslSCOREMATRIX_UTILITY2*/






/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/

#ifdef eslSCOREMATRIX_TESTDRIVE
#include <esl_dirichlet.h>

static void
utest_ReadWrite(ESL_ALPHABET *abc, ESL_SCOREMATRIX *S)
{
  char tmpfile[16]     = "esltmpXXXXXX";
  FILE            *fp  = NULL;
  ESL_SCOREMATRIX *S2  = NULL;
  ESL_FILEPARSER  *efp = NULL;
  
  if (esl_tmpfile_named(tmpfile, &fp)          != eslOK) esl_fatal("failed to open tmp file");
  if (esl_sco_Write(fp, S)                     != eslOK) esl_fatal("failed to write test matrix");
  fclose(fp);

  if (esl_fileparser_Open(tmpfile, NULL, &efp) != eslOK) esl_fatal("failed to open tmpfile containing BLOSUM62 matrix");
  if (esl_sco_Read(efp, abc, &S2)              != eslOK) esl_fatal("failed to read tmpfile containing BLOSUM62 matrix");
  if (esl_scorematrix_Compare(S, S2)           != eslOK) esl_fatal("the two test matrices aren't identical");
  
  remove(tmpfile); 
  esl_fileparser_Close(efp);
  esl_scorematrix_Destroy(S2);
  return;
}


static void
utest_ProbifyGivenBG(ESL_SCOREMATRIX *S0, ESL_DMATRIX *P0, double *wagpi, double lambda0)
{
  char *msg = "ProbifyGivenBG() unit test failed";
  ESL_DMATRIX     *P    = NULL;
  double           sum  = 0.0;
  double           lambda;
  int              a,b;

  if (esl_sco_ProbifyGivenBG(S0, wagpi, wagpi, &lambda, &P) != eslOK) esl_fatal(msg);

  if (esl_DCompare(lambda0, lambda, 1e-3)     != eslOK) esl_fatal("lambda is wrong");

  for (a = 0; a < 20; a++) 	/* you can't just call esl_dmx_Sum(P), because P includes */
    for (b = 0; b < 20; b++)    /* marginalized degeneracies */
      sum += P->mx[a][b];

  if (esl_DCompare(sum, 1.0, 1e-9)     != eslOK) esl_fatal("P doesn't sum to 1");

  for (a = 0; a < 20; a++)	/* for the same reason,  you can't dmatrix_Compare P and P0 */
    for (b = 0; b < 20; b++)
      if (esl_DCompare(P0->mx[a][b], P->mx[a][b], 1e-2) != eslOK) esl_fatal("P is wrong");

  esl_dmatrix_Destroy(P);
  return;
}
 

/* The scores->pij reverse engineering engine works with scores in doubles,
 * so we can separate effects of rounding to integers in standard
 * score matrices.
 */
static void 
utest_yualtschul(ESL_DMATRIX *P0, double *wagpi)
{
  char *msg = "reverse engineering engine test failed";
  ESL_DMATRIX     *S   = NULL;	/* original score matrix, in double form, not rounded to ints (calculated from P, fi, fj) */
  ESL_DMATRIX     *P   = NULL;	/* backcalculated P_ij joint probabilities */
  double          *fi  = NULL;	/* backcalculated f_i query composition */
  double          *fj  = NULL;	/* backcalculated f'_j target composition */
  double           lambda0;	/* true lambda */
  double           lambda;	/* backcalculated lambda */
  double           sum = 0.0;
  int              i,j;

  /* Allocations */
  if (( S  = esl_dmatrix_Create(20, 20))     == NULL)  esl_fatal(msg);
  if (( P  = esl_dmatrix_Create(20, 20))     == NULL)  esl_fatal(msg);
  if ((fi  = malloc(sizeof(double) * 20))    == NULL)  esl_fatal(msg);
  if ((fj  = malloc(sizeof(double) * 20))    == NULL)  esl_fatal(msg);

  /* Make a WAG-based score matrix in double-precision, without rounding to integers */
  lambda0 = 0.3;
  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      S->mx[i][j] = log(P0->mx[i][j] / (wagpi[i] * wagpi[j])) / lambda0;

  /* Reverse engineer it in double precision */
  if ( yualtschul_engine(S, P, fi, fj, &lambda) != eslOK) esl_fatal("reverse engineering engine failed");

  /* Validate the solution (expect more accuracy from this than from integer scores) */
  if (esl_DCompare(lambda0, lambda, 1e-4)      != eslOK) esl_fatal("failed to get right lambda");

  for (i = 0; i < 20; i++) 	/* you can't just call esl_dmx_Sum(P), because P includes */
    for (j = 0; j < 20; j++)    /* marginalized degeneracies */
      sum += P->mx[i][j];
  if (esl_DCompare(sum, 1.0, 1e-6) != eslOK) esl_fatal("reconstructed P doesn't sum to 1");

  for (i = 0; i < 20; i++)	/* for the same reason,  you can't dmatrix_Compare P and P0 */
    for (j = 0; j < 20; j++)
      if (esl_DCompare(P0->mx[i][j], P->mx[i][j], 1e-2) != eslOK) esl_fatal("failed to recover correct P_ij");
  for (i = 0; i < 20; i++) 
    {
      if (esl_DCompare(fi[i],    fj[i],  1e-6) != eslOK) esl_fatal("background fi, fj not the same");
      if (esl_DCompare(wagpi[i], fi[i],  1e-3) != eslOK) esl_fatal("failed to reconstruct WAG backgrounds");  
    }

  free(fj);
  free(fi);
  esl_dmatrix_Destroy(S);
  esl_dmatrix_Destroy(P);
  return;
}


/* utest_Probify()
 * This tests Probify on a matrix that was calculated from probabilities in the first
 * place. It verifies that the reconstructed Pij matrix matches the original Pij's
 * that the score matrix was built from.
 */
static void
utest_Probify(ESL_SCOREMATRIX *S0, ESL_DMATRIX *P0, double *wagpi, double lambda0)
{
  ESL_DMATRIX     *P  = NULL;
  double          *fi = NULL;
  double          *fj = NULL;
  double           lambda;	/* reconstructed lambda */
  double           sum = 0.0;
  int              i,j;

  if (esl_sco_Probify(S0, &P, &fi, &fj, &lambda) != eslOK) esl_fatal("reverse engineering failed");

  /* Validate the solution, gingerly (we expect significant error due to integer roundoff) */
  if (esl_DCompare(lambda0, lambda, 0.01)       != eslOK) esl_fatal("failed to get right lambda");
  for (i = 0; i < 20; i++) 	/* you can't just call esl_dmx_Sum(P), because P includes */
    for (j = 0; j < 20; j++)    /* marginalized degeneracies */
      sum += P->mx[i][j];
  if (esl_DCompare(sum, 1.0, 1e-6) != eslOK) esl_fatal("reconstructed P doesn't sum to 1");

  for (i = 0; i < 20; i++)	/* for the same reason,  you can't dmatrix_Compare P and P0 */
    for (j = 0; j < 20; j++)
      if (esl_DCompare(P0->mx[i][j], P->mx[i][j], 0.1) != eslOK) esl_fatal("failed to recover correct P_ij");
  free(fj);
  free(fi);
  esl_dmatrix_Destroy(P);
  return;
}

/* utest_ProbifyBLOSUM()
 * This tests Probify on a score matrix where the original Pij's are treated as
 * unknown. It verifies that if you create a new score matrix from the reconstructed
 * Pij's, you get the original score matrix back. BLOSUM62 makes a good example,
 * hence the name.
  */
static void
utest_ProbifyBLOSUM(ESL_SCOREMATRIX *BL62)
{
  char *msg = "failure in ProbifyBLOSUM() unit test";
  ESL_DMATRIX     *P  = NULL;
  double          *fi = NULL;
  double          *fj = NULL;
  double           lambda;	
  ESL_SCOREMATRIX *S2 = NULL;

  if (( S2 = esl_scorematrix_Clone(BL62))                  == NULL) esl_fatal(msg);
  if (esl_sco_Probify(BL62, &P, &fi, &fj, &lambda)        != eslOK) esl_fatal(msg);
  if (esl_scorematrix_SetFromProbs(S2, lambda, P, fi, fj) != eslOK) esl_fatal(msg);
  if (esl_scorematrix_CompareCanon(BL62, S2)              != eslOK) esl_fatal(msg);
  
  free(fj);
  free(fi);
  esl_scorematrix_Destroy(S2);
  esl_dmatrix_Destroy(P);
  return;
}

#endif /*eslSCOREMATRIX_TESTDRIVE*/


/*****************************************************************
 * 6. Test driver.
 *****************************************************************/
/* 
    gcc -g -Wall -I. -L. -o test -DeslSCOREMATRIX_TESTDRIVE esl_scorematrix.c -leasel -lm
    ./test
*/
#ifdef eslSCOREMATRIX_TESTDRIVE
#include "easel.h"
#include "esl_scorematrix.h"

int 
main(int argc, char **argv)
{
  ESL_ALPHABET    *abc = NULL;	/* amino acid alphabet */
  ESL_SCOREMATRIX *BL62= NULL;	/* BLOSUM62 matrix */
  ESL_SCOREMATRIX *S0  = NULL;	/* original score matrix (calculated from P, fi, fj) */
  ESL_DMATRIX     *P0  = NULL;	/* original P_ij joint probabilities */
  ESL_DMATRIX     *Q   = NULL;	/* WAG rate matrix */
  double           lambda0;	/* true lambda used to construct S */
  double           t;
  int              i,j;
  static double    wagpi[20];

  /* Allocations */
  if ((abc = esl_alphabet_Create(eslAMINO))      == NULL)  esl_fatal("allocation of alphabet failed");
  if ((BL62= esl_scorematrix_Create(abc))        == NULL)  esl_fatal("allocation of BLOSUM62 failed");
  if ((S0  = esl_scorematrix_Create(abc))        == NULL)  esl_fatal("allocation of scorematrix failed");
  if ((P0  = esl_dmatrix_Create(abc->K, abc->K)) == NULL)  esl_fatal("P allocation failed");
  if ((Q   = esl_dmatrix_Create(abc->K, abc->K)) == NULL)  esl_fatal("Q allocation failed");

  /* Make a BLOSUM matrix */
  if ( esl_scorematrix_SetBLOSUM62(BL62) != eslOK) esl_fatal("failed to set a BLOSUM matrix");

  /* Make a WAG-based score matrix with small lambda. */
  lambda0 = 0.00635;
  t    = 2.0;
  esl_scorematrix_SetWAG(S0, lambda0, t);
  esl_composition_WAG(wagpi);

  /* Redo some calculations to get the known probabilistic basis of that S */
  if ( esl_rmx_SetWAG(Q, wagpi)  != eslOK) esl_fatal("failed to set WAG");
  if ( esl_dmx_Exp(Q, t, P0)     != eslOK) esl_fatal("failed to exponentiate WAG");
  for (i = 0; i < 20; i++) 
    for (j = 0; j < 20; j++)
      P0->mx[i][j] *= wagpi[i];	/* P_ij = P(j|i) pi_i */

  /* The unit test battery
   */
  utest_ReadWrite(abc, BL62);
  utest_ReadWrite(abc, S0);
  utest_ProbifyGivenBG(S0, P0, wagpi, lambda0);
  utest_yualtschul(P0, wagpi);
  utest_Probify(S0, P0, wagpi, lambda0); 
  utest_ProbifyBLOSUM(BL62);

  esl_dmatrix_Destroy(Q);
  esl_dmatrix_Destroy(P0);
  esl_scorematrix_Destroy(BL62);
  esl_scorematrix_Destroy(S0);
  esl_alphabet_Destroy(abc);

  return 0;
}
#endif /*eslSCOREMATRIX_TESTDRIVE*/

/*****************************************************************
 * 7. Example program
 *****************************************************************/
#ifdef eslSCOREMATRIX_EXAMPLE
/*::cexcerpt::scorematrix_example::begin::*/
/*  gcc -g -Wall -I. -L. -o example -DeslSCOREMATRIX_EXAMPLE esl_scorematrix.c -leasel -lm
    ./example <score matrix file>
*/
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_fileparser.h"
#include "esl_dmatrix.h"
#include "esl_vectorops.h"
#include "esl_scorematrix.h"

int main(int argc, char **argv)
{
  char            *scorefile = argv[1];
  ESL_ALPHABET    *abc       = esl_alphabet_Create(eslAMINO);
  ESL_FILEPARSER  *efp       = NULL;
  ESL_SCOREMATRIX *S         = NULL;
  ESL_DMATRIX     *P         = NULL;
  double          *fi        = NULL;
  double          *fj        = NULL;
  double           lambda, D;
  
  /* Input an amino acid score matrix from a file. */
  if ( esl_fileparser_Open(scorefile, NULL, &efp) != eslOK) esl_fatal("failed to open score file %s", scorefile);
  if ( esl_sco_Read(efp, abc, &S)                 != eslOK) esl_fatal("failed to read matrix from %s", scorefile);
  esl_fileparser_Close(efp);

  /* Reverse engineer it to get implicit probabilistic model. */
  if ( esl_sco_Probify(S, &P, &fi, &fj, &lambda) != eslOK) esl_fatal("reverse engineering failed");

  /* Print some info, and the joint probabilities. */
  if (esl_scorematrix_IsSymmetric(S)) printf("Matrix is a standard symmetric matrix\n");
  else                                printf("Matrix is a nonstandard asymmetric matrix\n"); 
  printf("Lambda is %g\n\n", lambda);

  esl_sco_RelEntropy(S, fi, fj, lambda, &D);
  printf("Relative entropy = %.4f bits\n\n", D); 

  printf("Implicit joint probabilities are:\n");
  esl_dmatrix_Dump(stdout, P, abc->sym, abc->sym);
  printf("fi's are:\n");
  esl_vec_DDump(stdout, fi, S->K, abc->sym);
  printf("fj's are:\n");
  esl_vec_DDump(stdout, fj, S->K, abc->sym);

  esl_composition_BL62(fi);
  esl_sco_ProbifyGivenBG(S, fi, fi, &lambda, &P);
  esl_sco_RelEntropy(S, fi, fi, lambda, &D);
  printf("Using blosum62 background, lambda = %.4f, rel entropy = %.4f\n", lambda, D);

  free(fi);
  free(fj);
  esl_dmatrix_Destroy(P);
  esl_scorematrix_Destroy(S);
  esl_alphabet_Destroy(abc);
  return 0;
}
/*::cexcerpt::scorematrix_example::end::*/
#endif /*eslSCOREMATRIX_EXAMPLE*/


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/ 
