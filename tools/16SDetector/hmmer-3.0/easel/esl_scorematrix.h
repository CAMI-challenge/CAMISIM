/* Routines for manipulating sequence alignment score matrices.
 * 
 * SRE, Mon Apr  2 08:33:23 2007 [Janelia]
 * SVN $Id: esl_scorematrix.h 337 2009-05-12 02:13:02Z eddys $
 */
#ifndef ESL_SCOREMATRIX_INCLUDED
#define ESL_SCOREMATRIX_INCLUDED

#include <esl_alphabet.h>
#include <esl_fileparser.h>
#include <esl_dmatrix.h>
#include <esl_random.h>

/* ESL_SCOREMATRIX:
 * allocation is in one array in s[0].
 *
 * i,j can range from 0..Kp-1, including all characters valid in the alphabet.
 * Only values for 0..K-1 (canonical alphabet) are mandatory.
 */
typedef struct {
  int **s;			/* s[i][j] is the score of aligning residue i,j; i,j range 0..Kp-1 */
  int   K;			/* size of base alphabet (duplicate of S->abc_r->K) */
  int   Kp;			/* full size of s[][], including degeneracies (duplicate of S->abc_r->Kp) */

  /* bookkeeping for degenerate residues */
  char *isval;			/* array 0..Kp-1: which residues of alphabet have valid scores in S. */
  const ESL_ALPHABET *abc_r;	/* reference to the alphabet: includes K, Kp, and sym order */

  /* bookkeeping that lets us output exactly the residue order we read in a matrix file */
  int   nc;			/* number of residues with scores (inclusive of *, if present) */
  char *outorder;		/* NUL-terminated string 0..nc-1 giving order of residues in col/row labels   */

  char *name;			/* optional: name of score matrix; or NULL */
  char *path;			/* optional: full path to file that score matrix was read from; or NULL  */
} ESL_SCOREMATRIX;




/* 1. The ESL_SCOREMATRIX object. */
extern ESL_SCOREMATRIX *esl_scorematrix_Create(const ESL_ALPHABET *abc);
extern int              esl_scorematrix_SetIdentity(ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_SetBLOSUM62(ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_SetWAG(ESL_SCOREMATRIX *S, double lambda, double t);
extern int              esl_scorematrix_SetFromProbs(ESL_SCOREMATRIX *S, double lambda, const ESL_DMATRIX *P,
						     const double *fi, const double *fj);
extern int              esl_scorematrix_Copy(const ESL_SCOREMATRIX *src, ESL_SCOREMATRIX *dest);
extern ESL_SCOREMATRIX *esl_scorematrix_Clone(const ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_Compare(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2);
extern int              esl_scorematrix_CompareCanon(const ESL_SCOREMATRIX *S1, const ESL_SCOREMATRIX *S2);
extern int              esl_scorematrix_Max(const ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_Min(const ESL_SCOREMATRIX *S);
extern int              esl_scorematrix_IsSymmetric(const ESL_SCOREMATRIX *S);
extern void             esl_scorematrix_Destroy(ESL_SCOREMATRIX *S);

/* 2. Reading/writing score matrices. */
extern int  esl_sco_Read(ESL_FILEPARSER *efp, const ESL_ALPHABET *abc, ESL_SCOREMATRIX **ret_S);
extern int  esl_sco_Write(FILE *fp, const ESL_SCOREMATRIX *S);

/* 3. Interpreting score matrices probabilistically. */
extern int esl_sco_ProbifyGivenBG(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, 
				  double *opt_lambda, ESL_DMATRIX **opt_P);
extern int esl_sco_Probify(const ESL_SCOREMATRIX *S, ESL_DMATRIX **opt_P, 
			   double **opt_fi, double **opt_fj, double *opt_lambda);
extern int esl_sco_RelEntropy(const ESL_SCOREMATRIX *S, const double *fi, const double *fj, 
			      double lambda, double *ret_D);





#endif /*ESL_SCOREMATRIX_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/ 



