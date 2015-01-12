/* Shuffling or bootstrapping alignments.
 * 
 * Table of contents:
 *    1. Randomizing MSAs.
 *    2. Shuffling pairwise (QRNA) alignments.
 *    
 * SRE, Tue Jan 22 09:06:27 2008 [Market Street Cafe, Leesburg]
 * SVN $Id: esl_msashuffle.c 249 2008-04-24 19:19:50Z eddys $
 */
#include "esl_config.h"

#include <string.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif /*eslAUGMENT_ALPHABET*/
#include "esl_msa.h"
#include "esl_msashuffle.h"
#include "esl_random.h"


/*****************************************************************
 * 1. Randomizing MSAs
 *****************************************************************/ 

/* Function:  esl_msashuffle_Shuffle()
 * Synopsis:  Shuffle an alignment's columns.
 * Incept:    SRE, Tue Jan 22 10:10:57 2008 [Janelia]
 *
 * Purpose:   Returns a column-shuffled version of <msa> in <shuf>,
 *            using random generator <r>. Shuffling by columns
 *            preserves the \% identity of the original
 *            alignment. <msa> and <shuf> can be identical, to shuffle
 *            in place.
 *            
 *            The caller sets up the rest of the data (everything but
 *            the alignment itself) in <shuf> the way it wants,
 *            including sequence names, MSA name, and other
 *            annotation. The easy thing to do is to make <shuf>
 *            a copy of <msa>: the caller might create <shuf> by
 *            a call to <esl_msa_Clone()>.
 *            
 *            The alignments <msa> and <shuf> can both be in digital
 *            mode, or can both be in text mode; you cannot mix
 *            digital and text modes.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <msa>,<shuf> aren't in the same mode (digital vs. text).
 */
int
esl_msashuffle_Shuffle(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *shuf)
{
  int i, pos, alen;

  if (! (msa->flags & eslMSA_DIGITAL))
    {
      char c;
      if (shuf->flags & eslMSA_DIGITAL) ESL_EXCEPTION(eslEINVAL, "<shuf> must be in text mode if <msa> is");
      if (msa != shuf) {
	for (i = 0; i < msa->nseq; i++)
	  strcpy(shuf->aseq[i], msa->aseq[i]);
      }

      for (i = 0; i < msa->nseq; i++)
	shuf->aseq[i][msa->alen] = '\0';

      for (alen = msa->alen; alen > 1; alen--)
	{
	  pos = esl_rnd_Roll(r, alen);
	  for (i = 0; i < msa->nseq; i++)
	    {
	      c                     = msa->aseq[i][pos];
	      shuf->aseq[i][pos]    = shuf->aseq[i][alen-1];
	      shuf->aseq[i][alen-1] = c;
	    }
	}
    }
#ifdef eslAUGMENT_ALPHABET
  else 
    {
      ESL_DSQ x;
      if (! (shuf->flags & eslMSA_DIGITAL)) ESL_EXCEPTION(eslEINVAL, "<shuf> must be in digital mode if <msa> is");

      if (msa != shuf) {
	for (i = 0; i < msa->nseq; i++)
	  memcpy(shuf->ax[i], msa->ax[i], (msa->alen + 2) * sizeof(ESL_DSQ));
      }

      for (i = 0; i < msa->nseq; i++)
	shuf->ax[i][msa->alen+1] = eslDSQ_SENTINEL;

      for (alen = msa->alen; alen > 1; alen--)
	{
	  pos = esl_rnd_Roll(r, alen) + 1;
	  for (i = 0; i < msa->nseq; i++)
	    {
	      x                 = msa->ax[i][pos];
	      shuf->ax[i][pos]  = shuf->ax[i][alen];
	      shuf->ax[i][alen] = x;
	    }
	}
    }
#endif /*eslAUGMENT_ALPHABET*/

  return eslOK;
}

/* Function:  esl_msashuffle_Bootstrap()
 * Synopsis:  Bootstrap sample an MSA.
 * Incept:    SRE, Tue Jan 22 11:05:07 2008 [Janelia]
 *
 * Purpose:   Takes a bootstrap sample of <msa> (sample <alen> columns,
 *            with replacement) and puts it in <bootsample>, using
 *            random generator <r>. 
 *            
 *            The caller provides allocated space for <bootsample>.
 *            It must be different space than <msa>; you cannot take
 *            a bootstrap sample "in place". The caller sets up the
 *            rest of the data in <bootsample> (everything but the
 *            alignment itself) the way it wants, including sequence
 *            names, MSA name, and other annotation. The easy thing to
 *            do is to initialize <bootsample> by cloning <msa>.
 *
 *            The alignments <msa> and <bootsample> can both be in digital
 *            mode, or can both be in text mode; you cannot mix
 *            digital and text modes.
 *
 * Returns:   <eslOK> on success, and the alignment in <bootsample> is
 *            set to be a bootstrap resample of the alignment in <msa>.
 *
 * Throws:    <eslEINVAL> if <msa>,<bootsample> aren't in the same mode
 *            (digital vs. text).
 */
int 
esl_msashuffle_Bootstrap(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *bootsample)
{
  int i, pos, col;

  /* contract checks */
  if (  (msa->flags & eslMSA_DIGITAL) && ! (bootsample->flags & eslMSA_DIGITAL))
    ESL_EXCEPTION(eslEINVAL, "<msa> and <bootsample> must both be in digital or text mode");
  if (! (msa->flags & eslMSA_DIGITAL) &&   (bootsample->flags & eslMSA_DIGITAL))
    ESL_EXCEPTION(eslEINVAL, "<msa> and <bootsample> must both be in digital or text mode");

  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (pos = 0; pos < msa->alen; pos++)
	{
	  col = esl_rnd_Roll(r, msa->alen);
	  for (i = 0; i < msa->nseq; i++)
	    bootsample->aseq[i][pos] = msa->aseq[i][col];
	}

      for (i = 0; i < msa->nseq; i++)
	bootsample->aseq[i][msa->alen] = '\0';
    }
#ifdef eslAUGMENT_ALPHABET
  else
    {
      for (i = 0; i < msa->nseq; i++)
	bootsample->ax[i][0] = eslDSQ_SENTINEL;

      for (pos = 1; pos <= msa->alen; pos++)
	{
	  col = esl_rnd_Roll(r, msa->alen) + 1;
	  for (i = 0; i < msa->nseq; i++)
	    bootsample->ax[i][pos] = msa->ax[i][col];
	}

      for (i = 0; i < msa->nseq; i++)
	bootsample->ax[i][msa->alen+1] = eslDSQ_SENTINEL;
    }
#endif /*eslAUGMENT_ALPHABET*/

  return eslOK;
}

/*****************************************************************
 * 2. Shuffling pairwise (QRNA) alignments
 *****************************************************************/ 
#ifdef eslAUGMENT_ALPHABET
/* Function: esl_msashuffle_XQRNA()
 * Synopsis: Gap-preserving column shuffle of a digital pairwise alignment.
 * Incept:   SRE, Tue Jan 22 09:09:52 2008 [Market Street Cafe, Leesburg]
 *
 * Purpose:  Shuffle a digital pairwise alignment <x>,<y> while
 *           preserving the position of gaps, where both sequences are
 *           in digital alphabet <abc>, using the random number
 *           generator <r>. Return the shuffled alignment in <xs>,
 *           <ys>. Caller provides allocated space for <xs> and <ys>
 *           for at least the same length of <x>,<y>.
 *           
 *           Works by doing three separate
 *           shuffles, of (1) columns with residues in both
 *           <x> and <y>, (2) columns with residue in <x> and gap in <y>,
 *           and (3) columns with gap in <x> and residue in <y>.
 *           
 *           <xs>,<x> and <ys>,<y> may be identical: that is, to shuffle
 *           an alignment "in place", destroying the original
 *           alignment, just call <esl_msashuffle_XQRNA(r, abc, x,y,x,y)>.
 *
 * Returns:  <eslOK> on success, and the shuffled alignment is 
 *           returned in <xs>, <ys>.
 *           
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
esl_msashuffle_XQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, ESL_DSQ *x, ESL_DSQ *y, ESL_DSQ *xs, ESL_DSQ *ys)
{
  int  L;
  int *xycol = NULL;
  int *xcol  = NULL;
  int *ycol  = NULL;
  int  nxy, nx, ny;
  int  i;
  int  pos, c;
  char xsym, ysym;
  int  status;

  L = esl_abc_dsqlen(x);
  if (esl_abc_dsqlen(y) != L) ESL_XEXCEPTION(eslEINVAL, "sequences of different lengths in qrna shuffle");

  if (xs != x) esl_abc_dsqcpy(x, L, xs);
  if (ys != y) esl_abc_dsqcpy(y, L, ys);

  /* First, construct three arrays containing lists of the column positions
   * of the three types of columns. (If a column contains gaps in both x and y,
   * we've already simply copied it to the shuffled sequence.)
   */
  ESL_ALLOC(xycol, sizeof(int) * L);
  ESL_ALLOC(xcol,  sizeof(int) * L);
  ESL_ALLOC(ycol,  sizeof(int) * L);
  nxy = nx = ny = 0;

  for (i = 1; i <= L; i++)
    {
      if      (  esl_abc_XIsGap(abc, x[i]) &&   esl_abc_XIsGap(abc, y[i])) { continue; }
      else if (! esl_abc_XIsGap(abc, x[i]) && ! esl_abc_XIsGap(abc, y[i])) { xycol[nxy] = i; nxy++; }
      else if (  esl_abc_XIsGap(abc, x[i]))                                { ycol[ny] = i;   ny++;  }
      else if (  esl_abc_XIsGap(abc, y[i]))                                { xcol[nx] = i;   nx++;  }
    }

  /* Second, shuffle the sequences indirectly, via shuffling these arrays.
   * Yow, careful with those indices, and with order of the statements...
   */
  for (; nxy > 1; nxy--) {
    pos              = esl_rnd_Roll(r, nxy);
    xsym             = xs[xycol[pos]];   ysym             = ys[xycol[pos]];    c            = xycol[pos];   
    xs[xycol[pos]]   = xs[xycol[nxy-1]]; ys[xycol[pos]]   = ys[xycol[nxy-1]];  xycol[pos]   = xycol[nxy-1];
    xs[xycol[nxy-1]] = xsym;             ys[xycol[nxy-1]] = ysym;              xycol[pos]   = xycol[nxy-1];
  }
  for (; nx > 1; nx--) {
    pos            = esl_rnd_Roll(r, nx); 
    xsym           = xs[xcol[pos]];  ysym           = ys[xcol[pos]];  c          = xcol[pos];  
    xs[xcol[pos]]  = xs[xcol[nx-1]]; ys[xcol[pos]]  = ys[xcol[nx-1]]; xcol[pos]  = xcol[nx-1]; 
    xs[xcol[nx-1]] = xsym;           ys[xcol[nx-1]] = ysym;           xcol[nx-1] = c;          
  }
  for (; ny > 1; ny--) {
    pos            = esl_rnd_Roll(r, ny); 
    xsym           = xs[ycol[pos]];  ysym           = ys[ycol[pos]];  c          = ycol[pos]; 
    xs[ycol[pos]]  = xs[ycol[ny-1]]; ys[ycol[pos]]  = ys[ycol[ny-1]]; ycol[pos]  = ycol[ny-1];
    xs[ycol[ny-1]] = xsym;           ys[ycol[ny-1]] = ysym;           ycol[ny-1] = c;          
  }

  free(xycol); free(xcol); free(ycol);
  return eslOK;

 ERROR:
  if (xycol != NULL) free(xycol);
  if (xcol  != NULL) free(xcol);
  if (ycol  != NULL) free(ycol);
  return status;
}

/* Function: esl_msashuffle_CQRNA()
 * Synopsis: Gap-preserving column shuffle of a pairwise alignment.
 * Incept:   SRE, Tue Jan 22 08:45:34 2008 [Market Street Cafe, Leesburg]
 *
 * Purpose:  Shuffle a pairwise alignment <x>,<y> while preserving the
 *           position of gaps, using the random number generator <r>.
 *           Return the shuffled alignment in <xs>,
 *           <ys>. Caller provides allocated space for <xs> and <ys>.
 *           
 *           An alphabet <abc> must also be provided, solely for the
 *           definition of gap characters. Because Easel's default
 *           alphabets (DNA, RNA, and protein) all use the same
 *           definition of gap characters <-_.>, you can actually
 *           provide any alphabet here, and get the same results.
 *           (This may save having to determine the alphabet of input
 *           sequences.)
 *           
 *           Works by doing three separate
 *           shuffles, of (1) columns with residues in both
 *           <x> and <y>, (2) columns with residue in <x> and gap in <y>,
 *           and (3) columns with gap in <x> and residue in <y>.
 *           
 *           <xs>,<x> and <ys>,<y> may be identical: that is, to shuffle
 *           an alignment "in place", destroying the original
 *           alignment, just call <esl_msashuffle_CQRNA(r, abc, x,y,x,y)>.
 *
 * Returns:  <eslOK> on success, and the shuffled alignment is 
 *           returned in <xs>, <ys>.
 *           
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
esl_msashuffle_CQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, char *x, char *y, char *xs, char *ys)
{
  int  L;
  int *xycol = NULL;
  int *xcol  = NULL;
  int *ycol  = NULL;
  int  nxy, nx, ny;
  int  i;
  int  pos, c;
  char xsym, ysym;
  int  status;

  if (xs != x) strcpy(xs, x);
  if (ys != y) strcpy(ys, y);

  /* First, construct three arrays containing lists of the column positions
   * of the three types of columns. (If a column contains gaps in both x and y,
   * we've already simply copied it to the shuffled sequence.)
   */
  L = strlen(x);
  if (strlen(y) != L) ESL_XEXCEPTION(eslEINVAL, "sequences of different lengths in qrna shuffle");
  ESL_ALLOC(xycol, sizeof(int) * L);
  ESL_ALLOC(xcol,  sizeof(int) * L);
  ESL_ALLOC(ycol,  sizeof(int) * L);
  nxy = nx = ny = 0;

  for (i = 0; i < L; i++)
    {
      if      (  esl_abc_CIsGap(abc, x[i]) &&   esl_abc_CIsGap(abc, y[i])) { continue; }
      else if (! esl_abc_CIsGap(abc, x[i]) && ! esl_abc_CIsGap(abc, y[i])) { xycol[nxy] = i; nxy++; }
      else if (  esl_abc_CIsGap(abc, x[i]))                                { ycol[ny] = i;   ny++;  }
      else if (  esl_abc_CIsGap(abc, y[i]))                                { xcol[nx] = i;   nx++;  }
    }

  /* Second, shuffle the sequences indirectly, via shuffling these arrays.
   * Yow, careful with those indices, and with order of the statements...
   */
  for (; nxy > 1; nxy--) {
    pos              = esl_rnd_Roll(r, nxy);
    xsym             = xs[xycol[pos]];   ysym             = ys[xycol[pos]];    c            = xycol[pos];   
    xs[xycol[pos]]   = xs[xycol[nxy-1]]; ys[xycol[pos]]   = ys[xycol[nxy-1]];  xycol[pos]   = xycol[nxy-1];
    xs[xycol[nxy-1]] = xsym;             ys[xycol[nxy-1]] = ysym;              xycol[pos]   = xycol[nxy-1];
  }
  for (; nx > 1; nx--) {
    pos            = esl_rnd_Roll(r, nx); 
    xsym           = xs[xcol[pos]];  ysym           = ys[xcol[pos]];  c          = xcol[pos];  
    xs[xcol[pos]]  = xs[xcol[nx-1]]; ys[xcol[pos]]  = ys[xcol[nx-1]]; xcol[pos]  = xcol[nx-1]; 
    xs[xcol[nx-1]] = xsym;           ys[xcol[nx-1]] = ysym;           xcol[nx-1] = c;          
  }
  for (; ny > 1; ny--) {
    pos            = esl_rnd_Roll(r, ny); 
    xsym           = xs[ycol[pos]];  ysym           = ys[ycol[pos]];  c          = ycol[pos]; 
    xs[ycol[pos]]  = xs[ycol[ny-1]]; ys[ycol[pos]]  = ys[ycol[ny-1]]; ycol[pos]  = ycol[ny-1];
    xs[ycol[ny-1]] = xsym;           ys[ycol[ny-1]] = ysym;           ycol[ny-1] = c;          
  }

  free(xycol); free(xcol); free(ycol);
  return eslOK;

 ERROR:
  if (xycol != NULL) free(xycol);
  if (xcol  != NULL) free(xcol);
  if (ycol  != NULL) free(ycol);
  return status;
}
#endif /*eslAUGMENT_ALPHABET*/
