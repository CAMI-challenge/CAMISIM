/*::cexcerpt::header_example::begin::*/
/* Multiple sequence alignment file i/o.
 *    
 * Contents:   
 *    1. The <ESL_MSA> object
 *    2. The <ESL_MSAFILE> object
 *    3. Digital mode MSA's         (augmentation: alphabet)
 *    4. Random MSA database access (augmentation: ssi)
 *    5. General i/o API, for all alignment formats
 *    6. Miscellaneous functions for manipulating MSAs
 *    7. Stockholm (Pfam/Rfam) format
 *    8. A2M format
 *    9. PSIBLAST format
 *   10. SELEX format
 *   11. AFA (aligned FASTA) format
 *   12. Memory efficient routines for PFAM format
 *   13. Debugging/development routines
 *   14. Benchmark driver
 *   15. Unit tests
 *   16. Test driver
 *   17. Examples
 *   18. Copyright and license information
 *   
 * Augmentations:
 *   alphabet:  adds support for digital MSAs
 *   keyhash:   speeds up Stockholm file input
 *   ssi:       enables indexed random access in a file of many MSAs
 *
 * to do: SRE, Sat Jan  3 09:43:42 2009 (after selex parser added)
 * - SELEX parser is better in some respects than older Stockholm
 *    parser; stricter, better error detection, better modularity.  
 *    Generalize the SELEX parser routines and use them for Stockholm.
 * - Test files for SELEX parser are in esl_msa_testfiles/selex, with
 *    tabular summary list in 00MANIFEST. This is an experiment with
 *    writing tests that require lots of external files, such as
 *    format parsers. Write test driver routine that reads 00MANIFEST
 *    and runs esl_msa_Read() against these files, checking for proper
 *    return status, including errors.
 * - The selex parser's read_block() reads lines into memory and
 *    parses them later. afp->linenumber is thus no longer an
 *    accurate record of where a parse error occurs. read_xxx()
 *    format parsers now need to include line number in their 
 *    afp->errbuf[] message upon eslEFORMAT error. Stockholm parser
 *    doesn't do this. Make it so, and document in examples.
 * - Format autodetection doesn't work yet. Coordinate w/ how sqio
 *    does it, and implement. May require buffering input to make
 *    it work with .gz, pipes without rewinding a stream. Might be
 *    a good idea to generalize input buffering - perhaps making
 *    it part of ESL_FILEPARSER. 
 * - PSIBLAST, A2M format only supported on output, not input.
 *    Implement input parsers.
 * - SELEX format only supported on input, not output. 
 *    Implement output writer.
 * - More formats need to be parsed. Check on formats for current
 *    best MSA programs, such as MUSCLE, MAFFT; implement i/o.
 *    
 * SRE, Thu Jan 20 08:50:43 2005 [St. Louis]
 * SVN $Id: esl_msa.c 573 2010-03-27 15:13:52Z eddys $
 */
/*::cexcerpt::header_example::end::*/

/*::cexcerpt::include_example::begin::*/
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "easel.h"
#ifdef eslAUGMENT_KEYHASH
#include "esl_keyhash.h"
#endif
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"
#endif
#include "esl_msa.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"
/*::cexcerpt::include_example::end::*/



/******************************************************************************
 *# 1. The <ESL_MSA> object                                           
 *****************************************************************************/

/* create_mostly()
 * SRE, Sun Aug 27 16:40:00 2006 [Leesburg]
 *
 * This is the routine called by esl_msa_Create() and esl_msa_CreateDigital()
 * that does all allocation except the aseq/ax alignment data.
 * 
 * <nseq> may be the exact known # of seqs in an alignment; or <nseq>
 * may be an allocation block size (to be expanded by doubling, in
 * esl_msa_Expand(), as in:
 *     <if (msa->nseq == msa->sqalloc) esl_msa_Expand(msa);>
 * <nseq> should not be 0.
 *
 * <alen> may be the exact length of an alignment, in columns; or it
 * may be -1, which states that your parser will take responsibility
 * for expanding as needed as new input is read into a growing new
 * alignment.
 *
 * A created <msa> can only be <_Expand()>'ed if <alen> is -1.
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or -1     
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *            Note that msa->nseq is initialized to 0 here, even though space
 *            is allocated.
 *           
 * Throws:    <NULL> on allocation failure.          
 */
static ESL_MSA *
create_mostly(int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa     = NULL;
  int      i;

  ESL_ALLOC(msa, sizeof(ESL_MSA));
  msa->aseq    = NULL;
  msa->sqname  = NULL;
  msa->wgt     = NULL;
  msa->alen    = alen;		/* if -1, then we're growable. */
  msa->nseq    = 0;		/* our caller (text or digital allocation) sets this.  */
  msa->flags   = 0;

#ifdef eslAUGMENT_ALPHABET
  msa->abc     = NULL;
  msa->ax      = NULL;
#endif /*eslAUGMENT_ALPHABET*/

  msa->name    = NULL;
  msa->desc    = NULL;
  msa->acc     = NULL;
  msa->au      = NULL;
  msa->ss_cons = NULL;
  msa->sa_cons = NULL;
  msa->pp_cons = NULL;
  msa->rf      = NULL;
  msa->sqacc   = NULL;
  msa->sqdesc  = NULL;
  msa->ss      = NULL;
  msa->sa      = NULL;
  msa->pp      = NULL;
  for (i = 0; i < eslMSA_NCUTS; i++) {
    msa->cutoff[i] = 0.;
    msa->cutset[i] = FALSE;
  }
  msa->sqalloc = nseq;
  msa->sqlen   = NULL;
  msa->sslen   = NULL;
  msa->salen   = NULL;
  msa->pplen   = NULL;
  msa->lastidx = 0;

  /* Unparsed markup, including comments and Stockholm tags.
   * GS, GC, and GR Stockholm tags require keyhash augmentation
   */
  msa->comment        = NULL;
  msa->ncomment       = 0;
  msa->alloc_ncomment = 0;

  msa->gf_tag         = NULL;
  msa->gf             = NULL;
  msa->ngf            = 0;
  msa->alloc_ngf      = 0;

  msa->gs_tag         = NULL;
  msa->gs             = NULL;
  msa->ngs            = 0;

  msa->gc_tag         = NULL;
  msa->gc             = NULL;
  msa->ngc            = 0;

  msa->gr_tag         = NULL;
  msa->gr             = NULL;
  msa->ngr            = 0;

#ifdef eslAUGMENT_KEYHASH
  msa->index     = esl_keyhash_Create();
  msa->gs_idx    = NULL;
  msa->gc_idx    = NULL;
  msa->gr_idx    = NULL;
#endif /*eslAUGMENT_KEYHASH*/

#ifdef eslAUGMENT_SSI
  msa->offset    = 0;
#endif

  /* Allocation, round 2.
   */
  if(nseq > 0) { 
    ESL_ALLOC(msa->sqname, sizeof(char *) * nseq);
    ESL_ALLOC(msa->wgt,    sizeof(double) * nseq);
    ESL_ALLOC(msa->sqlen,  sizeof(int64_t)* nseq);
  }    
  /* Initialize at the second level.
   */
  for (i = 0; i < nseq; i++)
    {
      msa->sqname[i] = NULL;
      msa->sqlen[i]  = 0;
      msa->wgt[i]    = -1.0;	/* "unset so far" */
    }

  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}

/* get_seqidx()
 * 
 * Find the index of a given sequence <name> in an <msa>.
 * If caller has a good guess (for instance, the sequences are
 * coming in a previously seen order in a block of seqs or annotation),
 * the caller can pass this information in <guess>, or -1 if
 * it has no guess.
 * 
 * This function behaves differently depending on whether
 * keyhash augmentation is available or not. Without keyhashing,
 * the name is identified by bruteforce search of the names
 * in the <msa>. With keyhashing, we hash search, which should
 * improve performance for large alignments.
 * 
 * If the name does not already exist in the MSA, then it
 * is assumed to be a new sequence name that we need to store.
 * seqidx is set to msa->nseq, the MSA is Expand()'ed if necessary
 * to make room, the name is stored in msa->sqname[msa->nseq],
 * (and in the hash table, if we're keyhash augmented)
 * and msa->nseq is incremented.
 *
 * Returns:  <eslOK> on success, and the seqidx is 
 *           passed back via <ret_idx>. If <name> is new
 *           in the <msa>, the <name> is stored and the <msa> 
 *           may be internally reallocated if needed.
 *           
 * Throws:   <eslEMEM> if we try to add a name and allocation fails.
 *           <eslEINVAL> if we try to add a name to a non-growable MSA.
 */
static int
get_seqidx(ESL_MSA *msa, char *name, int guess, int *ret_idx)
{
  int seqidx;
  int status;

  *ret_idx = -1;

  /* can we guess? */
  if (guess >= 0 && 
      guess < msa->nseq && 
      strcmp(name, msa->sqname[guess]) == 0) 
    { *ret_idx = guess; return eslOK; }

  /* Else look it up - either brute force
   * or, if we're keyhash-augmented, by hashing.
   */
#ifdef eslAUGMENT_KEYHASH                  
  status = esl_key_Store(msa->index, name, &seqidx);
  if (status == eslEDUP) { *ret_idx = seqidx; return eslOK; }
  if (status != eslOK) return status; /* an error. */
#else
  for (seqidx = 0; seqidx < msa->nseq; seqidx++)
    if (strcmp(msa->sqname[seqidx], name) == 0) break;
  if (seqidx < msa->nseq) 
    { *ret_idx = seqidx; return eslOK; }
#endif

  /* If we reach here, then this is a new name that we're
   * adding.
   */
  if (seqidx >= msa->sqalloc &&  
     (status = esl_msa_Expand(msa)) != eslOK)
    return status; 
    
  status = esl_strdup(name, -1, &(msa->sqname[seqidx]));
  msa->nseq++;
  if (ret_idx != NULL) *ret_idx = seqidx;
  return status;
}


/* msa_get_rlen()
 *
 * Returns the raw (unaligned) length of sequence number <seqidx>
 * in <msa>. 
 */
static int64_t
msa_get_rlen(const ESL_MSA *msa, int seqidx)
{
  int64_t rlen = 0;
  int     pos;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL) rlen = esl_abc_dsqrlen(msa->abc, msa->ax[seqidx]);
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (pos = 0; pos < msa->alen; pos++)
	if (isalnum(msa->aseq[seqidx][pos])) rlen++;
    }
  return rlen;
}


/* set_seq_ss() 
 *
 * Set the secondary structure annotation for sequence number
 * <seqidx> in an alignment <msa> by copying the string <ss>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_ss(ESL_MSA *msa, int seqidx, const char *ss)
{
  int status;
  int i;

  if (msa->ss == NULL) 
    {
      ESL_ALLOC(msa->ss, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) msa->ss[i] = NULL;
    }
  if (msa->ss[seqidx] != NULL) free(msa->ss[seqidx]);
  return (esl_strdup(ss, -1, &(msa->ss[seqidx])));

 ERROR:
  return status;
}

/* set_seq_sa() 
 *
 * Set the surface accessibility annotation for sequence number
 * <seqidx> in an alignment <msa> by copying the string <sa>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_sa(ESL_MSA *msa, int seqidx, const char *sa)
{
  int status;
  int i;

  if (msa->sa == NULL) 
    {
      ESL_ALLOC(msa->sa, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) msa->sa[i] = NULL;
    }
  if (msa->sa[seqidx] != NULL) free(msa->sa[seqidx]);
  return (esl_strdup(sa, -1, &(msa->sa[seqidx])));

 ERROR:
  return status;
}

/* set_seq_pp() 
 *
 * Set the posterior probability annotation for sequence number
 * <seqidx> in an alignment <msa> by copying the string <pp>.
 *
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.
 */
static int
set_seq_pp(ESL_MSA *msa, int seqidx, const char *pp)
{
  int status;
  int i;

  if (msa->pp == NULL) 
    {
      ESL_ALLOC(msa->pp, sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) msa->pp[i] = NULL;
    }
  if (msa->pp[seqidx] != NULL) free(msa->pp[seqidx]);
  return (esl_strdup(pp, -1, &(msa->pp[seqidx])));

 ERROR:
  return status;
}




/* verify_parse()
 *
 * Last function called after a multiple alignment parser thinks it's
 * done. Checks that parse was successful; makes sure required
 * information is present; makes sure required information is
 * consistent. Some fields that are only use during parsing may be
 * freed (sqlen, for example), and some fields are finalized now
 * (<msa->alen> is set, for example). 
 * 
 * <errbuf> is a place to sprintf an informative message about the
 * reason for a parse error. The caller provides an <errbuf>
 * of at least 512 bytes.
 *
 * Returns:  <eslOK>, and errbuf is set to an empty string.
 *           
 * Throws:   <eslEFORMAT> if a problem is detected, and an
 *           informative message about the failure is in errbuf.
 */
static int
verify_parse(ESL_MSA *msa, char *errbuf)
{
  int idx;

  if (msa->nseq == 0) ESL_FAIL(eslEFORMAT, errbuf, "parse error: no alignment data found");

  /* set alen, until proven otherwise; we'll check that the other seqs
   * have the same length later.
   */
  msa->alen = msa->sqlen[0];

  /* We can rely on msa->sqname[] being valid for any index,
   * because of the way the line parsers always store any name
   * they add to the index.
   */
  for (idx = 0; idx < msa->nseq; idx++)
    {
#ifdef eslAUGMENT_ALPHABET
      if ((msa->flags & eslMSA_DIGITAL) &&  (msa->ax  == NULL || msa->ax[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 
#endif
      if (! (msa->flags & eslMSA_DIGITAL) && (msa->aseq == NULL || msa->aseq[idx] == NULL))
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: no sequence for %s",
		 msa->name != NULL ? msa->name : "", msa->sqname[idx]); 

      /* either all weights must be set, or none of them */
      if ((msa->flags & eslMSA_HASWGTS) && msa->wgt[idx] == -1.0)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: expected a weight for seq %s", 
		  msa->name != NULL ? msa->name : "", msa->sqname[idx]);

      /* all aseq must be same length. */
      if (msa->sqlen[idx] != msa->alen)
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: sequence %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->sqlen[idx], msa->alen);

      /* if individual SS is present, it must have length right too */
      if (msa->ss != NULL &&  msa->ss[idx] != NULL &&  msa->sslen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SS for %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->sslen[idx], msa->alen);

				/* if SA is present, must have length right */
      if (msa->sa != NULL && msa->sa[idx] != NULL && msa->salen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR SA for %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->salen[idx], msa->alen);

				/* if PP is present, must have length right */
      if (msa->pp != NULL && msa->pp[idx] != NULL && msa->pplen[idx] != msa->alen) 
	ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GR PP for %s: length %" PRId64 ", expected %" PRId64,
		 msa->name != NULL ? msa->name : "", msa->sqname[idx], msa->pplen[idx], msa->alen);
    }

  /* if cons SS is present, must have length right */
  if (msa->ss_cons != NULL && strlen(msa->ss_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SS_cons markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "",  strlen(msa->ss_cons), msa->alen);

  /* if cons SA is present, must have length right */
  if (msa->sa_cons != NULL && strlen(msa->sa_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC SA_cons markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "",  strlen(msa->sa_cons), msa->alen);

  /* if cons PP is present, must have length right */
  if (msa->pp_cons != NULL && strlen(msa->pp_cons) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC PP_cons markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "",  strlen(msa->pp_cons), msa->alen);

  /* if RF is present, must have length right */
  if (msa->rf != NULL && strlen(msa->rf) != msa->alen) 
    ESL_FAIL(eslEFORMAT, errbuf, "MSA %s parse error: GC RF markup: len %zd, expected %" PRId64,
	     msa->name != NULL ? msa->name : "", strlen(msa->rf), msa->alen);

  /* If no weights were set, set 'em all to 1.0 */
  if (!(msa->flags & eslMSA_HASWGTS))
    for (idx = 0; idx < msa->nseq; idx++)
      msa->wgt[idx] = 1.0;

  /* Clean up a little from the parser */
  if (msa->sqlen != NULL) { free(msa->sqlen); msa->sqlen = NULL; }
  if (msa->sslen != NULL) { free(msa->sslen); msa->sslen = NULL; }
  if (msa->salen != NULL) { free(msa->salen); msa->salen = NULL; }
  if (msa->pplen != NULL) { free(msa->pplen); msa->pplen = NULL; }
  return eslOK;
}


/* Function:  esl_msa_Create()
 * Synopsis:  Creates an <ESL_MSA> object.
 * Incept:    SRE, Sun Jan 23 08:25:26 2005 [St. Louis]
 *
 * Purpose:   Creates and initializes an <ESL_MSA> object, and returns a
 *            pointer to it. 
 *  
 *            If caller already knows the dimensions of the alignment,
 *            both <nseq> and <alen>, then <msa = esl_msa_Create(nseq,
 *            alen)> allocates the whole thing at once. The MSA's
 *            <nseq> and <alen> fields are set accordingly, and the
 *            caller doesn't have to worry about setting them; it can
 *            just fill in <aseq>.
 *            
 *            If caller doesn't know the dimensions of the alignment
 *            (for example, when parsing an alignment file), then
 *            <nseq> is taken to be an initial allocation size, and
 *            <alen> must be -1. <alen=-1> is used as a flag for a
 *            "growable" MSA. For example, the call <msa =
 *            esl_msa_Create(16, -1)>.  allocates internally for an
 *            initial block of 16 sequences, but without allocating
 *            any space for individual sequences.  This allocation can
 *            be expanded (by doubling) by calling <esl_msa_Expand()>.
 *            A created <msa> can only be <_Expand()>'ed if <alen> is
 *            -1.
 *            
 *            In a growable alignment, caller becomes responsible for
 *            memory allocation of each individual <aseq[i]>. Caller
 *            is also responsible for setting <nseq> and <alen> when
 *            it is done parsing and creating the new MSA. In
 *            particular, the <esl_msa_Destroy()> function relies on
 *            <nseq> to know how many individual sequences are
 *            allocated.
 *
 * Args:      <nseq> - number of sequences, or nseq allocation blocksize
 *            <alen> - length of alignment in columns, or -1      
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *           
 * Throws:    <NULL> on allocation failure.          
 *
 * Xref:      squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_Create(int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa;
  int      i;

  msa = create_mostly(nseq, alen); /* aseq is null upon successful return */
  if (msa == NULL) return NULL; /* already threw error in mostly_create, so percolate */

  ESL_ALLOC(msa->aseq,   sizeof(char *) * msa->sqalloc);
  for (i = 0; i < msa->sqalloc; i++)
    msa->aseq[i] = NULL;
  
  if (alen != -1) {
    for (i = 0; i < nseq; i++)
      {
	ESL_ALLOC(msa->aseq[i], sizeof(char) * (alen+1));
	msa->aseq[i][alen] = '\0'; /* caller might forget to null terminate; help the poor */
      }
    msa->nseq = nseq;
  }
  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}


/* Function:  esl_msa_Expand()
 * Synopsis:  Reallocate for more sequences.
 * Incept:    SRE, Sun Jan 23 08:26:30 2005 [St. Louis]
 *
 * Purpose:   Double the current sequence allocation in <msa>.
 *            Typically used when we're reading an alignment sequentially 
 *            from a file, so we don't know nseq 'til we're done.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on reallocation failure; <msa> is undamaged,
 *            and the caller may attempt to recover from the error.
 *            
 *            Throws <eslEINVAL> if <msa> is not growable: its <alen>
 *            field must be -1 to be growable.
 *
 * Xref:      squid's MSAExpand(), 1999.
 */
int
esl_msa_Expand(ESL_MSA *msa)
{
  int   status;
  int   old, new;		/* old & new allocation sizes (max # seqs) */
  void *p;			/* tmp ptr to realloc'ed memory */
  int   i,j;

  if (msa->alen != -1) 
    ESL_EXCEPTION(eslEINVAL, "that MSA is not growable");

  old = msa->sqalloc;
  new = 2*old;

  /* Normally either aseq (ascii) or ax (digitized) would be active, not both.
   * We could make sure that that's true, but that's checked elsewhere.           
   */
  if (msa->aseq != NULL) ESL_RALLOC(msa->aseq, p, sizeof(char *)    * new);
#ifdef eslAUGMENT_ALPHABET
  if (msa->ax   != NULL) ESL_RALLOC(msa->ax,   p, sizeof(ESL_DSQ *) * new);
#endif /*eslAUGMENT_ALPHABET*/

  ESL_RALLOC(msa->sqname, p, sizeof(char *) * new);
  ESL_RALLOC(msa->wgt,    p, sizeof(double) * new);
  ESL_RALLOC(msa->sqlen,  p, sizeof(int64_t)* new);

  if (msa->ss != NULL) 
    {
      ESL_RALLOC(msa->ss,    p, sizeof(char *)  * new);
      ESL_RALLOC(msa->sslen, p, sizeof(int64_t) * new);
    }
  
  if (msa->sa != NULL) 
    {
      ESL_RALLOC(msa->sa,    p, sizeof(char *)  * new);
      ESL_RALLOC(msa->salen, p, sizeof(int64_t) * new);
    }

  if (msa->pp != NULL) 
    {
      ESL_RALLOC(msa->pp,    p, sizeof(char *)  * new);
      ESL_RALLOC(msa->pplen, p, sizeof(int64_t) * new);
    }

  if (msa->sqacc != NULL)
    ESL_RALLOC(msa->sqacc,  p, sizeof(char *) * new);

  if (msa->sqdesc != NULL)
    ESL_RALLOC(msa->sqdesc, p, sizeof(char *) * new);

  for (i = old; i < new; i++)
    {
      if (msa->aseq != NULL) msa->aseq[i] = NULL;
#ifdef eslAUGMENT_ALPHABET
      if (msa->ax   != NULL) msa->ax[i]   = NULL;
#endif /*eslAUGMENT_ALPHABET*/
      msa->sqname[i] = NULL;
      msa->wgt[i]    = -1.0;	/* -1.0 means "unset so far" */
      msa->sqlen[i]  = 0;

      if (msa->ss != NULL) { msa->ss[i] = NULL; msa->sslen[i] = 0; }
      if (msa->sa != NULL) { msa->sa[i] = NULL; msa->salen[i] = 0; }
      if (msa->pp != NULL) { msa->pp[i] = NULL; msa->pplen[i] = 0; }

      if (msa->sqacc  != NULL) msa->sqacc[i]  = NULL;
      if (msa->sqdesc != NULL) msa->sqdesc[i] = NULL;
    }

  /* Reallocate and re-init for unparsed #=GS tags, if we have some.
   * gs is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gs != NULL)
    for (i = 0; i < msa->ngs; i++)
      {
	if (msa->gs[i] != NULL)
	  {
	    ESL_RALLOC(msa->gs[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gs[i][j] = NULL;
	  }
      }
  /* Reallocate and re-init for unparsed #=GR tags, if we have some.
   * gr is [0..ngs-1][0..nseq-1][], so we're reallocing the middle
   * set of pointers.
   */
  if (msa->gr != NULL)
    for (i = 0; i < msa->ngr; i++)
      {
	if (msa->gr[i] != NULL)
	  {
	    ESL_RALLOC(msa->gr[i], p, sizeof(char *) * new);
	    for (j = old; j < new; j++)
	      msa->gr[i][j] = NULL;
	  }
      }

  msa->sqalloc = new;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Copy()
 * Synopsis:  Copies an MSA.
 * Incept:    SRE, Tue Jan 22 15:30:32 2008 [Janelia]
 *
 * Purpose:   Makes a copy of <msa> in <new>. Caller has
 *            already allocated <new> to hold an MSA of
 *            at least <msa->nseq> sequences and <msa->alen>
 *            columns.
 *            
 * Note:      Because MSA's are not reusable, this function does a
 *            lot of internal allocation for optional fields, without
 *            checking <new> to see if space was already allocated. To
 *            reuse an MSA <new> and copy new data into it, we'll
 *            eventually need a <esl_msa_Reuse()> function, and/or
 *            recode this to reuse or free any already-allocated
 *            optional memory it encounters in <new>. Until then, 
 *            it's unlikely that <esl_msa_Copy()> is useful on its own;
 *            the caller would be expected to call <esl_msa_Clone()> 
 *            instead.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. In this case, <new>
 *            was only partially constructed, and should be treated
 *            as corrupt.
 */
int
esl_msa_Copy(const ESL_MSA *msa, ESL_MSA *new)
{
  int i, x, j;
  int status;

  /* aseq[0..nseq-1][0..alen-1] strings,
   * or ax[0..nseq-1][(0) 1..alen (alen+1)] digital seqs 
   * <new> must have one of them allocated already.
   */
  if (! (msa->flags & eslMSA_DIGITAL))
    for (i = 0; i < msa->nseq; i++)
      strcpy(new->aseq[i], msa->aseq[i]);
#ifdef eslAUGMENT_ALPHABET
  else
    {
      for (i = 0; i < msa->nseq; i++)
	memcpy(new->ax[i], msa->ax[i], (msa->alen+2) * sizeof(ESL_DSQ));
      new->abc = msa->abc;
    }
#endif
  
  for (i = 0; i < msa->nseq; i++) {
    esl_strdup(msa->sqname[i], -1, &(new->sqname[i]));
    new->wgt[i] = msa->wgt[i];
  }
  /* alen, nseq were already set by Create() */
  new->flags = msa->flags;

  esl_strdup(msa->name, -1, &(new->name));
  esl_strdup(msa->desc, -1, &(new->desc));
  esl_strdup(msa->acc,  -1, &(new->acc));
  esl_strdup(msa->au,    -1, &(new->au));
  esl_strdup(msa->ss_cons,   -1, &(new->ss_cons));
  esl_strdup(msa->sa_cons,   -1, &(new->sa_cons));
  esl_strdup(msa->pp_cons,   -1, &(new->pp_cons));
  esl_strdup(msa->rf,    -1, &(new->rf));

  if (msa->sqacc != NULL) {
    ESL_ALLOC(new->sqacc, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sqacc[i], -1, &(new->sqacc[i]));
  }
  if (msa->sqdesc != NULL) {
    ESL_ALLOC(new->sqdesc, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sqdesc[i], -1, &(new->sqdesc[i]));
  }
  if (msa->ss != NULL) {
    ESL_ALLOC(new->ss, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->ss[i], -1, &(new->ss[i]));
  }
  if (msa->sa != NULL) {
    ESL_ALLOC(new->sa, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->sa[i], -1, &(new->sa[i]));
  }
  if (msa->pp != NULL) {
    ESL_ALLOC(new->pp, sizeof(char **) * msa->nseq);
    for (i = 0; i < msa->nseq; i++)
      esl_strdup(msa->pp[i], -1, &(new->pp[i]));
  }
  
  for (x = 0; x < eslMSA_NCUTS; x++) {
    new->cutoff[x] = msa->cutoff[x];
    new->cutset[x] = msa->cutset[x];
  }

  if (msa->ncomment > 0) {
    ESL_ALLOC(new->comment, sizeof(char **) * msa->ncomment);
    new->ncomment       = msa->ncomment;
    new->alloc_ncomment = msa->ncomment;
    for (i = 0; i < msa->ncomment; i++)
      esl_strdup(msa->comment[i], -1, &(new->comment[i]));
  }

  if (msa->ngf > 0) {
    ESL_ALLOC(new->gf_tag, sizeof(char **) * msa->ngf);
    ESL_ALLOC(new->gf,     sizeof(char **) * msa->ngf);
    new->ngf       = msa->ngf;
    new->alloc_ngf = msa->ngf;
    for (i = 0; i < msa->ngf; i++) {
      esl_strdup(msa->gf_tag[i], -1, &(new->gf_tag[i]));
      esl_strdup(msa->gf[i],     -1, &(new->gf[i]));
    }
  }

  if (msa->ngs > 0) {
    ESL_ALLOC(new->gs_tag, sizeof(char **)  * msa->ngs);
    ESL_ALLOC(new->gs,     sizeof(char ***) * msa->ngs);
    new->ngs       = msa->ngs;
    for (i = 0; i < msa->ngs; i++) {
      ESL_ALLOC(new->gs[i], sizeof(char **) * msa->nseq);
      esl_strdup(msa->gs_tag[i], -1, &(new->gs_tag[i]));
      for (j = 0; j < msa->nseq; j++)
	esl_strdup(msa->gs[i][j],  -1, &(new->gs[i][j]));
    }
  }

  if (msa->ngc > 0) {
    ESL_ALLOC(new->gc_tag, sizeof(char **) * msa->ngc);
    ESL_ALLOC(new->gc,     sizeof(char **) * msa->ngc);
    new->ngc       = msa->ngc;
    for (i = 0; i < msa->ngc; i++) {
      esl_strdup(msa->gc_tag[i], -1, &(new->gc_tag[i]));
      esl_strdup(msa->gc[i],     -1, &(new->gc[i]));
    }
  }
  
  if (msa->ngr > 0) {
    ESL_ALLOC(new->gr_tag, sizeof(char **)  * msa->ngr);
    ESL_ALLOC(new->gr,     sizeof(char ***) * msa->ngr);
    new->ngr       = msa->ngr;
    for (i = 0; i < msa->ngr; i++) {
      ESL_ALLOC(new->gr[i], sizeof(char **) * msa->nseq);
      esl_strdup(msa->gr_tag[i], -1, &(new->gr_tag[i]));
      for (j = 0; j < msa->nseq; j++)
	esl_strdup(msa->gr[i][j],  -1, &(new->gr[i][j]));
    }
  }

#ifdef eslAUGMENT_KEYHASH
  esl_keyhash_Destroy(new->index);  new->index  = NULL;
  esl_keyhash_Destroy(new->gs_idx); new->gs_idx = NULL;
  esl_keyhash_Destroy(new->gc_idx); new->gc_idx = NULL;
  esl_keyhash_Destroy(new->gr_idx); new->gr_idx = NULL;

  if (msa->index  != NULL) new->index  = esl_keyhash_Clone(msa->index);
  if (msa->gs_idx != NULL) new->gs_idx = esl_keyhash_Clone(msa->gs_idx);
  if (msa->gc_idx != NULL) new->gc_idx = esl_keyhash_Clone(msa->gc_idx);
  if (msa->gr_idx != NULL) new->gr_idx = esl_keyhash_Clone(msa->gr_idx);
#endif

#ifdef eslAUGMENT_SSI
  new->offset = msa->offset;
#endif

  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Clone()
 * Synopsis:  Duplicates an MSA.
 * Incept:    SRE, Tue Jan 22 15:23:55 2008 [Janelia]
 *
 * Purpose:   Make a duplicate of <msa>, in newly 
 *            allocated space. 
 *
 * Returns:   a pointer to the newly allocated clone.
 *            Caller is responsible for free'ing it.
 *
 * Throws:    <NULL> on allocation error.
 */
ESL_MSA *
esl_msa_Clone(const ESL_MSA *msa)
{
  ESL_MSA *nw = NULL;
  int      status;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL) {
      if ((nw = esl_msa_CreateDigital(msa->abc, msa->nseq, msa->alen)) == NULL)  return NULL;
  } else
#endif
  if ((nw     = esl_msa_Create(msa->nseq, msa->alen)) == NULL)  return NULL;  

  if ((status = esl_msa_Copy(msa, nw) )               != eslOK) goto ERROR;
  return nw;

 ERROR:
  esl_msa_Destroy(nw);
  return NULL;
}


/* Function:  esl_msa_Destroy()
 * Synopsis:  Frees an <ESL_MSA>.
 * Incept:    SRE, Sun Jan 23 08:26:02 2005 [St. Louis]
 *
 * Purpose:   Destroys <msa>.
 *
 * Xref:      squid's MSADestroy().
 */
void
esl_msa_Destroy(ESL_MSA *msa)
{
  if (msa == NULL) return;

  if (msa->aseq != NULL) 
    esl_Free2D((void **) msa->aseq, msa->nseq);
#ifdef eslAUGMENT_ALPHABET
  if (msa->ax != NULL) 
    esl_Free2D((void **) msa->ax, msa->nseq);
#endif /*eslAUGMENT_ALPHABET*/

  esl_Free2D((void **) msa->sqname, msa->nseq);
  esl_Free2D((void **) msa->sqacc,  msa->nseq);
  esl_Free2D((void **) msa->sqdesc, msa->nseq);
  esl_Free2D((void **) msa->ss,     msa->nseq);
  esl_Free2D((void **) msa->sa,     msa->nseq);
  esl_Free2D((void **) msa->pp,     msa->nseq);

  if (msa->sqlen   != NULL) free(msa->sqlen);
  if (msa->wgt     != NULL) free(msa->wgt);

  if (msa->name    != NULL) free(msa->name);
  if (msa->desc    != NULL) free(msa->desc);
  if (msa->acc     != NULL) free(msa->acc);
  if (msa->au      != NULL) free(msa->au);
  if (msa->ss_cons != NULL) free(msa->ss_cons);
  if (msa->sa_cons != NULL) free(msa->sa_cons);
  if (msa->pp_cons != NULL) free(msa->pp_cons);
  if (msa->rf      != NULL) free(msa->rf);
  if (msa->sslen   != NULL) free(msa->sslen);
  if (msa->salen   != NULL) free(msa->salen);
  if (msa->pplen   != NULL) free(msa->pplen);  

  esl_Free2D((void **) msa->comment, msa->ncomment);
  esl_Free2D((void **) msa->gf_tag,  msa->ngf);
  esl_Free2D((void **) msa->gf,      msa->ngf);

  esl_Free2D((void **) msa->gs_tag,  msa->ngs);
  esl_Free3D((void ***)msa->gs,      msa->ngs, msa->nseq);
  esl_Free2D((void **) msa->gc_tag,  msa->ngc);
  esl_Free2D((void **) msa->gc,      msa->ngc);
  esl_Free2D((void **) msa->gr_tag,  msa->ngr);
  esl_Free3D((void ***)msa->gr,      msa->ngr, msa->nseq);

#ifdef eslAUGMENT_KEYHASH
  esl_keyhash_Destroy(msa->index);
  esl_keyhash_Destroy(msa->gs_idx);
  esl_keyhash_Destroy(msa->gc_idx);
  esl_keyhash_Destroy(msa->gr_idx);
#endif /* keyhash augmentation */  

  free(msa);
  return;
}


/* Function:  esl_msa_SetName()
 * Synopsis:  Set name of an MSA.
 * Incept:    SRE, Sat Feb 23 18:42:47 2008 [Casa de Gatos]
 *
 * Purpose:   Sets the name of the msa <msa> to <name>. 
 *
 *            <name> can be <NULL>, because the MSA name is an
 *            optional field; in which case any existing name in
 *            the <msa> is erased.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetName(ESL_MSA *msa, const char *name)
{
  int     status;

  if (msa->name != NULL) free(msa->name); 
  status = esl_strdup(name, -1, &(msa->name));
  return status;
}


/* Function:  esl_msa_SetDesc()
 * Synopsis:  Set the description line of an MSA.
 * Incept:    SRE, Sat Feb 23 18:47:06 2008 [Casa de Gatos]
 *
 * Purpose:   Sets the description line of the msa <msa> to <desc>. 
 *
 *            As a special case, <desc> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetDesc(ESL_MSA *msa, const char *desc)
{
  int     status;

  if (msa->desc != NULL) free(msa->desc);
  status = esl_strdup(desc, -1, &(msa->desc));
  return status;

}

/* Function:  esl_msa_SetAccession()
 * Synopsis:  Set the accession number of an MSA.
 * Incept:    SRE, Sat Feb 23 18:49:04 2008 [Casa de Gatos]
 *
 * Purpose:   Sets accession number of the msa <msa> to <acc>. 
 *
 *            As a special case, <acc> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetAccession(ESL_MSA *msa, const char *acc)
{
  int     status;

  if (msa->acc != NULL) free(msa->acc);
  status = esl_strdup(acc, -1, &(msa->acc));
  return status;
}


/* Function:  esl_msa_SetAuthor()
 * Synopsis:  Set the author string in an MSA.
 * Incept:    SRE, Wed Mar  4 10:41:21 2009 [Janelia]
 *
 * Purpose:   Sets the author string in <msa> to <author>.
 *            
 *            As a special case, <author> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetAuthor(ESL_MSA *msa, const char *author)
{
  int     status;

  if (msa->au != NULL) free(msa->au);
  status = esl_strdup(author, -1, &(msa->au));
  return status;
}


/* Function:  esl_msa_SetSeqName()
 * Synopsis:  Set an individual sequence name in an MSA.
 * Incept:    SRE, Wed Mar  4 10:56:28 2009 [Janelia]
 *
 * Purpose:   Set the name of sequence number <idx> in <msa>
 *            to <name>.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <name> is <NULL>;
 *            <eslEMEM> on allocation error.
 *
 * Note:      msa->sqname[] is not optional, so we may
 *            rely on it already being allocated for 
 *            i=0..sqalloc-1.
 */
int
esl_msa_SetSeqName(ESL_MSA *msa, int idx, const char *name)
{
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (name == NULL)         ESL_EXCEPTION(eslEINVAL, "seq names are mandatory; NULL is not a valid name");

  if (msa->sqname[idx] != NULL) free(msa->sqname[idx]);
  status = esl_strdup(name, -1, &(msa->sqname[idx]));
  return status;
}

/* Function:  esl_msa_SetSeqAccession()
 * Synopsis:  Sets individual sequence accession in an MSA.
 * Incept:    SRE, Wed Mar  4 11:03:26 2009 [Janelia]
 *
 * Purpose:   Set the accession of sequence number <idx> in <msa> to
 *            <acc>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int 
esl_msa_SetSeqAccession(ESL_MSA *msa, int idx, const char *acc)
{
  int     i;
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (acc == NULL) {
    if (msa->sqacc != NULL) { free(msa->sqacc[idx]); msa->sqacc[idx] = NULL; }
    return eslOK;
  }

  /* Allocate/initialize the optional sqacc array, if it's not already done: */
  if (msa->sqacc == NULL) {
    ESL_ALLOC(msa->sqacc, sizeof(char *) * msa->sqalloc);
    for (i = 0; i < msa->sqalloc; i++) msa->sqacc[i] = NULL;
  } 
  if (msa->sqacc[idx] != NULL) free(msa->sqacc[idx]);

  status = esl_strdup(acc, -1, &(msa->sqacc[idx]));
  return status;

 ERROR:
  return status;
}
  
/* Function:  esl_msa_SetSeqDescription()
 * Synopsis:  Sets individual sequence description in an MSA.
 * Incept:    SRE, Wed Mar  4 11:09:37 2009 [Janelia]
 *
 * Purpose:   Set the description of sequence number <idx> in <msa> to
 *            <desc>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
esl_msa_SetSeqDescription(ESL_MSA *msa, int idx, const char *desc)
{
  int     i;
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (desc == NULL) {
    if (msa->sqdesc != NULL) { free(msa->sqdesc[idx]); msa->sqdesc[idx] = NULL; }
    return eslOK;
  }

  /* Allocate/initialize the optional sqdesc array, if it's not already done: */
  if (msa->sqdesc == NULL) {
    ESL_ALLOC(msa->sqdesc, sizeof(char *) * msa->sqalloc);
    for (i = 0; i < msa->sqalloc; i++) msa->sqdesc[i] = NULL;
  } 
  if (msa->sqdesc[idx] != NULL) free(msa->sqdesc[idx]);

  status = esl_strdup(desc, -1, &(msa->sqdesc[idx]));
  return status;

 ERROR:
  return status;
}


/* Function:  esl_msa_FormatName()
 * Synopsis:  Format name of an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:33:34 2009 [Janelia]
 *
 * Purpose:   Sets the name of the msa <msa> using <name>, where 
 *            <name> is a <printf()>-style format with
 *            arguments; for example, <esl_msa_FormatName(msa, "random%d", i)>.
 *            
 *            <name> can be <NULL>, because the MSA name is an
 *            optional field; in which case any existing name in
 *            the <msa> is erased.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_msa_FormatName(ESL_MSA *msa, const char *name, ...)
{
  va_list ap;
  int     status;

  if (msa->name != NULL) free(msa->name); 
  if (name      == NULL) { msa->name = NULL; return eslOK; }

  va_start(ap, name);
  status = esl_vsprintf(&(msa->name), name, &ap);
  va_end(ap);
  return status;
}


/* Function:  esl_msa_FormatDesc()
 * Synopsis:  Format the description line of an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:34:25 2009 [Janelia]
 *
 * Purpose:   Format the description line of the msa <msa> using <desc>.
 *            where <desc> is a <printf()>-style format with
 *            arguments.
 *            For example, <esl_msa_FormatDesc(msa, "sample %d", i)>.
 *
 *            As a special case, <desc> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_msa_FormatDesc(ESL_MSA *msa, const char *desc, ...)
{
  va_list ap;
  int     status;

  if (msa->desc != NULL) free(msa->desc);
  va_start(ap, desc);
  status = esl_vsprintf(&(msa->desc), desc, &ap);
  va_end(ap);
  return status;

}

/* Function:  esl_msa_FormatAccession()
 * Synopsis:  Format the accession number of an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:35:24 2009 [Janelia].
 *
 * Purpose:   Sets accession number of the msa <msa> using <acc>, 
 *            where <acc> is a <printf()>-style format with arguments.
 *            For example, <esl_msa_FormatAccession(msa, "PF%06d", i)>.
 *
 *            As a special case, <acc> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_msa_FormatAccession(ESL_MSA *msa, const char *acc, ...)
{
  va_list ap;
  int     status;

  if (msa->acc != NULL) free(msa->acc);
  va_start(ap, acc);
  status = esl_vsprintf(&(msa->acc), acc, &ap);
  va_end(ap);
  return status;
}


/* Function:  esl_msa_FormatAuthor()
 * Synopsis:  Format the author string in an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:36:05 2009 [Janelia]
 *
 * Purpose:   Sets the author string in <msa>, using an <author> string
 *            and arguments in same format as <printf()> would take.
 *            
 *            As a special case, <author> may be <NULL>, to facilitate
 *            handling of optional annotation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_msa_FormatAuthor(ESL_MSA *msa, const char *author, ...)
{
  va_list ap;
  int     status;

  if (msa->au != NULL) free(msa->au);
  va_start(ap, author);
  status = esl_vsprintf(&(msa->au), author, &ap);
  va_end(ap);
  return status;
}


/* Function:  esl_msa_FormatSeqName()
 * Synopsis:  Formats an individual sequence name in an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:36:35 2009 [Janelia]
 *
 * Purpose:   Set the name of sequence number <idx> in <msa>
 *            to <name>, where <name> is a <printf()>
 *            style format and arguments.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <name> is <NULL>;
 *            <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 *
 * Note:      msa->sqname[] is not optional, so we may
 *            rely on it already being allocated for 
 *            i=0..sqalloc-1.
 */
int
esl_msa_FormatSeqName(ESL_MSA *msa, int idx, const char *name, ...)
{
  va_list ap;
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (name == NULL)         ESL_EXCEPTION(eslEINVAL, "seq names are mandatory; NULL is not a valid name");

  if (msa->sqname[idx] != NULL) free(msa->sqname[idx]);

  va_start(ap, name);
  status = esl_vsprintf(&(msa->sqname[idx]), name, &ap);
  va_end(ap);
  return status;
}

/* Function:  esl_msa_FormatSeqAccession()
 * Synopsis:  Format individual sequence accession in an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:37:08 2009 [Janelia]
 *
 * Purpose:   Set the accession of sequence number <idx> in <msa> to
 *            <acc>, where <acc> is a <printf()> style format and
 *            arguments.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int 
esl_msa_FormatSeqAccession(ESL_MSA *msa, int idx, const char *acc, ...)
{
  va_list ap;
  int     i;
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (acc == NULL) {
    if (msa->sqacc != NULL) { free(msa->sqacc[idx]); msa->sqacc[idx] = NULL; }
    return eslOK;
  }

  /* Allocate/initialize the optional sqacc array, if it's not already done: */
  if (msa->sqacc == NULL) {
    ESL_ALLOC(msa->sqacc, sizeof(char *) * msa->sqalloc);
    for (i = 0; i < msa->sqalloc; i++) msa->sqacc[i] = NULL;
  } 
  if (msa->sqacc[idx] != NULL) free(msa->sqacc[idx]);

  va_start(ap, acc);
  status = esl_vsprintf(&(msa->sqacc[idx]), acc, &ap);
  va_end(ap);
  return status;

 ERROR:
  return status;
}
  
/* Function:  esl_msa_FormatSeqDescription()
 * Synopsis:  Formats individual sequence description in an MSA, printf()-style.
 * Incept:    SRE, Fri Sep 11 11:37:35 2009 [Janelia]
 *
 * Purpose:   Set the description of sequence number <idx> in <msa> to
 *            <desc>, where <desc> may be a <printf()> style format and
 *            arguments.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_msa_FormatSeqDescription(ESL_MSA *msa, int idx, const char *desc, ...)
{
  va_list ap;
  int     i;
  int     status;

  if (idx  >= msa->sqalloc) ESL_EXCEPTION(eslEINVAL, "no such sequence %d (only %d allocated)", idx, msa->sqalloc);
  if (desc == NULL) {
    if (msa->sqdesc != NULL) { free(msa->sqdesc[idx]); msa->sqdesc[idx] = NULL; }
    return eslOK;
  }

  /* Allocate/initialize the optional sqdesc array, if it's not already done: */
  if (msa->sqdesc == NULL) {
    ESL_ALLOC(msa->sqdesc, sizeof(char *) * msa->sqalloc);
    for (i = 0; i < msa->sqalloc; i++) msa->sqdesc[i] = NULL;
  } 
  if (msa->sqdesc[idx] != NULL) free(msa->sqdesc[idx]);

  va_start(ap, desc);
  status = esl_vsprintf(&(msa->sqdesc[idx]), desc, &ap);
  va_end(ap);
  return status;

 ERROR:
  return status;
}
/*---------------------- end of ESL_MSA functions ---------------------------*/






/******************************************************************************
 *# 2. The <ESL_MSAFILE> object                                       
 *****************************************************************************/

/* msafile_open():
 * this is the routine that actually opens an ESL_MSAFILE;
 * esl_msafile_Open() and esl_msafile_OpenDigital() are wrappers around it.
 */
static int
msafile_open(const char *filename, int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  ESL_MSAFILE *afp     = NULL;
  char        *ssifile = NULL;
  char        *envfile = NULL;
  char        *cmd     = NULL;
  int          n       = strlen(filename);
  int          status;
  
  ESL_ALLOC(afp, sizeof(ESL_MSAFILE));
  afp->f          = NULL;
  afp->fname      = NULL;
  afp->linenumber = 0;
  afp->errbuf[0]  = '\0';
  afp->buf        = NULL;
  afp->buflen     = 0;
  afp->do_gzip    = FALSE;
  afp->do_stdin   = FALSE;
  afp->format     = eslMSAFILE_UNKNOWN;	
  afp->do_digital = FALSE;
#ifdef eslAUGMENT_ALPHABET
  afp->abc        = NULL;	        
#endif
#ifdef eslAUGMENT_SSI		
  afp->ssi        = NULL;	         
#endif  
  afp->msa_cache  = NULL;

  if (strcmp(filename, "-") == 0)
    {
      afp->f         = stdin;
      afp->do_stdin  = TRUE; 
      if ((status = esl_strdup("[STDIN]", -1, &(afp->fname))) != eslOK) goto ERROR;
    }
#ifdef HAVE_POPEN
  /* popen(), pclose() aren't portable to non-POSIX systems; 
   * disable this section in strict ANSI C mode.
   */
  /* tricky: if n= length of a string s, then
   * s+n-i repositions pointer s at the last i chars
   * of the string.
   */
  else if (n > 3 && strcmp(filename+n-3, ".gz") == 0)
    {
      /* Note that popen() will return "successfully"
       * if file doesn't exist, because gzip works fine
       * and prints an error! So we have to check for
       * existence of file ourself.
       */
      if (! esl_FileExists(filename))	      { status = eslENOTFOUND; goto ERROR; }
      ESL_ALLOC(cmd, sizeof(char) * (n+1+strlen("gzip -dc ")));
      sprintf(cmd, "gzip -dc %s", filename);
      if ((afp->f = popen(cmd, "r")) == NULL) { status = eslENOTFOUND; goto ERROR; }
      if ((status = esl_strdup(filename, n, &(afp->fname))) != eslOK)  goto ERROR;
      afp->do_gzip  = TRUE;
    }
#endif /*HAVE_POPEN*/
  else	/* Normal file open or env file open: set ssifile */
    {
      /* When we open a file, it may be either in the current
       * directory, or in the directory indicated by the env
       * argument - and we construct an SSI filename accordingly.
       * (Whether or not we're SSI augmented, in fact, for simplicity.)
       */
      if ((afp->f = fopen(filename, "r")) != NULL)
	{
	  if (esl_strdup(filename, n, &ssifile)                      != eslOK) goto ERROR;
	  if (esl_strcat(&ssifile, n, ".ssi", 4)                     != eslOK) goto ERROR;
	  if ((status = esl_strdup(filename, n, &(afp->fname)))      != eslOK) goto ERROR;
	}
      else if (esl_FileEnvOpen(filename, env, &(afp->f), &envfile) == eslOK)
	{
	  if (esl_strdup(envfile, n, &ssifile)                      != eslOK) goto ERROR;
	  if (esl_strcat(&ssifile, n, ".ssi", 4)                    != eslOK) goto ERROR;
	  if ((status = esl_strdup(envfile, n, &(afp->fname)))      != eslOK) goto ERROR;
	}
      else 
	{ status = eslENOTFOUND; goto ERROR;}
    }

#ifdef eslAUGMENT_SSI
  /* If augmented by SSI indexing:
   * Open the SSI index file. If it doesn't exist, or
   * it's corrupt, or some error happens, afp->ssi stays NULL.
   * We should warn, probably, or provide some way for caller to 
   * to know that we've opened the index successfully or not.
   */
  status = esl_ssi_Open(ssifile, &(afp->ssi)); 
#endif

  /* Invoke autodetection if we haven't already been told what
   * to expect.
   */
  if (format == eslMSAFILE_UNKNOWN)
    {
      if (afp->do_stdin == TRUE || afp->do_gzip)
	ESL_XEXCEPTION(eslEINVAL, "Can't autodetect alignment file fmt in stdin, gzip pipe");
      if (esl_msa_GuessFileFormat(afp) != eslOK)
	{ status = eslEFORMAT; goto ERROR; }
    }
  else 
    afp->format     = format;

  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  if (cmd     != NULL) free(cmd);
  *ret_msafp = afp;
  return eslOK;

 ERROR:
  if (envfile != NULL) free(envfile);
  if (ssifile != NULL) free(ssifile);
  if (cmd     != NULL) free(cmd);
  if (afp     != NULL) esl_msafile_Close(afp); 
  *ret_msafp = NULL;
  return status;
}


/* Function: esl_msafile_Open()
 * Synopsis: Open an MSA file for input.
 * Date:     SRE, Sun Jan 23 08:30:33 2005 [St. Louis]
 *
 * Purpose:  Open an alignment database file <filename> and prepare for
 *           reading one alignment, or sequentially in the case of 
 *           multiple MSA databases (e.g. Stockholm format); returns
 *           the opened file pointer in <ret_msafp>.
 *          
 *           There are one or two special cases for <filename>. If
 *           <filename> is "-", then the alignment is read from
 *           <stdin>. If <filename> ends in ".gz", then the file is
 *           assumed to be compressed by gzip, and it is opened as a
 *           pipe from <gunzip -dc>. (Auto-decompression of gzipp'ed files
 *           is only available on POSIX-compliant systems w/ popen(), when 
 *           <HAVE_POPEN> is defined at compile-time.)
 *          
 *           If <env> is non-NULL, then we look for <filename> in
 *           one or more directories in a colon-delimited list
 *           that is the value of the environment variable <env>.
 *           For example, if we had 
 *              <setenv HMMERDB /nfs/db/Pfam:/nfs/db/Rfam> 
 *           in the environment, a profile HMM application
 *           might pass "HMMERDB" as <env>.
 *          
 *          The file is asserted to be in format <fmt>, which is
 *          either a known format like <eslMSAFILE_STOCKHOLM>, or
 *          <eslMSAFILE_UNKNOWN>; if <fmt> is <eslMSAFILE_UNKNOWN>,
 *          then format autodetection is invoked.
 *
 * Returns:  <eslOK> on success, and <ret_msafp> is set to point at
 *           an open <ESL_MSAFILE>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           Returns <eslENOTFOUND> if <filename> cannot be opened,
 *           or <eslEFORMAT> if autodetection is attempted and 
 *           format cannot be determined.
 *           
 * Throws:   <eslEMEM> on allocation failure.
 *           <eslEINVAL> if format autodetection is attempted on 
 *           stdin or a gunzip pipe.
 *
 * Xref:     squid's MSAFileOpen(), 1999.
 * 
 * Note      Implemented as a wrapper around msafile_open(), because
 *           esl_msafile_OpenDigital() shares almost all the same code.
 */
int
esl_msafile_Open(const char *filename, int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  return msafile_open(filename, format, env, ret_msafp);
}

/* Function:  esl_msafile_Close()
 * Synopsis:  Closes an open MSA file.
 * Incept:    SRE, Sun Jan 23 08:18:39 2005 [St. Louis]
 *
 * Purpose:   Close an open <ESL_MSAFILE>.
 *
 * Xref:      squid's MSAFileClose().
 */
void
esl_msafile_Close(ESL_MSAFILE *afp)
{
  if (afp == NULL) return;

#ifdef HAVE_POPEN /* gzip functionality */
  if (afp->do_gzip && afp->f != NULL)    pclose(afp->f);
#endif
  if (!afp->do_gzip && ! afp->do_stdin && afp->f != NULL) fclose(afp->f);
  if (afp->fname != NULL) free(afp->fname);
  if (afp->buf  != NULL)  free(afp->buf);
#ifdef eslAUGMENT_SSI
  if (afp->ssi  != NULL)  esl_ssi_Close(afp->ssi); 
#endif /* eslAUGMENT_SSI*/
  if (afp->msa_cache != NULL) esl_msa_Destroy(afp->msa_cache);
  free(afp);
}
/*-------------------- end of ESL_MSAFILE functions -------------------------*/






/******************************************************************************
 *# 3. Digital mode MSA's (augmentation: alphabet)
 *****************************************************************************/
#ifdef eslAUGMENT_ALPHABET
/* Function:  esl_msa_GuessAlphabet()
 * Synopsis:  Guess alphabet of MSA.
 * Incept:    SRE, Fri May 18 09:55:08 2007 [Janelia]
 *
 * Purpose:   Guess whether the sequences in the <msa> are
 *            <eslDNA>, <eslRNA>, or <eslAMINO>, and return
 *            that guess in <*ret_type>.
 *            
 *            The determination is made based on the classifications
 *            of the individual sequences in the alignment. At least
 *            one sequence must contain ten residues or more to be
 *            classified. If one or more sequences is called
 *            <eslAMINO> and one or more is called <eslDNA>/<eslRNA>,
 *            the alignment's alphabet is considered to be
 *            indeterminate (<eslUNKNOWN>). If some sequences are
 *            <eslDNA> and some are <eslRNA>, the alignment is called
 *            <eslDNA>; this should cause no problems, because Easel
 *            reads U as a synonym for T in DNA sequence anyway.
 *            
 *            Tested on Pfam 21.0 and Rfam 7.0, this routine correctly
 *            classified all 8957 Pfam alignments as protein, and 503
 *            Rfam alignments as RNA (both seed and full alignments).
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set
 *            to <eslDNA>, <eslRNA>, or <eslAMINO>. 
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to
 *            <eslUNKNOWN> if the alphabet cannot be reliably guessed.
 *
 * Xref:      J1/62
 */
int
esl_msa_GuessAlphabet(const ESL_MSA *msa, int *ret_type)
{
  int64_t namino   = 0,
          ndna     = 0,
          nrna     = 0,
          nunknown = 0;
  int     type;
  int     i,x;
  int64_t j,n;
  int64_t ct[26];

  if (msa->flags & eslMSA_DIGITAL) { *ret_type = msa->abc->type; return eslOK; }

  *ret_type = eslUNKNOWN;

  /* On wide alignments, we're better off looking at individual sequence
   * classifications. We don't want to end up calling the whole alignment
   * indeterminate just because a few sequences have degenerate residue
   * codes.
   */
  for (i = 0; i < msa->nseq; i++) 
    {
      for (x = 0; x < 26; x++) ct[x] = 0;
      for (n = 0, j = 0; j < msa->alen; j++) {
	x = toupper(msa->aseq[i][j]) - 'A';
	if (x < 0 || x > 26) continue;
	ct[x]++;
	n++;
	if (n > 10000) break;	/* ought to know by now */
      }
      esl_abc_GuessAlphabet(ct, &type);

      switch (type) {
      case eslAMINO:   namino++; break;
      case eslDNA:     ndna++;   break;
      case eslRNA:     nrna++;   break;
      default:         nunknown++; 
      }
    }
  if      (namino    > 0 && (ndna+nrna)   == 0) *ret_type = eslAMINO;
  else if (ndna      > 0 && (nrna+namino) == 0) *ret_type = eslDNA;
  else if (nrna      > 0 && (ndna+namino) == 0) *ret_type = eslRNA;
  else if (ndna+nrna > 0 && namino        == 0) *ret_type = eslDNA;

  /* On narrow alignments, no single sequence may be long enough to 
   * be classified, but we can determine alphabet from composition
   * of the complete alignment. Of course, degenerate residue codes in
   * a DNA alignment will still screw us.
   */
  if (*ret_type == eslUNKNOWN)
    {

      n = 0;
      for (x = 0; x < 26; x++) ct[x] = 0;
      for (i = 0; i < msa->nseq; i++) {
	for (j = 0; j < msa->alen; j++) {
	  x = toupper(msa->aseq[i][j]) - 'A';
	  if (x < 0 || x > 26) continue;
	  ct[x]++;
	  n++;
	  if (n > 10000) break;	/* ought to know by now */
	}
	if (n > 10000) break;	
      }
      esl_abc_GuessAlphabet(ct, ret_type);
    }

  if (*ret_type == eslUNKNOWN) return eslEAMBIGUOUS;
  else                         return eslOK;
}


/* Function:  esl_msa_CreateDigital()
 * Synopsis:  Create a digital <ESL_MSA>.
 * Incept:    SRE, Sun Aug 27 16:49:58 2006 [Leesburg]
 *
 * Purpose:   Same as <esl_msa_Create()>, except the returned MSA is configured
 *            for a digital alignment using internal alphabet <abc>, instead of 
 *            a text alignment.
 *   
 *            Internally, this means the <ax> field is allocated instead of
 *            the <aseq> field, and the <eslMSA_DIGITAL> flag is raised.
 *
 * Args:     <nseq> - number of sequences, or nseq allocation blocksize
 *           <alen> - length of alignment in columns, or -1
 *
 * Returns:   pointer to new MSA object, w/ all values initialized.
 *            Note that <msa->nseq> is initialized to 0, even though space
 *            is allocated.
 *           
 * Throws:    NULL on allocation failure.          
 *
 * Xref:      squid's MSAAlloc()
 */
ESL_MSA *
esl_msa_CreateDigital(const ESL_ALPHABET *abc, int nseq, int64_t alen)
{
  int      status;
  ESL_MSA *msa;
  int      i;

  msa = create_mostly(nseq, alen); /* aseq is null upon successful return */
  if (msa == NULL) return NULL; /* already threw error in mostly_create, so percolate */

  ESL_ALLOC(msa->ax,   sizeof(ESL_DSQ *) * msa->sqalloc); 
  for (i = 0; i < msa->sqalloc; i++)
    msa->ax[i] = NULL;

  if (alen != -1)
    {
      for (i = 0; i < nseq; i++) {
	ESL_ALLOC(msa->ax[i], sizeof(ESL_DSQ) * (alen+2));
	msa->ax[i][0] = msa->ax[i][alen+1] = eslDSQ_SENTINEL; /* help the poor */
      }
      msa->nseq = nseq;
    }

  msa->abc    = (ESL_ALPHABET *) abc; /* this cast away from const-ness is deliberate & safe. */
  msa->flags |= eslMSA_DIGITAL;
  return msa;

 ERROR:
  esl_msa_Destroy(msa);
  return NULL;
}

/* Function:  esl_msa_Digitize()
 * Synopsis:  Digitizes an msa, converting it from text mode.
 * Incept:    SRE, Sat Aug 26 17:33:08 2006 [AA 5302 to Dulles]
 *
 * Purpose:   Given an alignment <msa> in text mode, convert it to
 *            digital mode, using alphabet <abc>.
 *            
 *            Internally, the <ax> digital alignment field is filled,
 *            the <aseq> text alignment field is destroyed and free'd,
 *            a copy of the alphabet pointer is kept in the msa's
 *            <abc> reference, and the <eslMSA_DIGITAL> flag is raised
 *            in <flags>.
 *            
 *            Because <esl_msa_Digitize()> may be called on
 *            unvalidated user data, <errbuf> may be passed, for
 *            capturing an informative error message. For example, in
 *            reading alignments from files, invalid characters in the
 *            alignment are caught at the digitization step.
 *            
 * Args:      abc    - digital alphabet
 *            msa    - multiple alignment to digitize
 *            errbuf - optional: error message buffer, or <NULL>
 *
 * Returns:   <eslOK> on success;
 *            <eslEINVAL> if one or more sequences contain invalid characters
 *            that can't be digitized. If this happens, the <msa> is returned
 *            unaltered - left in text mode, with <aseq> as it was. (This is
 *            a normal error, because <msa->aseq> may be user input that we 
 *            haven't validated yet.)
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, state of <msa> may be 
 *            wedged, and it should only be destroyed, not used.
 */
int
esl_msa_Digitize(const ESL_ALPHABET *abc, ESL_MSA *msa, char *errbuf)
{
  char errbuf2[eslERRBUFSIZE];
  int  i;
  int  status;

  /* Contract checks */
  if (msa->aseq == NULL)           ESL_EXCEPTION(eslEINVAL, "msa has no text alignment");
  if (msa->ax   != NULL)           ESL_EXCEPTION(eslEINVAL, "msa already has digital alignment");
  if (msa->flags & eslMSA_DIGITAL) ESL_EXCEPTION(eslEINVAL, "msa is flagged as digital");

  /* Validate before we convert. Then we can leave the <aseq> untouched if
   * any of the sequences contain invalid characters.
   */
  for (i = 0; i < msa->nseq; i++)
    if (esl_abc_ValidateSeq(abc, msa->aseq[i], msa->alen, errbuf2) != eslOK) 
      ESL_FAIL(eslEINVAL, errbuf, "%s: %s", msa->sqname[i], errbuf2);

  /* Convert, sequence-by-sequence, free'ing aseq as we go.  */
  ESL_ALLOC(msa->ax, msa->sqalloc * sizeof(ESL_DSQ *));
  for (i = 0; i < msa->nseq; i++)
    {
      ESL_ALLOC(msa->ax[i], (msa->alen+2) * sizeof(ESL_DSQ));
      status = esl_abc_Digitize(abc, msa->aseq[i], msa->ax[i]);
      if (status != eslOK) goto ERROR;
      free(msa->aseq[i]);
    }    
  for (; i < msa->sqalloc; i++) 
    msa->ax[i] = NULL;
  free(msa->aseq);
  msa->aseq = NULL;

  msa->abc   =  (ESL_ALPHABET *) abc; /* convince compiler that removing const-ness is safe */
  msa->flags |= eslMSA_DIGITAL;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_msa_Textize()
 * Synopsis:  Convert a digital msa to text mode.
 * Incept:    SRE, Sat Aug 26 18:14:30 2006 [AA 5302 to Dulles]
 *
 * Purpose:   Given an alignment <msa> in digital mode, convert it
 *            to text mode.
 *            
 *            Internally, the <aseq> text alignment field is filled, the
 *            <ax> digital alignment field is destroyed and free'd, the
 *            msa's <abc> digital alphabet reference is nullified, and 
 *            the <eslMSA_DIGITAL> flag is dropped in <flags>.
 *            
 * Args:      msa   - multiple alignment to convert to text
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslECORRUPT> if one or more of the digitized alignment strings
 *            contain invalid characters.
 */
int
esl_msa_Textize(ESL_MSA *msa)
{
  int status;
  int i;

  /* Contract checks
   */
  if (msa->ax   == NULL)               ESL_EXCEPTION(eslEINVAL, "msa has no digital alignment");
  if (msa->aseq != NULL)               ESL_EXCEPTION(eslEINVAL, "msa already has text alignment");
  if (! (msa->flags & eslMSA_DIGITAL)) ESL_EXCEPTION(eslEINVAL, "msa is not flagged as digital");
  if (msa->abc  == NULL)               ESL_EXCEPTION(eslEINVAL, "msa has no digital alphabet");

  /* Convert, sequence-by-sequence, free'ing ax as we go.
   */
  ESL_ALLOC(msa->aseq, msa->sqalloc * sizeof(char *));
  for (i = 0; i < msa->nseq; i++)
    {
      ESL_ALLOC(msa->aseq[i], (msa->alen+1) * sizeof(char));
      status = esl_abc_Textize(msa->abc, msa->ax[i], msa->alen, msa->aseq[i]);
      if (status != eslOK) goto ERROR;
      free(msa->ax[i]);
    }
  for (; i < msa->sqalloc; i++)
    msa->aseq[i] = NULL;
  free(msa->ax);
  msa->ax = NULL;
  
  msa->abc    = NULL;      	 /* nullify reference (caller still owns real abc) */
  msa->flags &= ~eslMSA_DIGITAL; /* drop the flag */
  return eslOK;

 ERROR:
  return status;
}



/* Function:  esl_msafile_GuessAlphabet()
 * Synopsis:  Guess what kind of sequences the alignment file contains.
 * Incept:    SRE, Wed May 16 10:45:37 2007 [Janelia]
 *
 * Purpose:   Guess the alphabet of the sequences in the open
 *            <ESL_MSAFILE> -- <eslDNA>, <eslRNA>, or <eslAMINO> --
 *            based on the composition of the next MSA in the
 *            file. Usually this would be the first MSA, because we
 *            would call <esl_msafile_GuessAlphabet()> immediately
 *            after opening a new MSA file.
 *            
 * Returns:   Returns <eslOK> on success, and <*ret_type> is set
 *            to <eslDNA>, <eslRNA>, or <eslAMINO>. 
 *            
 *            Returns <eslEAMBIGUOUS> and sets <*ret_type> to
 *            <eslUNKNOWN> if the first alignment in the file contains
 *            no more than ten residues total, or if its alphabet
 *            cannot be reliably guessed (it contains IUPAC degeneracy
 *            codes, but no amino acid specific residues).
 * 
 *            Returns <eslEFORMAT> if a parse error is encountered
 *            in trying to read the alignment file. <msafp->errbuf>
 *            is set to a useful error message if this occurs; 
 *            <*ret_type> is set to <eslUNKNOWN>.
 *
 *            Returns <eslENODATA> if the file is empty and no
 *            alignment was found; <*ret_type> is set to <eslUNKNOWN>.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      J1/62.
 */
int
esl_msafile_GuessAlphabet(ESL_MSAFILE *msafp, int *ret_type)
{
  int      status;
  ESL_MSA *msa = NULL;
  
  /* If the msafp is digital mode, we already know the type (so why
   * are we being called?) If do_digital is TRUE, msafp->abc is
   * non-NULL: 
   */
  if (msafp->abc != NULL) { *ret_type = msafp->abc->type; return eslOK; } /* that was easy */

  /* If there's already an MSA cached, we've already called
   * GuessAlphabet(); don't read another one, or we'll overwrite the
   * first.
   */
  if (msafp->msa_cache != NULL) return esl_msa_GuessAlphabet(msafp->msa_cache, ret_type);

  /* Read first alignment, collect its residue composition for
   * passing off to esl_abc_GuessAlphabet()
   */
  status = esl_msa_Read(msafp, &msa);
  if      (status == eslEOF)     return eslENODATA;
  else if (status != eslOK)      return status;

  /* Store the msa in the MSAFILE's cache.  This is because if we're
   * reading from stdin or a gzip pipe, we'll have trouble rewinding
   * to prepare for the first msa_Read() call.
   */
  msafp->msa_cache = msa;

  /* And over to msa_GuessAlphabet() for the decision.
   */
  return esl_msa_GuessAlphabet(msa, ret_type);
}



/* Function:  esl_msafile_OpenDigital()
 * Synopsis:  Open an msa file for digital input.
 * Incept:    SRE, Sun Aug 27 17:40:33 2006 [Leesburg]
 *
 * Purpose:   Same as <esl_msafile_Open()>, except the alignment file
 *            will be read into a digitized internal representation,
 *            using internal alphabet <abc>, rather than the default
 *            internal ASCII text representation.
 *            
 * Args:      abc      - pointer to internal alphabet
 *            filename - name of alignment data file to open;
 *                       if "*.gz", attempt to read through <gunzip -dc> using <popen()>;
 *                       or "-" for stdin 
 *            format   - file format code (e.g. <eslMSAFILE_STOCKHOLM>);
 *                       or <eslMSAFILE_UNKNOWN> to invoke format autodetection.
 *            env      - NULL, or the name of an environment variable from which
 *                       to retrieve a colon-delimited directory list to search
 *                       for <filename> in. (e.g. "HMMERDB")
 *            ret_msafp - RETURN: open MSAFILE.
 *
 * Returns:  <eslOK> on success, and <ret_msafp> is set to point at
 *           an open <ESL_MSAFILE>. Caller frees this file pointer with
 *           <esl_msafile_Close()>.
 *           
 *           <eslENOTFOUND> if <filename> cannot be opened;
 *           <eslEFORMAT> if autodetection is attempted and format
 *           cannot be determined.
 *           
 * Throws:   <eslEMEM> on allocation failure.
 *           <eslEINVAL> if format autodetection is attempted on 
 *           stdin or a gunzip pipe.
 */
int
esl_msafile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, 
			int format, const char *env, ESL_MSAFILE **ret_msafp)
{
  ESL_MSAFILE *msafp;
  int          status;

  if ((status = msafile_open(filename, format, env, &msafp)) != eslOK) return status;
  esl_msafile_SetDigital(msafp, abc);

  *ret_msafp = msafp;
  return eslOK;
}



/* Function:  esl_msafile_SetDigital()
 * Synopsis:  Set an open <ESL_MSAFILE> to read in digital mode.
 * Incept:    SRE, Wed May 16 10:40:24 2007 [Janelia]
 *
 * Purpose:   Given an open <ESL_MSAFILE>, set it so that all subsequent
 *            calls to <esl_msa_Read()> will read multiple alignments
 *            in digital mode instead of text mode, using alphabet
 *            <abc>.
 *            
 *            <esl_msafile_Open(); esl_msafile_SetDigital()> is
 *            equivalent to <esl_msafile_OpenDigital()>. The two-step
 *            version is useful if you don't already know the alphabet
 *            type for your msa file, and you need to call
 *            <esl_msafile_GuessAlphabet()> after opening the file but
 *            before setting its digital alphabet.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_msafile_SetDigital(ESL_MSAFILE *msafp, const ESL_ALPHABET *abc)
{
  msafp->abc        = abc;
  msafp->do_digital = TRUE;
  return eslOK;
}
#endif /* eslAUGMENT_ALPHABET */
/*---------------------- end of digital MSA functions -----------------------*/





/******************************************************************************
 *# 4. Random MSA database access (augmentation: ssi)
 *****************************************************************************/
#ifdef eslAUGMENT_SSI
/* Function:  esl_msafile_PositionByKey()
 * Synopsis:  Use SSI to reposition file to start of named MSA.
 * Incept:    SRE, Mon May 28 11:04:59 2007 [Janelia]
 *
 * Purpose:   Reposition <afp> so that the next MSA we read
 *            will be the one named (or accessioned) <key>.
 *
 * Returns:   <eslOK> on success, and the file <afp> is repositioned
 *            such that the next <esl_msafile_Read()> call will read the
 *            alignment named <key>.
 *            
 *            Returns <eslENOTFOUND> if <key> isn't found in the index
 *            for <afp>. 
 *            
 *            Returns <eslEFORMAT> if something goes wrong trying to
 *            read the index, indicating some sort of file format
 *            problem in the SSI file.
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINVAL> if there's no open SSI index;
 *            <eslESYS> if an <fseek()> fails.
 *            In all these cases, the state of the <afp> is uncertain
 *            and may be corrupt; the application should not continue
 *            to use it.
 */
int
esl_msafile_PositionByKey(ESL_MSAFILE *afp, const char *key)
{
  uint16_t fh;
  off_t    offset;
  int      status;

  if (afp->ssi == NULL) ESL_EXCEPTION(eslEINVAL, "Need an open SSI index to call esl_msafile_PositionByKey()");
  if ((status = esl_ssi_FindName(afp->ssi, key, &fh, &offset, NULL, NULL)) != eslOK) return status;
  if (fseeko(afp->f, offset, SEEK_SET) != 0)    ESL_EXCEPTION(eslESYS, "fseek failed");
  
  /* If the <afp> had an MSA cached, we will probably have to discard
   * it, unless by chance it's exactly the MSA we're looking for.
   */
  if (afp->msa_cache != NULL)
    {
      if ( (afp->msa_cache->name == NULL || strcmp(afp->msa_cache->name, key) != 0) &&
	   (afp->msa_cache->acc  == NULL || strcmp(afp->msa_cache->acc,  key) != 0))
	{
	  esl_msa_Destroy(afp->msa_cache);
	  afp->msa_cache = NULL;
	}
    }

  /* The linenumber gets messed up after a file positioning. Best we can do
   * is to reset it to zero.
   */
  afp->linenumber = 0; 
  return eslOK;
}
#endif /*eslAUGMENT_SSI*/
/*------------- end of functions added by SSI augmentation -------------------*/




/******************************************************************************
 *# 5. General i/o API, all alignment formats                                 
 *****************************************************************************/
static int write_stockholm(FILE *fp, const ESL_MSA *msa);
static int write_pfam     (FILE *fp, const ESL_MSA *msa);
static int write_a2m      (FILE *fp,       ESL_MSA *msa);
static int write_psiblast (FILE *fp,       ESL_MSA *msa);
static int write_afa      (FILE *fp,       ESL_MSA *msa);

static int read_stockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa);
static int read_selex    (ESL_MSAFILE *afp, ESL_MSA **ret_msa);
static int read_afa      (ESL_MSAFILE *afp, ESL_MSA **ret_msa);

/* Function:  esl_msa_Read()
 * Synopsis:  Read next MSA from a file.
 * Incept:    SRE, Fri Jan 28 08:10:49 2005 [St. Louis]
 *
 * Purpose:   Reads the next MSA from an open MSA file <afp>,
 *            and returns it via <ret_msa>. 
 *
 * Returns:   <eslOK> on success, and <ret_msa> points at the
 *            new MSA object.
 *            
 *            Returns <eslEOF> if there are no more alignments in the file.
 *            
 *            Returns <eslEFORMAT> if there is a parse error, and <afp->errbuf>
 *            is set to an informative message.
 *            
 *            <eslEINVAL> if we're trying to read a digital alignment,
 *            but one or more residues are seen in the file that
 *            aren't valid in our alphabet.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Read(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA *msa;
  int      status;

  *ret_msa = NULL;

  /* If we've just used GuessAlphabet(), we have an MSA already read
   * and stored in the MSAFILE's cache. Just return it, after worrying
   * about whether it's supposed to be in digital or text mode. (It
   * should always be in text mode, and maybe in need of Digitize(),
   * given how GuessAlphabet works now; but this is coded for more
   * generality in case we use the MSA cache some other way in the
   * future.)
   */
  if (afp->msa_cache != NULL) 
    {
#ifdef eslAUGMENT_ALPHABET
      if      (afp->do_digital   && !(afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Digitize(afp->abc, afp->msa_cache, afp->errbuf)) != eslOK) return status; 
      }
      else if (! afp->do_digital && (afp->msa_cache->flags & eslMSA_DIGITAL)) {
	if ((status = esl_msa_Textize(afp->msa_cache)) != eslOK) return status;
      }
#endif

      *ret_msa         = afp->msa_cache;
      afp->msa_cache = NULL;
      return eslOK;
    }

  /* Otherwise, read the next MSA from the file.
   */      
  switch (afp->format) {
  case eslMSAFILE_STOCKHOLM: status = read_stockholm(afp, &msa); break;
  case eslMSAFILE_PFAM:      status = read_stockholm(afp, &msa); break;
  case eslMSAFILE_A2M:       ESL_FAIL(eslEFORMAT, afp->errbuf, "A2M format input parser not implemented yet.");
  case eslMSAFILE_PSIBLAST:  ESL_FAIL(eslEFORMAT, afp->errbuf, "PSIBLAST format input parser not implemented yet.");
  case eslMSAFILE_SELEX:     status = read_selex    (afp, &msa); break;
  case eslMSAFILE_AFA:       status = read_afa      (afp, &msa); break;
  default:                   ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
  }

  *ret_msa = msa;
  return status;
}

/* Function:  esl_msa_Write()
 * Synopsis:  Write an MSA to a file.
 * Incept:    SRE, Fri Jan 28 09:29:28 2005 [St. Louis]
 *
 * Purpose:   Writes an alignment <msa> to an open stream <fp>,
 *            in format specified by <fmt>.
 *            
 *            In general, the <msa> is unchanged; however, in certain
 *            formats and under certain conditions, modifications may
 *            be made. For example, writing an alignment in A2M format
 *            will alter the alignment data (marking missing data
 *            symbols on heuristically defined sequence fragments) and
 *            create an <\#=RF> annotation line, if an <msa->rf>
 *            annotation line isn't already present in the <msa>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal error.
 */
int
esl_msa_Write(FILE *fp, ESL_MSA *msa, int fmt)
{
  int status;
  switch (fmt) {
  case eslMSAFILE_STOCKHOLM: status = write_stockholm(fp, msa); break;
  case eslMSAFILE_PFAM:      status = write_pfam(fp, msa);      break;
  case eslMSAFILE_A2M:       status = write_a2m(fp, msa);       break;
  case eslMSAFILE_PSIBLAST:  status = write_psiblast(fp, msa);  break;
  case eslMSAFILE_SELEX:     ESL_EXCEPTION(eslEUNIMPLEMENTED, "selex format writing isn't implemented yet");
  case eslMSAFILE_AFA:       status = write_afa(fp, msa);       break;
  default: ESL_EXCEPTION(eslEINCONCEIVABLE, "no such format");
  } 
  return status;
}


/* Function:  esl_msa_EncodeFormat()
 * Synopsis:  Convert text string to an MSA file format code.
 * Incept:    SRE, Fri Oct 24 13:21:08 2008 [Janelia]
 *
 * Purpose:   Given a text string, match it case-insensitively
 *            against a list of possible formats, and return the
 *            appropriate MSA file format code. For example,
 *            <esl_msa_EncodeFormat("Stockholm")> returns
 *            <eslMSAFILE_STOCKHOLM>.
 *            
 *            If the format is unrecognized, return
 *            <eslMSAFILE_UNKNOWN>.
 *            
 * Note:      Keep in sync with <esl_sqio_EncodeFormat()>, 
 *            which decodes all possible sequence file formats,
 *            both unaligned and aligned.           
 */
int
esl_msa_EncodeFormat(char *fmtstring)
{
  if (strcasecmp(fmtstring, "stockholm") == 0) return eslMSAFILE_STOCKHOLM;
  if (strcasecmp(fmtstring, "pfam")      == 0) return eslMSAFILE_PFAM;
  if (strcasecmp(fmtstring, "a2m")       == 0) return eslMSAFILE_A2M;
  if (strcasecmp(fmtstring, "psiblast")  == 0) return eslMSAFILE_PSIBLAST;
  if (strcasecmp(fmtstring, "selex")     == 0) return eslMSAFILE_SELEX;
  if (strcasecmp(fmtstring, "afa")       == 0) return eslMSAFILE_AFA;
  return eslMSAFILE_UNKNOWN;
}


/* Function:  esl_msa_DecodeFormat()
 * Synopsis:  Convert internal file format code to text string.
 * Incept:    SRE, Fri May 18 11:59:58 2007 [Janelia]
 *
 * Purpose:   Given an internal file format code <fmt> 
 *            (<eslMSAFILE_STOCKHOLM>, for example), returns
 *            a string suitable for printing ("Stockholm",
 *            for example).
 *            
 * Returns:   a pointer to a static description string.
 * 
 * Throws:    If code isn't valid, throws an <eslEINVAL> exception 
 *            internally, and returns <NULL>.
 *            
 * Note:      Keep in sync with <esl_sqio_DecodeFormat()>.
 */
char *
esl_msa_DecodeFormat(int fmt)
{
  switch (fmt) {
  case eslMSAFILE_UNKNOWN:   return "unknown";
  case eslMSAFILE_STOCKHOLM: return "Stockholm";
  case eslMSAFILE_PFAM:      return "Pfam";
  case eslMSAFILE_A2M:       return "UCSC A2M";
  case eslMSAFILE_PSIBLAST:  return "PSI-BLAST";
  case eslMSAFILE_SELEX:     return "SELEX";
  case eslMSAFILE_AFA:       return "aligned FASTA";
  default:                   break;
  }
  esl_exception(eslEINVAL, __FILE__, __LINE__, "no such msa format code %d\n", fmt);
  return NULL;
}


/* Function:  esl_msa_GuessFileFormat()
 * Synopsis:  Determine the format of an open MSA file.
 * Incept:    SRE, Fri Jan 28 07:29:00 2005 [St. Louis]
 *
 * Purpose:   Attempts to determine the format of an open alignment file
 *            <afp>, for which <afp->format> is <eslMSAFILE_UNKNOWN>. 
 *            If successful, sets <afp->format>.
 *            
 *            Currently a placeholder: it always guesses Stockholm!
 *
 * Returns:   <eslOK> on success, and sets <afp->format>. 
 *            <eslEFORMAT> if format can't be determined.
 *
 * Xref:      squid's MSAFileFormat()
 */
int
esl_msa_GuessFileFormat(ESL_MSAFILE *afp)
{
  /* Placeholder: FIXME: autodetection code goes here.
   */
  afp->format = eslMSAFILE_STOCKHOLM;
  return eslOK;
}
/*-------------------- end of general i/o functions -------------------------*/




/*****************************************************************
 *# 6. Miscellaneous functions for manipulating MSAs
 *****************************************************************/

/* Function:  esl_msa_ReasonableRF()
 * Synopsis:  Determine a reasonable #=RF line marking "consensus" columns.
 * Incept:    SRE, Wed Sep  3 10:42:05 2008 [Janelia]
 *
 * Purpose:   Define an <rfline> for the multiple alignment <msa> that
 *            marks consensus columns with an 'x', and non-consensus 
 *            columns with a '.'.
 *            
 *            Consensus columns are defined as columns with fractional
 *            occupancy of $\geq$ <symfrac> in residues. For example,
 *            if <symfrac> is 0.7, columns containing $\geq$ 70\%
 *            residues are assigned as 'x' in the <rfline>, roughly
 *            speaking. "Roughly speaking", because the fractional
 *            occupancy is in fact calculated as a weighted frequency
 *            using sequence weights in <msa->wgt>, and because
 *            missing data symbols are ignored in order to be able to
 *            deal with sequence fragments. 
 *            
 *            The greater <symfrac> is, the more stringent the
 *            definition, and the fewer columns will be defined as
 *            consensus. <symfrac=0> will define all columns as
 *            consensus. <symfrac=1> will only define a column as
 *            consensus if it contains no gap characters at all.
 *            
 *            If the caller wants to designate any sequences as
 *            fragments, it must convert all leading and trailing gaps
 *            to the missing data symbol '~'.
 *
 *            For text mode alignments, any alphanumeric character is
 *            considered to be a residue, and any non-alphanumeric
 *            character is considered to be a gap.
 *            
 *            The <rfline> is a NUL-terminated string, indexed
 *            <0..alen-1>.
 *
 *            The <rfline> result can be <msa->rf>, if the caller
 *            wants to set the <msa's> own RF line; or it can be any
 *            alternative storage provided by the caller. In either
 *            case, the caller must provide allocated space for at
 *            least <msa->alen+1> chars.
 *            
 * Args:      msa      - MSA to define a consensus RF line for
 *            symfrac  - threshold for defining consensus columns
 *            rfline   - RESULT: string containing x for consensus, . for not
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      HMMER p7_Fastmodelmaker() uses an essentially identical
 *            calculation to define model architecture, and could be
 *            rewritten now to use this function. 
 *            
 *            A2M format alignment output uses this to define
 *            consensus columns when #=RF annotation isn't available.
 */
int
esl_msa_ReasonableRF(ESL_MSA *msa, double symfrac, char *rfline)
{
  int    apos;
  int    idx;
  double r;
  double totwgt;
  
#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      for (apos = 1; apos <= msa->alen; apos++) 
	{  
	  r = totwgt = 0.;
	  for (idx = 0; idx < msa->nseq; idx++) 
	    {
	      if       (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
	      else if  (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                     totwgt += msa->wgt[idx]; }
	      else if  (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos])) continue;
	    }
	  if (r > 0. && r / totwgt >= symfrac) msa->rf[apos-1] = 'x';
	  else                                 msa->rf[apos-1] = '.';
	}
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (apos = 0; apos < msa->alen; apos++) 
	{  
	  r = totwgt = 0.;
	  for (idx = 0; idx < msa->nseq; idx++) 
	    {
	      if    (isalpha(msa->aseq[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
	      else                                                        totwgt += msa->wgt[idx];
	    }
	  if (r > 0. && r / totwgt >= symfrac) msa->rf[apos] = 'x';
	  else                                 msa->rf[apos] = '.';
	}
    }

  msa->rf[msa->alen] = '\0';
  return eslOK;
}


/* Function:  esl_msa_MarkFragments()
 * Synopsis:  Heuristically define seq fragments in an alignment.
 * Incept:    SRE, Wed Sep  3 11:49:25 2008 [Janelia]
 *
 * Purpose:   Use a heuristic to define sequence fragments (as opposed
 *            to "full length" sequences in alignment <msa>.
 *            
 *            The rule is that if the sequence has a raw (unaligned)
 *            length of less than <fragthresh> times the alignment
 *            length in columns, the sequence is defined as a fragment.
 *            
 *            For each fragment, all leading and trailing gap symbols
 *            (all gaps before the first residue and after the last
 *            residue) are converted to missing data symbols
 *            (typically '~', but nonstandard digital alphabets may
 *            have defined another character).
 *            
 *            As a special case, if <fragthresh> is negative, then all
 *            sequences are defined as fragments.
 *
 * Args:      msa        - alignment in which to define and mark seq fragments 
 *            fragthresh - define frags if rlen < fragthresh * alen;
 *                         or if fragthresh < 0, all seqs are marked as frags.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_msa_MarkFragments(ESL_MSA *msa, double fragthresh)
{
  int    i;
  int    pos;

  for (i = 0; i < msa->nseq; i++)
    if (fragthresh < 0.0 || msa_get_rlen(msa, i) < fragthresh * msa->alen)
      {  
#ifdef eslAUGMENT_ALPHABET
	if (msa->flags & eslMSA_DIGITAL) {
	  for (pos = 1; pos <= msa->alen; pos++) {
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[i][pos])) break;
	    msa->ax[i][pos] = esl_abc_XGetMissing(msa->abc);
	  }
	  for (pos = msa->alen; pos >= 1; pos--) {	  
	    if (esl_abc_XIsResidue(msa->abc, msa->ax[i][pos])) break;
	    msa->ax[i][pos] = esl_abc_XGetMissing(msa->abc);
	  }
	}
#endif
	if (! (msa->flags & eslMSA_DIGITAL)) 
	  {
	    for (pos = 0; pos < msa->alen; pos++) {
	      if (isalnum(msa->aseq[i][pos])) break;
	      msa->aseq[i][pos] = '~';
	    }
	    for (pos = msa->alen-1; pos >= 0; pos--) {	  
	      if (isalnum(msa->aseq[i][pos])) break;
	      msa->aseq[i][pos] = '~';
	    }
	  }
      }
  return eslOK;
}


/* Function:  esl_msa_SequenceSubset()
 * Synopsis:  Select subset of sequences into a smaller MSA.
 * Incept:    SRE, Wed Apr 13 10:05:44 2005 [St. Louis]
 *
 * Purpose:   Given an array <useme> (0..nseq-1) of TRUE/FALSE flags for each
 *            sequence in an alignment <msa>; create a new alignment containing
 *            only those seqs which are flagged <useme=TRUE>. Return a pointer
 *            to this newly allocated alignment through <ret_new>. Caller is
 *            responsible for freeing it.
 *            
 *            The smaller alignment might now contain columns
 *            consisting entirely of gaps or missing data, depending
 *            on what sequence subset was extracted. The caller may
 *            want to immediately call <esl_msa_MinimGaps()> on the
 *            new alignment to clean this up.
 *
 *            Unparsed GS and GR Stockholm annotation that is presumably still
 *            valid is transferred to the new alignment. Unparsed GC, GF, and
 *            comments that are potentially invalidated by taking the subset
 *            of sequences are not transferred to the new MSA.
 *            
 *            Weights are transferred exactly. If they need to be
 *            renormalized to some new total weight (such as the new,
 *            smaller total sequence number), the caller must do that.
 *            
 *            <msa> may be in text mode or digital mode. The new MSA
 *            in <ret_new> will have the same mode.
 *
 * Returns:   <eslOK> on success, and <ret_new> is set to point at a new
 *            (smaller) alignment.
 *
 * Throws:    <eslEINVAL> if the subset has no sequences in it;
 *            <eslEMEM> on allocation error.
 *
 * Xref:      squid's MSASmallerAlignment(), 1999.
 */
int
esl_msa_SequenceSubset(const ESL_MSA *msa, const int *useme, ESL_MSA **ret_new)
{
  ESL_MSA *new = NULL;
  int  nnew;			/* number of seqs in the new MSA */
  int  oidx, nidx;		/* old, new indices */
  int  i;
  int  status;
  
  *ret_new = NULL;

  nnew = 0; 
  for (oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx]) nnew++;
  if (nnew == 0) ESL_EXCEPTION(eslEINVAL, "No sequences selected");

  /* Note that the Create() calls allocate exact space for the sequences,
   * so we will strcpy()/memcpy() into them below.
   */
#ifdef eslAUGMENT_ALPHABET
  if ((msa->flags & eslMSA_DIGITAL) &&
      (new = esl_msa_CreateDigital(msa->abc, nnew, msa->alen)) == NULL)
    {status = eslEMEM; goto ERROR; }
#endif
  if (! (msa->flags & eslMSA_DIGITAL) &&
      (new = esl_msa_Create(nnew, msa->alen)) == NULL) 
    {status = eslEMEM; goto ERROR; }
  if (new == NULL) 
    {status = eslEMEM; goto ERROR; }
  

  /* Copy the old to the new */
  for (nidx = 0, oidx = 0; oidx < msa->nseq; oidx++)
    if (useme[oidx])
      {
#ifdef eslAUGMENT_ALPHABET
	if (msa->flags & eslMSA_DIGITAL)
	  memcpy(new->ax[nidx], msa->ax[oidx], sizeof(ESL_DSQ) * (msa->alen+2));
#endif
	if (! (msa->flags & eslMSA_DIGITAL))
	  strcpy(new->aseq[nidx], msa->aseq[oidx]);
	if ((status = esl_strdup(msa->sqname[oidx], -1, &(new->sqname[nidx])))    != eslOK) goto ERROR;

	new->wgt[nidx] = msa->wgt[oidx];
      
	if (msa->sqacc != NULL && msa->sqacc[oidx] != NULL) {
	  if ((status = esl_msa_SetSeqAccession(new, nidx, msa->sqacc[oidx])) != eslOK) goto ERROR;
	}
	if (msa->sqdesc != NULL && msa->sqdesc[oidx] != NULL) {
	  if ((status = esl_msa_SetSeqDescription(new, nidx, msa->sqdesc[oidx])) != eslOK) goto ERROR;
	}
	if (msa->ss != NULL && msa->ss[oidx] != NULL) {
	  if ((status = set_seq_ss(new, nidx, msa->ss[oidx])) != eslOK) goto ERROR;
	}
	if (msa->sa != NULL && msa->sa[oidx] != NULL) {
	  if ((status = set_seq_sa(new, nidx, msa->sa[oidx])) != eslOK) goto ERROR;
	}
	if (msa->pp != NULL && msa->pp[oidx] != NULL) {
	  if ((status = set_seq_pp(new, nidx, msa->pp[oidx])) != eslOK) goto ERROR;
	}
	/* unparsed annotation */
	for(i = 0; i < msa->ngs; i++) {
	  if(msa->gs[i] != NULL) 
	    if ((status = esl_msa_AddGS(new, msa->gs_tag[i], nidx, msa->gs[i][oidx])) != eslOK) goto ERROR;
	}
	for(i = 0; i < msa->ngr; i++) {
	  if(msa->gr[i] != NULL) 
	    if ((status = esl_msa_AppendGR(new, msa->gr_tag[i], nidx, msa->gr[i][oidx])) != eslOK) goto ERROR;
	}

	nidx++;
      }

  new->flags = msa->flags;

  if ((status = esl_strdup(msa->name,           -1, &(new->name)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->desc,           -1, &(new->desc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->acc,            -1, &(new->acc)))     != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->au,             -1, &(new->au)))      != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->ss_cons, msa->alen, &(new->ss_cons))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->sa_cons, msa->alen, &(new->sa_cons))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->pp_cons, msa->alen, &(new->pp_cons))) != eslOK) goto ERROR;
  if ((status = esl_strdup(msa->rf,      msa->alen, &(new->rf)))      != eslOK) goto ERROR;

  for (i = 0; i < eslMSA_NCUTS; i++) {
    new->cutoff[i] = msa->cutoff[i];
    new->cutset[i] = msa->cutset[i];
  }
  
  new->nseq  = nnew;
  new->sqalloc = nnew;

  /* Since we have a fully constructed MSA, we don't need the
   * aux info used by parsers.
   */
  if (new->sqlen != NULL) { free(new->sqlen);  new->sqlen = NULL; }
  if (new->sslen != NULL) { free(new->sslen);  new->sslen = NULL; }
  if (new->salen != NULL) { free(new->salen);  new->salen = NULL; }
  if (new->pplen != NULL) { free(new->pplen);  new->pplen = NULL; }
  new->lastidx = -1;

  *ret_new = new;
  return eslOK;

 ERROR:
  if (new != NULL) esl_msa_Destroy(new);
  *ret_new = NULL;
  return status;
}

/* remove_broken_basepairs_from_ss_string()
 * 
 * Given an array <useme> (0..alen-1) of TRUE/FALSE flags, remove
 * any basepair from an SS string that is between alignment
 * columns (i,j) for which either <useme[i-1]> or <useme[j-1]> is FALSE.
 * Helper function for remove_broken_basepairs_from_msa(). 
 * 
 * The input SS string will be overwritten. If it was not in 
 * full WUSS format when pass in, it will be upon exit. 
 * Note that that means if there's residues in the input ss
 * that correspond to gaps in an aligned sequence or RF 
 * annotation, they will not be treated as gaps in the 
 * returned SS. For example, a gap may become a '-' character,
 * a '_' character, or a ':' character. I'm not sure how
 * to deal with this in a better way. We could demand an
 * aligned sequence to use to de-gap the SS string, but 
 * that would require disallowing any gap to be involved
 * in a basepair, which I'm not sure is something we want
 * to forbid.
 * 
 * If the original SS is inconsistent it's left untouched and
 * non-eslOK is returned as listed below.
 *
 * Returns:   <eslOK> on success.
 *            <eslESYNTAX> if SS string 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted 
 *            SS, this should be impossible.
 */
static int
remove_broken_basepairs_from_ss_string(char *ss, char *errbuf, int len, const int *useme)
{
  int64_t  apos;                 /* alignment position */
  int     *ct = NULL;	         /* 0..alen-1 base pair partners array for current sequence */
  char    *ss_nopseudo = NULL;   /* no-pseudoknot version of structure */
  int      status;

  ESL_ALLOC(ct,          sizeof(int)  * (len+1));
  ESL_ALLOC(ss_nopseudo, sizeof(char) * (len+1));

  esl_wuss_nopseudo(ss, ss_nopseudo);
  if ((status = esl_wuss2ct(ss_nopseudo, len, ct)) != eslOK) 
    ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");
  for (apos = 1; apos <= len; apos++) { 
    if (!(useme[apos-1])) { 
      if (ct[apos] != 0) ct[ct[apos]] = 0;
      ct[apos] = 0;
    }
  }
  /* All broken bps removed from ct, convert to WUSS SS string and overwrite SS */
  if ((status = esl_ct2wuss(ct, len, ss)) != eslOK) 
    ESL_FAIL(status, errbuf, "Error converting de-knotted bp ct array to WUSS notation.");
  
  free(ss_nopseudo);
  free(ct);
  return eslOK;

 ERROR: 
  if (ct          != NULL) free(ct);
  if (ss_nopseudo != NULL) free(ss_nopseudo);
  return status; 
}  

/* remove_broken_basepairs_from_msa()
 * 
 * Given an array <useme> (0..alen-1) of TRUE/FALSE flags, remove
 * any basepair from SS_cons and individual SS annotation in alignment
 * columns (i,j) for which either <useme[i-1]> or <useme[j-1]> is FALSE.
 * Called automatically from esl_msa_ColumnSubset() with same <useme>. 
 * 
 * If the original structure data is inconsistent it's left untouched.
 *
 * Returns:   <eslOK> on success.
 *            <eslESYNTAX> if WUSS string for SS_cons or msa->ss 
 *            following esl_wuss_nopseudo() is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted 
 *            SS, this should be impossible
 */
static int
remove_broken_basepairs_from_msa(ESL_MSA *msa, char *errbuf, const int *useme)
{
  int status;
  int  i;

  if (msa->ss_cons != NULL) { 
    if((status = remove_broken_basepairs_from_ss_string(msa->ss_cons, errbuf, msa->alen, useme)) != eslOK) return status; 
  }
  /* per-seq SS annotation */
  if (msa->ss != NULL) { 
    for(i = 0; i < msa->nseq; i++) { 
      if (msa->ss[i] != NULL) { 
	if((status = remove_broken_basepairs_from_ss_string(msa->ss[i], errbuf, msa->alen, useme)) != eslOK) return status; 
      }
    }
  }
  return eslOK;
}  

/* Function:  esl_msa_ColumnSubset()
 * Synopsis:  Remove a selected subset of columns from the MSA
 *
 * Incept:    SRE, Sun Feb 27 10:05:07 2005
 *            From squid's MSAShorterAlignment(), 1999
 * 
 * Purpose:   Given an array <useme> (0..alen-1) of TRUE/FALSE flags,
 *            where TRUE means "keep this column in the new alignment"; 
 *            remove all columns annotated as FALSE in the <useme> 
 *            array. This is done in-place on the MSA, so the MSA is 
 *            modified: <msa->alen> is reduced, <msa->aseq> is shrunk 
 *            (or <msa->ax>, in the case of a digital mode alignment), 
 *            and all associated per-residue or per-column annotation
 *            is shrunk.
 * 
 * Returns:   <eslOK> on success.
 *            Possibilities from <remove_broken_basepairs_from_msa()> call:
 *            <eslESYNTAX> if WUSS string for <SS_cons> or <msa->ss>
 *            following <esl_wuss_nopseudo()> is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 */
int
esl_msa_ColumnSubset(ESL_MSA *msa, char *errbuf, const int *useme)
{
  int     status;
  int64_t opos;			/* position in original alignment */
  int64_t npos;			/* position in new alignment      */
  int     idx;			/* sequence index */
  int     i;			/* markup index */

  /* Remove any basepairs from SS_cons and individual sequence SS
   * for aln columns i,j for which useme[i-1] or useme[j-1] are FALSE 
   */
  if((status = remove_broken_basepairs_from_msa(msa, errbuf, useme)) != eslOK) return status;

  /* Since we're minimizing, we can overwrite in place, within the msa
   * we've already got. 
   * opos runs all the way to msa->alen to include (and move) the \0
   * string terminators (or sentinel bytes, in the case of digital mode)
   */
  for (opos = 0, npos = 0; opos <= msa->alen; opos++)
    {
      if (opos < msa->alen && useme[opos] == FALSE) continue;

      if (npos != opos)	/* small optimization */
	{
	  /* The alignment, and per-residue annotations */
	  for (idx = 0; idx < msa->nseq; idx++)
	    {
#ifdef eslAUGMENT_ALPHABET
	      if (msa->flags & eslMSA_DIGITAL) /* watch off-by-one in dsq indexing */
		msa->ax[idx][npos+1] = msa->ax[idx][opos+1];
	      else
		msa->aseq[idx][npos] = msa->aseq[idx][opos];
#else
	      msa->aseq[idx][npos] = msa->aseq[idx][opos];
#endif /*eslAUGMENT_ALPHABET*/
	      if (msa->ss != NULL && msa->ss[idx] != NULL) msa->ss[idx][npos] = msa->ss[idx][opos];
	      if (msa->sa != NULL && msa->sa[idx] != NULL) msa->sa[idx][npos] = msa->sa[idx][opos];
	      if (msa->pp != NULL && msa->pp[idx] != NULL) msa->pp[idx][npos] = msa->pp[idx][opos];
	      for (i = 0; i < msa->ngr; i++)
		if (msa->gr[i][idx] != NULL)
		  msa->gr[i][idx][npos] = msa->gr[i][idx][opos];
	    }	  
	  /* The per-column annotations */
	  if (msa->ss_cons != NULL) msa->ss_cons[npos] = msa->ss_cons[opos];
	  if (msa->sa_cons != NULL) msa->sa_cons[npos] = msa->sa_cons[opos];
	  if (msa->pp_cons != NULL) msa->pp_cons[npos] = msa->pp_cons[opos];
	  if (msa->rf      != NULL) msa->rf[npos]      = msa->rf[opos];
	  for (i = 0; i < msa->ngc; i++)
	    msa->gc[i][npos] = msa->gc[i][opos];
	}
      npos++;
    }
  msa->alen = npos-1;	/* -1 because npos includes NUL terminators */
  return eslOK;
}

/* Function:  esl_msa_MinimGaps()
 * Synopsis:  Remove columns containing all gap symbols.
 * Incept:    SRE, Sun Feb 27 11:03:42 2005 [St. Louis]
 *
 * Purpose:   Remove all columns in the multiple alignment <msa>
 *            that consist entirely of gaps or missing data.
 *            
 *            For a text mode alignment, <gaps> is a string defining
 *            the gap characters, such as <"-_.">. For a digital mode
 *            alignment, <gaps> may be passed as <NULL>, because the
 *            internal alphabet already knows what the gap and missing
 *            data characters are.
 *            
 *            <msa> is changed in-place to a narrower alignment
 *            containing fewer columns. All per-residue and per-column
 *            annotation is altered appropriately for the columns that
 *            remain in the new alignment.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation failure.
 *            Possibilities from <esl_msa_ColumnSubset()> call:
 *            <eslESYNTAX> if WUSS string for <SS_cons> or <msa->ss>
 *            following <esl_wuss_nopseudo()> is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 *
 * Xref:      squid's MSAMingap().
 */
int
esl_msa_MinimGaps(ESL_MSA *msa, char *errbuf, const char *gaps)
{
  int    *useme = NULL;	/* array of TRUE/FALSE flags for which cols to keep */
  int64_t apos;		/* column index   */
  int     idx;		/* sequence index */
  int     status;

  ESL_ALLOC(useme, sizeof(int) * (msa->alen+1)); /* +1 is just to deal w/ alen=0 special case */

#ifdef eslAUGMENT_ALPHABET	   /* digital mode case */
  if (msa->flags & eslMSA_DIGITAL) /* be careful of off-by-one: useme is 0..L-1 indexed */
    {
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (! esl_abc_XIsGap    (msa->abc, msa->ax[idx][apos]) &&
		! esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      break;
	  if (idx == msa->nseq) useme[apos-1] = FALSE; else useme[apos-1] = TRUE;
	}
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL)) /* text mode case */
    {
      for (apos = 0; apos < msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (strchr(gaps, msa->aseq[idx][apos]) == NULL)
	      break;
	  if (idx == msa->nseq) useme[apos] = FALSE; else useme[apos] = TRUE;
	}
    }

  if((status = esl_msa_ColumnSubset(msa, errbuf, useme)) != eslOK) return status;
  free(useme);
  return eslOK;

 ERROR:
  if (useme != NULL) free(useme);
  return status;
}

/* Function:  esl_msa_NoGaps()
 * Synopsis:  Remove columns containing any gap symbol.
 * Incept:    SRE, Sun Feb 27 10:17:58 2005 [St. Louis]
 *
 * Purpose:   Remove all columns in the multiple alignment <msa> that
 *            contain any gaps or missing data, such that the modified
 *            MSA consists only of ungapped columns (a solid block of
 *            residues). 
 *            
 *            This is useful for filtering alignments prior to
 *            phylogenetic analysis using programs that can't deal
 *            with gaps.
 *            
 *            For a text mode alignment, <gaps> is a string defining
 *            the gap characters, such as <"-_.">. For a digital mode
 *            alignment, <gaps> may be passed as <NULL>, because the
 *            internal alphabet already knows what the gap and
 *            missing data characters are.
 *    
 *            <msa> is changed in-place to a narrower alignment
 *            containing fewer columns. All per-residue and per-column
 *            annotation is altered appropriately for the columns that
 *            remain in the new alignment.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            Possibilities from <esl_msa_ColumnSubset()> call:
 *            <eslESYNTAX> if WUSS string for <SS_cons> or <msa->ss>
 *            following <esl_wuss_nopseudo()> is inconsistent.
 *            <eslEINVAL> if a derived ct array implies a pknotted SS.
 *
 * Xref:      squid's MSANogap().
 */
int
esl_msa_NoGaps(ESL_MSA *msa, char *errbuf, const char *gaps)
{
  int    *useme = NULL;	/* array of TRUE/FALSE flags for which cols to keep */
  int64_t apos;		/* column index */
  int     idx;		/* sequence index */
  int     status;

  ESL_ALLOC(useme, sizeof(int) * (msa->alen+1)); /* +1 is only to deal with alen=0 special case */

#ifdef eslAUGMENT_ALPHABET	   /* digital mode case */
  if (msa->flags & eslMSA_DIGITAL) /* be careful of off-by-one: useme is 0..L-1 indexed */
    {
      for (apos = 1; apos <= msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (esl_abc_XIsGap    (msa->abc, msa->ax[idx][apos]) ||
		esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
	      break;
	  if (idx == msa->nseq) useme[apos-1] = TRUE; else useme[apos-1] = FALSE;
	}
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL)) /* text mode case */
    {
      for (apos = 0; apos < msa->alen; apos++)
	{
	  for (idx = 0; idx < msa->nseq; idx++)
	    if (strchr(gaps, msa->aseq[idx][apos]) != NULL)
	      break;
	  if (idx == msa->nseq) useme[apos] = TRUE; else useme[apos] = FALSE;
	}
    }

  esl_msa_ColumnSubset(msa, errbuf, useme);
  free(useme);
  return eslOK;

 ERROR:
  if (useme != NULL) free(useme);
  return status;
}


/* Function:  esl_msa_SymConvert()
 * Synopsis:  Global search/replace of symbols in an MSA.
 * Incept:    SRE, Sun Feb 27 11:20:41 2005 [St. Louis]
 *
 * Purpose:   In the aligned sequences in a text-mode <msa>, convert any
 *            residue in the string <oldsyms> to its counterpart (at the same
 *            position) in string <newsyms>.
 * 
 *            To convert DNA to RNA, <oldsyms> could be "Tt" and
 *            <newsyms> could be "Uu". To convert IUPAC symbols to
 *            N's, <oldsyms> could be "RYMKSWHBVDrymkswhbvd" and
 *            <newsyms> could be "NNNNNNNNNNnnnnnnnnnn". 
 *            
 *            As a special case, if <newsyms> consists of a single
 *            character, then any character in the <oldsyms> is 
 *            converted to this character. 
 *            
 *            Thus, <newsyms> must either be of the same length as
 *            <oldsyms>, or of length 1. Anything else will cause
 *            undefined behavior (and probably segfault). 
 *            
 *            The conversion is done in-place, so the <msa> is
 *            modified.
 *            
 *            This is a poor man's hack for processing text mode MSAs
 *            into a more consistent text alphabet. It is unnecessary
 *            for digital mode MSAs, which are already in a standard
 *            internal alphabet. Calling <esl_msa_SymConvert()> on a
 *            digital mode alignment throws an <eslEINVAL> error.
 *            
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEINVAL> if <msa> is in digital mode, or if the <oldsyms>
 *            and <newsyms> strings aren't valid together.
 */
int
esl_msa_SymConvert(ESL_MSA *msa, const char *oldsyms, const char *newsyms)
{
  int64_t apos;			/* column index */
  int     idx;			/* sequence index */
  char   *sptr;
  int     special;

  if (msa->flags & eslMSA_DIGITAL)
    ESL_EXCEPTION(eslEINVAL, "can't SymConvert on digital mode alignment");
  if ((strlen(oldsyms) != strlen(newsyms)) && strlen(newsyms) != 1)
    ESL_EXCEPTION(eslEINVAL, "invalid newsyms/oldsyms pair");

  special = (strlen(newsyms) == 1 ? TRUE : FALSE);

  for (apos = 0; apos < msa->alen; apos++)
    for (idx = 0; idx < msa->nseq; idx++)
      if ((sptr = strchr(oldsyms, msa->aseq[idx][apos])) != NULL)
	msa->aseq[idx][apos] = (special ? *newsyms : newsyms[sptr-oldsyms]);
  return eslOK;
}

/* Function:  esl_msa_AddComment()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an (unparsed) comment line to the MSA structure, 
 *            allocating as necessary.
 *
 * Args:      msa - a multiple alignment
 *            s   - comment line to add
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddComment(ESL_MSA *msa, char *s)
{
  void *p;
  int   status;

  /* If this is our first recorded comment, we need to allocate;
   * and if we've filled available space, we need to reallocate.
   */
  if (msa->comment == NULL) {
    ESL_ALLOC(msa->comment, sizeof(char *) * 16);
    msa->alloc_ncomment = 16;
  }
  if (msa->ncomment == msa->alloc_ncomment) {
    ESL_RALLOC(msa->comment, p, sizeof(char *) * msa->alloc_ncomment * 2);
    msa->alloc_ncomment *= 2;
  }
  if ((status = esl_strdup(s, -1, &(msa->comment[msa->ncomment]))) != eslOK) goto ERROR;
  msa->ncomment++;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_msa_AddGF()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GF+ markup line to the MSA, 
 *            allocating as necessary. <tag> is the GF markup 
 *            tag; <value> is the text associated w/ that tag.
 *
 * Args:      msa - a multiple alignment
 *            tag - markup tag 
 *            value - markup text
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddGF(ESL_MSA *msa, char *tag, char *value)
{  
  void *p;
  int   n;
  int   status;

  /* If this is our first recorded unparsed #=GF line, we need to allocate().
   */
  if (msa->gf_tag == NULL) {
    ESL_ALLOC(msa->gf_tag, sizeof(char *) * 16);
    ESL_ALLOC(msa->gf,     sizeof(char *) * 16);
    msa->alloc_ngf = 16;
  }
  /* or if we're out of room for new GF's, reallocate by doubling
   */
  if (msa->ngf == msa->alloc_ngf) {
    n = msa->alloc_ngf * 2;
    ESL_RALLOC(msa->gf_tag, p, sizeof(char *) * n);
    ESL_RALLOC(msa->gf,     p, sizeof(char *) * n);
    msa->alloc_ngf = n;
  }

  if ((status = esl_strdup(tag,  -1,  &(msa->gf_tag[msa->ngf]))) != eslOK) goto ERROR;
  if ((status = esl_strdup(value, -1, &(msa->gf[msa->ngf])))     != eslOK) goto ERROR;
  msa->ngf++;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  esl_msa_AddGS()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GS+ markup line to the MSA, 
 *            allocating as necessary. It's possible that we 
 *            could get more than one of the same type of GS 
 *            tag per sequence; for example, "DR PDB;" structure 
 *            links in Pfam.  Hack: handle these by appending to 
 *            the string, in a \verb+\n+ separated fashion.
 *
 * Args:      msa    - multiple alignment structure
 *            tag    - markup tag (e.g. "AC")
 *            sqidx  - index of sequence to assoc markup with (0..nseq-1)
 *            value  - markup (e.g. "P00666")
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AddGS(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int   tagidx;
  int   i;
  int   status;

  /* first GS tag? init&allocate  */
  if (msa->gs_tag == NULL)	
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gs_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gs_tag, sizeof(char *));  /* one at a time. */
      ESL_ALLOC(msa->gs,     sizeof(char **));
      ESL_ALLOC(msa->gs[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++)
	msa->gs[0][i] = NULL;
    }
  else 
    {
      /* Get a tagidx for this GS tag.
       * tagidx < ngs; we already saw this tag;
       * tagidx == ngs; this is a new one.
       */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gs_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
#else
      for (tagidx = 0; tagidx < msa->ngs; tagidx++)
	if (strcmp(msa->gs_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocation (in blocks of 1) */
      if (tagidx == msa->ngs ) 
	{
	  ESL_RALLOC(msa->gs_tag, p, (msa->ngs+1) * sizeof(char *));
	  ESL_RALLOC(msa->gs,     p, (msa->ngs+1) * sizeof(char **));
	  ESL_ALLOC(msa->gs[msa->ngs], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gs[msa->ngs][i] = NULL;
	}
    }

  /* Store the tag, if it's new.
   */
  if (tagidx == msa->ngs) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gs_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngs++;
    }
  
  /* Store the annotation on the sequence.
   * If seq is unannotated, dup the value; if
   * seq already has a GS annotation, cat a \n, then cat the value.
   */
  if (msa->gs[tagidx][sqidx] == NULL)
    {
      if ((status = esl_strdup(value, -1, &(msa->gs[tagidx][sqidx]))) != eslOK) goto ERROR;
    }
  else 
    {			
      int n1,n2;
      n1 = strlen(msa->gs[tagidx][sqidx]);
      n2 = strlen(value);
      ESL_RALLOC(msa->gs[tagidx][sqidx], p, sizeof(char) * (n1+n2+2));
      msa->gs[tagidx][sqidx][n1] = '\n';
      strcpy(msa->gs[tagidx][sqidx]+n1+1, value);
    }
  return eslOK;

 ERROR:
  return status;
} 

/* Function:  esl_msa_AppendGC()
 * Incept:    SRE, Tue Jun  1 17:37:21 1999 [St. Louis]
 *
 * Purpose:   Add an unparsed \verb+#=GC+ markup line to the MSA 
 *            structure, allocating as necessary. When called 
 *            multiple times for the same tag, appends value 
 *            strings together -- used when parsing multiblock 
 *            alignment files, for example.
 *
 * Args:      msa   - multiple alignment structure
 *            tag   - markup tag (e.g. "CS")
 *            value - markup, one char per aligned column      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AppendGC(ESL_MSA *msa, char *tag, char *value)
{
  int   tagidx;
  int   status;
  void *p;

  /* Is this an unparsed tag name that we recognize?
   * If not, handle adding it to index, and reallocating
   * as needed.
   */
  if (msa->gc_tag == NULL)	/* first tag? init&allocate  */
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gc_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gc_idx, tag, &tagidx);      
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gc_tag, sizeof(char **));
      ESL_ALLOC(msa->gc,     sizeof(char **));
      msa->gc[0]  = NULL;
    }
  else
    {			/* new tag? */
      /* get tagidx for this GC tag. existing tag: <ngc; new: == ngc. */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gc_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) goto ERROR;
#else
      for (tagidx = 0; tagidx < msa->ngc; tagidx++)
	if (strcmp(msa->gc_tag[tagidx], tag) == 0) break;
#endif
      /* Reallocate, in block of one tag at a time
       */
      if (tagidx == msa->ngc)
	{
	  ESL_RALLOC(msa->gc_tag, p, (msa->ngc+1) * sizeof(char **));
	  ESL_RALLOC(msa->gc,     p, (msa->ngc+1) * sizeof(char **));
	  msa->gc[tagidx] = NULL;
	}
    }
  /* new tag? store it.
   */
  if (tagidx == msa->ngc) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gc_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngc++;
    }
  return (esl_strcat(&(msa->gc[tagidx]), -1, value, -1));

 ERROR:
  return status;
}

/* Function:  esl_msa_AppendGR()
 * Incept:    SRE, Thu Jun  3 06:34:38 1999 [Madison]
 *
 * Purpose:   Add an unparsed \verb+#=GR+ markup line to the MSA structure, 
 *            allocating as necessary.
 *              
 *            When called multiple times for the same tag, appends 
 *            value strings together -- used when parsing multiblock 
 *            alignment files, for example.
 *
 * Args:      msa    - multiple alignment structure
 *            tag    - markup tag (e.g. "SS")
 *            sqidx  - index of seq to assoc markup with (0..nseq-1)
 *            value  - markup, one char per aligned column      
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_msa_AppendGR(ESL_MSA *msa, char *tag, int sqidx, char *value)
{
  void *p;
  int tagidx;
  int i;
  int status;

  if (msa->gr_tag == NULL)	/* first tag? init&allocate  */
    {
#ifdef eslAUGMENT_KEYHASH
      msa->gr_idx = esl_keyhash_Create();
      status = esl_key_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
      ESL_DASSERT1((tagidx == 0));
#else
      tagidx = 0;
#endif
      ESL_ALLOC(msa->gr_tag, sizeof(char *));
      ESL_ALLOC(msa->gr,     sizeof(char **));
      ESL_ALLOC(msa->gr[0],  sizeof(char *) * msa->sqalloc);
      for (i = 0; i < msa->sqalloc; i++) 
	msa->gr[0][i] = NULL;
    }
  else 
    {
      /* get tagidx for this GR tag. existing<ngr; new=ngr.
       */
#ifdef eslAUGMENT_KEYHASH
      status = esl_key_Store(msa->gr_idx, tag, &tagidx);
      if (status != eslOK && status != eslEDUP) return status;
#else
      for (tagidx = 0; tagidx < msa->ngr; tagidx++)
	if (strcmp(msa->gr_tag[tagidx], tag) == 0) break;
#endif
      /* if a new tag, realloc for it */      
      if (tagidx == msa->ngr)
	{ 
	  ESL_RALLOC(msa->gr_tag, p, (msa->ngr+1) * sizeof(char *));
	  ESL_RALLOC(msa->gr,     p, (msa->ngr+1) * sizeof(char **));
	  ESL_ALLOC(msa->gr[msa->ngr], sizeof(char *) * msa->sqalloc);
	  for (i = 0; i < msa->sqalloc; i++) 
	    msa->gr[msa->ngr][i] = NULL;
	}
    }

  if (tagidx == msa->ngr) 
    {
      if ((status = esl_strdup(tag, -1, &(msa->gr_tag[tagidx]))) != eslOK) goto ERROR;
      msa->ngr++;
    }
  return (esl_strcat(&(msa->gr[tagidx][sqidx]), -1, value, -1));

 ERROR:
  return status;
}


/* Function:  esl_msa_Checksum()
 * Synopsis:  Calculate a checksum for an MSA.
 * Incept:    SRE, Tue Sep 16 13:23:34 2008 [Janelia]
 *
 * Purpose:   Calculates a 32-bit checksum for <msa>.
 * 
 *            Only the alignment data are considered, not the sequence
 *            names or other annotation. For text mode alignments, the
 *            checksum is case sensitive.
 *            
 *            This is used as a quick way to try to verify that a
 *            given alignment is identical to an expected one; for
 *            example, when HMMER is mapping new sequence alignments
 *            onto exactly the same seed alignment an HMM was built
 *            from.
 *
 * Returns:   <eslOK> on success.
 *
 * Xref:      The checksum is a modified version of Jenkin's hash;
 *            see <esl_keyhash> for the original and citations.
 */
int
esl_msa_Checksum(const ESL_MSA *msa, uint32_t *ret_checksum)
{
  uint32_t val = 0;
  int      i,pos;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      for (i = 0; i < msa->nseq; i++)
	for (pos = 1; pos <= msa->alen; pos++)
	  {
	    val += msa->ax[i][pos];
	    val += (val << 10);
	    val ^= (val >>  6);
	  }
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      for (i = 0; i < msa->nseq; i++)
	for (pos = 0; pos < msa->alen; pos++)
	  {
	    val += msa->aseq[i][pos];
	    val += (val << 10);
	    val ^= (val >>  6);
	  }
    }
  val += (val <<  3);
  val ^= (val >> 11);
  val += (val << 15);

  *ret_checksum = val;
  return eslOK;
}
/*-------------------- end of misc MSA functions ----------------------*/

/*****************************************************************
 * 7. Stockholm (Pfam/Rfam) format
 *****************************************************************/

/* msafile_getline():
 * load the next line of <afp> into <afp->buf>. 
 * Returns eslOK on success, eslEOF on normal eof.
 * Throws eslEMEM on alloc failure.
 */
static int
msafile_getline(ESL_MSAFILE *afp)
{
  int status;
  status = esl_fgets(&(afp->buf), &(afp->buflen), afp->f);
  afp->linenumber++;
  return status;
}

/* maxwidth()
 * Return the length of the longest string in 
 * an array of strings.
 */
static int64_t
maxwidth(char **s, int n)
{
  int64_t max,len;
  int     i; 
  
  max = 0;
  for (i = 0; i < n; i++)
    if (s[i] != NULL)
      {
	len = strlen(s[i]);
	if (len > max) max = len;
      }
  return max;
}
static int
is_blankline(char *s)
{
  for (; *s != '\0'; s++)
    if (! isspace((int) *s)) return FALSE;
  return TRUE;
}

/* Format of a GF line:
 *    #=GF <tag> <text>
 * Returns eslOK on success; eslEFORMAT on parse failure.
 * Throws eslEMEM on allocation failure.
 */
static int
parse_gf(ESL_MSA *msa, char *buf)
{
  char *gf;
  char *tag;
  char *text;
  char *tok;
  char *s;
  int   n;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gf)  != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag) != eslOK) return eslEFORMAT;

  /* text might be empty; watch out for this. (for example, a blank #=GF CC line) */
  status = esl_strtok_adv(&s, "\n\r",    &text, &n, NULL);
  if      (status == eslOK) { while (*text && (*text == ' ' || *text == '\t')) text++; }
  else if (status == eslEOL){ text = NULL; n = 0; } 
  else return eslEFORMAT;

  if      (strcmp(tag, "ID") == 0) status = esl_strdup(text, n, &(msa->name));
  else if (strcmp(tag, "AC") == 0) status = esl_strdup(text, n, &(msa->acc));
  else if (strcmp(tag, "DE") == 0) status = esl_strdup(text, n, &(msa->desc));
  else if (strcmp(tag, "AU") == 0) status = esl_strdup(text, n, &(msa->au));
  else if (strcmp(tag, "GA") == 0) 
    {				/* Pfam has GA1, GA2. Rfam just has GA1. */
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_GA1] = atof(tok);
      msa->cutset[eslMSA_GA1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_GA2] = atof(tok);
	  msa->cutset[eslMSA_GA2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "NC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_NC1] = atof(tok);
      msa->cutset[eslMSA_NC1] = TRUE;
      if ((esl_strtok(&s, " \t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_NC2] = atof(tok);
	  msa->cutset[eslMSA_NC2] = TRUE;
	}
      status = eslOK;
    }
  else if (strcmp(tag, "TC") == 0) 
    {
      s = text;
      if ((esl_strtok(&s, " \t\n\r", &tok)) != eslOK) 
	return eslEFORMAT;
      msa->cutoff[eslMSA_TC1] = atof(tok);
      msa->cutset[eslMSA_TC1] = TRUE;
      if ((esl_strtok(&s, "\t\n\r", &tok)) == eslOK) 
	{
	  msa->cutoff[eslMSA_TC2] = atof(tok);
	  msa->cutset[eslMSA_TC2] = TRUE;
	}
      status = eslOK;
    }
  else 				/* an unparsed #=GF: */
    status = esl_msa_AddGF(msa, tag, text);

  return status;
}


/* Format of a GS line:
 *    #=GS <seqname> <tag> <text>
 * Return <eslOK> on success; <eslEFORMAT> on parse error.
 * Throws <eslEMEM> on allocation error (trying to grow for a new
 *        name; <eslEINVAL> if we try to grow an ungrowable MSA.
 */
static int
parse_gs(ESL_MSA *msa, char *buf)
{
  char *gs;
  char *seqname;
  char *tag;
  char *text; 
  int   seqidx;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok(&s, " \t\n\r", &gs)      != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &seqname) != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, " \t\n\r", &tag)     != eslOK) return eslEFORMAT;
  if (esl_strtok(&s, "\n\r",    &text)    != eslOK) return eslEFORMAT;
  while (*text && (*text == ' ' || *text == '\t')) text++;
  
  /* GS usually follows another GS; guess lastidx+1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "WT") == 0)
    {
      msa->wgt[seqidx] = atof(text);
      msa->flags      |= eslMSA_HASWGTS;
      status           = eslOK;
    }
  else if (strcmp(tag, "AC") == 0)
    status = esl_msa_SetSeqAccession(msa, seqidx, text);
  else if (strcmp(tag, "DE") == 0)
    status = esl_msa_SetSeqDescription(msa, seqidx, text);
  else				
    status = esl_msa_AddGS(msa, tag, seqidx, text);

  return status;
}


/* parse_gc():
 * Format of a GC line:
 *    #=GC <tag> <aligned text>
 */
static int 
parse_gc(ESL_MSA *msa, char *buf)
{
  char *gc;
  char *tag;
  char *text; 
  char *s;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &gc)               != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &tag)              != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT;
  
  if      (strcmp(tag, "SS_cons") == 0)  status = esl_strcat(&(msa->ss_cons), -1, text, len);
  else if (strcmp(tag, "SA_cons") == 0)  status = esl_strcat(&(msa->sa_cons), -1, text, len);
  else if (strcmp(tag, "PP_cons") == 0)  status = esl_strcat(&(msa->pp_cons), -1, text, len);
  else if (strcmp(tag, "RF")      == 0)  status = esl_strcat(&(msa->rf),      -1, text, len);
  else                                   status = esl_msa_AppendGC(msa, tag, text);

  return status;
}

/* parse_gr():
 * Format of a GR line:
 *    #=GR <seqname> <featurename> <text>
 */
static int
parse_gr(ESL_MSA *msa, char *buf)
{
  char *gr;
  char *seqname;
  char *tag;
  char *text;
  int   seqidx;
  int   len;
  int   j;
  char *s;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &gr)               != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &seqname)          != eslOK) return eslEFORMAT;
  if (esl_strtok    (&s, " \t\n\r", &tag)              != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT;

  /* GR usually follows sequence it refers to; guess msa->lastidx */
  status = get_seqidx(msa, seqname, msa->lastidx, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

  if (strcmp(tag, "SS") == 0) 
    {
      if (msa->ss == NULL)
	{
	  ESL_ALLOC(msa->ss,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC(msa->sslen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++)
	    {
	      msa->ss[j]    = NULL;
	      msa->sslen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->ss[seqidx]), msa->sslen[seqidx], text, len);
      msa->sslen[seqidx] += len;
    }
  else if (strcmp(tag, "SA") == 0)
    {
      if (msa->sa == NULL)
	{
	  ESL_ALLOC(msa->sa,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC(msa->salen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->sa[j]    = NULL;
	      msa->salen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->sa[seqidx]), msa->salen[seqidx], text, len);
      msa->salen[seqidx] += len;
    }
  else if (strcmp(tag, "PP") == 0)
    {
      if (msa->pp == NULL)
	{
	  ESL_ALLOC(msa->pp,    sizeof(char *) * msa->sqalloc);
	  ESL_ALLOC(msa->pplen, sizeof(int64_t)* msa->sqalloc);
	  for (j = 0; j < msa->sqalloc; j++) 
	    {
	      msa->pp[j]    = NULL;
	      msa->pplen[j] = 0;
	    }
	}
      status = esl_strcat(&(msa->pp[seqidx]), msa->pplen[seqidx], text, len);
      msa->pplen[seqidx] += len;
    }
  else 
    status = esl_msa_AppendGR(msa, tag, seqidx, text);
  return status;

 ERROR:
  return status;
}


/* parse_comment():
 * comments are simply stored verbatim, not parsed
 */
static int
parse_comment(ESL_MSA *msa, char *buf)
{
  char *s;
  char *comment;

  s = buf + 1;			               /* skip leading '#' */
  if (*s == '\n' || *s == '\r') { *s = '\0'; comment = s; }  /* deal with blank comment */
  else if (esl_strtok(&s, "\n\r", &comment)!= eslOK) return eslEFORMAT;
  return (esl_msa_AddComment(msa, comment));
}

/* parse_sequence():
 * Format of line is:
 *     <name>  <aligned text>
 * 
 * On digital sequence, returns <eslEINVAL> if any of the residues can't be digitized.
 */
static int
parse_sequence(ESL_MSA *msa, char *buf)
{
  char *s;
  char *seqname;
  char *text;
  int   seqidx;
  int   len;
  int   status;

  s = buf;
  if (esl_strtok    (&s, " \t\n\r", &seqname)          != eslOK) return eslEFORMAT;
  if (esl_strtok_adv(&s, " \t\n\r", &text, &len, NULL) != eslOK) return eslEFORMAT; 
  
  /* seq usually follows another seq; guess msa->lastidx +1 */
  status = get_seqidx(msa, seqname, msa->lastidx+1, &seqidx);
  if (status != eslOK) return status;
  msa->lastidx = seqidx;

#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), text, len);
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
      msa->sqlen[seqidx] += len;
    }

  return status;
}

/* read_stockholm():
 * SRE, Sun Jan 23 08:33:32 2005 [St. Louis]
 *
 * Purpose:   Parse the next alignment from an open Stockholm format alignment
 *            file <afp>, leaving the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            Returns <eslEOF> if there are no more alignments in <afp>,
 *            and <ret_msa> is set to NULL.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem, and <ret_msa> is
 *            set to NULL. 
 *
 *            Returns <eslEINVAL> if we're trying to read a digital alignment,
 *            and an invalid residue is found that can't be digitized.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      squid's ReadStockholm(), 1999.
 */
static int
read_stockholm(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa = NULL;
  char      *s;
  int        status;
  int        status2;
#ifdef eslAUGMENT_SSI
  off_t      offset;
#endif

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
#ifdef eslAUGMENT_SSI
    offset = ftello(afp->f);
#endif
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR; /* includes EOF  */
  } while (is_blankline(afp->buf));

  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);

#ifdef eslAUGMENT_SSI
  msa->offset = offset;
#endif

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */

      if (*s == '#') {

	if      (strncmp(s, "#=GF", 4) == 0)
	  {
	    if ((status = parse_gf(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GF line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    if ((status = parse_gs(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GS line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if  ((status = parse_gc(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GC line", afp->linenumber);
	  }

	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if ((status = parse_gr(msa, s)) != eslOK)
	      ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GR line", afp->linenumber);
	  }

	else if ((status = parse_comment(msa, s)) != eslOK)
	  ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad comment line", afp->linenumber);
      } 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n' || *s == '\r')     continue;
      else if ((status = parse_sequence(msa, s)) != eslOK)
	ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad sequence line", afp->linenumber);
    }
  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): didn't find // at end of alignment", afp->linenumber);
  
  /* Stockholm fmt is complex, so give the newly parsed MSA a good
   * going-over, and finalize the fields of the MSA data structure.
   * verify_parse will fill in errbuf if it sees a problem.
   */
  if (verify_parse(msa, afp->errbuf) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}

/* actually_write_stockholm()
 * SRE, Fri May 21 17:39:22 1999 [St. Louis]
 *
 * Write an alignment in Stockholm format to an open file. This is the
 * function that actually does the work. The API's WriteStockholm()
 * and WriteStockholmOneBlock() are wrappers.
 *
 * Args:     fp    - file that's open for writing
 *           msa   - alignment to write        
 *           cpl   - characters to write per line in alignment block
 *
 * Returns:  eslOK on success.
 * 
 * Throws:   eslEMEM on allocation failure.
 */
static int
actually_write_stockholm(FILE *fp, const ESL_MSA *msa, int cpl)
{
  int  i, j;
  int  maxname;		/* maximum name length     */
  int  maxgf;		/* max #=GF tag length     */
  int  maxgc;		/* max #=GC tag length     */
  int  maxgr; 		/* max #=GR tag length     */
  int  margin;        	/* total left margin width */
  int  gslen;		/* length of a #=GS tag    */
  char *buf = NULL;
  int  currpos;
  char *s, *tok;
  int  acpl;            /* actual number of character per line */
  int  status;
  
  /* Figure out how much space we need for name + markup
   * to keep the alignment in register. Required by Stockholm
   * spec, even though our Stockholm parser doesn't care (Erik's does).
   *
   * The left margin of an alignment block can be composed of:
   * 
   * <seqname>                      max length: maxname + 1
   * #=GC <gc_tag>                  max length: 4 + 1 + maxgc + 1
   * #=GR <seqname> <gr_tag>        max length: 4 + 1 + maxname + 1 + maxgr + 1
   * 
   * <margin> is the max of these. It is the total length of the
   * left margin that we need to leave, inclusive of the last space.
   * 
   * Then when we output, we do:
   * name:  <leftmargin-1>
   * gc:    #=GC <leftmargin-6>
   * gr:    #=GR <maxname> <leftmargin-maxname-7>
   *
   * xref STL9/p17
   */
  maxname = maxwidth(msa->sqname, msa->nseq);
  
  maxgf   = maxwidth(msa->gf_tag, msa->ngf);
  if (maxgf < 2) maxgf = 2;

  maxgc   = maxwidth(msa->gc_tag, msa->ngc);
  if (msa->rf      !=NULL && maxgc < 2) maxgc = 2;
  if (msa->ss_cons !=NULL && maxgc < 7) maxgc = 7;
  if (msa->sa_cons !=NULL && maxgc < 7) maxgc = 7;
  if (msa->pp_cons !=NULL && maxgc < 7) maxgc = 7;

  maxgr   = maxwidth(msa->gr_tag, msa->ngr);
  if (msa->ss != NULL && maxgr < 2) maxgr = 2;
  if (msa->sa != NULL && maxgr < 2) maxgr = 2;
  if (msa->pp != NULL && maxgr < 2) maxgr = 2;

  margin = maxname + 1;
  if (maxgc > 0 && maxgc+6 > margin) margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  
  /* Allocate a tmp buffer to hold sequence chunks in
   */
  ESL_ALLOC(buf, sizeof(char) * (cpl+1));

  /* Magic Stockholm header
   */
  fprintf(fp, "# STOCKHOLM 1.0\n");

  /* Free text comments
   */
  for (i = 0;  i < msa->ncomment; i++)
    fprintf(fp, "#%s\n", msa->comment[i]);
  if (msa->ncomment > 0) fprintf(fp, "\n");

  /* GF section: per-file annotation
   */
  if (msa->name != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "ID", msa->name);
  if (msa->acc  != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AC", msa->acc);
  if (msa->desc != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "DE", msa->desc);
  if (msa->au   != NULL) fprintf(fp, "#=GF %-*s %s\n", maxgf, "AU", msa->au);
  
  /* Thresholds are hacky. Pfam has two. Rfam has one.
   */
  if      (msa->cutset[eslMSA_GA1] && msa->cutset[eslMSA_GA2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1], msa->cutoff[eslMSA_GA2]);
  else if (msa->cutset[eslMSA_GA1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "GA", msa->cutoff[eslMSA_GA1]);

  if      (msa->cutset[eslMSA_NC1] && msa->cutset[eslMSA_NC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1], msa->cutoff[eslMSA_NC2]);
  else if (msa->cutset[eslMSA_NC1])
    fprintf(fp, "#=GF %-*s %.1f\n",
	    maxgf, "NC", msa->cutoff[eslMSA_NC1]);

  if      (msa->cutset[eslMSA_TC1] && msa->cutset[eslMSA_TC2])
    fprintf(fp, "#=GF %-*s %.1f %.1f\n",
	    maxgf, "TC", msa->cutoff[eslMSA_TC1], msa->cutoff[eslMSA_TC2]);
  else if (msa->cutset[eslMSA_TC1])
    fprintf(fp, "#=GF %-*s %.1f\n", 
	    maxgf, "TC", msa->cutoff[eslMSA_TC1]);

  for (i = 0; i < msa->ngf; i++)
    fprintf(fp, "#=GF %-*s %s\n", maxgf, msa->gf_tag[i], msa->gf[i]); 
  fprintf(fp, "\n");


  /* GS section: per-sequence annotation
   */
  if (msa->flags & eslMSA_HASWGTS) 
    {
      for (i = 0; i < msa->nseq; i++) 
	fprintf(fp, "#=GS %-*s WT %.2f\n", 
		maxname, msa->sqname[i], msa->wgt[i]);		
      fprintf(fp, "\n");
    }

  if (msa->sqacc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqacc[i] != NULL)
	  fprintf(fp, "#=GS %-*s AC %s\n", 
		  maxname, msa->sqname[i], msa->sqacc[i]);
      fprintf(fp, "\n");
    }

  if (msa->sqdesc != NULL) 
    {
      for (i = 0; i < msa->nseq; i++) 
	if (msa->sqdesc[i] != NULL)
	  fprintf(fp, "#=GS %-*s DE %s\n", 
		  maxname, msa->sqname[i], msa->sqdesc[i]);
      fprintf(fp, "\n");
    }

  for (i = 0; i < msa->ngs; i++)
    {
      /* Multiannotated GS tags are possible; for example, 
       *     #=GS foo DR PDB; 1xxx;
       *     #=GS foo DR PDB; 2yyy;
       * These are stored, for example, as:
       *     msa->gs[0][0] = "PDB; 1xxx;\nPDB; 2yyy;"
       * and must be decomposed.
       */
      gslen = strlen(msa->gs_tag[i]);
      for (j = 0; j < msa->nseq; j++)
	if (msa->gs[i][j] != NULL)
	  {
	    s = msa->gs[i][j];
	    while (esl_strtok(&s, "\n", &tok) == eslOK)
	      fprintf(fp, "#=GS %-*s %-*s %s\n", 
		      maxname, msa->sqname[j],
		      gslen,   msa->gs_tag[i], 
		      tok);
	  }
      fprintf(fp, "\n");
    }

  /* Alignment section:
   * contains aligned sequence, #=GR annotation, and #=GC annotation
   */
  for (currpos = 0; currpos < msa->alen; currpos += cpl)
    {
      acpl = (msa->alen - currpos > cpl)? cpl : msa->alen - currpos;

      if (currpos > 0) fprintf(fp, "\n");
      for (i = 0; i < msa->nseq; i++)
	{
#ifdef eslAUGMENT_ALPHABET
	  if (msa->flags & eslMSA_DIGITAL)
	    esl_abc_TextizeN(msa->abc, msa->ax[i] + currpos + 1, acpl, buf);
	  else
	    strncpy(buf, msa->aseq[i] + currpos, acpl);
#else
	  strncpy(buf, msa->aseq[i] + currpos, acpl);
#endif
	  
	  buf[acpl] = '\0';	      
	  fprintf(fp, "%-*s %s\n", 
		  margin-1, msa->sqname[i], buf);

	  if (msa->ss != NULL && msa->ss[i] != NULL) {
	    strncpy(buf, msa->ss[i] + currpos, acpl);
	    buf[acpl] = '\0';	 
	    fprintf(fp, "#=GR %-*s %-*s %s\n", 
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SS",
		    buf);
	  }
	  if (msa->sa != NULL && msa->sa[i] != NULL) {
	    strncpy(buf, msa->sa[i] + currpos, acpl);
	    buf[acpl] = '\0';
	    fprintf(fp, "#=GR %-*s %-*s %s\n",
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "SA",
		    buf);
	  }
	  if (msa->pp != NULL && msa->pp[i] != NULL) {
	    strncpy(buf, msa->pp[i] + currpos, acpl);
	    buf[acpl] = '\0';
	    fprintf(fp, "#=GR %-*s %-*s %s\n",
		    maxname,          msa->sqname[i],
		    margin-maxname-7, "PP",
		    buf);
	  }
	  for (j = 0; j < msa->ngr; j++)
	    if (msa->gr[j][i] != NULL) {
	      strncpy(buf, msa->gr[j][i] + currpos, acpl);
	      buf[acpl] = '\0';
	      fprintf(fp, "#=GR %-*s %-*s %s\n", 
		      maxname,          msa->sqname[i],
		      margin-maxname-7, msa->gr_tag[j],
		      buf);
	    }
	}

      if (msa->ss_cons != NULL) {
	strncpy(buf, msa->ss_cons + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SS_cons", buf);
      }
      if (msa->sa_cons != NULL) {
	strncpy(buf, msa->sa_cons + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "SA_cons", buf);
      }
      if (msa->pp_cons != NULL) {
	strncpy(buf, msa->pp_cons + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "PP_cons", buf);
      }

      if (msa->rf != NULL) {
	strncpy(buf, msa->rf + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, "RF", buf);
      }
      for (j = 0; j < msa->ngc; j++) {
	strncpy(buf, msa->gc[j] + currpos, acpl);
	buf[acpl] = '\0';
	fprintf(fp, "#=GC %-*s %s\n", margin-6, msa->gc_tag[j], buf);
      }
    }
  fprintf(fp, "//\n");
  free(buf);
  return eslOK;

 ERROR:
  if (buf != NULL) free(buf);
  return status;
}


/* write_stockholm():
 * SRE, Fri Jan 28 09:24:02 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in multiblock format, with
 *            50 residues per line.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholm(), 1999.
 */
static int
write_stockholm(FILE *fp, const ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, 50)); /* 50 char per block */
}

/* write_pfam():
 * SRE, Fri Jan 28 09:25:42 2005 [St. Louis]
 *
 * Purpose:   Write an alignment <msa> in Stockholm format 
 *            to a stream <fp>, in single block (Pfam) format.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's WriteStockholmOneBlock(), 1999.
 */
static int
write_pfam(FILE *fp, const ESL_MSA *msa)
{
  return (actually_write_stockholm(fp, msa, msa->alen)); /* one big block */
}
/*---------------- end, stockholm/pfam format -------------------*/

/*****************************************************************
 * 8. A2M format
 *****************************************************************/

/* write_a2m()
 * SRE, Wed Sep  3 09:44:36 2008 [Janelia]
 *
 * Purpose:   Write alignment <msa> in dotless UCSC A2M format to a
 *            stream <fp>.
 *            
 *            The <msa> should have a valid reference line <msa->rf>,
 *            with alphanumeric characters marking consensus (match)
 *            columns, and non-alphanumeric characters marking
 *            nonconsensus (insert) columns. 
 *            
 *            As a fallback, if the <msa> does not have a reference
 *            line, a default one is created. In this case,
 *            <esl_msa_Fragmentize()> will be called with
 *            <fragthresh=0.5> to heuristically define sequence
 *            fragments in the alignment, and <esl_msa_ReasonableRF()>
 *            will be called with <symfrac=0.5> to mark consensus
 *            columns.  This will modify the alignment data
 *            (converting leading/trailing gap symbols on "fragments"
 *            to '~' missing data symbols) in addition to adding
 *            <#=RF> annotation to it. If caller doesn't want the
 *            alignment data to be modified, it must provide its own
 *            <msa->rf> line, or it must avoid writing alignments in
 *            A2M format.
 *            
 *            In "dotless" A2M format, gap characters in insert
 *            columns are omitted; therefore sequences can be of
 *            different lengths, but each sequence has the same number
 *            of consensus columns (residue or -).
 *            
 *            A2M format cannot represent missing data symbols
 *            (Easel's ~). Any missing data symbols are converted to
 *            gaps.
 *            
 *            A2M format cannot represent pyrrolysine residues in
 *            amino acid sequences, because it treats 'O' symbols
 *            specially, as indicating a position at which a
 *            free-insertion module (FIM) should be created. Any 'O'
 *            in the <msa> is converted to 'X' (the unknown residue
 *            symbol for amino acid sequences).
 *            
 *            Because of the 'O' issue, A2M format should not be used
 *            for nonstandard/custom Easel alphabets that include 'O'
 *            as a symbol. For <msa>'s in digital mode, where digital
 *            alphabet information is available, <eslEINVAL> is
 *            returned if the <msa>'s has a nonstandard digital
 *            alphabet that uses 'O' as a symbol. For <msa>'s in text
 *            mode, because text-mode sequences are treated literally,
 *            'O' is always converted to 'X', which may not be what
 *            the caller intends.
 *
 * Args:      fp  - open output stream
 *            msa - MSA to write       
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> if the <msa> is in digital mode and uses a nonstandard
 *            alphabet that expects 'O' to be a legal symbol.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      http://www.soe.ucsc.edu/compbio/a2m-desc.html
 */
static int
write_a2m(FILE *fp, ESL_MSA *msa)
{
  int     i;
  int64_t pos;
  int     bpos;
  char    buf[61];
  int     is_consensus;
  int     is_residue;
  int     sym;
  int     do_dotless = TRUE;
  int     status;
  
#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    if (strchr(msa->abc->sym, 'O') != NULL) return eslEINVAL;
#endif

  /* if <msa> lacks an RF line, make a default one */
  if (msa->rf == NULL) 
    {
      ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
      if ((status = esl_msa_MarkFragments(msa, 0.5))          !=  eslOK) return status;
      if ((status = esl_msa_ReasonableRF (msa, 0.5, msa->rf)) !=  eslOK) return status;
    }


  for (i = 0; i < msa->nseq; i++)
    {
      /* Construct the description line */
      fprintf(fp, ">%s", msa->sqname[i]);
      if (msa->sqacc  != NULL && msa->sqacc[i]  != NULL) fprintf(fp, " %s", msa->sqacc[i]);
      if (msa->sqdesc != NULL && msa->sqdesc[i] != NULL) fprintf(fp, " %s", msa->sqdesc[i]);
      fputc('\n', fp);

#ifdef eslAUGMENT_ALPHABET
      if (msa->flags & eslMSA_DIGITAL)
	{
	  pos = 0;
	  while (pos < msa->alen)
	    {
	      for (bpos = 0; pos < msa->alen && bpos < 60; pos++)
		{
		  sym          = msa->abc->sym[msa->ax[i][pos+1]]; /* note off-by-one in digitized aseq: 1..alen */
		  is_consensus = (isalnum(msa->rf[pos]) ? TRUE : FALSE);
		  is_residue   = esl_abc_XIsResidue(msa->abc, msa->ax[i][pos+1]);

		  if (msa->abc->type == eslAMINO && sym == 'O') 
		    sym = esl_abc_XGetUnknown(msa->abc); /* watch out: O means "insert a FIM" in a2m format, not pyrrolysine */

		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos++] = toupper(sym);
		      else            buf[bpos++] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if      (is_residue)   buf[bpos++] = tolower(sym);
		      else if (! do_dotless) buf[bpos++] = '.';          /* includes conversion of missing data to . */
		    }
		}
	      if (bpos) 
		{
		  buf[bpos] = '\0';
		  fprintf(fp, "%s\n", buf);	      
		}
	    }
	}
#endif
      if (! (msa->flags & eslMSA_DIGITAL))
	{
	  pos = 0;
	  while (pos < msa->alen)
	    {
	      for (bpos = 0; pos < msa->alen && bpos < 60; pos++)
		{
		  sym          = msa->aseq[i][pos];
		  is_consensus = (isalnum(msa->rf[pos])) ? TRUE : FALSE;
		  is_residue   = isalpha(msa->aseq[i][pos]);

		  if (sym == 'O') sym = 'X';

		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos++] = toupper(sym);
		      else            buf[bpos++] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if      (is_residue)   buf[bpos++] = tolower(sym);
		      else if (! do_dotless) buf[bpos++] = '.';          /* includes conversion of missing data to . */
		    }
		}
	      if (bpos) 
		{
		  buf[bpos] = '\0';
		  fprintf(fp, "%s\n", buf);	      
		}
	    } 
	}
    } /* end, loop over sequences in the MSA */

  return eslOK;

 ERROR:
  return status;
}
/*---------------------- end, A2M format ------------------------*/




/*****************************************************************
 * 9. PSIBLAST format
 *****************************************************************/

/* write_psiblast()
 * SRE, Wed Sep  3 13:05:10 2008 [Janelia]
 *
 * Purpose:   Write alignment <msa> in NCBI PSI-BLAST format to 
 *            stream <fp>.
 *            
 *            The <msa> should have a valid reference line <msa->rf>,
 *            with alphanumeric characters marking consensus (match)
 *            columns, and non-alphanumeric characters marking
 *            nonconsensus (insert) columns. 
 *            
 *            As a fallback, if the <msa> does not have a reference
 *            line, a default one is created. In this case,
 *            <esl_msa_Fragmentize()> will be called with
 *            <fragthresh=0.5> to heuristically define sequence
 *            fragments in the alignment, and <esl_msa_ReasonableRF()>
 *            will be called with <symfrac=0.5> to mark consensus
 *            columns.  This will modify the alignment data
 *            (converting leading/trailing gap symbols on "fragments"
 *            to '~' missing data symbols) in addition to adding
 *            <#=RF> annotation to it. If caller doesn't want the
 *            alignment data to be modified, it must provide its own
 *            <msa->rf> line, or it must avoid writing alignments in
 *            PSI-BLAST format.
 *            
 *            PSI-BLAST format allows only one symbol ('-') for gaps,
 *            and cannot represent missing data symbols (Easel's
 *            '~'). Any missing data symbols are converted to gaps.
 *
 * Args:      fp  - open output stream
 *            msa - MSA to write       
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
write_psiblast(FILE *fp, ESL_MSA *msa)
{
  int     cpl = 60;
  char    buf[61];
  int     i;
  int     sym;
  int     pos;
  int     bpos;
  int     acpl;
  int     maxnamewidth;
  int     is_consensus;
  int     is_residue;
  int     status;

  maxnamewidth = (int) maxwidth(msa->sqname, msa->nseq);

  /* if <msa> lacks an RF line, make a default one */
  if (msa->rf == NULL) 
    {
      ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen+1));
      if ((status = esl_msa_MarkFragments(msa, 0.5))          !=  eslOK) return status;
      if ((status = esl_msa_ReasonableRF (msa, 0.5, msa->rf)) !=  eslOK) return status;
    }

  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      for (i = 0; i < msa->nseq; i++)
	{
	  acpl =  (msa->alen - pos > cpl)? cpl : msa->alen - pos;

#ifdef eslAUGMENT_ALPHABET
	  if (msa->flags & eslMSA_DIGITAL)
	    {
	      for (bpos = 0; bpos < acpl; bpos++)
		{
		  sym          = msa->abc->sym[msa->ax[i][pos + bpos + 1]];
		  is_consensus = (isalnum(msa->rf[pos + bpos + 1]) ? TRUE : FALSE);
		  is_residue   = esl_abc_XIsResidue(msa->abc, msa->ax[i][pos+bpos+1]);
				      
		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos] = toupper(sym);
		      else            buf[bpos] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if (is_residue) buf[bpos] = tolower(sym);
		      else            buf[bpos] = '-';     /* includes conversion of missing data to - */
		    }
		}
	    }
#endif
	  if (! (msa->flags & eslMSA_DIGITAL))
	    {
	      for (bpos = 0; bpos < acpl; bpos++)
		{
		  sym          = msa->aseq[i][pos + bpos];
		  is_consensus = (isalnum(msa->rf[pos + bpos]) ? TRUE : FALSE);
		  is_residue   = isalnum(sym);

		  if (is_consensus) 
		    {
		      if (is_residue) buf[bpos] = toupper(sym);
		      else            buf[bpos] = '-';          /* includes conversion of missing data to - */
		    }
		  else
		    {
		      if (is_residue) buf[bpos] = tolower(sym);
		      else            buf[bpos] = '-';     /* includes conversion of missing data to - */
		    }
		}
	    }
	  buf[acpl] = '\0';	      
	  fprintf(fp, "%-*s  %s\n", maxnamewidth, msa->sqname[i], buf);
	}  /* end loop over sequences */

      if (pos + cpl < msa->alen) fputc('\n', fp);
    } /* end loop over alignment blocks */
  return eslOK;
  
 ERROR:
  return status;
}
/*------------------------- end, psiblast format -----------------------------*/



/*****************************************************************
 * 10. SELEX format
 *****************************************************************/
#define eslMSA_LINE_SQ 1
#define eslMSA_LINE_RF 2
#define eslMSA_LINE_CS 3
#define eslMSA_LINE_SS 4
#define eslMSA_LINE_SA 5

static int read_block(ESL_MSAFILE *afp, char ***line_p, int **llen_p, int **lpos_p, int **rpos_p, int *lalloc_p, int *nlines_p, int *ret_starti);
static int first_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA **ret_msa, int **ret_ltype);
static int other_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA      *msa, int      *ltype);
static int append_selex_block(ESL_MSA *msa, char **line, int *ltype, int *lpos, int *rpos, int nlines);

/* read_selex()
 * Read an alignment in SELEX format.
 * SRE, Mon Dec 29 10:19:32 2008 [Pamplona]
 *
 * Purpose:  Parse an alignment from an open SELEX format alignment
 *           file <afp>, returning the alignment in <ret_msa>.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *
 *            Returns <eslEFORMAT> if parse fails because of a file
 *            format problem. 
 *            Returns <eslEOF> if no alignment is found in the file.
 *            Returns <eslEINVAL> if we're trying to read a digital
 *            alignment, and an invalid residue is found that 
 *            can't be digitized.
 *
 *            On all normal error conditions, <afp->errbuf> contains
 *            an informative error message for the user, and the 
 *            <*ret_msa> is <NULL>. The error message looks like
 *            "parse failed (line 156): too many #=SS lines for seq"
 *            The caller can prefix with filename if it likes.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
static int
read_selex(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA  *msa     = NULL;
  char    **line    = NULL;
  int      *ltype   = NULL;
  int      *llen    = NULL;
  int      *lpos    = NULL;
  int      *rpos    = NULL;
  int       lalloc  = 0;
  int       nlines  = 0;
  int       nblocks = 0;
  int       starti;
  int       i, apos;
  int       status;

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* For each alignment block: */
  while ( (status = read_block(afp, &line, &llen, &lpos, &rpos, &lalloc, &nlines, &starti)) == eslOK)
    { /* now line[0..nlines-1] are data lines; llen[0..nlines-1] are max line lengths exc. \0 */
      /* lpos[0..nlines-1] are 0; and rpos[0..nlines-1] are idx of last nonwhitespace char on lines */
      nblocks++;

      if (nblocks == 1) status = first_selex_block(afp->errbuf, starti, line, lpos, rpos, nlines, &msa, &ltype);
      else              status = other_selex_block(afp->errbuf, starti, line, lpos, rpos, nlines,  msa,  ltype);
      if (status != eslOK) goto ERROR;

      if ((status = append_selex_block(msa, line, ltype, lpos, rpos, nlines)) != eslOK) goto ERROR;
    }
  if (status != eslEOF || nblocks == 0) goto ERROR; 

#ifdef eslAUGMENT_SSI
  msa->offset = 0; /* SELEX files are single MSA only; offset is always 0. */
#endif

  /* SELEX format allows ' ' as gaps, but easel doesn't */
  if (msa->rf != NULL)  
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->rf[apos] == ' ') msa->rf[apos] = '.';
  if (msa->ss_cons != NULL)  
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->ss_cons[apos] == ' ') msa->ss_cons[apos] = '.';
  if (msa->ss != NULL) 
    for (i = 0; i < msa->nseq; i++)
      if (msa->ss[i] != NULL) 
	for (apos = 0; apos < msa->alen; apos++)
	  if (msa->ss[i][apos] == ' ') msa->ss[i][apos] = '.';
  if (msa->sa != NULL) 
    for (i = 0; i < msa->nseq; i++)
      if (msa->sa[i] != NULL) 
	for (apos = 0; apos < msa->alen; apos++)
	  if (msa->sa[i][apos] == ' ') msa->sa[i][apos] = '.';
  for (i = 0; i < msa->nseq; i++)
    for (apos = 0; apos < msa->alen; apos++)
      if (msa->aseq[i][apos] == ' ') msa->aseq[i][apos] = '.';

#ifdef eslAUGMENT_ALPHABET 
  if (afp->do_digital) status = esl_msa_Digitize(afp->abc, msa, afp->errbuf);
#endif  
    
  *ret_msa = msa;
  free(ltype);
  return eslOK;
  
 ERROR:
  if (msa     != NULL) esl_msa_Destroy(msa);
  if (ltype   != NULL) free(ltype); 
  return status;
}

/* read_block()
 * Read one block of alignment data into memory.
 *
 * If we're trying to read the *first* block in an alignment, then at start:
 *   No lines of the alignment file have been read into <afp->buf> yet. 
 *   *line_p, *llen_p, *lpos_p, *rpos_p are all NULL
 *   *lalloc_p is 0
 *   *nlines_p is 0
 * On success, returns <eslOK> (even if last line is end of file), and:
 *   afp->buf either contains a blank line (immediately after block end), or <afp> is at EOF.
 *   *nlines_p points to the number of lines stored
 *   *line_p points to line[0..nlines-1][] array of \0 terminated strings
 *   *llen_p points to llen[0..nlines-1] array of string allocations in chars, not including \0
 *   *lpos_p points to lpos[0..nlines-1] array, all initialized to 0
 *   *rpos_p points to rpos[0..nlines-1] array, all initialized to idx of last nonwhitespace char on line
 *   *lalloc_p is >= *nlines
 * If file is empty (no data), returns <eslEOF>.
 * If an allocation fails at any point, throw <eslEMEM>. 
 *
 * If we are trying to read a *subsequent* block in the file, at start:
 *   <afp->buf> is where a previous read left it: on a blank line, or <afp> is at EOF.
 *   *line_p points to previous line[][] array; we'll reuse it, reallocating if needed.
 *   *llen_p points to previous llen[] array; ditto
 *   *lpos_p points to lpos[] array: all 0, same as first block
 *   *rpos_p points to rpos[] array; idx of last nonwhitespace char on line
 *   *lalloc_p is >0, and is what the previous read_block call reported
 *   *nlines_p points to the number of lines we saw in the *first* block.
 * On success, returns <eslOK> as above.
 * If the number of lines seen in the block doesn't match the expected number, return <eslEFORMAT>.
 * If no more data remain in the file, return <eslEOF>.
 * If an allocation fails at any point, throw <eslEMEM>. 
 *   
 * Memory for line[][] and llen[] are entirely managed here - not by caller. They are
 * initially allocated on the first block (*lalloc_p == 0); reallocated
 * as needed; and free'd when we reach <EOF> or an error is detected.
 */
static int
read_block(ESL_MSAFILE *afp, char ***line_p, int **llen_p, int **lpos_p, int **rpos_p, int *lalloc_p, int *nlines_p, int *ret_starti)
{
  void  *tmpp;
  char **line;
  int   *llen;
  int   *lpos;
  int   *rpos;
  char  *s;
  int    lalloc;
  int    blen;
  int    nlines = 0;
  int    i;
  int    starti;
  int    status;
  
  afp->errbuf[0] = '\0';
  if (*lalloc_p == 0) 		/* first block? allocate. */
    {
      ESL_ALLOC(line, sizeof(char *) * 16);
      ESL_ALLOC(llen, sizeof(int)    * 16);
      ESL_ALLOC(lpos, sizeof(int)    * 16);
      ESL_ALLOC(rpos, sizeof(int)    * 16);
      for (i = 0; i < 16; i++) 
	{ line[i] = NULL; llen[i] = 0; lpos[i] = 0; rpos[i] = 0; }
      lalloc = 16;
    }
  else 				/* second or later block? reuse existing arrays   */
    {
      line   = *line_p;
      llen   = *llen_p;
      lpos   = *lpos_p;
      rpos   = *rpos_p;
      lalloc = *lalloc_p;
    }
  
  /* Advance 'til afp->buf contains first line of the block. */
  do { status = msafile_getline(afp); }
  while (status == eslOK && 
	 (is_blankline(afp->buf) || 
	  (*afp->buf == '#' && (strncmp(afp->buf, "#=", 2) != 0))));
  if      (status == eslEOF && *lalloc_p == 0) ESL_XFAIL(eslEOF, afp->errbuf, "parse failed: no alignment data found");
  else if (status != eslOK)                    goto ERROR; /* includes true (normal) EOF, EMEM */

  starti =  afp->linenumber;

  do {
    if (nlines == lalloc) 
      {
	ESL_RALLOC(line, tmpp, sizeof(char *) * lalloc * 2);
	ESL_RALLOC(llen, tmpp, sizeof(int)    * lalloc * 2);
	ESL_RALLOC(lpos, tmpp, sizeof(int)    * lalloc * 2);
	ESL_RALLOC(rpos, tmpp, sizeof(int)    * lalloc * 2);
	for (i = lalloc; i < lalloc*2; i++) { line[i] = NULL; llen[i] = 0; lpos[i] = 0; rpos[i] = 0; }
	lalloc*=2;
      }

    blen = strlen(afp->buf);
    if (blen > llen[nlines]) ESL_RALLOC(line[nlines], tmpp, sizeof(char) * (blen+1)); /* +1 for \0 */
    strcpy(line[nlines], afp->buf);
    llen[nlines] = blen;

    /* rpos is most efficiently determined in read_block() rather than elsewhere, 
     *  because we know blen here; saves a strlen() elsewhere.
     */
    for (s = line[nlines]+blen-1; isspace(*s); s--) ;
    rpos[nlines] = s - line[nlines];

    nlines++;

    do { status = msafile_getline(afp); }
    while (status == eslOK && (*afp->buf == '#' && strncmp(afp->buf, "#=", 2) != 0)); /* skip comments */
  } while (status == eslOK && ! is_blankline(afp->buf));

  if (status != eslOK && status != eslEOF) goto ERROR; /* EMEM */
  if (*lalloc_p != 0 && *nlines_p != nlines) 
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): expected %d lines in block, saw %d", afp->linenumber, *nlines_p, nlines);

  *line_p     = line;
  *llen_p     = llen; 
  *lpos_p     = lpos;
  *rpos_p     = rpos;
  *lalloc_p   = lalloc;
  *nlines_p   = nlines;
  *ret_starti = starti;
  return eslOK; /* an EOF is turned into OK the first time we see it: so last block read gets dealt with */

 ERROR:  /* includes final EOF, when we try to read another block and fail.  */
  if (line != NULL) { 
    for (i = 0; i < lalloc; i++) 
      if (line[i] != NULL) free(line[i]); 
    free(line);
  }
  if (llen   != NULL) free(llen);
  if (lpos   != NULL) free(lpos);
  if (rpos   != NULL) free(rpos);
  *line_p     = NULL;
  *llen_p     = NULL;
  *lpos_p     = NULL;
  *rpos_p     = NULL;
  *lalloc_p   = 0;
  *nlines_p   = 0;
  *ret_starti = 0;
  return status;	
}

/* First block: determine and store line types, in ltype[0..nlines-1].
 * From that, we know the number of sequences, nseq.
 * From that, we can allocate a new MSA object for <nseq> sequences.
 * Then parse we store all the sequence names in msa->sqname[].
 * This gives us information we will use to validate subsequent blocks,
 * making sure they contain exactly the same line order.
 *
 * We also set lpos[], rpos[] here to the position of the leftmost,
 * rightmost non-whitespace sequence residue character.
 * In the special case of lines with all whitespace data (which SELEX 
 * format allows!), set both lpos[] = -1; we'll catch this
 * as a special case when we need to.
 * 
 * <msa> and <ltype> are allocated here, and must be free'd by caller.
 */
static int
first_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA **ret_msa, int **ret_ltype)
{
  ESL_MSA *msa    = NULL;
  int     *ltype  = NULL;
  int      nseq   = 0;
  int      nrf, ncs, nss, nsa;
  int      has_ss, has_sa;
  int      li, i;
  char    *s, *tok;
  int      n;
  int      status;

  if (errbuf != NULL) errbuf[0] = '\0';

  /* Determine ltype[]; count sequences */
  ESL_ALLOC(ltype, sizeof(int) * nlines);
  nrf = ncs = nss = nsa = 0;
  has_ss = has_sa = FALSE;
  for (nseq = 0, li = 0; li < nlines; li++)
    {
      if      (strncmp(line[li], "#=RF", 4) == 0) { ltype[li] = eslMSA_LINE_RF; nrf++; }
      else if (strncmp(line[li], "#=CS", 4) == 0) { ltype[li] = eslMSA_LINE_CS; ncs++; }
      else if (strncmp(line[li], "#=SS", 4) == 0) { ltype[li] = eslMSA_LINE_SS; nss++;  has_ss = TRUE; }
      else if (strncmp(line[li], "#=SA", 4) == 0) { ltype[li] = eslMSA_LINE_SA; nsa++;  has_sa = TRUE; }
      else                                        { ltype[li] = eslMSA_LINE_SQ; nseq++; nss = nsa = 0; }
      if (nss > 0 && nseq==0) ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SS must follow a sequence", li+starti);
      if (nsa > 0 && nseq==0) ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SA must follow a sequence", li+starti);
      if (nrf > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=RF lines for block", li+starti);
      if (ncs > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=CS lines for block", li+starti);
      if (nss > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=SS lines for seq",   li+starti);
      if (nsa > 1)            ESL_XFAIL(eslEFORMAT, errbuf, "parse failed (line %d): too many #=SA lines for seq",   li+starti);
    }

  /* Allocate the MSA, now that we know nseq */
  if ((msa = esl_msa_Create(nseq, -1)) == NULL) { status = eslEMEM; goto ERROR; } 
  if (has_ss) 
    {
      ESL_ALLOC(msa->ss, sizeof(char *) * nseq); 
      for (i = 0; i < nseq; i++) msa->ss[i] = NULL;
    }
  if (has_sa) 
    {
      ESL_ALLOC(msa->sa, sizeof(char *) * nseq);
      for (i = 0; i < nseq; i++) msa->sa[i] = NULL;
    }
  msa->nseq = nseq;
  msa->alen = 0;
  /* msa->aseq[], msa->sqname[], msa->ss[], msa->sa[] arrays are all ready (all [i] are NULL) */
  
  for (i = 0, li = 0; li < nlines; li++)
    if (ltype[li] == eslMSA_LINE_SQ)
      {
	s = line[li];
	if (esl_strtok_adv(&s, " \t\n\r", &tok, &n, NULL)     != eslOK) ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");
	if ((status = esl_strdup(tok, -1, &(msa->sqname[i]))) != eslOK) goto ERROR;

	while (*s && isspace(*s)) s++;   /* advance s to first residue */
	lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	i++;
      }
    else
      {
	for (s = line[li]; *s && !isspace(*s); s++) ; /* advance s past #=XX tag       */
	for (           ;  *s &&  isspace(*s); s++) ; /* advance s to first residue    */
	lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
      }
  *ret_msa   = msa;
  *ret_ltype = ltype;
  return eslOK;

 ERROR:
  if (msa   != NULL) esl_msa_Destroy(msa);
  if (ltype != NULL) free(ltype);
  *ret_msa   = NULL;
  *ret_ltype = NULL;
  return status;
}


/* Subsequent blocks: 
 * validate that lines are coming in same order as first block (including sqname);
 * set lpos[] as in first_selex_block().
 */
static int
other_selex_block(char *errbuf, int starti, char **line, int *lpos, int *rpos, int nlines, ESL_MSA *msa, int *ltype)
{
  char *s, *tok;
  int   i, li;

  /* Compare order of line types. */
  for (li = 0; li < nlines; li++)
    {
      if      (strncmp(line[li], "#=RF", 4) == 0) { if (ltype[li] != eslMSA_LINE_RF) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=RF line isn't in expected order", li+starti); }
      else if (strncmp(line[li], "#=CS", 4) == 0) { if (ltype[li] != eslMSA_LINE_CS) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=CS line isn't in experted order", li+starti); }
      else if (strncmp(line[li], "#=SS", 4) == 0) { if (ltype[li] != eslMSA_LINE_SS) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SS line isn't in expected order", li+starti); }
      else if (strncmp(line[li], "#=SA", 4) == 0) { if (ltype[li] != eslMSA_LINE_SA) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): #=SA line isn't in expected order", li+starti); }
      else                                        { if (ltype[li] != eslMSA_LINE_SQ) ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): seq line isn't in expected order",  li+starti); }
    }

  /* Compare order of sequence names, and set lpos[]. */
  for (i = 0, li = 0; li < nlines; li++)
    {
      if (ltype[li] == eslMSA_LINE_SQ)
	{
	  s = line[li];
	  if (esl_strtok(&s, " \t\n\r", &tok) != eslOK) ESL_EXCEPTION(eslEINCONCEIVABLE, "can't happen");
	  if (strcmp(tok, msa->sqname[i])     != 0)     ESL_FAIL(eslEFORMAT, errbuf, "parse failed (line %d): expected seq %s, saw %s",  li+starti, msa->sqname[i], tok);
	  
	  while (*s && isspace(*s)) s++;              /* advance s to first residue */
	  lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	  i++;
	}
      else
	{
	  for (s = line[li]; *s && !isspace(*s); s++) ; /* advance s past #=XX tag    */
	  for (            ; *s &&  isspace(*s); s++) ; /* advance s to first residue */
	  lpos[li] = ((*s == '\0') ? -1 : s-line[li]);
	}
    }
  return eslOK;
}


static int 
append_selex_block(ESL_MSA *msa, char **line, int *ltype, int *lpos, int *rpos, int nlines)
{
  void *tmpp;
  char *s;
  int   leftmost, rightmost;
  int   li, i;
  int   nleft, ntext, nright;
  int   nadd;
  int   pos;
  int   status;

  /* Determine rightmost, leftmost columns for data */
  /* Watch out for special case of empty data lines: lpos= -1 flag */
  /* Watch out for extra special case where *no* line on block has data! */
  leftmost  = INT_MAX;
  rightmost = -1;
  for (li = 0; li < nlines; li++) {
    leftmost  = (lpos[li] == -1) ? leftmost  : ESL_MIN(leftmost,  lpos[li]);
    rightmost = (lpos[li] == -1) ? rightmost : ESL_MAX(rightmost, rpos[li]);
  }
  if (rightmost == -1) return eslOK; /* extra special case: no data in block at all! */
  nadd = rightmost - leftmost + 1; /* block width in aligned columns */

  for (i = 0, li = 0; li < nlines; li++)
    {
      nleft  = ((lpos[li] != -1) ? lpos[li] - leftmost     : nadd); /* watch special case of all whitespace on data line, lpos>rpos */
      ntext  = ((lpos[li] != -1) ? rpos[li] - lpos[li] + 1 : 0);
      nright = ((lpos[li] != -1) ? rightmost - rpos[li]    : 0);

      if      (ltype[li] == eslMSA_LINE_SQ) { ESL_RALLOC(msa->aseq[i], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->aseq[i]; i++; }
      else if (ltype[li] == eslMSA_LINE_RF) { ESL_RALLOC(msa->rf,      tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->rf;      }
      else if (ltype[li] == eslMSA_LINE_CS) { ESL_RALLOC(msa->ss_cons, tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->ss_cons; }
      else if (ltype[li] == eslMSA_LINE_SS) { ESL_RALLOC(msa->ss[i-1], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->ss[i-1]; }
      else if (ltype[li] == eslMSA_LINE_SA) { ESL_RALLOC(msa->sa[i-1], tmpp, sizeof(char) * (msa->alen + nadd + 1)); s = msa->sa[i-1]; }

      for (pos = msa->alen; pos < msa->alen+nleft; pos++) s[pos] = ' ';
      if (ntext > 0) memcpy(s+msa->alen+nleft, line[li]+lpos[li], sizeof(char) * ntext);
      for (pos = msa->alen+nleft+ntext; pos < msa->alen+nadd; pos++) s[pos] = ' ';
      s[pos] = '\0';
    }
  msa->alen += nadd;
  return eslOK;

 ERROR:
  return status;
}
/*------------------- end, selex format -------------------------*/



/*****************************************************************
 * 11. AFA (aligned FASTA) format
 *****************************************************************/

/* write_afa()
 * EPN, Mon Nov  2 15:55:10 2009
 *
 * Purpose:   Write alignment <msa> in aligned FASTA format to a 
 *            stream <fp>.
 *            
 *            If <msa> is in text mode, residues and gaps are written
 *            as they exist in the data structure. If <msa> is
 *            digitized, all residues are written in uppercase, all
 *            gaps as -.
 *            
 * Args:      fp     - open output stream
 *            msa    - MSA to write       
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
static int
write_afa(FILE *fp, ESL_MSA *msa)
{
  int     i;
  int64_t pos;
  char    buf[61];
  int     acpl;       /* actual number of character per line */
  
  for (i = 0; i < msa->nseq; i++)
    {
      /* Construct the description line */
      fprintf(fp, ">%s", msa->sqname[i]);
      if (msa->sqacc  != NULL && msa->sqacc[i]  != NULL) fprintf(fp, " %s", msa->sqacc[i]);
      if (msa->sqdesc != NULL && msa->sqdesc[i] != NULL) fprintf(fp, " %s", msa->sqdesc[i]);
      fputc('\n', fp);

      pos = 0;
      while (pos < msa->alen)
	{
	  acpl = (msa->alen - pos > 60)? 60 : msa->alen - pos;
#ifdef eslAUGMENT_ALPHABET
	  if (msa->flags & eslMSA_DIGITAL)
	    esl_abc_TextizeN(msa->abc, msa->ax[i] + pos + 1, acpl, buf);
	  else 
#endif
	    strncpy(buf, msa->aseq[i] + pos, acpl);

	  buf[acpl] = '\0';
	  fprintf(fp, "%s\n", buf);	      
	  pos += 60;
	}
    } /* end, loop over sequences in the MSA */

  return eslOK;
}


/* read_afa()
 * EPN, Mon Nov  2 17:24:40 2009
 *
 * Purpose:   Parse the one-and-only alignment from an open AFA (aligned
 *            fasta) format alignment file <afp>, leaving the
 *            alignment in <ret_msa>.
 *
 *            The current implementation reads the file one line at a
 *            time. Blank lines are skipped. Lines with '>' as the
 *            first non-whitespace character begin a new sequence,
 *            first word is sequence name, remainder of line is
 *            sequence description. All other lines are sequence lines
 *            currently processed one whitespace-delimited token at a
 *            time (to permit whitespace in the file).  A possibly
 *            more efficient route would be to read each complete
 *            sequence at a time (since AFA is not interleaved) and
 *            write it to the msa in a single step.
 *            
 *            Starting with the second sequence, all sequence lengths
 *            are confirmed to be identical to the length of the first
 *            sequence. If any are not, afp->errbuf is filled, <ret_msa>
 *            is set to NULL and <eslEINVAL> Is returned.
 *
 * Returns:   <eslOK> on success, and the alignment is in <ret_msa>.
 *            If no sequences exist, return <eslEOF> and <ret_msa>
 *            is set to NULL.
 * 
 *            <eslEFORMAT> if parse fails because of a file format
 *            problem, in which case afp->errbuf is set to contain a
 *            formatted message that indicates the cause of the
 *            problem, and <ret_msa> is set to NULL.
 *            
 *            Returns <eslEINVAL> if we're trying to read a digital
 *            alignment, and an invalid residue is found that can't be
 *            digitized.
 * 
 * Throws:    <eslEMEM> on allocation error.
 *
 */
static int
read_afa(ESL_MSAFILE *afp, ESL_MSA **ret_msa)
{
  ESL_MSA   *msa = NULL;
  char      *s;
  int        status;
  int        status2;
  int        seqidx;
  char      *seqname;
  char      *desc;
  char      *text;
  int        len, i;
  char       errbuf2[eslERRBUFSIZE];

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Initialize allocation of the MSA:
   * make it growable, by giving it an initial blocksize of
   * 16 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }

#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }


#ifdef eslAUGMENT_SSI
  /* EPN: not sure if this is appropriate/necessary, we assume only one alignment in AFA files */
  msa->offset = ftello(afp->f);
#endif

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++; /* skip leading whitespace */

      if (*s == '\n' || *s == '\r') continue; /* skip blank lines */
      if (*s == '>') { /* header line */
	/* if nec, make room for the new seq */
	if (msa->nseq >= msa->sqalloc && (status = esl_msa_Expand(msa)) != eslOK) return status; 

	/* store the name (space delimited) */
	s++; /* move past the '>' */
	seqidx = msa->nseq;
	msa->nseq++;
	if (esl_strtok(&s, " \t\n\r", &seqname) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, problem reading name of sequence %d at line %d\n", seqidx+1, afp->linenumber);
	status = esl_strdup(seqname, -1, &(msa->sqname[seqidx]));

	status = esl_strtok(&s, "\n\r", &desc);
	if     (status == eslOK) status = esl_msa_SetSeqDescription(msa, seqidx, desc);
	else if(status != eslEOL) ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, problem reading description of sequence %d at line %d\n", seqidx, afp->linenumber);
	/* else, no description */

	if((seqidx > 1) && (msa->sqlen[(seqidx-1)] != msa->sqlen[0])) { /* make sure the aligned seq we just read is the same length as the first seq we read */
	  ESL_XFAIL(eslEFORMAT, afp->errbuf, "sequence %d length (%" PRId64 ") is not equal to the expected length (%" PRId64 ") (the length of first seq in file)", seqidx, msa->sqlen[(seqidx-1)], msa->sqlen[0]);
	}
      }
      else {  /* not a '>' */
	if(msa->nseq == 0) { /* shouldn't happen, we haven't yet seen a '>' */
	  ESL_XFAIL(eslEFORMAT, afp->errbuf, "AFA MSA parse error, first non-whitespace character is not a '>' at line %d\n", afp->linenumber);	  
	}
	/* A sequence line: it doesn't begin with, but may contain, whitespace (' ' or '\t').
	 * We add whitespace-delimited tokens one at a time to the aseq (or ax).
	 * (Note: if we're digitized I think we could use a single call to esl_abc_dsqcat()
	 *  instead of splitting into tokens, which may be slightly more efficient).
	 */
	while(esl_strtok_adv(&s, " \t\n", &text, &len, NULL) == eslOK) 
	  { 
#ifdef eslAUGMENT_ALPHABET
	    if (msa->flags & eslMSA_DIGITAL)
	      {
		if((status = esl_abc_dsqcat(msa->abc, &(msa->ax[seqidx]), &(msa->sqlen[seqidx]), text, len)) != eslOK) {
		  /* invalid char(s), get informative error message */
		  if (esl_abc_ValidateSeq(msa->abc, text, len, afp->errbuf) != eslOK) 
		    ESL_XFAIL(eslEFORMAT, errbuf2, "%s (line %d): %s", msa->sqname[i], afp->linenumber, afp->errbuf);
		}
	      }
#endif
	  if (! (msa->flags & eslMSA_DIGITAL))
	    {
	      status = esl_strcat(&(msa->aseq[seqidx]), msa->sqlen[seqidx], text, len);
	      msa->sqlen[seqidx] += len;
	    } 
	  }
      }
    }

  /* check the length of the final sequence */
  if((msa->nseq > 1) && (msa->sqlen[seqidx] != msa->sqlen[0])) { /* make sure the aligned seq we just read is the same length as the first seq we read */
    ESL_XFAIL(eslEINVAL, afp->errbuf, "sequence %d length (%" PRId64 ") is not equal to the expected length (%" PRId64 ") (the length of first seq in file)", seqidx+1, msa->sqlen[seqidx], msa->sqlen[0]);
  }

  if (status2 == eslEMEM) ESL_XFAIL(status, afp->errbuf, "out of memory");
  if (status2 != eslEOF)  ESL_XFAIL(status, afp->errbuf, "unexpected error reading AFA alignment");

  /* Verify the msa */
  if (verify_parse(msa, afp->errbuf) != eslOK) { status = eslEFORMAT; goto ERROR; } 

  /* if alignment is empty set <ret_msa> to NULL and return eslEOF, (verification still works in this case) */
  if (msa->nseq == 0) { status = eslEOF; goto ERROR; }

  if (ret_msa != NULL) *ret_msa = msa; else esl_msa_Destroy(msa);
  return eslOK;

 ERROR:
  if (msa != NULL)      esl_msa_Destroy(msa);
  if (ret_msa != NULL) *ret_msa = NULL;
  return status;
}

/*---------------------- end, AFA format ------------------------*/

/******************************************************************************
 * 12. Memory efficient routines for PFAM format
 *****************************************************************************/

/* get_pp_idx
 *                   
 * Given a #=GR PP or #=GC PP_cons character, return the appropriate index
 * in a pp_ct[] vector. 
 * '0' return 0;
 * '1' return 1;
 * '2' return 2;
 * '3' return 3;
 * '4' return 4;
 * '5' return 5;
 * '6' return 6;
 * '7' return 7;
 * '8' return 8;
 * '9' return 9;
 * '*' return 10;
 * gap return 11;
 * 
 * Anything else (including missing or nonresidue) return -1;
 *
 * This mapping of PP chars to return values should probably be 
 * stored in some internal map structure somewhere, instead of 
 * only existing in this function as used by esl_msa_ReadNonSeqInfoPfam().
 */
int
get_pp_idx(ESL_ALPHABET *abc, char ppchar)
{
  if(esl_abc_CIsGap(abc, ppchar)) return 11;
  if(ppchar == '*')               return 10;
  if(ppchar == '9')               return 9;
  if(ppchar == '8')               return 8;
  if(ppchar == '7')               return 7;
  if(ppchar == '6')               return 6;
  if(ppchar == '5')               return 5;
  if(ppchar == '4')               return 4;
  if(ppchar == '3')               return 3;
  if(ppchar == '2')               return 2;
  if(ppchar == '1')               return 1;
  if(ppchar == '0')               return 0;
  return -1;
}

/* Function: esl_msa_ReadNonSeqInfoPfam()
 * Synopsis: Read non-sequence information in next Pfam formatted MSA.
 * Incept:   EPN, Sat Dec  5 07:56:42 2009
 *
 * Purpose:  Read the next alignment from an open Stockholm Pfam
 *           (non-interleaved, one line per seq) format alignment file
 *           <afp> and store all non per-sequence information
 *           (comments, \verb+#=GF+, \verb+#=GC+) in a new msa.
 *
 *           Note: this function could work on regular interleaved
 *           (non-Pfam) Stockholm if either (a) we didn't optionally
 *           return nseq or (b) we stored all seqnames in a keyhash,
 *           so we knew how many unique sequences there are.
 *
 *           We can't be as rigorous about validating the input as the
 *           other read functions that store the full alignment. Here,
 *           we only verify that there is only one line for the first
 *           sequence read and that's it (verifying that all sequences
 *           are only one line would require storing and looking up
 *           all sequence names which we could do at a cost).
 *
 *           Many optional return values (<opt_*>) make this function
 *           flexible and able to accomodate the diverse needs of the
 *           memory efficient enabled easel miniapps that use it
 *           (esl-alimerge, esl-alimask, esl-ssdraw). For any that are
 *           unwanted, pass <NULL>.
 *
 * Args:     afp           - open alignment file pointer
 *           abc           - alphabet to use, only nec and used if one 
 *                           of the opt_*_ct arrays is non-NULL
 *           known_alen    - known length of the alignment, -1 if unknown
 *                           must not be -1, if known_rf != NULL or
 *                           known_ss_cons != NULL.
 *           known_rf      - known RF annot. (msa->rf) for this alignment, 
 *                           might be known from prev call of this func,
 *                           for example. NULL if unknown.
 *           known_ss_cons - the known SS_cons annotation (msa->ss_cons) 
 *                           for this alignment, NULL if unknown.
 *           ret_msa       - RETURN: msa with comments, #=GC, #=GF annotation 
 *                           but no sequence info (nor #=GS,#=GR) 
 *                           pass NULL if not wanted
 *           opt_nseq      - optRETURN: number of sequences in msa 
 *           opt_alen      - optRETURN: length of first aligned sequence 
 *           opt_ngs       - optRETURN: number of #=GS lines in alignment 
 *           opt_maxname   - optRETURN: maximum seqname length 
 *           opt_maxgf     - optRETURN: maximum GF tag length
 *           opt_maxgc     - optRETURN: maximum GC tag length 
 *           opt_maxgr     - optRETURN: maximum GR tag length 
 *           opt_abc_ct    - optRETURN: [0..apos..alen-1][0..abc->K] 
 *                           per position count of each symbol in abc over all seqs
 *           opt_pp_ct     - optRETURN: [0..apos..alen-1][0..11], 
 *                           per position count of PPs over all seqs, 
 *                           [11] is gaps, [10] is '*', [0-9] are '0'-'9'
 *           opt_bp_ct     - optRETURN: [0..apos..alen-1][0..abc->K-1][0..abc->K-1]
 *                           per position count of each possible basepair 
 *                           in alignment, for pair apos:apos2, where 
 *                           apos < apos2 and apos:apos2 form a basepair 
 *                           in <known_ss_cons>. If non-NULL, <known_ss_cons> 
 *                           must be non-NULL.
 *           opt_srfpos_ct - optRETURN: [0..rfpos..rflen-1] per nongap RF 
 *                           position count of first nongap RF residue in 
 *                           each sequence, ex: opt_spos_ct[100] = x means x
 *                           seqs have their first nongap RF residue at rfpos 100
 *           opt_erfpos_ct - optRETURN: [0..rfpos..rflen-1] same as opt_srfpos_ct,
 *                           except for final nongap RF residue instead of first
 * 
 * Returns:  <eslOK> on success.  Returns <eslEOF> if there are no more
 *           alignments in <afp>, and <ret_msa> is set to NULL and
 *           <opt_*> are set to 0.
 *           <eslEFORMAT> if parse fails because of a file format
 *           problem, in which case <afp->errbuf> is set to contain a
 *           formatted message that indicates the cause of the
 *           problem, <ret_msa> is set to <NULL> and <opt_*> are set 
 *           to 0.
 *
 * Throws:    <eslEMEM> on allocation error.
 * 
 * Xref:      ~nawrockie/notebook/9_1206_esl_msa_mem_efficient/
 */
int
esl_msa_ReadNonSeqInfoPfam(ESL_MSAFILE *afp, ESL_ALPHABET *abc, int64_t known_alen, char *known_rf, char *known_ss_cons, ESL_MSA **ret_msa, 
			   int *opt_nseq, int64_t *opt_alen, int *opt_ngs, int *opt_maxname, int *opt_maxgf, int *opt_maxgc, int *opt_maxgr, 
			   double ***opt_abc_ct, int ***opt_pp_ct, double ****opt_bp_ct, int **opt_srfpos_ct, int **opt_erfpos_ct)
{
  char      *s;                    /* pointer to current character in afp */
  int        status;               /* easel status code */
  int        status2;              /* another easel status code */
  ESL_MSA   *msa = NULL;           /* the msa we're creating */
  int        nseq = 0;             /* number of sequences read */
  int64_t    alen = -1;            /* length of the alignment */
  int        ngs = 0;              /* number of #=GS lines read */
  int        maxname = 0;          /* max length seq name */
  int        maxgf = 0;            /* max length GF tag */
  int        maxgc = 0;            /* max length GC tag */
  int        maxgr = 0;            /* max length GR tag */
  char      *seqname;              /* ptr to a sequence name */
  int        namelen;              /* length of a sequence name */
  char      *first_seqname = NULL; /* name of first sequence read */
  char      *gc, *gr;              /* for storing "#=GC", "#=GR", temporarily */
  char      *tag;                  /* a GC or GR tag */
  int        taglen;               /* length of a tag */
  char      *text;                 /* text string */
  int        textlen;              /* length of text string */
  int        i, x;                 /* counters */
  int        j;                    /* position for a right half of a bp */
  int        apos;                 /* counter over alignment positions */
  int        rfpos;                /* counter over nongap RF positions */
  double   **abc_ct = NULL;        /* [0..alen-1][0..abc->K] per position count of each residue in abc and gaps over all seqs */
  double  ***bp_ct = NULL;         /* [0..alen-1][0..abc->Kp][0..abc->Kp], count of each possible base pair at each position, over all sequences, missing and nonresidues are *not counted* 
                                       base pairs are indexed by 'i' for a base pair between positions i and j, where i < j. */
  int        nppvals = 12;         /* '0'-'9' = 0-9, '*' = 10, gap = '11' */
  int      **pp_ct = NULL;         /* [0..alen-1][0..nppvals-1] per position count of each possible PP char over all seqs */
  int        ppidx;                /* index for 2nd dim of pp_ct array */
  int       *srfpos_ct = NULL;     /* [0..rflen-1] number of seqs that start (have first RF residue) at each RF position */
  int       *erfpos_ct = NULL;     /* [0..rflen-1] number of seqs that end   (have final RF residue) at each RF position */
  ESL_DSQ   *tmp_dsq = NULL;       /* temporary digitized sequence, only used if opt_abc_ct != NULL */
  int        rflen = 0;            /* non-gap RF residues, only used if known_rf != NULL */
  int       *a2rf_map = NULL;      /* [0..apos..known_alen-1] = rfpos, nongap RF position apos maps to, 
				    * -1 if apos is not a nongap RF position */
  int       *ct = NULL; 	   /* 0..known_alen-1 base pair partners array for known_ss_cons */
  char      *ss_nopseudo = NULL;   /* no-pseudoknot version of known_ss_cons */

  if(afp->format   != eslMSAFILE_PFAM) ESL_EXCEPTION(eslEINCONCEIVABLE, "only non-interleaved (1 line /seq, Pfam) Stockholm formatted files can be read in small memory mode");
  if(opt_abc_ct    != NULL && abc == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, abc == NULL, opt_abc_ct  != NULL");
  if(opt_pp_ct     != NULL && abc == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, abc == NULL, opt_pp_ct   != NULL");
  if(opt_bp_ct     != NULL && abc == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, abc == NULL, opt_bp_ct != NULL");
  if(opt_srfpos_ct != NULL && abc == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, abc == NULL, opt_srfpos_ct != NULL");
  if(opt_erfpos_ct != NULL && abc == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, abc == NULL, opt_erfpos_ct != NULL");
  if(opt_srfpos_ct != NULL && known_rf == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, known_rf == NULL, opt_srfpos_ct != NULL");
  if(opt_erfpos_ct != NULL && known_rf == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, known_rf == NULL, opt_erfpos_ct != NULL");
  if(opt_bp_ct     != NULL && known_ss_cons == NULL) ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, known_ss_cons == NULL, opt_bp_ct != NULL");
  if(known_rf      != NULL && known_alen == -1)      ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, known_rf != NULL, known_alen == -1");
  if(known_ss_cons != NULL && known_alen == -1)      ESL_FAIL(eslEINVAL, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam() contract violation, known_ss_cons != NULL, known_alen == -1");

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Preliminaries */
  /* if opt_srfpos_ct != NULL || opt_erfpos_ct != NULL, allocate them and determine RF positions */
  /* count non-gap RF columns */
  if(opt_srfpos_ct != NULL || opt_erfpos_ct != NULL) { 
    ESL_ALLOC(a2rf_map, sizeof(int) * known_alen);
    esl_vec_ISet(a2rf_map, known_alen, -1);
    for(apos = 0; apos < known_alen; apos++) { 
      if((! esl_abc_CIsGap(abc, known_rf[apos])) && 
	 (! esl_abc_CIsMissing(abc, known_rf[apos])) && 
	 (! esl_abc_CIsNonresidue(abc, known_rf[apos])))
	{ 
	  a2rf_map[apos] = rflen;
	  rflen++;
	}
    }
    ESL_ALLOC(srfpos_ct, sizeof(int) * rflen); 
    ESL_ALLOC(erfpos_ct, sizeof(int) * rflen);
    esl_vec_ISet(srfpos_ct, rflen, 0); 
    esl_vec_ISet(erfpos_ct, rflen, 0);   
  }
  /* if bp_ct != NULL, determine the ct array from the known_ss_cons, and allocate the bp_ct */
  if(opt_bp_ct != NULL) { /* contract enforces that if this is true, known_ss_cons != NULL and known_alen != -1 */
    /* get ct array which defines the consensus base pairs */
    ESL_ALLOC(ct,  sizeof(int)  * (known_alen+1));
    ESL_ALLOC(ss_nopseudo, sizeof(char) * (known_alen+1));
    esl_wuss_nopseudo(known_ss_cons, ss_nopseudo);
    if ((status = esl_wuss2ct(ss_nopseudo, known_alen, ct)) != eslOK) ESL_FAIL(status, afp->errbuf, "esl_msa_ReadNonSeqInfoPfam(), consensus structure string is inconsistent.");
    ESL_ALLOC(bp_ct,  sizeof(double **) * known_alen); 
    for(apos = 0; apos < known_alen; apos++) { 
      /* careful ct is indexed 1..alen, not 0..alen-1 */
      if(ct[(apos+1)] > (apos+1)) { /* apos+1 is an 'i' in an i:j pair, where i < j */
	ESL_ALLOC(bp_ct[apos], sizeof(double *) * (abc->Kp));
	for(x = 0; x < abc->Kp; x++) { 
	  ESL_ALLOC(bp_ct[apos][x], sizeof(double) * (abc->Kp));
	  esl_vec_DSet(bp_ct[apos][x], abc->Kp, 0.);
	}
      }
      else { /* apos+1 is not an 'i' in an i:j pair, where i < j, set to NULL */
	bp_ct[apos] = NULL;
      }
    }
  }
  /* end of preliminaries */

  /* Initialize allocation of the MSA:
   * We won't store any sequence information, so initial blocksize is
   * 0 seqs of 0 length.
   */
#ifdef eslAUGMENT_ALPHABET
  if (afp->do_digital == TRUE && (msa = esl_msa_CreateDigital(afp->abc, 16, -1))  == NULL) 
    { status = eslEMEM; goto ERROR; }
#endif
  if (afp->do_digital == FALSE && (msa = esl_msa_Create(16, -1))  == NULL)
    { status = eslEMEM; goto ERROR; }
  if (msa == NULL)    
    { status = eslEMEM; goto ERROR; }
  
  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR; /* includes EOF  */
  } while (is_blankline(afp->buf));
  
  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);
  
  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */
      
      if (*s == '#') {
	if(strncmp(s, "#=GF", 4) == 0)
	  {
	    if (ret_msa != NULL) { 
	      if ((status = parse_gf(msa, s)) != eslOK)
		ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GF line", afp->linenumber);
	    }
	    else if (opt_maxgf != NULL) { /* we need to parse out GF tag len to see if it is > maxgf */
	      s = afp->buf;
	      if (esl_strtok    (&s, " \t\n\r", &gc)                  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &tag,  &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
	      maxgf = ESL_MAX(maxgf, taglen); 
	    }
	  }

	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if  (ret_msa != NULL) { 
	      if  ((status = parse_gc(msa, s)) != eslOK)
		ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad #=GC line", afp->linenumber);
	    }
	    else if (opt_maxgc != NULL) { /* we need to parse out GC tag len to see if it is > maxgc */
	      s = afp->buf;
	      if (esl_strtok    (&s, " \t\n\r", &gc)                  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &tag,  &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	      maxgc = ESL_MAX(maxgc, taglen); 
	    }
	  }
	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    ngs++; 
	  }
	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if(opt_maxgr != NULL || opt_pp_ct != NULL) { 
	      s = afp->buf;
	      if (esl_strtok    (&s, " \t\n\r", &gr)                      != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &tag,      &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      maxgr = ESL_MAX(maxgr, taglen); 
	      if(opt_pp_ct != NULL) { 
		if (strncmp(tag, "PP", 2) == 0) { 
		  if (esl_strtok_adv(&s, " \t\n\r", &text, &textlen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR PP line", afp->linenumber);
		  /* verify, or set alignment length */
		  if(alen == -1) { /* first aligned text line, need to allocate pp_ct, and possibly abc_ct, srfpos_ct, erfpos_ct */
		    alen = textlen;
		    if(known_alen != -1 && known_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): known alen (%" PRId64 " passed in) != actual alen (%d)", afp->linenumber, known_alen, textlen);
		    ESL_ALLOC(pp_ct, sizeof(int *) * alen);
		    for(apos = 0; apos < alen; apos++) { 
		      ESL_ALLOC(pp_ct[apos], sizeof(int) * nppvals);
		      esl_vec_ISet(pp_ct[apos], nppvals, 0);
		    }
		    if(opt_abc_ct != NULL || opt_bp_ct != NULL) { 
		      ESL_ALLOC(tmp_dsq, (alen+2) * sizeof(ESL_DSQ));
		    }
		    if(opt_abc_ct != NULL) { 
		      ESL_ALLOC(abc_ct, sizeof(double *) * alen); 
		      for(apos = 0; apos < alen; apos++) { 
			ESL_ALLOC(abc_ct[apos], sizeof(double) * (abc->K+1));
			esl_vec_DSet(abc_ct[apos], (abc->K+1), 0.);
		      }
		    }
		  }
		  else if(alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR PP line, len %d, expected %" PRId64, afp->linenumber, textlen, alen);
		  for(apos = 0; apos < alen; apos++) { /* update appropriate PP count */
		    if((ppidx = get_pp_idx(abc, text[apos])) == -1) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR PP char: %c", afp->linenumber, text[apos]);
		    pp_ct[apos][ppidx]++;
		  }
		}
	      }
	    }
	  }
	else if (ret_msa != NULL && ((status = parse_comment(msa, s)) != eslOK)) { 
	  ESL_XFAIL(status, afp->errbuf, "parse failed (line %d): bad comment line", afp->linenumber);
	}
      } /* end of 'if (*s == '#')' */ 
      else if (strncmp(s, "//",   2) == 0)   break; /* normal way out */
      else if (*s == '\n' || *s == '\r')     continue;
      else { /* sequence line */
	if(opt_maxname != NULL || opt_alen != NULL || opt_abc_ct != NULL || opt_srfpos_ct != NULL || opt_erfpos_ct != NULL) { /* we need to parse out the seqname */
	  s = afp->buf;
	  if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad sequence line", afp->linenumber);
	  maxname = ESL_MAX(maxname, namelen);
	  if (opt_alen != NULL || opt_abc_ct != NULL || opt_srfpos_ct != NULL || opt_erfpos_ct != NULL) { /* we need to parse out the seq */
	    if (esl_strtok_adv(&s, " \t\n\r", &text, &textlen, NULL)  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad sequence line", afp->linenumber);
	    /* if first aligned seq read, store it's name, else see if it is an additional line of first aseq */
	    if(nseq == 0) { 
	      if ((status = esl_strdup(seqname, -1, &(first_seqname))) != eslOK) goto ERROR; 
	    }
	    else if(strcmp(first_seqname, seqname) == 0) { ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): two seqs with same name. Alignment may be in interleaved Stockholm. Reformat to Pfam with esl-reformat.", afp->linenumber); }
	    if(alen == -1) { /* first aligned text line, need to allocate pp_ct, and possibly abc_ct, srfpos_ct, erfpos_ct */
	      alen = textlen;
	      if(known_alen != -1 && known_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): known alen (%" PRId64 " passed in) != actual alen (%d)", afp->linenumber, known_alen, textlen);
	      if(opt_abc_ct != NULL || opt_bp_ct != NULL) { 
		ESL_ALLOC(tmp_dsq, (alen+2) * sizeof(ESL_DSQ));
	      }
	      if(opt_abc_ct != NULL) { 
		ESL_ALLOC(abc_ct, sizeof(double *) * alen); 
		for(apos = 0; apos < alen; apos++) { 
		  ESL_ALLOC(abc_ct[apos], sizeof(double) * (abc->K+1));
		  esl_vec_DSet(abc_ct[apos], (abc->K+1), 0.);
		}
	      }
	      if(opt_pp_ct != NULL) { 
		ESL_ALLOC(pp_ct, sizeof(int *) * alen);
		for(apos = 0; apos < alen; apos++) { 
		  ESL_ALLOC(pp_ct[apos], sizeof(int) * nppvals);
		  esl_vec_ISet(pp_ct[apos], nppvals, 0);
		}
	      }
	    }
	    else if(alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad aligned seq line, len %d, expected %" PRId64, afp->linenumber, textlen, alen);
	    if(opt_abc_ct != NULL || opt_bp_ct != NULL) { 
	      /* update appropriate abc and/or bp count. first, digitize the text */
	      if((status = esl_abc_Digitize(abc, text, tmp_dsq)) != eslOK) ESL_XFAIL(status, afp->errbuf, "small mem parse failed (line %d): problem digitizing sequence", afp->linenumber);
	    }
	    if(opt_abc_ct != NULL) { 
	      for(apos = 0; apos < alen; apos++) { /* update appropriate abc count, careful, tmp_dsq ranges from 1..alen (not 0..alen-1) */
		if((status = esl_abc_DCount(abc, abc_ct[apos], tmp_dsq[apos+1], 1.0)) != eslOK) ESL_XFAIL(status, afp->errbuf, "small mem parse failed (line %d): problem counting residue %d", afp->linenumber, apos+1);
	      }
	    }	    
	    if(opt_bp_ct != NULL) { 
	      for(apos = 0; apos < alen; apos++) { /* update appropriate abc count, careful, tmp_dsq ranges from 1..alen (not 0..alen-1) */
		if(bp_ct[apos] != NULL) { /* our flag for whether position (apos+1) is an 'i' in an i:j pair where i < j */
		  j = ct[apos+1] - 1; /* ct is indexed 1..alen */
		  bp_ct[apos][tmp_dsq[(apos+1)]][tmp_dsq[(j+1)]]++;
		}
	      }
	    }
	    if(opt_srfpos_ct != NULL) { 
	      for(apos = 0; apos < alen; apos++) { /* find first non-gap RF position */
		if((! esl_abc_XIsGap(abc, tmp_dsq[apos+1])) && ((rfpos = a2rf_map[apos]) != -1)) { 
		  srfpos_ct[rfpos]++; 
		  break;
		}
	      }
	    }
	    if(opt_erfpos_ct != NULL) { /* find final non-gap RF position */
	      for(apos = alen-1; apos >= 0; apos--) { 
		if((! esl_abc_XIsGap(abc, tmp_dsq[apos+1])) && ((rfpos = a2rf_map[apos]) != -1)) { 
		  erfpos_ct[rfpos]++;
		  break;
		}
	      }
	    }
	  }
	}
	nseq++; 
      }
    }

  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): didn't find // at end of alignment", afp->linenumber);

  /* if we're returning maxgc and an msa, determine maxgc, which we didn't do above b/c we parsed GC lines with parse_gc()
   * If msa != NULL, we already know maxgc */
  if(ret_msa != NULL && opt_maxgc != NULL) { 
    for(i = 0; i < msa->ngc; i++) maxgc = ESL_MAX(maxgc, strlen(msa->gc_tag[i])); 
    if (msa->rf      != NULL)     maxgc = ESL_MAX(maxgc, 2); /* 2 == strlen("RF") */
    if (msa->ss_cons != NULL)     maxgc = ESL_MAX(maxgc, 7); /* 7 == strlen("SS_cons") */
    if (msa->sa_cons != NULL)     maxgc = ESL_MAX(maxgc, 7); /* 7 == strlen("SA_cons") */
    if (msa->pp_cons != NULL)     maxgc = ESL_MAX(maxgc, 7); /* 7 == strlen("PP_cons") */
  }
  /* same as for maxgc, but now for maxgf */
  if(ret_msa != NULL && opt_maxgf != NULL) { 
    for(i = 0; i < msa->ngf; i++) maxgf = ESL_MAX(maxgf, strlen(msa->gf[i])); 
    if (msa->name != NULL) maxgf = ESL_MAX(maxgf, 2); /* 2 == strlen("ID") */
    if (msa->desc != NULL) maxgf = ESL_MAX(maxgf, 2); /* 2 == strlen("DE") */
    if (msa->acc  != NULL) maxgf = ESL_MAX(maxgf, 2); /* 2 == strlen("AC") */
    if (msa->au   != NULL) maxgf = ESL_MAX(maxgf, 2); /* 2 == strlen("AU") */
    if (msa->rf   != NULL) maxgf = ESL_MAX(maxgf, 2); /* 2 == strlen("RF") */
  }

  /* Note that we don't verify the parse, b/c we didn't read any sequence data, a verify_parse() would fail */


  if (first_seqname) free(first_seqname);
  if (tmp_dsq)       free(tmp_dsq);
  if (ct)            free(ct);
  if (ss_nopseudo)   free(ss_nopseudo);
  if (a2rf_map)      free(a2rf_map);

  if (ret_msa)       *ret_msa       = msa;       else if (msa)      esl_msa_Destroy(msa);
  if (opt_nseq)      *opt_nseq      = nseq; 
  if (opt_alen)      *opt_alen      = alen;
  if (opt_ngs)       *opt_ngs       = ngs;
  if (opt_maxname)   *opt_maxname   = maxname;
  if (opt_maxgf)     *opt_maxgf     = maxgf;
  if (opt_maxgc)     *opt_maxgc     = maxgc;
  if (opt_maxgr)     *opt_maxgr     = maxgr;
  if (opt_abc_ct)    *opt_abc_ct    = abc_ct;    else if (abc_ct)    esl_Free2D((void **) abc_ct, alen);
  if (opt_pp_ct)     *opt_pp_ct     = pp_ct;     else if (pp_ct)     esl_Free2D((void **) pp_ct, alen);
  if (opt_bp_ct)     *opt_bp_ct     = bp_ct;     else if (bp_ct)     esl_Free3D((void ***) bp_ct, known_alen, abc->Kp);
  if (opt_srfpos_ct) *opt_srfpos_ct = srfpos_ct; else if (srfpos_ct) free(srfpos_ct);
  if (opt_erfpos_ct) *opt_erfpos_ct = erfpos_ct; else if (erfpos_ct) free(erfpos_ct);
  return eslOK;

 ERROR:
  if (first_seqname)  free(first_seqname);
  if (tmp_dsq)        free(tmp_dsq);
  if (ct)             free(ct);
  if (ss_nopseudo)    free(ss_nopseudo);
  if (a2rf_map)       free(a2rf_map);

  if (msa)            esl_msa_Destroy(msa);
  if (pp_ct)          esl_Free2D((void **)  pp_ct,  alen);
  if (abc_ct)         esl_Free2D((void **)  abc_ct, alen);
  if (bp_ct)          esl_Free3D((void ***) bp_ct,  known_alen, abc->Kp);
  if (srfpos_ct)      free(srfpos_ct);
  if (erfpos_ct)      free(erfpos_ct);

  if (ret_msa)       *ret_msa       = NULL;
  if (opt_nseq)      *opt_nseq      = 0;
  if (opt_alen)      *opt_alen      = 0;
  if (opt_ngs)       *opt_ngs       = 0;
  if (opt_maxname)   *opt_maxname   = 0;
  if (opt_maxgf)     *opt_maxgf     = 0;
  if (opt_maxgc)     *opt_maxgc     = 0;
  if (opt_maxgr)     *opt_maxgr     = 0;
  if (opt_pp_ct)     *opt_pp_ct     = NULL;
  if (opt_abc_ct)    *opt_abc_ct    = NULL;
  if (opt_bp_ct)     *opt_bp_ct     = NULL;
  if (opt_srfpos_ct) *opt_srfpos_ct = NULL;
  if (opt_erfpos_ct) *opt_erfpos_ct = NULL;
  return status;
}

/* gapize_string
 *                   
 * Given a string, create a new one that is a copy of it, 
 * but with gaps added before each position (apos) as specified 
 * by ngapA[0..apos..len]. <gapchar> specifies the gap character
 * to add.
 * 
 * ngapA[0]    - number of gaps to add before first posn
 * ngapA[apos] - number of gaps to add before posn apos
 * ngapA[len]  - number of gaps to add after  final posn
 * 
 * ret_str is allocated here.
 *
 * Returns eslOK on success.
 *         eslEMEM on memory error.
 */
int 
gapize_string(char *src_str, int64_t src_len, int64_t dst_len, int *ngapA, char gapchar, char **ret_dst_str)
{
  int status;
  int src_apos = 0;
  int dst_apos  = 0;
  int i;
  char *dst_str;

  ESL_ALLOC(dst_str, sizeof(char) * (dst_len+1));
  dst_str[dst_len] = '\0';

  /* add gaps before first position */
  for(i = 0; i < ngapA[0]; i++) dst_str[dst_apos++] = gapchar;

  /* add gaps after every position */
  for(src_apos = 0; src_apos < src_len; src_apos++) { 
    dst_str[dst_apos++] = src_str[src_apos];
    for(i = 0; i < ngapA[(src_apos+1)]; i++) dst_str[dst_apos++] = gapchar;
  }

  *ret_dst_str = dst_str;
  return eslOK;

 ERROR: 
  return eslEMEM;
}

/* shrink_string
 *                   
 * Remove some characters of a string in place.
 * If useme[0..pos..len-1] == FALSE, remove position pos.
 * 
 */
void
shrink_string(char *str, const int *useme, int len)
{
  int opos, npos;

  for (opos = 0, npos = 0; opos < len; opos++) { 
      if (useme[opos] == FALSE) continue;
      str[npos++] = str[opos];
  }
  str[npos] = '\0';
  return;
}


/* determine_spacelen
 *                   
 * Determine number of consecutive ' ' characters 
 * in the string pointed to by s.
 */
int
determine_spacelen(char *s)
{
  int spacelen = 0;
  while (*s == ' ') { spacelen++; s++; } 
  return spacelen;
}

#ifdef eslAUGMENT_KEYHASH
/* Function: esl_msa_RegurgitatePfam()
 * Synopsis: Read and write next Pfam formatted MSA without storing it.
 * Incept:   EPN, Sun Dec  6 11:25:34 2009
 *
 * Purpose:  Read and immediately write each line of a MSA after
 *           optionally modifying aligned data by either adding all
 *           gap columns or removing some columns. Do this without
 *           creating an msa object, so memory usage is minimized.
 *
 *           The alignment file <afp> must be in Pfam format (1
 *           line/seq non-interleaved Stockholm). The <do_*> arguments
 *           specify which parts of the alignment to write.  <useme>
 *           specifies which positions to keep in aligned strings, if
 *           NULL then all are kept. <add2me> specifies how many gap
 *           characters to add after each aligned position, if NULL
 *           then none are added. Only one of <useme> and <add2me> 
 *           can be non-NULL. 
 * 
 *           If the <keepme> keyhash is non-NULL, it specifies the
 *           names of sequences (and affiliated annotation) to
 *           output. All others will be ignored. If <keepme> is NULL,
 *           all sequences will be regurgitated.
 *
 *           <maxname>, <maxgf>, <maxgc> and <maxgr> specify the max
 *           length sequence name, GF tag, GC tag, and GR tag, and can
 *           be provided by a caller that knows their values, e.g. as
 *           revealed by a previous call to <esl_msa_ReadNonSeqPfam()>.
 *           If any are -1, the caller didn't know the value, and the
 *           spacing in the alignment file we read in will be
 *           preserved. An example of useful non -1 values is if we're
 *           merging multiple alignments into a single alignment, and
 *           the spacing of any given alignment should change when all
 *           alignments are considered (like what esl-alimerge does).
 *
 *           We can't be as rigorous about validating the input as the
 *           other read functions that store the full alignment. Here,
 *           we verify that there is only one line for the first
 *           sequence read (verifying that all sequences are only one
 *           line would require storing and looking up all sequence
 *           names which we could do at a cost), that each aligned
 *           data line (aseq, GC, GR) are all the same length <alen>
 *           which should equal <exp_alen> unless <exp_alen> is -1,
 *           indicating the caller doesn't know what alen should be.
 *           If useme or add2me is non-<NULL>, <exp_alen> must not be -1.
 *           No validation of aligned sequence characters being part
 *           of an alphabet is done, though we could, at a cost in
 *           time.
 *
 * Args:     afp         - open alignment file pointer
 *           ofp         - output file pointer
 *           maxname     - maximum length of a sequence name (-1 if unknown) 
 *           maxgf       - maximum length of a GF tag (-1 if unknown) 
 *           maxgc       - maximum length of a GC tag (-1 if unknown) 
 *           maxgr       - maximum length of a GR tag (-1 if unknown) 
 *           do_header   - TRUE to write magic Stockholm header at top to ofp 
 *           do_trailer  - TRUE to write '//' at end to ofp
 *           do_blanks   - TRUE to regurgitate blank lines, FALSE not to
 *           do_comments - TRUE to write comments to ofp
 *           do_gf       - TRUE to write #=GF annotation to ofp
 *           do_gs       - TRUE to write #=GS annotation to ofp
 *           do_gc       - TRUE to write #=GC annotation to ofp
 *           do_gr       - TRUE to write #=GR annotation to ofp
 *           do_aseq     - TRUE to write aligned sequences to ofp
 *           seqs2regurg - keyhash of names of the sequences to write, all others
 *                         will not be written. Associated annotation (#=GS, #=GR) 
 *                         will be written for these sequences only. Must be NULL
 *                         if seqs2skip is non-NULL (enforced by contract).
 *                         If both are NULL all seqs are written.
 *           seqs2skip   - keyhash of names of the sequences to skip (not write), 
 *                         all others will be written. Associated annotation (#=GS, #=GR) 
 *                         will not be written for these sequences. Must be NULL
 *                         if seqs2regurg is NULL (enforced by contract).
 *                         If both are NULL all seqs are written.
 *           useme       - [0..apos..exp_alen-1] TRUE to include position apos in output of 
 *                         aligned data (GC,GR,aseq), FALSE to remove it, can be NULL
 *           add2me      - [0..apos..exp_alen-1] number of all gaps to add after each
 *                         position of aligned data (GC,GR,aseq), can be NULL
 *           exp_alen    - expected alignment length, -1 if unknown, which
 *                         is okay as long as useme == add2me == NULL
 *           gapchar     - gap character, only relevant if add2me != NULL
 * 
 * Returns:   <eslOK> on success. 
 *            Returns <eslEOF> if there are no more alignments in <afp>.
 *            <eslEFORMAT> if parse fails because of a file format problem,
 *            in which case afp->errbuf is set to contain a formatted message 
 *            that indicates the cause of the problem.
 *
 * Throws:    <eslEMEM> on allocation error.
 * 
 * Xref:      ~nawrockie/notebook/9_1206_esl_msa_mem_efficient/
 */
int
esl_msa_RegurgitatePfam(ESL_MSAFILE *afp, FILE *ofp, int maxname, int maxgf, int maxgc, int maxgr, 
			int do_header, int do_trailer, int do_blanks, int do_comments, int do_gf, 
			int do_gs, int do_gc, int do_gr, int do_aseq, 
			ESL_KEYHASH *seqs2regurg, ESL_KEYHASH *seqs2skip, int *useme, int *add2me, int exp_alen, char gapchar)
{
  char      *s = NULL;
  int        status;
  int        status2;
  int        nseq = 0;
  char      *seqname = NULL;
  char      *first_seqname = NULL;
  char      *text = NULL;
  char      *gapped_text = NULL;
  char      *tag = NULL;
  char      *gf = NULL;
  char      *gc = NULL;
  char      *gs = NULL;
  char      *gr = NULL;
  int       curmargin, curmargin2, namelen, spacelen, spacelen2, textlen, taglen;
  int       gaps2addlen;
  int       margin;           /* width of left hand side margin */
  int       flushpoint = 100; /* number of lines read at which to flush ofp */

  /* contract check */
  if(ofp == NULL) ESL_EXCEPTION(eslEINCONCEIVABLE, "ofp is NULL");
  if((afp->format != eslMSAFILE_STOCKHOLM) && (afp->format != eslMSAFILE_PFAM)) {
    ESL_EXCEPTION(eslEINCONCEIVABLE, "only non-interleaved (1 line /seq, Pfam) Stockholm formatted files can be read in small memory mode");
  }
  if(add2me != NULL && useme != NULL) { 
    ESL_EXCEPTION(eslEINCONCEIVABLE, "add2me and useme both non-NULL");
  }
  if((add2me != NULL || useme != NULL) && exp_alen == -1) { 
    ESL_EXCEPTION(eslEINCONCEIVABLE, "exp_alen == -1, but add2me or useme non-NULL");
  }
  if(seqs2regurg != NULL && seqs2skip != NULL) {
    ESL_EXCEPTION(eslEINVAL, "seqs2regurg and seqs2skip both non-NULL, only one may be");
  }

  gaps2addlen = (add2me == NULL) ? 0 : esl_vec_ISum(add2me, (exp_alen+1));

  margin = -1;
  if (maxgf != -1 && maxgf < 2) maxgf = 2;
  if (maxname != -1)                         margin = maxname+1; 
  if (maxgc > 0 && maxgc+6 > margin)         margin = maxgc+6;
  if (maxgr > 0 && maxname+maxgr+7 > margin) margin = maxname+maxgr+7; 
  /* if margin is still -1, we'll use the same spacing we read in from the file */

  if (feof(afp->f))  { status = eslEOF; goto ERROR; }
  afp->errbuf[0] = '\0';

  /* Check the magic Stockholm header line.
   * We have to skip blank lines here, else we perceive
   * trailing blank lines in a file as a format error when
   * reading in multi-record mode.
   */
  do {
    if ((status = msafile_getline(afp)) != eslOK) goto ERROR; /* includes EOF  */
  } while (is_blankline(afp->buf));
  
  if (strncmp(afp->buf, "# STOCKHOLM 1.", 14) != 0)
    ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): missing \"# STOCKHOLM\" header", afp->linenumber);
  if(do_header) fprintf(ofp, "%s", afp->buf);

  /* Read the alignment file one line at a time.
   */
  while ((status2 = msafile_getline(afp)) == eslOK) 
    {
      if(afp->linenumber % flushpoint == 0) fflush(ofp);
      s = afp->buf;
      while (*s == ' ' || *s == '\t') s++;  /* skip leading whitespace */
      
      if (*s == '#') {
	
	if      (strncmp(s, "#=GF", 4) == 0)
	  {
	    if (do_gf) { 
	      if(maxgf == -1) { /* just print line as is */
		fprintf(ofp, "%s", afp->buf); 
	      }
	      else { /* parse line into temporary strings, then print it out with correct formatting */
		s = afp->buf;
		if (esl_strtok(&s, " \t\n\r", &gf)   != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		if (esl_strtok(&s, " \t\n\r", &tag)  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		if (esl_strtok(&s, "\n\r",    &text) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		fprintf(ofp, "#=GF %-*s %s\n", maxgf, tag, text);
	      }
	    }
	  }
	else if (strncmp(s, "#=GC", 4) == 0)
	  {
	    if(do_gc) { 
	      /* parse line into temporary strings */
	      s = afp->buf;
	      if (esl_strtok    (&s, " \t\n\r", &gc)                  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &tag,  &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);
	      spacelen = determine_spacelen(s);
	      if (esl_strtok_adv(&s, " \t\n\r", &text, &textlen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line", afp->linenumber);

	      /* verify alignment length */
	      if(exp_alen == -1) exp_alen = textlen;
	      else if(exp_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GC line, len %d, expected %d", afp->linenumber, textlen, exp_alen);

	      /* determine margin length of sequence name field for output formatting */
	      curmargin = (margin == -1) ? taglen + spacelen : margin - 6; 

	      /* output, after optionally removing some characters (if useme != NULL) or adding gaps (if add2me != NULL) (contract enforces only one can be non-null) */
	      if(useme  != NULL) { 
		/* if this is a GC SS_cons line, remove broken basepairs first */
		if(strncmp(tag, "SS_cons", 7) == 0) {
		  if((status = remove_broken_basepairs_from_ss_string(text, afp->errbuf, textlen, useme)) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR SS line", afp->linenumber);
		}
		shrink_string(text, useme, exp_alen); /* this is done in place on text */
	      }
	      if(add2me != NULL) { 
		if((status = gapize_string(text, textlen, textlen + gaps2addlen, add2me, gapchar, &gapped_text)) != eslOK) goto ERROR; 
		fprintf(ofp, "#=GC %-*s %s\n", curmargin, tag, gapped_text);
		free(gapped_text);
	      }
	      else { 
		fprintf(ofp, "#=GC %-*s %s\n", curmargin, tag, text);
	      }
	    }
	  }
	else if (strncmp(s, "#=GS", 4) == 0)
	  {
	    /* we don't validate the sequence exists, this would require storing all seqnames */
	    if (do_gs) { 
	      if(maxname == -1 && seqs2regurg == NULL) { /* just print line as is */
		fprintf(ofp, "%s", afp->buf); 
	      }
	      else { /* parse line into temporary strings, then print it out with correct formatting */
		if(seqs2regurg == NULL || (status = esl_key_Lookup(seqs2regurg, seqname, NULL)) == eslOK) 
		  { /* this if() will evaluate as TRUE if seqs2regurg is NULL, or the seqname exists in seqs2regurg, else it will return FALSE */
		    s = afp->buf;
		    if (esl_strtok(&s, " \t\n\r", &gs)   != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		    if (esl_strtok(&s, " \t\n\r", &tag)  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		    if (esl_strtok(&s, "\n\r",    &text) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GF line", afp->linenumber);
		    fprintf(ofp, "#=GS %-*s %s\n", maxname, tag, text);
		  }
	      }
	    }
	  }
	else if (strncmp(s, "#=GR", 4) == 0)
	  {
	    if(do_gr) { 
	      /* parse line into temporary strings */
	      s = afp->buf;
	      if (esl_strtok    (&s, " \t\n\r", &gr)                      != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      spacelen = determine_spacelen(s);
	      if (esl_strtok_adv(&s, " \t\n\r", &tag,      &taglen, NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);
	      spacelen2 = determine_spacelen(s);
	      if (esl_strtok_adv(&s, " \t\n\r", &text, &textlen, NULL)   != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line", afp->linenumber);

	      /* verify alignment length */
	      if(exp_alen == -1) exp_alen = textlen;
	      else if(exp_alen != textlen) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR line, len %d, expected %d", afp->linenumber, textlen, exp_alen);

	      /* determine length of fields for output formatting */
	      curmargin  = (maxname == -1) ? namelen + spacelen : maxname; 
	      curmargin2 = (maxname == -1) ? taglen + spacelen2 : margin - maxname - 7;

	      /* if seqs2regurg != NULL, skip this sequence unless its in the seqs2regurg hash */
	      if(seqs2regurg == NULL || (status = esl_key_Lookup(seqs2regurg, seqname, NULL)) == eslOK) 
		{ /* this if() will evaluate as TRUE if seqs2regurg is NULL, or the seqname exists in seqs2regurg, else it will return FALSE */
		  
		  /* output, after optionally removing some characters (if useme != NULL) or adding gaps (if add2me != NULL) (contract enforces only one can be non-null) */
		  if(useme  != NULL) { 
		    /* if this is a GR SS line, remove broken basepairs first */
		    if(strncmp(tag, "SS", 2) == 0) {
		      if((status = remove_broken_basepairs_from_ss_string(text, afp->errbuf, textlen, useme)) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad #=GR SS line", afp->linenumber);
		    }
		    shrink_string(text, useme, exp_alen); /* this is done in place on text */
		  }
		  if(add2me != NULL) { 
		    if((status = gapize_string(text, textlen, textlen + gaps2addlen, add2me, gapchar, &gapped_text)) != eslOK) goto ERROR; 
		    fprintf(ofp, "#=GR %-*s %-*s %s\n", curmargin, seqname, curmargin2, tag, gapped_text);
		    free(gapped_text);
		  }
		  else { 
		    fprintf(ofp, "#=GR %-*s %-*s %s\n", curmargin, seqname, curmargin2, tag, text);
		  }
		}
	    }
	  }
	else if (do_comments) fprintf(ofp, "%s", afp->buf); /* print comment line, if desired */
      } /* end of 'if (*s == '#')' */ 
      else if (strncmp(s, "//",   2) == 0)   { if(do_trailer) fprintf(ofp, "%s", afp->buf); break; /* normal way out */ }
      else if (*s == '\n' || *s == '\r')     { if(do_blanks)  { fprintf(ofp, "%s", afp->buf); } continue; } 
      else { /* sequence line */
	if(do_aseq) { 
	  /* parse line into temporary strings */
	  s = afp->buf;
	  if (esl_strtok_adv(&s, " \t\n\r", &seqname, &namelen,  NULL) != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad sequence line", afp->linenumber);
	  spacelen = determine_spacelen(s);
	  if (esl_strtok_adv(&s, " \t\n\r", &text,    &textlen, NULL)  != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad sequence line", afp->linenumber);

	  /* verify alignment length */
	  if((exp_alen != -1) && (exp_alen != textlen)) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): bad seq line, len %d, expected %d", afp->linenumber, textlen, exp_alen);

	  /* determine length of sequence name field for output formatting */
	  curmargin = (margin == -1) ? namelen + spacelen : margin-1; 

	  /* make sure we haven't just read a second line of the first sequence in file (we must be in Pfam 1 line/seq file) */
	  if(nseq == 0) { if ((status = esl_strdup(seqname, -1, &(first_seqname))) != eslOK) goto ERROR; }
	  else if(strcmp(first_seqname, seqname) == 0) { ESL_XFAIL(eslEFORMAT, afp->errbuf, "parse failed (line %d): two seqs named %s. Alignment appears to be in Stockholm format. Reformat to Pfam with esl-reformat.", afp->linenumber, seqname); }
	  nseq++;

	  /* if seqs2regurg != NULL, skip this sequence unless its in the seqs2regurg hash */
	  /* seq usually follows another seq; guess msa->lastidx +1 */
	  if(seqs2regurg == NULL || (status = esl_key_Lookup(seqs2regurg, seqname, NULL)) == eslOK) 
	    { 
	      /* this if() will evaluate as TRUE if seqs2regurg is NULL, or the seqname exists in seqs2regurg, else it will return FALSE */

	      /* output, after optionally removing some characters (if useme != NULL) or adding gaps (if add2me != NULL) (contract enforces only one can be non-null) */
	      if(useme  != NULL) shrink_string(text, useme, exp_alen); /* this is done in place on text */
	      if(add2me != NULL) { 
		if((status = gapize_string(text, textlen, textlen + gaps2addlen, add2me, gapchar, &gapped_text)) != eslOK) goto ERROR; 
		fprintf(ofp, "%-*s %s\n", curmargin, seqname, gapped_text);
		free(gapped_text);
	      }
	      else { 
		fprintf(ofp, "%-*s %s\n", curmargin, seqname, text);
	      }
	    }
	}
      }
    }

  /* If we saw a normal // end, we would've successfully read a line,
   * so when we get here, status (from the line read) should be eslOK.
   */ 
  if (status2 != eslOK) ESL_XFAIL(eslEFORMAT, afp->errbuf, "small mem parse failed (line %d): didn't find // at end of alignment", afp->linenumber);
  if (first_seqname != NULL) free(first_seqname);
  return eslOK;

 ERROR:
  return status;
}
#endif

/*---------------- end, memory efficient Pfam routines  -------------------*/

/******************************************************************************
 * 13. Debugging/development routines.
 *****************************************************************************/

/* Function:  esl_msa_CreateFromString()
 * Synopsis:  Creates a small <ESL_MSA> from a test case string.
 * Incept:    SRE, Sat Nov 11 12:09:04 2006 [Janelia]
 *
 * Purpose:   A convenience for making small test cases in the test
 *            suites: given the contents of a complete multiple
 *            sequence alignment file as a single string <s> in
 *            alignment format <fmt>, convert it to an <ESL_MSA>.
 *            
 *            For example, 
 *            {\small\begin{verbatim}
 *            esl_msa_CreateFromString("# STOCKHOLM 1.0\n\nseq1 AAAAA\nseq2 AAAAA\n//\n", 
 *                                     eslMSAFILE_STOCKHOLM)
 *            \end{verbatim}}
 *            creates an ungapped alignment of two AAAAA sequences.
 *
 * Returns:   a pointer to the new <ESL_MSA> on success.
 *
 * Throws:    <NULL> if it fails to obtain, open, or read the temporary file
 *            that it puts the string <s> in.
 */
ESL_MSA *
esl_msa_CreateFromString(const char *s, int fmt)
{
  char         tmpfile[16] = "esltmpXXXXXX";
  FILE        *fp          = NULL;
  ESL_MSAFILE *mfp         = NULL;
  ESL_MSA     *msa         = NULL;

  if (esl_tmpfile_named(tmpfile, &fp)            != eslOK) goto ERROR;
  fprintf(fp, "%s", s);
  fclose(fp); 
  fp = NULL;
  if (esl_msafile_Open(tmpfile, fmt, NULL, &mfp) != eslOK) goto ERROR;
  if (esl_msa_Read(mfp, &msa)                    != eslOK) goto ERROR;

  esl_msafile_Close(mfp);
  remove(tmpfile);
  return msa;

 ERROR:
  if (fp  != NULL) fclose(fp);
  if (mfp != NULL) esl_msafile_Close(mfp);
  if (strcmp(tmpfile, "esltmpXXXXXX") != 0) remove(tmpfile);
  if (msa != NULL) esl_msa_Destroy(msa);                        
  return NULL;
}

/* Function:  esl_msa_Compare()
 * Synopsis:  Compare two MSAs for equality.
 * Incept:    SRE, Wed Jun 13 10:40:05 2007 [Janelia]
 *
 * Purpose:   Returns <eslOK> if the mandatory and optional contents
 *            of MSAs <a1> and <a2> are identical; otherwise return
 *            <eslFAIL>.
 *            
 *            Only mandatory and parsed optional information is
 *            compared. Unparsed Stockholm markup is not compared.
 */
int
esl_msa_Compare(ESL_MSA *a1, ESL_MSA *a2)
{
  if (esl_msa_CompareMandatory(a1, a2) != eslOK) return eslFAIL;
  if (esl_msa_CompareOptional(a1, a2)  != eslOK) return eslFAIL;
  return eslOK;
}

/* Function:  esl_msa_CompareMandatory()
 * Synopsis:  Compare mandatory subset of MSA contents.
 * Incept:    SRE, Wed Jun 13 09:42:56 2007 [Janelia]
 *
 * Purpose:   Compare mandatory contents of two MSAs, <a1> and <a2>.
 *            This comprises <aseq> (or <ax>, for a digital alignment);
 *            <sqname>, <wgt>, <alen>, <nseq>, and <flags>.
 *
 * Returns:   <eslOK> if the MSAs are identical; 
 *            <eslFAIL> if they are not.
 */
int
esl_msa_CompareMandatory(ESL_MSA *a1, ESL_MSA *a2)
{
  int i;

  if (a1->nseq  != a2->nseq)  return eslFAIL;
  if (a1->alen  != a2->alen)  return eslFAIL;
  if (a1->flags != a2->flags) return eslFAIL;

  for (i = 0; i < a1->nseq; i++)
    {
      if (strcmp(a1->sqname[i], a2->sqname[i])        != 0)     return eslFAIL;
      if (esl_DCompare(a1->wgt[i], a2->wgt[i], 0.001) != eslOK) return eslFAIL;
#ifdef eslAUGMENT_ALPHABET
      if ((a1->flags & eslMSA_DIGITAL) &&
	  memcmp(a1->ax[i], a2->ax[i], sizeof(ESL_DSQ) * (a1->alen+2)) != 0) 
	return eslFAIL;
#endif
      if (! (a1->flags & eslMSA_DIGITAL) && strcmp(a1->aseq[i], a2->aseq[i]) != 0) return eslFAIL;
    }
  return eslOK;
}

/* Function:  esl_msa_CompareOptional()
 * Synopsis:  Compare optional subset of MSA contents.
 * Incept:    SRE, Wed Jun 13 09:52:48 2007 [Janelia]
 *
 * Purpose:   Compare optional contents of two MSAs, <a1> and <a2>.
 *
 * Returns:   <eslOK> if the MSAs are identical; 
 *            <eslFAIL> if they are not.
 */
int
esl_msa_CompareOptional(ESL_MSA *a1, ESL_MSA *a2)
{
  int i;

  if (esl_CCompare(a1->name,    a2->name)    != eslOK) return eslFAIL;
  if (esl_CCompare(a1->desc,    a2->desc)    != eslOK) return eslFAIL;
  if (esl_CCompare(a1->acc,     a2->acc)     != eslOK) return eslFAIL;
  if (esl_CCompare(a1->au,      a2->au)      != eslOK) return eslFAIL;
  if (esl_CCompare(a1->ss_cons, a2->ss_cons) != eslOK) return eslFAIL;
  if (esl_CCompare(a1->sa_cons, a2->sa_cons) != eslOK) return eslFAIL;
  if (esl_CCompare(a1->pp_cons, a2->pp_cons) != eslOK) return eslFAIL;
  if (esl_CCompare(a1->rf,      a2->rf)      != eslOK) return eslFAIL;
  
  if (a1->sqacc != NULL && a2->sqacc != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sqacc[i], a2->sqacc[i]) != eslOK) return eslFAIL;
  } else if (a1->sqacc != NULL || a2->sqacc != NULL) return eslFAIL;

  if (a1->sqdesc != NULL && a2->sqdesc != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sqdesc[i], a2->sqdesc[i]) != eslOK) return eslFAIL;
  } else if (a1->sqdesc != NULL || a2->sqdesc != NULL) return eslFAIL;

  if (a1->ss != NULL && a2->ss != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->ss[i], a2->ss[i]) != eslOK) return eslFAIL;
  } else if (a1->ss != NULL || a2->ss != NULL) return eslFAIL;

  if (a1->sa != NULL && a2->sa != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->sa[i], a2->sa[i]) != eslOK) return eslFAIL;
  } else if (a1->sa != NULL || a2->sa != NULL) return eslFAIL;

  if (a1->pp != NULL && a2->pp != NULL) {
    for (i = 0; i < a1->nseq; i++) if (esl_CCompare(a1->pp[i], a2->pp[i]) != eslOK) return eslFAIL;
  } else if (a1->pp != NULL || a2->pp != NULL) return eslFAIL;
  
  for (i = 0; i < eslMSA_NCUTS; i++)
    {
      if (a1->cutset[i] && a2->cutset[i]) {
	if (esl_FCompare(a1->cutoff[i], a2->cutoff[i], 0.01) != eslOK) return eslFAIL;
      } else if (a1->cutset[i] || a2->cutset[i]) return eslFAIL;
    }
  return eslOK;
}
/*---------------- end of debugging/development routines  -------------------*/



/******************************************************************************
 * 13. Benchmark driver.
 *****************************************************************************/
#ifdef eslMSA_BENCHMARK
/* gcc -O2 -o msa_benchmark -I. -L. -DeslMSA_BENCHMARK esl_msa.c -leasel -lm
 * ./benchmark Pfam
 */
#include <stdio.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS options[] = {
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-v",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "be verbose: show info as each msa is read",        0 },
  { 0,0,0,0,0,0,0,0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "benchmarking speed of MSA reading";

int
main(int argc, char **argv)
{
  ESL_GETOPTS   *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  ESL_STOPWATCH *w         = esl_stopwatch_Create();
  char          *filename  = esl_opt_GetArg(go, 1);
  int            fmt       = eslMSAFILE_UNKNOWN;
  ESL_MSAFILE   *afp       = NULL;
  ESL_MSA       *msa       = NULL;
  int            status;
  int            nali;
  int            alphatype;

  /* Open the alignment file in text mode */
  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND)     esl_fatal("Alignment file %s doesn't exist or is not readable\n", filename);
  else if (status == eslEFORMAT)  esl_fatal("Couldn't determine format of alignment %s\n", filename);
  else if (status != eslOK)       esl_fatal("Alignment file open failed with error %d\n", status);

  /* Loop over all the alignments. */
  esl_stopwatch_Start(w);
  nali = 0;
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;

      if (esl_opt_GetBoolean(go, "-v")) 
	{
	  esl_msa_GuessAlphabet(msa, &alphatype);
	  printf("%5d %15s %6d %5d %7s\n", 
		 nali, msa->name, msa->nseq, msa->alen, esl_abc_DecodeType(alphatype));
	}

      esl_msa_Destroy(msa);
    }
  if (status == eslEFORMAT)
    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nOffending line is:\n%s\n", 
	      afp->linenumber, afp->fname, afp->errbuf, afp->buf);
  else if (status != eslEOF)
    esl_fatal("Alignment file read failed with error code %d\n", status);

  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "CPU Time: ");

  esl_stopwatch_Destroy(w);
  esl_getopts_Destroy(go);
  esl_msafile_Close(afp);
  exit(0);
}
#endif /*eslMSA_BENCHMARK*/
/*---------------------- end of benchmark driver ----------------------*/



/******************************************************************************
 * 15. Unit tests
 *****************************************************************************/
#ifdef eslMSA_TESTDRIVE

/* write_known_msa()
 * Write a known MSA to a tmpfile in Stockholm format.
 */
static void
write_known_msa(FILE *ofp)
{
  fprintf(ofp, "# STOCKHOLM 1.0\n");
  fprintf(ofp, "seq1 --ACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "seq2 aaACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "seq3 aaACDEFGHIK~LMNPQRS-TVWY\n");
  fprintf(ofp, "\n");
  fprintf(ofp, "seq1 ACDEFGHIKLMNPQRSTVWY~~~\n");
  fprintf(ofp, "seq2 ACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "seq3 ACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "//\n");
  return;
}
  
/* write_known_pfam_msa()
 * Write a known MSA to a tmpfile in Pfam Stockholm format.
 */
static void
write_known_pfam_msa(FILE *ofp)
{
  fprintf(ofp, "# STOCKHOLM 1.0\n");
  fprintf(ofp, "#=GF ID pfam-test\n");
  fprintf(ofp, "#=GS seq2 DE seq2 is interesting\n");
  fprintf(ofp, "seq1    --ACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWY~~~\n");
  fprintf(ofp, "seq2    aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "seq3    aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy\n");
  fprintf(ofp, "#=GC RF ..xxxxxxxxx.xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...\n");
  fprintf(ofp, "//\n");
  return;
}

/* compare_to_known() 
 * SRE, Thu Sep  7 09:52:07 2006 [Janelia]
 * Spotcheck an ESL_MSA to make sure it matches the test known alignment.
 */
static void
compare_to_known(ESL_MSA *msa)
{
  if (msa->alen != 47)                     esl_fatal("bad alen");
  if (msa->nseq != 3)                      esl_fatal("bad nseq");
  if (strcmp(msa->sqname[1], "seq2") != 0) esl_fatal("bad sqname");
#ifdef eslAUGMENT_ALPHABET
  if (msa->flags & eslMSA_DIGITAL)
    {
      if (! esl_abc_XIsGap(msa->abc, msa->ax[0][2]))      esl_fatal("no gap where expected");
      if (! esl_abc_XIsMissing(msa->abc, msa->ax[0][47])) esl_fatal("no missing-data symbol where expected");
      if (msa->ax[1][1]  != 0)                            esl_fatal("spotcheck on ax failed"); /* 0=A */
      if (msa->ax[1][47] != 19)                           esl_fatal("spotcheck on ax failed"); /*19=Y */
    }
#endif
  if (! (msa->flags & eslMSA_DIGITAL))
    {
      if (strcasecmp(msa->aseq[0], "--ACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWY~~~") != 0) esl_fatal("aseq 0 is bad");
      if (strcasecmp(msa->aseq[1], "aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy") != 0) esl_fatal("aseq 1 is bad");
      if (strcasecmp(msa->aseq[2], "aaACDEFGHIK~LMNPQRS-TVWYACDEFGHIKLMNPQRSTVWYyyy") != 0) esl_fatal("aseq 2 is bad");
    }
  return;
}

/* msa_compare()
 * SRE, Fri Sep  8 08:13:38 2006 [Janelia]
 * 
 * Compares two MSAs; returns eslOK if they appear to be the same, eslFAIL if not.
 * (Not worth putting in external API. Only useful for testing purposes.)
 * Not a complete comparison: just checks mode, sequence names, and aligned data.
 * MSAs have to be in same mode (text vs. digital).
 */
static void
msa_compare(ESL_MSA *m1, ESL_MSA *m2)
{
  int i;

  if (m1->nseq  != m2->nseq   ||
      m1->alen  != m2->alen   ||
      m1->flags != m2->flags)
    esl_fatal("msa1, msa2 differ in nseq, alen, or flags");

  for (i = 0; i < m1->nseq; i++)
    {
      if (strcmp(m1->sqname[i], m2->sqname[i]) != 0) esl_fatal("msa1, msa2 sqnames differ for %d", i);
#ifdef eslAUGMENT_ALPHABET
      if ((m1->flags & eslMSA_DIGITAL) && 
	  memcmp(m1->ax[i], m2->ax[i], sizeof(ESL_DSQ) * (m1->alen+2)) != 0) 
	esl_fatal("msa1, msa2 digital sequences differ for %d", i);
#endif
      if (! (m1->flags & eslMSA_DIGITAL) && 
	  strcmp(m1->aseq[i], m2->aseq[i]) != 0) 
	esl_fatal("msa1, msa2 sequences differ for %d", i);
    }
  return;
}

/* Unit tests for every function in the exposed API
 */
static void
utest_Create(void)
{
  ESL_MSA *msa = NULL;

  msa = esl_msa_Create(16, -1);	  /* nseq blocksize 16, growable */
  esl_msa_Destroy(msa);
  msa = esl_msa_Create(16, 100);  /* nseq=16, alen=100, not growable */
  esl_msa_Destroy(msa);

  return;
}

static void
utest_Destroy(void)
{
  ESL_MSA *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;
#endif

  msa = esl_msa_Create(16, -1);	
  esl_msa_Destroy(msa);	 	  /* normal usage */

#ifdef eslAUGMENT_ALPHABET
  abc = esl_alphabet_Create(eslRNA);
  msa = esl_msa_CreateDigital(abc, 16, 100);	
  esl_msa_Destroy(msa);	 	  /* normal usage, digital mode */
  esl_alphabet_Destroy(abc);
#endif

  esl_msa_Destroy(NULL);	  /* should tolerate NULL argument */
  return;
}

static void
utest_Expand(void)
{
  ESL_MSA *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;
#endif

  msa = esl_msa_Create(16, -1);                	    /* growable */
  if (esl_msa_Expand(msa) != eslOK) esl_fatal("Expand failed"); /* expand by 2x in nseq */
  esl_msa_Destroy(msa);

  msa = esl_msa_Create(16, 100);                        /* not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal("Expand should have failed but didn't"); /* should fail w/ EINVAL*/
#endif
  esl_msa_Destroy(msa);
  
#ifdef eslAUGMENT_ALPHABET
  abc = esl_alphabet_Create(eslDNA);
  msa = esl_msa_CreateDigital(abc, 16, -1);               /* growable */
  if (esl_msa_Expand(msa) != eslOK) esl_fatal("Expand failed"); /* expand by 2x in nseq */
  esl_msa_Destroy(msa);

  msa = esl_msa_CreateDigital(abc, 16, 100);                 /* not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal("Expand should have failed but didn't"); /* should fail w/ EINVAL*/
#endif /* eslTEST_THROWING*/
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}

static void
utest_Open(char *tmpfile)	/* filename must be in /tmp */
{
  char        *msg      = "Open() unit test failed";
  ESL_MSAFILE *msafp    = NULL;
  int          status;
  
  status = esl_msafile_Open(tmpfile, eslMSAFILE_UNKNOWN, NULL, &msafp); 
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);

  status = esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &msafp);
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);

#ifdef HAVE_PUTENV
  putenv("ESLTEST=./");
  esl_FileTail(tmpfile, FALSE, &filename);
  status = esl_msafile_Open(filename, eslMSAFILE_STOCKHOLM, "ESLTEST", &msafp);
  if (status != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);
  free(filename);
#endif

  return;
}

static void
utest_Close(char *filename)
{
  ESL_MSAFILE *msafp    = NULL;
  int status;

  status = esl_msafile_Open(filename, eslMSAFILE_UNKNOWN, NULL, &msafp); 
  if (status != eslOK) esl_fatal("Close() unit test failed");
  esl_msafile_Close(msafp);
  esl_msafile_Close(NULL);	/* should work */
  return;
}

#ifdef eslAUGMENT_ALPHABET
static void
utest_CreateDigital(ESL_ALPHABET *abc)
{
  char    *msg = "CreateDigital() unit test failure";
  ESL_MSA *msa = NULL;

  msa = esl_msa_CreateDigital(abc, 16, -1);	  /* nseq blocksize 16, growable */
  if (! (msa->flags & eslMSA_DIGITAL)) esl_fatal(msg);
  if (msa->ax   == NULL)               esl_fatal(msg);
  if (msa->aseq != NULL)               esl_fatal(msg);
  if (esl_msa_Expand(msa) != eslOK)    esl_fatal(msg);
  esl_msa_Destroy(msa);

  msa = esl_msa_CreateDigital(abc, 16, 100);  /* nseq=16, alen=100, not growable */
#ifdef eslTEST_THROWING
  if (esl_msa_Expand(msa) != eslEINVAL) esl_fatal(msg); /* shouldn't grow */
#endif
  esl_msa_Destroy(msa);

  return;
}
#endif /*eslAUGMENT_ALPHABET*/

#ifdef eslAUGMENT_ALPHABET
static void
utest_Digitize(ESL_ALPHABET *abc, char *filename)
{
  char        *msg = "Digitize() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
  int c, i, pos;

  /* Get ourselves a copy of the known alignment that we can muck with */
  if (esl_msafile_Open(filename, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK)  esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                       esl_fatal(msg);
  esl_msafile_Close(mfp);
  
  /* Deliberately corrupt it with inval character in the middle */
  i   = msa->nseq / 2;
  pos = msa->alen / 2;
  c   = msa->aseq[i][pos];
  msa->aseq[i][pos] = '%';
  if (esl_msa_Digitize(abc, msa, NULL) != eslEINVAL) esl_fatal(msg); /* should detect corruption as normal error */
  msa->aseq[i][pos] = c;	                               /* restore original         */
  compare_to_known(msa);
  if (esl_msa_Digitize(abc, msa, NULL) != eslOK)     esl_fatal(msg); /* should be fine now       */
  compare_to_known(msa);

  esl_msa_Destroy(msa);
  return;
}
#endif /*eslAUGMENT_ALPHABET*/


#ifdef eslAUGMENT_ALPHABET
static void
utest_Textize(ESL_ALPHABET *abc, char *filename)
{
  char        *msg = "Textize() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;

  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_UNKNOWN, NULL, &mfp) != eslOK)  esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)   esl_fatal(msg);
  if (esl_msa_Textize(msa)    != eslOK)   esl_fatal(msg);
  compare_to_known(msa);

  esl_msafile_Close(mfp);
  esl_msa_Destroy(msa);
  return;
}
#endif /*eslAUGMENT_ALPHABET*/

#ifdef eslAUGMENT_ALPHABET
static void
utest_OpenDigital(ESL_ALPHABET *abc, char *filename)  /* filename must be in /tmp */
{
  char        *msg   = "OpenDigital() unit test failure";
  ESL_MSAFILE *msafp = NULL;
  
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_UNKNOWN,   NULL, &msafp) != eslOK) esl_fatal(msg);  esl_msafile_Close(msafp);
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_STOCKHOLM, NULL, &msafp) != eslOK) esl_fatal(msg);  esl_msafile_Close(msafp);
#ifdef HAVE_PUTENV
  putenv("ESLTEST=./");
  esl_FileTail(tmpfile, FALSE, &filename);
  if (esl_msafile_OpenDigital(abc, filename, eslMSAFILE_STOCKHOLM, "ESLTEST", &msafp) != eslOK) esl_fatal(msg);
  esl_msafile_Close(msafp);
  free(filename);
#endif
  return;
}
#endif /*eslAUGMENT_ALPHABET*/

static void
utest_Read(char *filename)
{
  char        *msg = "Read() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;

  if (esl_msafile_Open(filename, eslMSAFILE_UNKNOWN, NULL, &mfp) != eslOK)  esl_fatal(msg);  
  if (esl_msa_Read(mfp, &msa) != eslOK)  esl_fatal(msg);
  compare_to_known(msa);
  esl_msa_Destroy(msa);

  if (esl_msa_Read(mfp, &msa) != eslEOF) esl_fatal(msg);
  if (msa != NULL)                       esl_fatal(msg);

  esl_msafile_Close(mfp);
  return;
}

static void
utest_Write(ESL_MSA *msa1)
{
  char        *msg  = "Write() unit test failure";
  ESL_MSAFILE *mfp  = NULL;
  ESL_MSA     *msa2 = NULL;
  FILE        *fp   = NULL;
  int      i;
  int      formats[] = { eslMSAFILE_STOCKHOLM, eslMSAFILE_PFAM, eslMSAFILE_AFA, -1 }; /* -1=sentinel */
  char     template[16] = "esltmpXXXXXX";
  char     tmpfile[16];

  for (i = 0; formats[i] != -1; i++)
    {
      strcpy(tmpfile, template);
      if (esl_tmpfile_named(tmpfile, &fp) != eslOK)     esl_fatal(msg);
      if (esl_msa_Write(fp, msa1, formats[i]) != eslOK) esl_fatal(msg);
      fclose(fp);
  
#ifdef eslAUGMENT_ALPHABET
      if ((msa1->flags & eslMSA_DIGITAL) &&
	  esl_msafile_OpenDigital(msa1->abc, tmpfile, formats[i], NULL, &mfp) != eslOK)	esl_fatal(msg);
#endif
      if (! (msa1->flags & eslMSA_DIGITAL) &&
	  esl_msafile_Open(tmpfile, formats[i], NULL, &mfp) != eslOK) esl_fatal(msg);

      if (esl_msa_Read(mfp, &msa2) != eslOK) esl_fatal(msg);
      msa_compare(msa1, msa2);
      
      esl_msafile_Close(mfp);
      esl_msa_Destroy(msa2);
      remove(tmpfile);      
    }      
  return;
}

static void
utest_GuessFileFormat(void)
{
  /* SRE: To be filled in. Currently, esl_msa_GuessFileFormat() is a placeholder that
   * always guesses Stockholm
   */
  return;
}


static void
utest_SequenceSubset(ESL_MSA *m1)
{
  char    *msg   = "SequenceSubset() unit test failure";
  ESL_MSA *m2    = NULL;
  int     *useme = NULL;
  int      i,j;
  int      n2;

  /* Make every other sequence (1,3..) get excluded from the subset */
  useme = malloc(m1->nseq * sizeof(int));
  for (i = 0, n2 = 0; i < m1->nseq; i++)
    if (i%2 == 0) { useme[i] = TRUE; n2++; }
    else          useme[i] = FALSE;

  if (esl_msa_SequenceSubset(m1, useme, &m2) != eslOK) esl_fatal(msg);
  if (m2->nseq != n2) esl_fatal(msg);
  
  for (i = 0, j = 0; i < m1->nseq; i++)
    {
      if (useme[i])
	{
	  if (strcmp(m1->sqname[i], m2->sqname[j]) != 0) esl_fatal(msg);
	  if (! (m1->flags & eslMSA_DIGITAL) && (strcmp(m1->aseq[i],   m2->aseq[j])  != 0)) esl_fatal(msg);
#ifdef eslAUGMENT_ALPHABET
	  if (  (m1->flags & eslMSA_DIGITAL) && memcmp(m1->ax[i], m2->ax[j], sizeof(ESL_DSQ) * (m1->alen+2)) != 0) esl_fatal(msg);
#endif
	  j++;
	}
    }  
  esl_msa_Destroy(m2);
  free(useme);
  return;
}

static void
utest_MinimGaps(char *tmpfile)
{
  char        *msg = "MinimGaps() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_MinimGaps(msa, NULL, "-~") != eslOK) esl_fatal(msg);
  if (msa->alen        != 45)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (msa->aseq[0][11] != 'L') esl_fatal(msg); /* L shifted from column 13->12 */
  if (msa->aseq[0][18] != 'T') esl_fatal(msg); /* T shifted from column 21->19 */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_MinimGaps(msa, NULL, NULL) != eslOK) esl_fatal(msg);
  if (msa->alen        != 45)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (esl_msa_Textize(msa) != eslOK) esl_fatal(msg);
  if (msa->aseq[0][11] != 'L') esl_fatal(msg); /* L shifted from column 13->12 */
  if (msa->aseq[0][18] != 'T') esl_fatal(msg); /* T shifted from column 21->19 */
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}  

static void
utest_NoGaps(char *tmpfile)
{
  char        *msg = "NoGaps() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_NoGaps(msa, NULL, "-~") != eslOK) esl_fatal(msg);
  if (msa->alen        != 40)  esl_fatal(msg); /* orig =47, w/ 7 columns with gaps */
  if (msa->aseq[0][9]  != 'L') esl_fatal(msg); /* L shifted from column 13->10  */
  if (msa->aseq[0][16] != 'T') esl_fatal(msg); /* T shifted from column 21->17 */
  if (msa->aseq[0][39] != 'Y') esl_fatal(msg); /* Y shifted from column 47->40 */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
  if (esl_msa_NoGaps(msa, NULL, NULL) != eslOK) esl_fatal(msg);
  if (msa->alen        != 40)  esl_fatal(msg); /* orig =47, with one all - column and one all ~ column */
  if (esl_msa_Textize(msa) != eslOK) esl_fatal(msg);
  if (msa->aseq[0][9]  != 'L') esl_fatal(msg); /* L shifted from column 13->10  */
  if (msa->aseq[0][16] != 'T') esl_fatal(msg); /* T shifted from column 21->17 */
  if (msa->aseq[0][39] != 'Y') esl_fatal(msg); /* Y shifted from column 47->40 */
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}  

static void
utest_SymConvert(char *tmpfile)
{
  char        *msg = "SymConvert() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc = NULL;
#endif

  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK)                                     esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* many->one version */
  if (esl_msa_SymConvert(msa, "VWY", "-")   != eslOK) esl_fatal(msg); /* 6 columns convert to all-gap: now 8/47 */
  if (esl_msa_MinimGaps(msa, NULL, "-~")    != eslOK) esl_fatal(msg); /* now we're 39 columns long */
  if (msa->alen                             != 39)    esl_fatal(msg);

  /* many->many version */
  if (esl_msa_SymConvert(msa, "DEF", "VWY") != eslOK) esl_fatal(msg);
  if (msa->aseq[0][4]                       != 'V')   esl_fatal(msg);
  if (msa->aseq[0][5]                       != 'W')   esl_fatal(msg);
  if (msa->aseq[0][23]                      != 'Y')   esl_fatal(msg); /* F in orig col 29; -5; converted to Y */

  /* bad calls */
#ifdef eslTEST_THROWING
  if (esl_msa_SymConvert(msa, "XXX", "XX")  != eslEINVAL) esl_fatal(msg); /* check for clean fail on mismatched args */
#endif
  esl_msa_Destroy(msa);
  
#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);
#ifdef eslTEST_THROWING
  if (esl_msa_SymConvert(msa, "Tt", "Uu") != eslEINVAL) esl_fatal(msg); /* must cleanly fail on digital mode msa */
#endif
  esl_msa_Destroy(msa);
  esl_alphabet_Destroy(abc);
#endif
  return;
}

/* Exercise a boundary case: zero length MSA (alen=0) */
/* Given an input *digital* MSA as a starting point, we clone it, 
 * column subset it to zero length, then make sure that 
 * various MSA functions operate correctly on it;
 * then we textize it and test it in text mode; then we 
 * digitize it again, and throw it away.
 * (The input <msa> is unchanged.)
 */
static void
utest_ZeroLengthMSA(const char *tmpfile)
{
  char    *msg      = "zero length msa unit test failed";
  ESL_MSAFILE *mfp  = NULL;
  ESL_MSA *z1       = NULL;
  ESL_MSA *z2       = NULL;
  ESL_MSA *z3       = NULL;
  int     *useme    = NULL;
  int      nuseme   = 0;
  int      i;
  char     errbuf[eslERRBUFSIZE];


  /* Read a text mode alignment from the tmpfile */
  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &z1) != eslOK)                                      esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* make an alen=0 text alignment by column subsetting */
  nuseme = ESL_MAX(z1->alen, z1->nseq);
  if ((useme = malloc(sizeof(int) * nuseme)) == NULL)  esl_fatal(msg);
  for (i = 0; i < z1->alen; i++) useme[i] = 0;
  if (esl_msa_ColumnSubset(z1, errbuf, useme) != eslOK) esl_fatal(msg);

  /* These should all no-op if alen=0*/
  if (esl_msa_MinimGaps(z1, NULL, "-")!= eslOK) esl_fatal(msg);
  if (esl_msa_NoGaps(z1, NULL, "-")   != eslOK) esl_fatal(msg);
  if (esl_msa_SymConvert(z1,"RY","NN")!= eslOK) esl_fatal(msg);
  
  /* test sequence subsetting by removing the first sequence */
  for (i = 1; i < z1->nseq; i++) useme[i] = 1;  
  if (esl_msa_SequenceSubset(z1, useme, &z2) != eslOK) esl_fatal(msg);
  esl_msa_Destroy(z1);
  /* keep z2; we'll compare it to z3 in the end */
      
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET *abc;

  /* Now read the same alignment, in digital mode */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal(msg);
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &z1) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* Now make an alen=0 alignment in digital mode */
  for (i = 0; i < z1->alen; i++) useme[i] = 0;
  if (esl_msa_ColumnSubset(z1, errbuf, useme) != eslOK) esl_fatal(msg);

  /* again these should all no-op if alen=0*/
  if (esl_msa_MinimGaps(z1, NULL, NULL) != eslOK) esl_fatal(msg);
  if (esl_msa_NoGaps(z1, NULL, NULL)    != eslOK) esl_fatal(msg);
  /* SymConvert throws EINVAL on a digital mode alignment */

  /* test sequence subsetting by removing the first sequence */
  for (i = 1; i < z1->nseq; i++) useme[i] = 1;  
  if (esl_msa_SequenceSubset(z1, useme, &z3) != eslOK) esl_fatal(msg);
  esl_msa_Destroy(z1);

  if ((z1 = esl_msa_Clone(z3))        == NULL)  esl_fatal(msg); /* z1 is now alen=0, digital */
  if (esl_msa_Textize(z3)             != eslOK) esl_fatal(msg); /* convert z3 back to text mode */
  if (esl_msa_Compare(z2, z3)         != eslOK) esl_fatal(msg); /* compare in text mode */
  if (esl_msa_Digitize(abc, z2, NULL) != eslOK) esl_fatal(msg); /* now z2 is digital */
  if (esl_msa_Compare(z1, z2)         != eslOK) esl_fatal(msg); /* compare digital mode z1,z2 */

  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(z1);
  esl_msa_Destroy(z3);
#endif /*eslAUGMENT_ALPHABET*/

  esl_msa_Destroy(z2);
  free(useme);
}

static void
utest_ReadNonSeqInfoPfam(char *filename)
{
  char        *msg = "ReadNonSeqInfo() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  ESL_MSA     *msa = NULL;
  int          nseq = 0;
  int64_t      alen = 0;
  int          ngs = 0;
  int          maxname = 0;
  int          maxgf = 0;
  int          maxgc = 0;
  int          maxgr = 0;

  if (esl_msafile_Open(filename, eslMSAFILE_PFAM, NULL, &mfp) != eslOK) esl_fatal(msg);  /* don't autodetect, assert pfam, ReadNonSeqInfo() requires it */
  if (esl_msa_ReadNonSeqInfoPfam(mfp, NULL, -1, NULL, NULL, &msa, &nseq, &alen, &ngs, &maxname, &maxgf, &maxgc, &maxgr, NULL, NULL, NULL, NULL, NULL) != eslOK)  esl_fatal(msg);

  if (msa->nseq != 0)  esl_fatal("bad msa->nseq");
  if (msa->alen != -1) esl_fatal("bad msa->alen");
  if (nseq      != 3)  esl_fatal("bad nseq");
  if (alen      != 47) esl_fatal("bad alen");
  if (ngs       != 1)  esl_fatal("bad ngs");
  if (maxname   != 4)  esl_fatal("bad maxname");
  if (maxgf     != 2)  esl_fatal("bad maxgf");
  if (maxgc     != 2)  esl_fatal("bad maxgc");
  if (maxgr     != 0)  esl_fatal("bad maxgr");
  esl_msa_Destroy(msa);

  if (esl_msa_ReadNonSeqInfoPfam(mfp, NULL, -1, NULL, NULL, &msa, &nseq, &alen, &ngs, &maxname, &maxgf, &maxgc, &maxgr, NULL, NULL, NULL, NULL, NULL) != eslEOF) esl_fatal(msg);
  if (msa  != NULL) esl_fatal(msg);
  if (nseq != 0 || alen != 0 || ngs != 0 || maxname != 0 || maxgf != 0 || maxgc != 0 || maxgr != 0) esl_fatal("bad nseq");

  esl_msafile_Close(mfp);
  return;
}

static void
utest_RegurgitatePfam(char *filename)
{
  char        *msg = "RegurgitatePfam() unit test failure";
  ESL_MSAFILE *mfp = NULL;
  char         tmpfile[16] = "esltmpXXXXXX";
  FILE        *fp = NULL;
  ESL_MSA     *msa1 = NULL;
  ESL_MSA     *msa2 = NULL;

  /* regurgitate msa in filename to tmpfile (an msa structure will not be created) */
  if (esl_msafile_Open(filename, eslMSAFILE_PFAM, NULL, &mfp) != eslOK) esl_fatal(msg);  /* don't autodetect, assert pfam, ReadNonSeqInfo() requires it */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal(msg);
  if (esl_msa_RegurgitatePfam(mfp, fp, 
			      -1, -1, -1, -1, /* maxname, maxgf, maxgc, maxgr unknown: output msa formatting will match input msa formatting */
			      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, /* do_header, do_trailer, do_blanks, do_comments, do_gf, do_gs, do_gc, do_gr, do_aseq: print all components */
			      NULL, /* seqs2regurg: if non-NULL specifies which sequences to keep in output */
			      NULL, /* seqs2skip:   if non-NULL specifies which sequences to skip in output */
			      NULL, /* useme:  if non-NULL specifies which columns to keep in output */
			      NULL, /* add2me: if non-NULL specifies how many gap columns to add in output */
			      -1,   /* expected alignment length, unknown (must not be if useme != NULL or add2me != NULL */
			      '.')  /* gapchar, irrelevant since add2me is NULL */
      != eslOK) esl_fatal(msg);
  fclose(fp);
  esl_msafile_Close(mfp);

  /* read in msa from filename as msa1 */
  if (esl_msafile_Open(filename, eslMSAFILE_PFAM, NULL, &mfp) != eslOK) esl_fatal(msg); 
  if (esl_msa_Read(mfp, &msa1) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);

  /* read in msa from tmpfile as msa2 */
  if (esl_msafile_Open(tmpfile, eslMSAFILE_PFAM, NULL, &mfp) != eslOK) esl_fatal(msg);
  if (esl_msa_Read(mfp, &msa2) != eslOK) esl_fatal(msg);
  esl_msafile_Close(mfp);

  msa_compare(msa1, msa2);

  esl_msa_Destroy(msa1);
  esl_msa_Destroy(msa2);
  remove(tmpfile);

  return;
}
#endif /* eslMSA_TESTDRIVE */
/*------------------------ end of unit tests --------------------------------*/



/*****************************************************************************
 * 16. Test driver
 *****************************************************************************/
#ifdef eslMSA_TESTDRIVE
/* 
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_KEYHASH esl_msa.c esl_keyhash.c easel.c -lm
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_ALPHABET esl_msa.c esl_alphabet.c easel.c -lm
 * gcc -g -Wall -o msa_utest -I. -DeslMSA_TESTDRIVE -DAUGMENT_SSI esl_msa.c esl_ssi.c easel.c -lm
 * gcc -g -Wall -o msa_utest -L. -I. -DeslMSA_TESTDRIVE esl_msa.c -leasel -lm
 * gcc -g -Wall -o msa_utest -L. -I. -DeslTEST_THROWING -DeslMSA_TESTDRIVE esl_msa.c -leasel -lm
 * ./msa_utest
 */
#include <stdlib.h>
#include <stdio.h>

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_KEYHASH
#include "esl_keyhash.h"
#endif
#ifdef eslAUGMENT_RANDOM
#include "esl_random.h"
#endif
#ifdef eslAUGMENT_SSI
#include "esl_ssi.h"
#endif
#include "esl_msa.h"


int
main(int argc, char **argv)
{
  ESL_MSAFILE    *mfp  = NULL;
  ESL_MSA        *msa  = NULL;
  FILE           *fp   = NULL;
  char            tmpfile[16]  = "esltmpXXXXXX"; /* tmpfile template */
  char            tmpfile2[16] = "esltmpXXXXXX"; /* tmpfile template */
#ifdef eslAUGMENT_ALPHABET
  ESL_ALPHABET   *abc  = NULL;
#endif

#ifdef eslTEST_THROWING
  esl_exception_SetHandler(&esl_nonfatal_handler);
#endif

  /* Create a known Stockholm test alignment in a tempfile.
   */
  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal("failed to create tmpfile");
  write_known_msa(fp);
  fclose(fp);

  /* Read it back in for use in tests.
   */
  if (esl_msafile_Open(tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) esl_fatal("Failed to open MSA tmp file");
  if (esl_msa_Read(mfp, &msa)                                     != eslOK) esl_fatal("Failed to read MSA tmp file");
  esl_msafile_Close(mfp);

  /* Unit tests
   */
  utest_Create();
  utest_Destroy();
  utest_Expand();
  utest_Open(tmpfile);
  utest_Close(tmpfile);
  utest_Read(tmpfile);
  utest_Write(msa);
  utest_GuessFileFormat();
  utest_SequenceSubset(msa);
  utest_MinimGaps(tmpfile);
  utest_NoGaps(tmpfile);
  utest_SymConvert(tmpfile);
  utest_ZeroLengthMSA(tmpfile);	/* this tests in digital mode too if eslAUGMENT_ALPHABET */
  esl_msa_Destroy(msa);

#ifdef eslAUGMENT_ALPHABET
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) 
    esl_fatal("alphabet creation failed");
  if (esl_msafile_OpenDigital(abc, tmpfile, eslMSAFILE_STOCKHOLM, NULL, &mfp) != eslOK) 
    esl_fatal("MSA digital open failed");
  if (esl_msa_Read(mfp, &msa) != eslOK) 
    esl_fatal("MSA digital read failed");
  esl_msafile_Close(mfp);

  utest_CreateDigital(abc);
  utest_Digitize(abc, tmpfile);
  utest_Textize(abc, tmpfile);
  utest_OpenDigital(abc, tmpfile);
  utest_Write(msa);

  esl_alphabet_Destroy(abc);
  esl_msa_Destroy(msa);
#endif

  /* Unit tests of memory efficient functions
   */
  if (esl_tmpfile_named(tmpfile2, &fp) != eslOK) esl_fatal("failed to create tmpfile2");
  write_known_pfam_msa(fp);
  fclose(fp);

  utest_ReadNonSeqInfoPfam(tmpfile2);
  utest_RegurgitatePfam(tmpfile2);

  remove(tmpfile);
  remove(tmpfile2);
  exit(0);	/* success  */
}
#endif /*eslMSA_TESTDRIVE*/
/*-------------------- end of test driver ---------------------*/




/******************************************************************************
 * 17. Examples.
 *****************************************************************************/
/* The examples are also useful as i/o speed benchmarks. */


/* Example 1: 
 *   time ./example SSU_rRNA_5 > /dev/null
 *     without keyhash: [345.118u 31.564s 10:45.60 58.3%  SRE, Tue Sep  5 11:52:41 2006]
 *     with keyhash:    [33.353u 1.681s 0:35.04 99.9% SRE, Tue Sep  5 11:55:00 2006]
 */
#ifdef eslMSA_EXAMPLE
/*::cexcerpt::msa_example::begin::*/
/* An example of reading an MSA in text mode, and handling any returned errors.
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE esl_msa.c easel.c 
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE -DeslAUGMENT_KEYHASH esl_msa.c esl_keyhash.c easel.c
   gcc -g -Wall -o example -I. -L. -DeslMSA_EXAMPLE esl_msa.c -leasel -lm
   ./example <MSA file>
 */
#include <stdio.h>
#include "easel.h"
#include "esl_msa.h"

int
main(int argc, char **argv)
{
  char        *filename = argv[1];
  int          fmt      = eslMSAFILE_SELEX;
  ESL_MSAFILE *afp      = NULL;
  ESL_MSA     *msa      = NULL;
  int          nali     = 0;
  int          status;


  status = esl_msafile_Open(filename, fmt, NULL, &afp);
  if (status == eslENOTFOUND)    esl_fatal("Alignment file %s isn't readable\n", filename);
  else if (status == eslEFORMAT) esl_fatal("Couldn't determine format of %s\n",  filename);
  else if (status != eslOK)      esl_fatal("Alignment file open failed (error %d)\n", status);

  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, msa->nseq, (int) msa->alen);
      esl_msa_Write(stdout, msa, eslMSAFILE_PFAM);
      esl_msa_Destroy(msa);
    }
  if      (status == eslEFORMAT) esl_fatal("alignment file %s: %s\n", afp->fname, afp->errbuf);
  else if (status != eslEOF)     esl_fatal("alignment file %s: read failed, error %d\n", afp->fname, status);
  else if (nali == 0)            esl_fatal("alignment file %s: %s\n", afp->fname, afp->errbuf);

  esl_msafile_Close(afp);
  exit(0);
}
/*::cexcerpt::msa_example::end::*/
#endif /*eslMSA_EXAMPLE*/




#ifdef eslAUGMENT_ALPHABET
#ifdef eslMSA_EXAMPLE2
/*::cexcerpt::msa_example2::begin::*/
/* An example of reading an MSA in digital mode, after guessing
 * the alphabet by looking at the first alignment.
 *
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE2 -DeslAUGMENT_ALPHABET esl_msa.c esl_alphabet.c easel.c 
   gcc -g -Wall -o example -I. -DeslMSA_EXAMPLE2 -DeslAUGMENT_ALPHABET -DeslAUGMENT_KEYHASH\
       esl_msa.c esl_keyhash.c esl_alphabet.c easel.c
   gcc -g -Wall -o example -I. -L. -DeslMSA_EXAMPLE2 esl_msa.c -leasel -lm
   ./example <MSA file>
 */
#include <stdio.h>
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",    0 },
  { "--dna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use DNA alphabet",                        0 },
  { "--rna",     eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use RNA alphabet",                        0 },
  { "--amino",   eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "use protein alphabet",                    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <msafile>";
static char banner[] = "example of digital MSA reading using the msa module";

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go        = esl_getopts_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char         *msafile   = esl_opt_GetArg(go, 1);
  int           fmt       = eslMSAFILE_UNKNOWN;
  int           alphatype = eslUNKNOWN;
  ESL_ALPHABET *abc       = NULL;
  ESL_MSAFILE  *afp       = NULL;
  ESL_MSA      *msa       = NULL;
  int           nali      = 0;
  int           status;

  /* First you open the msa file in normal text mode */
  status = esl_msafile_Open(msafile, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s isn't readable", msafile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of %s",  msafile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed (error code %d)", status);

  /* Now you can set or guess the alphabet type - this looks at the first alignment */
  if      (esl_opt_GetBoolean(go, "--rna"))   alphatype = eslRNA;
  else if (esl_opt_GetBoolean(go, "--dna"))   alphatype = eslDNA;
  else if (esl_opt_GetBoolean(go, "--amino")) alphatype = eslAMINO;
  else {
    status = esl_msafile_GuessAlphabet(afp, &alphatype);
    if      (status == eslEAMBIGUOUS) esl_fatal("Couldn't guess alphabet from first alignment in %s", msafile);
    else if (status == eslEFORMAT)    esl_fatal("Alignment file parse error, line %d of file %s:\n%s\nBad line is: %s\n",
						afp->linenumber, afp->fname, afp->errbuf, afp->buf);
    else if (status == eslENODATA)    esl_fatal("Alignment file %s contains no data?", msafile);
    else if (status != eslOK)         esl_fatal("Failed to guess alphabet (error code %d)\n", status);
  }
    
  /* Now you know how to create the alphabet */
  abc = esl_alphabet_Create(alphatype);

  /* Then you set the msafile into digital mode */
  esl_msafile_SetDigital(afp, abc);

  /* Now the MSA's that you read are digital data in msa->ax, not text in msa->aseq */
  while ((status = esl_msa_Read(afp, &msa)) == eslOK)
    {
      nali++;
      printf("alignment %5d: %15s: %6d seqs, %5d columns\n", 
	     nali, msa->name, msa->nseq, msa->alen);
      esl_msa_Write(stdout, msa, eslMSAFILE_STOCKHOLM);
      esl_msa_Destroy(msa);
    }
  if      (status == eslEFORMAT) esl_fatal("alignment file %s: %s\n", afp->fname, afp->errbuf);
  else if (status != eslEOF)     esl_fatal("alignment file %s: read failed, error %d\n", afp->fname, status);
  else if (nali == 0)            esl_fatal("alignment file %s: %s\n", afp->fname, afp->errbuf);

  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  exit(0);
}
#endif /*eslMSA_EXAMPLE2*/
#endif /*eslAUGMENT_ALPHABET*/
/*::cexcerpt::msa_example2::end::*/
/*------------------------ end of examples -----------------------*/

 


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
