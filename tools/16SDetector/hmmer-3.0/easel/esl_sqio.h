/* Unaligned sequence file i/o.
 * 
 * SVN $Id: esl_sqio.h 509 2010-02-07 22:56:55Z eddys $
 */
#ifndef ESL_SQIO_INCLUDED
#define ESL_SQIO_INCLUDED

#include <stdio.h>

#include "esl_sqio_ascii.h"
#ifdef eslAUGMENT_NCBI
#include "esl_sqio_ncbi.h"
#endif
#include "esl_sq.h"

#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif
#ifdef eslAUGMENT_MSA
#include "esl_msa.h"
#endif

/* ESL_SQDATA:
 * Data for different sequence formats.
 */
typedef union {
  ESL_SQASCII_DATA ascii;
#ifdef eslAUGMENT_NCBI
  ESL_SQNCBI_DATA  ncbi;
#endif
} ESL_SQDATA;

/* ESL_SQFILE:
 * An open sequence file for reading.
 */
typedef struct esl_sqio_s {
  char *filename;	      /* Name of file (for diagnostics)           */

  /* In digital mode, we have an alphabet ptr                             */
  int   do_digital;	      /* TRUE if we're reading in digital mode    */
#if defined(eslAUGMENT_ALPHABET)  
  const ESL_ALPHABET *abc;
#else
  void               *abc;
#endif

  /* Format-specific configuration                                        */
  int     format;	      /* Format code of this file                 */
  ESL_DSQ inmap[128];	      /* an input map, 0..127                     */

  /* function pointers to format specific routines                        */
  int   (*position)        (struct esl_sqio_s *sqfp, off_t offset);
  void  (*close)           (struct esl_sqio_s *sqfp);

  int   (*set_digital)     (struct esl_sqio_s *sqfp, const ESL_ALPHABET *abc);
  int   (*guess_alphabet)  (struct esl_sqio_s *sqfp, int *ret_type);

  int   (*read)            (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_info)       (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_seq)        (struct esl_sqio_s *sqfp, ESL_SQ *sq);
  int   (*read_window)     (struct esl_sqio_s *sqfp, int C, int W, ESL_SQ *sq);
  int   (*echo)            (struct esl_sqio_s *sqfp, const ESL_SQ *sq, FILE *ofp);

  int   (*read_block)      (struct esl_sqio_s *sqfp, ESL_SQ_BLOCK *sqBlock);

#ifdef eslAUGMENT_SSI
  int   (*open_ssi)        (struct esl_sqio_s *sqfp, const char *ssifile_hint);
  int   (*pos_by_key)      (struct esl_sqio_s *sqfp, const char *key);
  int   (*pos_by_number)   (struct esl_sqio_s *sqfp, int which);

  int   (*fetch)           (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_info)      (struct esl_sqio_s *sqfp, const char *key, ESL_SQ *sq);
  int   (*fetch_subseq)    (struct esl_sqio_s *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

  int   (*is_rewindable)   (const struct esl_sqio_s *sqfp);
  const char *(*get_error) (const struct esl_sqio_s *sqfp);

  ESL_SQDATA data;            /* format specific data                     */
} ESL_SQFILE;

/* Unaligned file format codes
 * These codes are coordinated with the msa module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqio, for unaligned formats
 *   - >100  is reserved for msa, for aligned formats
 */
#define eslSQFILE_UNKNOWN 0
#define eslSQFILE_FASTA   1
#define eslSQFILE_EMBL    2	/* EMBL/Swissprot/TrEMBL */
#define eslSQFILE_GENBANK 3	/* Genbank */
#define eslSQFILE_DDBJ    4	/* DDBJ (currently passed to Genbank parser */
#define eslSQFILE_UNIPROT 5     /* Uniprot (passed to EMBL parser) */
#define eslSQFILE_NCBI    6     /* NCBI (blast db) */


/* eslREADBUFSIZE is the fixed size of a block to bring in at one time,
 * in character-based (fread()) parsers (like the FASTA parser).
 */
#define eslREADBUFSIZE  4096

extern int  esl_sqfile_Open(const char *seqfile, int fmt, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_Position(ESL_SQFILE *sqfp, off_t offset);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);

#ifdef eslAUGMENT_ALPHABET
extern int  esl_sqfile_OpenDigital(const ESL_ALPHABET *abc, const char *filename, int format, const char *env, ESL_SQFILE **ret_sqfp);
extern int  esl_sqfile_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
extern int  esl_sqfile_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type);
#endif

const char  *esl_sqfile_GetErrorBuf(const ESL_SQFILE *sqfp);
extern int   esl_sqfile_IsRewindable(const ESL_SQFILE *sqfp);
extern int   esl_sqio_Ignore(ESL_SQFILE *sqfp, const char *ignoredchars);
extern int   esl_sqio_AcceptAs(ESL_SQFILE *sqfp, char *xchars, char readas);
extern int   esl_sqio_EncodeFormat(char *fmtstring);
extern char *esl_sqio_DecodeFormat(int fmt);
extern int   esl_sqio_IsAlignment(int fmt);

extern int   esl_sqio_Read      (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadInfo  (ESL_SQFILE *sqfp, ESL_SQ *sq);
extern int   esl_sqio_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
extern int   esl_sqio_Echo      (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);

extern int   esl_sqio_ReadBlock (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock);

#ifdef eslAUGMENT_SSI
extern int   esl_sqfile_OpenSSI         (ESL_SQFILE *sqfp, const char *ssifile_hint);
extern int   esl_sqfile_PositionByKey   (ESL_SQFILE *sqfp, const char *key);
extern int   esl_sqfile_PositionByNumber(ESL_SQFILE *sqfp, int which);

extern int   esl_sqio_Fetch      (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchInfo  (ESL_SQFILE *sqfp, const char *key, ESL_SQ *sq);
extern int   esl_sqio_FetchSubseq(ESL_SQFILE *sqfp, const char *source, int64_t start, int64_t end, ESL_SQ *sq);
#endif

extern int   esl_sqio_Write(FILE *fp, ESL_SQ *s, int format, int update);

#endif /*!ESL_SQIO_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
