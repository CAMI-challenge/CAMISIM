/* Unaligned ncbi sequence file i/o.
 * 
 * SVN $Id: esl_sqio_ncbi.h 361 2009-06-30 00:40:48Z farrarm $
 */
#ifndef ESL_SQIO_NCBI_INCLUDED
#define ESL_SQIO_NCBI_INCLUDED

#include <stdio.h>
#include "esl_sq.h"
#include "esl_sqio.h"

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif


/* forward declaration */
struct esl_sqio_s;

/* ESL_SQNCBI:
 * An open sequence file for reading.
 */
typedef struct esl_sqncbi_s {
  FILE      *fppin;                /* Open .pin file ptr                       */
  FILE      *fpphr;                /* Open .phr file ptr                       */
  FILE      *fppsq;                /* Open .psq file ptr                       */
  char       errbuf[eslERRBUFSIZE];/* parse error mesg.  Size must match msa.h */

  char      *title;                /* database title                           */
  int        version;              /* database version                         */
  char      *timestamp;            /* time stamp of database creation          */

  uint32_t   num_seq;              /* number of sequences in the database      */
  uint64_t   total_res;            /* total number of residues                 */
  uint32_t   max_seq;              /* longest sequence in the database         */

  uint32_t   hdr_off;              /* disk offset in .pin to header index      */
  uint32_t   seq_off;              /* disk offset to .pin to sequence index    */
  uint32_t   amb_off;              /* disk offset to .pin to ambiguous index   */
  
  int        index;                /* current sequence index in the database   */

  uint32_t   cur_indexes;          /* start of indexes currently loaded        */
  uint32_t  *hdr_indexes;          /* block of header indexes from .pin        */
  uint32_t  *seq_indexes;          /* block of header indexes from .pin        */
  uint32_t  *amb_indexes;          /* block of header indexes from .pin        */

  /* information for the current header */
  unsigned char *hdr_buf;          /* buffer for holding unparsed header       */
  unsigned char *hdr_ptr;          /* current parser position                  */
  int            hdr_alloced;      /* size of the allocated buffer             */
  int            hdr_size;         /* size of the current header               */
  uint32_t       hdr_fpos;         /* offset into the .phr file                */

  /* information on the current sequence */
  uint32_t       seq_apos;         /* position of ambiguity table              */
  uint32_t       seq_alen;         /* size of ambiguity table                  */
  uint32_t       seq_cpos;         /* current position in ambiguity table      */
  int32_t        seq_L;            /* true sequence length                     */

  /* alphabet used to convert ncbi to hmmer to ascii */
  int            alphatype;        /* amino or dna                             */
  char          *alphasym;         /* string of residues                       */

} ESL_SQNCBI_DATA;


extern int  esl_sqncbi_Open(char *seqfile, int format, struct esl_sqio_s *sqfp);


#endif /*!ESL_SQIO_NCBI_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
