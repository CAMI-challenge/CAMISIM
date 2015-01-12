/* Unaligned ncbi sequence file i/o.
 * 
 * Contents:
 *    1. An <ESL_SQFILE> object, in text mode.
 *    2. An <ESL_SQFILE> object, in digital mode. [with <alphabet>]
 *    3. Miscellaneous routines.
 *    4. Sequence reading (sequential).
 *    5. Parsing routines
 *    6. Copyright and license.
 * 
 * MSF, Mon Dec 10, 2009
 * SVN $Id$
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_ENDIAN_H
#include <endian.h>
#endif

#include "easel.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"	/* alphabet aug adds digital sequences */
#endif 
#include "esl_sqio.h"
#include "esl_sq.h"

#ifndef htobe32
#ifdef  WORDS_BIGENDIAN
#define htobe32(x) (x)
#else
#define htobe32(x) \
     ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) |		      \
      (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24))
#endif
#endif

/* format specific routines */
static int   sqncbi_Position       (ESL_SQFILE *sqfp, off_t offset);
static void  sqncbi_Close          (ESL_SQFILE *sqfp);
static int   sqncbi_SetDigital     (ESL_SQFILE *sqfp, const ESL_ALPHABET *abc);
static int   sqncbi_GuessAlphabet  (ESL_SQFILE *sqfp, int *ret_type);
static int   sqncbi_Read           (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadInfo       (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadSequence   (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int   sqncbi_ReadWindow     (ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq);
static int   sqncbi_ReadBlock      (ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock);
static int   sqncbi_Echo           (ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp);

static int   sqncbi_IsRewindable   (const ESL_SQFILE *sqfp);
static const char *sqncbi_GetError (const ESL_SQFILE *sqfp);

/* common routines for processing ncbi database */
static int  get_offsets         (ESL_SQNCBI_DATA *ncbi, int inx, off_t *hdr, off_t *seq, off_t *amb);
static int  read_amino          (ESL_SQFILE *sqfp, ESL_SQ *sq);
static int  read_dna            (ESL_SQFILE *sqfp, ESL_SQ *sq, off_t amb);
static int  read_nres_amino     (ESL_SQFILE *sqfp, ESL_SQ *sq, int len, uint64_t *nres);
static int  read_nres_dna       (ESL_SQFILE *sqfp, ESL_SQ *sq, int len, uint64_t *nres);

static int  inmap_ncbi          (ESL_SQFILE *sqfp);
static int  inmap_ncbi_amino    (ESL_SQFILE *sqfp);
static int  inmap_ncbi_dna      (ESL_SQFILE *sqfp);

static int  sqncbi_OpenAmino(ESL_SQNCBI_DATA *ncbi, char *filename);
static int  sqncbi_OpenDna  (ESL_SQNCBI_DATA *ncbi, char *filename);

/* parsing routines */
static int  parse_header              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_def_line            (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_seq_id              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_textseq_id          (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_object_id           (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_dbtag               (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_patent_seq_id       (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_id_pat              (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_pdb_seq_id          (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_date_std            (ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq);
static int  parse_string              (ESL_SQNCBI_DATA *ncbi, int max, char **str);
static int  parse_integer             (ESL_SQNCBI_DATA *ncbi, int *value);
static int  ignore_sequence_of_integer(ESL_SQNCBI_DATA *ncbi);

#define INDEX_TABLE_SIZE      1024
#define INIT_HDR_BUFFER_SIZE  2048

#define NCBI_VERSION_4             4
#define NCBI_DNA_DB                0
#define NCBI_AMINO_DB              1

/* set the max residue count to 1 meg when reading a block */
#define MAX_RESIDUE_COUNT (1024 * 1024)

/*****************************************************************
 *# 1. An <ESL_SQFILE> object, in text mode.
 *****************************************************************/ 

/* Function:  esl_sqncbi_Open()
 * Synopsis:  Open a sequence file for reading.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Open a sequence file <filename> for reading. 
 *            The opened <ESL_SQFILE> is returned through <ret_sqfp>.
 * 
 *            The .pin, .phr and .psq files are required for the
 *            open function to succeed.  Only protien version 4
 *            databases are currently supported.
 *            
 * Returns:   <eslOK> on success, and <*ret_sqfp> points to a new
 *            open <ESL_SQFILE>. Caller deallocates this object with
 *            <esl_sqfile_Close()>. 
 *            
 *            Returns <eslENOTFOUND> if <filename> can't be found or
 *            opened.  Returns <eslEFORMAT> if the file is empty, or
 *            if autodetection is attempted and the format can't be
 *            determined.  On any error condition, <*ret_sqfp> is
 *            returned NULL.
 *             
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_sqncbi_Open(char *filename, int format, ESL_SQFILE *sqfp)
{
  int         status = eslOK;	/* return status from an ESL call */

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  /* before we go any further, make sure we can handle the format */
  if (format != eslSQFILE_NCBI && format != eslSQFILE_UNKNOWN) return eslENOTFOUND;

  ncbi->fppin        = NULL;
  ncbi->fpphr        = NULL;
  ncbi->fppsq        = NULL;

  ncbi->title        = NULL;
  ncbi->timestamp    = NULL;

  ncbi->index        = -1;

  ncbi->hdr_off      = -1;
  ncbi->seq_off      = -1;
  ncbi->amb_off      = -1;

  ncbi->cur_indexes  = -1;
  ncbi->hdr_indexes  = NULL;
  ncbi->seq_indexes  = NULL;
  ncbi->amb_indexes  = NULL;

  ncbi->hdr_buf      = NULL;

  ncbi->amb_off      = 0;

  ncbi->alphatype    = eslUNKNOWN;
  ncbi->alphasym     = NULL;

  if ((status = sqncbi_OpenAmino(ncbi, filename)) != eslOK) {
    sqncbi_Close(sqfp);
    if ((status = sqncbi_OpenDna(ncbi, filename)) != eslOK) goto ERROR;
  }

  sqfp->format = eslSQFILE_NCBI;
  if ((status = inmap_ncbi(sqfp)) != eslOK) goto ERROR;

  /* initialize the function pointers for the ncbi routines */
  sqfp->position          = &sqncbi_Position;
  sqfp->close             = &sqncbi_Close;

  sqfp->set_digital       = &sqncbi_SetDigital;
  sqfp->guess_alphabet    = &sqncbi_GuessAlphabet;

  sqfp->is_rewindable     = &sqncbi_IsRewindable;

  sqfp->read              = &sqncbi_Read;
  sqfp->read_info         = &sqncbi_ReadInfo;
  sqfp->read_seq          = &sqncbi_ReadSequence;
  sqfp->read_window       = &sqncbi_ReadWindow;
  sqfp->echo              = &sqncbi_Echo;

  sqfp->read_block        = &sqncbi_ReadBlock;

  sqfp->get_error         = &sqncbi_GetError;

  return eslOK;

 ERROR:
  sqncbi_Close(sqfp); 
  return status;
}


static int
sqncbi_OpenAmino(ESL_SQNCBI_DATA *ncbi, char *filename)
{
  int         status = eslOK;	/* return status from an ESL call */
  int         len;

  uint32_t    info[4];
  char       *name = NULL;

  len = strlen(filename);
  ESL_ALLOC(name, sizeof(char) * (len+5));
  strcpy(name, filename);

  /* Check the current working directory first. */
  strcpy(name+len, ".pin");
  if ((ncbi->fppin = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".phr");
  if ((ncbi->fpphr = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".psq");
  if ((ncbi->fppsq = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }

  /* make sure we are looking at a version 4 protien db.
   * the values are stored in big endian, so we will just
   * against the values in big endian format
   */

  if (fread(&info[0], sizeof(uint32_t), 3, ncbi->fppin) != 3) status = eslFAIL;
  if (htobe32(info[0]) != NCBI_VERSION_4)                     status = eslEFORMAT;
  if (htobe32(info[1]) != NCBI_AMINO_DB)                      status = eslEUNIMPLEMENTED;

  if (status != eslOK) goto ERROR;
  ncbi->version = htobe32(info[0]);
  ncbi->alphatype = eslAMINO;
  ncbi->index = 0;

  /* read the database title */
  len = htobe32(info[2]);
  ESL_ALLOC(ncbi->title, sizeof(char) * (len + 1));
  if (fread(ncbi->title, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->title[len] = 0;

  /* read the database time stamp */
  if (fread(&info[0], sizeof(uint32_t), 1, ncbi->fppin) != 1) { status = eslFAIL; goto ERROR; }
  len = htobe32(info[0]);
  ESL_ALLOC(ncbi->timestamp, sizeof(char) * (len + 1));
  if (fread(ncbi->timestamp, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->timestamp[len] = 0;

  /* read in database stats */
  if (fread(&info[0], sizeof(uint32_t), 4, ncbi->fppin) != 4) { status = eslFAIL; goto ERROR; }
  ncbi->num_seq   = htobe32(info[0]);
  ncbi->total_res = *(uint64_t *)(info+1);
  ncbi->max_seq   = htobe32(info[3]);

  /* save the offsets to the index tables */
  ncbi->hdr_off = ftell(ncbi->fppin);
  ncbi->seq_off = ncbi->hdr_off + sizeof(uint32_t) * (ncbi->num_seq + 1);

  /* allocate buffers used in parsing the database files */
  ESL_ALLOC(ncbi->hdr_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);
  ESL_ALLOC(ncbi->seq_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);

  ncbi->hdr_alloced = INIT_HDR_BUFFER_SIZE;
  ESL_ALLOC(ncbi->hdr_buf, sizeof(char) * INIT_HDR_BUFFER_SIZE);

  /* skip the first sentinal byte in the .psq file */
  fgetc(ncbi->fppsq);

  if (name != NULL) free(name);

  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  return status;
}


static int
sqncbi_OpenDna(ESL_SQNCBI_DATA *ncbi, char *filename)
{
  int         status = eslOK;	/* return status from an ESL call */
  int         len;

  uint32_t    info[4];
  char       *name = NULL;

  len = strlen(filename);
  ESL_ALLOC(name, sizeof(char) * (len+5));
  strcpy(name, filename);

  /* Check the current working directory first. */
  strcpy(name+len, ".nin");
  if ((ncbi->fppin = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".nhr");
  if ((ncbi->fpphr = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }
  strcpy(name+len, ".nsq");
  if ((ncbi->fppsq = fopen(name, "r")) == NULL) {
    status = eslENOTFOUND; 
    goto ERROR;
  }

  /* make sure we are looking at a version 4 dna db. */

  if (fread(&info[0], sizeof(uint32_t), 3, ncbi->fppin) != 3) status = eslFAIL;
  if (htobe32(info[0]) != NCBI_VERSION_4)                     status = eslEFORMAT;
  if (htobe32(info[1]) != NCBI_DNA_DB)                        status = eslEUNIMPLEMENTED;

  if (status != eslOK) goto ERROR;
  ncbi->version = htobe32(info[0]);
  ncbi->alphatype = eslDNA;
  ncbi->index = 0;

  /* read the database title */
  len = htobe32(info[2]);
  ESL_ALLOC(ncbi->title, sizeof(char) * (len + 1));
  if (fread(ncbi->title, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->title[len] = 0;

  /* read the database time stamp */
  if (fread(&info[0], sizeof(uint32_t), 1, ncbi->fppin) != 1) { status = eslFAIL; goto ERROR; }
  len = htobe32(info[0]);
  ESL_ALLOC(ncbi->timestamp, sizeof(char) * (len + 1));
  if (fread(ncbi->timestamp, sizeof(char), len, ncbi->fppin) != len) { status = eslFAIL; goto ERROR; }
  ncbi->timestamp[len] = 0;

  /* read in database stats */
  if (fread(&info[0], sizeof(uint32_t), 4, ncbi->fppin) != 4) { status = eslFAIL; goto ERROR; }
  ncbi->num_seq   = htobe32(info[0]);
  ncbi->total_res = *(uint64_t *)(info+1);
  ncbi->max_seq   = htobe32(info[3]);

  /* save the offsets to the index tables */
  ncbi->hdr_off = ftell(ncbi->fppin);
  ncbi->seq_off = ncbi->hdr_off + sizeof(uint32_t) * (ncbi->num_seq + 1);
  ncbi->amb_off = ncbi->seq_off + sizeof(uint32_t) * (ncbi->num_seq + 1);

  /* allocate buffers used in parsing the database files */
  ESL_ALLOC(ncbi->hdr_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);
  ESL_ALLOC(ncbi->seq_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);
  ESL_ALLOC(ncbi->amb_indexes, sizeof(uint32_t) * INDEX_TABLE_SIZE);

  ncbi->hdr_alloced = INIT_HDR_BUFFER_SIZE;
  ESL_ALLOC(ncbi->hdr_buf, sizeof(char) * INIT_HDR_BUFFER_SIZE);

  /* skip the first sentinal byte in the .nsq file */
  fgetc(ncbi->fppsq);

  if (name != NULL) free(name);

  return eslOK;

 ERROR:
  if (name != NULL) free(name);
  return status;
}


/* Function:  sqncbi_Position()
 * Synopsis:  Reposition an open sequence file to an offset.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reposition an open <sqfp> to offset <offset>.
 *            <offset> for the ncbi db format specified the sequence
 *            index, not file offset.  Both the sequence and header
 *            files are repositioned.
 *            
 * Returns:   <eslOK>     on success;
 *
 * Throws:    <eslESYS> if the fseeko() or fread() call fails.
 *            On errors, the state of <sqfp> is indeterminate, and
 *            it should not be used again.
 */
static int
sqncbi_Position(ESL_SQFILE *sqfp, off_t offset)
{
  off_t    hdr_start;
  off_t    seq_start;

  int      status;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if ((status = get_offsets(ncbi, offset, &hdr_start, &seq_start, NULL)) != eslOK) return status;

  if (fseek(ncbi->fpphr, hdr_start, SEEK_SET) != 0) return eslESYS;
  if (fseek(ncbi->fppsq, seq_start, SEEK_SET) != 0) return eslESYS;

  ncbi->index = offset;

  return eslOK;
}


/* Function:  sqncbi_Close()
 * Synopsis:  Close a sequence file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Closes an open <sqfp>.
 *
 * Returns:   (void).
 */
static void
sqncbi_Close(ESL_SQFILE *sqfp)
{
  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->title != NULL)       free(ncbi->title);
  if (ncbi->timestamp != NULL)   free(ncbi->timestamp);

  if (ncbi->hdr_buf != NULL)     free(ncbi->hdr_buf);

  if (ncbi->hdr_indexes != NULL) free(ncbi->hdr_indexes);
  if (ncbi->seq_indexes != NULL) free(ncbi->seq_indexes);
  if (ncbi->amb_indexes != NULL) free(ncbi->amb_indexes);

  if (ncbi->alphasym != NULL)    free(ncbi->alphasym);

  if (ncbi->fppin != NULL) fclose(ncbi->fppin);
  if (ncbi->fpphr != NULL) fclose(ncbi->fpphr);
  if (ncbi->fppsq != NULL) fclose(ncbi->fppsq);

  ncbi->fppin        = NULL;
  ncbi->fpphr        = NULL;
  ncbi->fppsq        = NULL;

  ncbi->title        = NULL;
  ncbi->timestamp    = NULL;

  ncbi->index        = -1;

  ncbi->hdr_off      = -1;
  ncbi->seq_off      = -1;
  ncbi->amb_off      = -1;

  ncbi->cur_indexes  = -1;
  ncbi->hdr_indexes  = NULL;
  ncbi->seq_indexes  = NULL;
  ncbi->amb_indexes  = NULL;

  ncbi->hdr_buf      = NULL;

  ncbi->alphatype    = eslUNKNOWN;
  ncbi->alphasym     = NULL;

  return;
}
/*------------------- SQNCBI open/close -----------------------*/


/*****************************************************************
 *# 2. An <ESL_SQFILE> object, in digital mode [with <alphabet>]
 *****************************************************************/
#ifdef eslAUGMENT_ALPHABET

/* Function:  sqncbi_SetDigital()
 * Synopsis:  Set an open <ESL_SQFILE> to read in digital mode.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Given an <ESL_SQFILE> that's already been opened,
 *            configure it to expect subsequent input to conform
 *            to the digital alphabet <abc>.
 *            
 *            Calling <esl_sqfile_Open(); esl_sqfile_SetDigital()> is
 *            equivalent to <esl_sqfile_OpenDigital()>. The two-step
 *            version is useful when you need a
 *            <esl_sqfile_GuessAlphabet()> call in between, guessing
 *            the file's alphabet in text mode before you set it to
 *            digital mode.
 *
 * Returns:   <eslOK> on success.
 */
static int
sqncbi_SetDigital(ESL_SQFILE *sqfp, const ESL_ALPHABET *abc)
{
  return eslOK;
}

/* Function:  sqncbi_GuessAlphabet()
 * Synopsis:  Guess the alphabet of an open <ESL_SQFILE>.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   The only ncbi db format supported is protein.
 *
 * Returns:   <eslOK> on success, and <*ret_type> is set to <eslAMINO>.
 */
static int
sqncbi_GuessAlphabet(ESL_SQFILE *sqfp, int *ret_type)
{
  *ret_type = sqfp->data.ncbi.alphatype;
  return eslOK;
}
#endif /*eslAUGMENT_ALPHABET*/
/*-------------- end, digital mode SQNCBI -------------------*/




/*****************************************************************
 *# 3. Miscellaneous routines 
 *****************************************************************/ 

/* Function:  sqncbi_IsRewindable()
 * Synopsis:  Return <TRUE> if <sqfp> can be rewound.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Returns <TRUE> if <sqfp> can be rewound (positioned 
 *            to an offset of zero), in order to read it a second
 *            time.
 */
static int
sqncbi_IsRewindable(const ESL_SQFILE *sqfp)
{
  return TRUE;
}

/* Function:  sqncbi_GetError()
 * Synopsis:  Returns error buffer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Return a pointer to the error buffer.
 */
static const char *
sqncbi_GetError(const ESL_SQFILE *sqfp)
{
  return sqfp->data.ncbi.errbuf;
}




/*****************************************************************
 *# 4. Sequence reading (sequential)
 *****************************************************************/ 

/* Function:  sqncbi_Read()
 * Synopsis:  Read the next sequence from a file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads the next sequence from open sequence file <sqfp> into 
 *            <sq>. Caller provides an allocated and initialized <s>, which
 *            will be internally reallocated if its space is insufficient.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character; the line number that the parse
 *            error occurs on is in <sqfp->linenumber>, and an informative
 *            error message is placed in <sqfp->errbuf>. 
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
static int
sqncbi_Read(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  off_t   amb;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start, &amb)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end, NULL)) != eslOK) return status;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->hoff = hdr_end;
  sq->eoff = seq_end;

  if (ncbi->alphatype == eslAMINO) 
    status = read_amino(sqfp, sq);
  else                             
    status = read_dna(sqfp, sq, amb);
  if (status != eslOK) return status;

  /* read and parse the ncbi header */
  ncbi->hdr_fpos = hdr_start;
  ncbi->hdr_size = hdr_end - hdr_start;
  if ((status = parse_header(ncbi, sq)) != eslOK) return status;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;
}


/* Function:  sqncbi_ReadInfo()
 * Synopsis:  Read sequence info, but not the sequence itself.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Read the next sequence from open sequence file <sqfp>,
 *            but don't store the sequence (or secondary structure).
 *            Upon successful return, <s> holds all the available 
 *            information about the sequence -- its name, accession,
 *            description, and overall length <sq->L>. 
 *            
 *            This is useful for indexing sequence files, where
 *            individual sequences might be ginormous, and we'd rather
 *            avoid reading complete seqs into memory.
 *
 * Returns:   <eslOK> on success.
 */
static int
sqncbi_ReadInfo(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start, NULL)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end, NULL)) != eslOK) return status;

  /* advance the sequence file pointer to point to the next sequence in case
   * the user calls a Read after a ReadInfo.  this will garuntee that the
   * header and sequence match up for the Read.
   */
  if (fseek(ncbi->fppsq, seq_end, SEEK_SET) != 0) return eslEFORMAT;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->hoff = hdr_end;
  sq->eoff = seq_end;

  /* figure out the sequence length */
  sq->L = -1;

  /* read and parse the ncbi header */
  ncbi->hdr_fpos = hdr_start;
  ncbi->hdr_size = hdr_end - hdr_start;
  if ((status = parse_header(ncbi, sq)) != eslOK) return status;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;
}


/* Function:  sqncbi_ReadSequence()
 * Synopsis:  Read the sequence, not the sequence header.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Read the next sequence from open sequence file <sqfp>,
 *            but not the header information.  Upon successful return,
 *            <s> holds all the sequence.
 *            
 *            This is useful reading binary formats and delaying the
 *            over heads of reading the sequence name until needed by
 *            the report generator.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 *
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
sqncbi_ReadSequence(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     status;

  off_t   hdr_start;
  off_t   hdr_end;
  off_t   seq_start;
  off_t   seq_end;

  off_t   amb;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start, &amb)) != eslOK) return status;
  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end, NULL)) != eslOK) return status;

  /* Disk offset bookkeeping */
  sq->idx  = ncbi->index;
  sq->roff = hdr_start;
  sq->doff = seq_start;
  sq->hoff = hdr_end;
  sq->eoff = seq_end;

  if (ncbi->alphatype == eslAMINO) 
    status = read_amino(sqfp, sq);
  else                             
    status = read_dna(sqfp, sq, amb);
  if (status != eslOK) return status;

  /* advance the header file pointer to point to the next header in case
   * the user calls a Read after a ReadSequence.  this will garuntee that
   * the header and sequence match up for the Read.
   */
  if (fseek(ncbi->fpphr, hdr_end, SEEK_SET) != 0) return eslEFORMAT;

  /* update the sequence index */
  ++ncbi->index;

  return eslOK;
}


/* Function:  sqncbi_ReadWindow()
 * Synopsis:  Read next window of sequence.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Read a next window of <W> residues from open file <sqfp>,
 *            keeping <C> residues from the previous window as
 *            context, and keeping previous annotation in the <sq>
 *            as before. 
 *            
 *            If this is the first window of a new sequence record,
 *            <C> is ignored (there's no previous context yet), and
 *            the annotation fields of the <sq> (name, accession, and
 *            description) are initialized by reading the sequence
 *            record's header. This is the only time the annotation
 *            fields are initialized.
 *            
 *            On return, <sq->dsq[]> contains the window and its
 *            context; residues <1..sq->C> are the previous context,
 *            and residues <sq->C+1..sq->n> are the new window.  The
 *            start and end coordinates of the whole <dsq[1..n]>
 *            (including context) in the original source sequence are
 *            <sq->start..sq->end>. (Or, for text mode sequences,
 *            <sq->seq[0..sq->C-1,sq->C..sq->n-1]>, while <start> and
 *            <end> coords are still <1..L>.)
 *
 *            When a sequence record is completed and no more data
 *            remain, <eslEOD> is returned, with an ``info'' <sq>
 *            structure (containing the annotation and the total
 *            sequence length <L>, but no sequence). (The total
 *            sequence length <L> is unknown in <sq> until this
 *            <eslEOD> return.)
 *            
 *            The caller may then do one of two things before calling
 *            <esl_sq_ReadWindow()> again; it can reset the sequence
 *            with <esl_sq_Reuse()> to continue reading the next
 *            sequence in the file, or it can set a negative <W> as a
 *            signal to read windows from the reverse complement
 *            (Crick) strand. Reverse complement reading only works
 *            for nucleic acid sequence. 
 *            
 *            If you read the reverse complement strand, you must read
 *            the whole thing, calling <esl_sqio_ReadWindow()> with
 *            negative <W> windows until <eslEOD> is returned again
 *            with an empty (info-only) <sq> structure. When that
 *            <EOD> is reached, the <sqfp> is repositioned at the
 *            start of the next sequence record; the caller should now
 *            <Reuse()> the <sq>, and the next <esl_sqio_ReadWindow()>
 *            call must have a positive <W>, corresponding to starting
 *            to read the Watson strand of the next sequence.
 *
 *            Note that the <ReadWindow()> interface is designed for
 *            an idiom of sequential reading of complete sequences in
 *            overlapping windows, possibly on both strands; if you
 *            want more freedom to move around in the sequence
 *            grabbing windows in another order, you can use the
 *            <FetchSubseq()> interface.
 *
 *            Reading the reverse complement strand requires file
 *            repositioning, so it will not work on non-repositionable
 *            streams like gzipped files or a stdin pipe. Moreover,
 *            for reverse complement input to be efficient, the
 *            sequence file should have consistent line lengths, 
 *            suitable for SSI's fast subsequence indexing.
 *            
 * Returns:   <eslOK> on success; <sq> now contains next window of
 *            sequence, with at least 1 new residue. The number
 *            of new residues is <sq->W>; <sq->C> residues are 
 *            saved from the previous window. Caller may now
 *            process residues <sq->dsq[sq->C+1]..sq->dsq[sq->n]>.
 *            
 *            <eslEOD> if no new residues were read for this sequence
 *            and strand, and <sq> now contains an empty info-only
 *            structure (annotation and <L> are valid). Before calling
 *            <esl_sqio_ReadWindow()> again, caller will either want
 *            to make <W> negative (to start reading the Crick strand
 *            of the current sequence), or it will want to reset the
 *            <sq> (with <esl_sq_Reuse()>) to go on the next sequence.
 *            
 *            <eslEOF> if we've already returned <eslEOD> before to
 *            signal the end of the previous seq record, and moreover,
 *            there's no more sequence records in the file.
 *            
 *            <eslEINVAL> if an invalid residue is found in the
 *            sequence, or if you attempt to take the reverse
 *            complement of a sequence that can't be reverse
 *            complemented.
 *
 * Throws:    <eslESYNTAX> if you try to read a reverse window before
 *            you've read forward strand.
 *            
 *            <eslECORRUPT> if something goes awry internally in the
 *            coordinate system.
 *            
 *            <eslEMEM> on allocation error.
 */
static int
sqncbi_ReadWindow(ESL_SQFILE *sqfp, int C, int W, ESL_SQ *sq)
{
  uint64_t  nres;
  int       status;

  off_t     hdr_start;
  off_t     hdr_end;
  off_t     seq_start;
  off_t     seq_end;

  off_t     amb;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  /* Negative W indicates reverse complement direction */
  if (W < 0)	
    {
      if (sq->L == -1) ESL_EXCEPTION(eslESYNTAX, "Can't read reverse complement until you've read forward strand");

      /* update the sequence index */
      if ((status = sqncbi_Position(sqfp, sq->idx)) != eslOK) 
	ESL_FAIL(eslEINVAL, ncbi->errbuf, "Unexpected error positioning datbase to sequence %d", ncbi->index+1);

      if (sq->end == 1) 
	{ /* last end == 1 means last window was the final one on reverse strand,
	   * so we're EOD; jump back to last forward position. 
	   */
	  sq->start      = 0;
	  sq->end        = 0;
	  sq->C          = 0;
	  sq->W          = 0;
	  sq->n          = 0;
	  /* sq->L stays as it is */
	  if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL;
	  else                 sq->seq[0] = '\0';

	  /* update the sequence index */
	  if ((status = sqncbi_Position(sqfp, ncbi->index+1)) != eslOK) 
	    ESL_FAIL(eslEINVAL, ncbi->errbuf, "Unexpected error positioning datbase to sequence %d", ncbi->index+1);

	  return eslEOD;
	}

      /* If s == 0, we haven't read any reverse windows yet; 
       * init reading from sq->L
       */
      W = -W;
      if (sq->start == 0)	
	{
	  sq->start        = ESL_MAX(1, (sq->L - W + 1)); 
	  sq->end          = sq->start;
	  sq->C            = 0;
	}
      else
	{ /* Else, we're continuing to next window; prv was <end>..<start> */
	  sq->C     = ESL_MIN(C, sq->L - sq->end + 1);  /* based on prev window's end */
	  sq->start = ESL_MAX(1, (sq->end - W));
	  W         = sq->end - sq->start + sq->C;
	  sq->end   = sq->start;
	}

      /* grab the subseq and rev comp it */
      if ((status = esl_sq_GrowTo(sq, W)) != eslOK) return status;
      sq->n = 0;
      if (ncbi->alphatype == eslAMINO) status = read_nres_amino(sqfp, sq, W, &nres);
      else                             status = read_nres_dna(sqfp, sq, W, &nres);
      
      if (status != eslOK || nres != W) {
	ESL_EXCEPTION(eslECORRUPT, "Failed to extract %d..%d", sq->start, sq->end);
      } else {
	sq->end        = sq->start + nres - 1;	  
	sq->W          = nres - sq->C;	  
      }

      status = esl_sq_ReverseComplement(sq);
      if      (status    == eslEINVAL) ESL_FAIL(eslEINVAL, ncbi->errbuf, "can't reverse complement that seq - it's not DNA/RNA");
      else if (status    != eslOK)     return status;

      return eslOK;
    } 

  /* Else, we're reading the forward strand */
  else 
    { /* sq->start == 0 means we haven't read any windows on this sequence yet...
       * it's a new record, and we need to initialize with the header and
       * the first window. This is the only case that we're allowed to return
       * EOF from.
       */
      if (sq->start == 0)
	{
	  if (ncbi->index >= ncbi->num_seq) return eslEOF;

	  /* get the sequence and header offsets */
	  if ((status = get_offsets(ncbi, ncbi->index, &hdr_start, &seq_start, &amb)) != eslOK) return status;
	  if ((status = get_offsets(ncbi, ncbi->index + 1, &hdr_end, &seq_end, NULL)) != eslOK) return status;

	  /* Disk offset bookkeeping */
	  sq->idx  = ncbi->index;
	  sq->roff = hdr_start;
	  sq->doff = seq_start;
	  sq->hoff = hdr_end;
	  sq->eoff = seq_end;

	  /* sequence bookkeeping */
	  if (ncbi->alphatype == eslAMINO) {
	    ncbi->seq_apos = 0;
	    ncbi->seq_alen = 0;
	  } else {
	    ncbi->seq_apos = amb;
	    ncbi->seq_alen = seq_end - amb;
	  }
	  ncbi->seq_cpos = -1;
	  ncbi->seq_L    = -1;

	  /* read and parse the ncbi header */
	  ncbi->hdr_fpos = hdr_start;
	  ncbi->hdr_size = hdr_end - hdr_start;
	  if ((status = parse_header(ncbi, sq)) != eslOK) return status;

	  sq->start    = 1;
	  sq->C        = 0;	/* no context in first window                   */
	  sq->L        = -1;	/* won't be known 'til EOD.                     */
	  ncbi->seq_L  = -1;	/* init to 0, so we can count residues as we go */
	  esl_sq_SetSource(sq, sq->name);
	}
      else
	{ /* else we're reading a window other than first; slide context over. */
	  sq->C = ESL_MIN(C, sq->n);
	  if (sq->seq != NULL) memmove(sq->seq,   sq->seq + sq->n - sq->C,     sq->C);
	  else                 memmove(sq->dsq+1, sq->dsq + sq->n - sq->C + 1, sq->C);
	  sq->start = sq->end - sq->C + 1;
	  sq->n = C;
	}      

      if ((status = esl_sq_GrowTo(sq, C+W)) != eslOK)                return status; /* EMEM    */
      if (ncbi->alphatype == eslAMINO) status = read_nres_amino(sqfp, sq, W, &nres);
      else                             status = read_nres_dna(sqfp, sq, W, &nres);

      if (status == eslEOD)	
	{
	  sq->start  = 0;
	  sq->end    = 0;
	  sq->C      = 0;
	  sq->W      = 0;
	  sq->n      = 0;

	  if (sq->dsq != NULL) sq->dsq[1] = eslDSQ_SENTINEL; /* erase the saved context */
	  else                 sq->seq[0] = '\0';

	  /* update the sequence index */
	  if ((status = sqncbi_Position(sqfp, ncbi->index+1)) != eslOK) 
	    ESL_FAIL(eslEINVAL, ncbi->errbuf, "Unexpected error positioning datbase to sequence %d", ncbi->index+1);

	  return eslEOD;
	}
      else if (status == eslOK)
	{ /* Forward strand is still in progress. <= W residues were read. Return eslOK. */
	  sq->end        = sq->start + sq->C + nres - 1;	  
	  sq->W          = nres;	  
	  return eslOK;
	}
      else return status;	/* EFORMAT,EMEM */
    }
  /*NOTREACHED*/
  return eslOK;
}

/* Function:  sqncbi_ReadBlock()
 * Synopsis:  Read the next block of sequences from a file.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads a block of sequences from open sequence file <sqfp> into 
 *            <sqBlock>.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <sqBlock>.
 * 
 *            Returns <eslEOF> when there is no sequence left in the
 *            file (including first attempt to read an empty file).
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character;
 *
 * Throws:    <eslEMEM> on allocation failure;
 *            <eslEINCONCEIVABLE> on internal error.
 */
static int
sqncbi_ReadBlock(ESL_SQFILE *sqfp, ESL_SQ_BLOCK *sqBlock)
{
  int     i;
  int     size = 0;
  int     status = eslOK;

  sqBlock->count = 0;
  for (i = 0; i < sqBlock->listSize && size < MAX_RESIDUE_COUNT; ++i)
    {
      status = sqncbi_Read(sqfp, sqBlock->list + i);
      if (status != eslOK) break;
      size += sqBlock->list[i].n;
      ++sqBlock->count;
    }

  /* EOF will be returned only in the case were no sequences were read */
  if (status == eslEOF && i > 0) status = eslOK;

  return status;
}

/* Function:  sqncbi_Echo()
 * Synopsis:  Echo a sequence's record onto output stream.
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Returns:   <eslEUNIMPLEMENTED>.
 */
static int
sqncbi_Echo(ESL_SQFILE *sqfp, const ESL_SQ *sq, FILE *ofp)
{
  ESL_EXCEPTION(eslEINVAL, "can't Echo() a sequence from NCBI database");
  return eslEUNIMPLEMENTED;
}
/*------------------ end, sequential sequence input -------------*/



/* Function:  get_offsets()
 * Synopsis:  Return the header and sequence offsets
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   For sequence <inx> reads the offsets in the sequence
 *            and header files.  If <hdr> or <seq> are non-null, the
 *            offsets for the .phr and .psq files will be set.
 *
 *            Before reading the offsets from the file, check if the
 *            current offset is cached.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there is an error reading the index file.
 */
static int
get_offsets(ESL_SQNCBI_DATA *ncbi, int inx, off_t *hdr, off_t *seq, off_t *amb)
{
  int        cnt;

  uint32_t   offset;
  uint32_t   start;
  uint32_t   end;

  if (inx < 0 || inx > ncbi->num_seq) return eslEINVAL;

  start = ncbi->cur_indexes;
  end   = start + INDEX_TABLE_SIZE - 1;

  if (ncbi->cur_indexes == -1 || inx < start || inx > end) {

    /* when calculating the count be sure to take into account the fact that the
     * index tables contain one index more that the number of sequences and this
     * last index is used to point to the end of the last header and sequences.
     */
    cnt = ncbi->num_seq - inx + 1;
    cnt = (cnt > INDEX_TABLE_SIZE) ? INDEX_TABLE_SIZE : cnt;

    offset = ncbi->hdr_off + (sizeof(uint32_t) * inx);
    if (fseek(ncbi->fppin, offset, SEEK_SET) != 0) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error seeking header index %d\n", offset);
    }
    if (fread(ncbi->hdr_indexes, sizeof(uint32_t), cnt, ncbi->fppin) != cnt) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error reading header index %d(%d)\n", offset, cnt);
    }

    offset = ncbi->seq_off + (sizeof(uint32_t) * inx);
    if (fseek(ncbi->fppin, offset, SEEK_SET) != 0) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error seeking sequence index %d\n", offset);
    }
    if (fread(ncbi->seq_indexes, sizeof(uint32_t), cnt, ncbi->fppin) != cnt) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error reading sequence index %d(%d)\n", offset, cnt);
    }

    if (ncbi->alphatype == eslDNA) {
      offset = ncbi->amb_off + (sizeof(uint32_t) * inx);
      if (fseek(ncbi->fppin, offset, SEEK_SET) != 0) {
	ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error seeking ambiguity index %d\n", offset);
      }
      if (fread(ncbi->amb_indexes, sizeof(uint32_t), cnt, ncbi->fppin) != cnt) {
	ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Error reading ambiguity index %d(%d)\n", offset, cnt);
      }
    }

    ncbi->cur_indexes = inx;
  }

  inx -= ncbi->cur_indexes;
  if (hdr != NULL) *hdr = htobe32(ncbi->hdr_indexes[inx]);
  if (seq != NULL) *seq = htobe32(ncbi->seq_indexes[inx]);

  if (amb != NULL && ncbi->alphatype == eslDNA) {
    *amb = htobe32(ncbi->amb_indexes[inx]);
  }

  return eslOK;
}

/* Function:  read_amino()
 * Synopsis:  Read in the amino sequence
 * Incept:    MSF, Wed Jan 27, 2010 [Janelia]
 *
 * Purpose:   Read and translate the amino acid sequence.
 */
static int
read_amino(ESL_SQFILE *sqfp, ESL_SQ *sq)
{
  int     inx;
  int     size;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  size = sq->eoff - sq->doff;

  /* figure out the sequence length */
  if (esl_sq_GrowTo(sq, size) != eslOK) return eslEMEM;

  /* figure out if the sequence is in digital mode or not */
  if (sq->dsq != NULL) {
    ESL_DSQ *ptr = sq->dsq + 1;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      ++ptr;
    }
    *ptr = eslDSQ_SENTINEL;
  } else {
    char *ptr = sq->seq;
    if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
    for (inx = 0; inx < size - 1; ++inx) {
      *ptr = sqfp->inmap[(int) *ptr];
      *ptr = ncbi->alphasym[(int) *ptr];
      ++ptr;
    }
    *ptr = '\0';
  }

  sq->start = 1;
  sq->end   = size - 1;
  sq->C     = 0;
  sq->W     = size - 1;
  sq->L     = size - 1;
  sq->n     = size - 1;

  return eslOK;
}


/* Function:  read_dna()
 * Synopsis:  Read in the dna sequence
 * Incept:    MSF, Wed Jan 27, 2010 [Janelia]
 *
 * Purpose:   Read and translate the dna sequence.
 */
static int
read_dna(ESL_SQFILE *sqfp, ESL_SQ *sq, off_t amb)
{
  int     inx;
  int     cnt;
  int     off;
  int     size;
  int     text;
  int     status;

  int     remainder;
  int     length;
  int     ssize;
  int     n;

  char   *ptr;
  void   *t;

  unsigned char c;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  /* calculate the max size of the sequence.  It is most likely
   * a bit smaller, but we need to read the sequence in first
   * before the real size can be figured out.
   */
  size = sq->eoff - sq->doff;
  if (ncbi->hdr_alloced < size) {
    while (ncbi->hdr_alloced < size) ncbi->hdr_alloced += ncbi->hdr_alloced;
    ESL_RALLOC(ncbi->hdr_buf, t, sizeof(char) * ncbi->hdr_alloced);
  }
  if (fread(ncbi->hdr_buf, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;

  ssize     = amb - sq->doff - 1;
  remainder = *(ncbi->hdr_buf + ssize) & 0x03;
  length    = ssize * 4 + remainder;

  /* figure out the sequence length */
  if (esl_sq_GrowTo(sq, length) != eslOK) return eslEMEM;

  /* figure out if the sequence is in digital mode or not */
  if (sq->dsq != NULL) {
    text = FALSE;
    ptr = (char *)sq->dsq + 1;
  } else {
    text = TRUE;
    ptr = sq->seq;
  }

  for (inx = 0; inx < ssize; ++inx) {
    c = ncbi->hdr_buf[inx];
    n = 1 << ((c >> 6) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 4) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 2) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 0) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  /* handle the remainder */
  c = ncbi->hdr_buf[inx];
  for (inx = 0; inx < remainder; ++inx) {
    n = 1 << ((c >> (6 - inx * 2)) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  *ptr = (text) ? '\0' : eslDSQ_SENTINEL;

  /* skip past the count and start processing the abmiguity table */
  ssize = amb - sq->doff + 4;
  ptr = (text) ? sq->seq : (char *)sq->dsq + 1;

  while (ssize < size) {
    /* get the ambiguity character */
    n = ((ncbi->hdr_buf[ssize] >> 4) & 0x0f);
    c = sqfp->inmap[n];
    if (text) c = ncbi->alphasym[(int) c];

    /* get the repeat count */
    cnt = (ncbi->hdr_buf[ssize] & 0x0f) + 1;

    /* get the offset */
    off = ncbi->hdr_buf[ssize+1];
    off = (off << 8) | ncbi->hdr_buf[ssize+2];
    off = (off << 8) | ncbi->hdr_buf[ssize+3];

    for (inx = 0; inx < cnt; ++inx) ptr[off+inx] = c;

    ssize += 4;
  }

  sq->start = 1;
  sq->end   = length;
  sq->C     = 0;
  sq->W     = length;
  sq->L     = length;
  sq->n     = length;

  return eslOK;

 ERROR:
  return eslEMEM;
}


/* Function:  read_nres_amino()
 * Synopsis:  Read in the amino sequence
 * Incept:    MSF, Wed Jan 27, 2010 [Janelia]
 *
 * Purpose:   Read and translate the dna sequence.
 */
static int
read_nres_amino(ESL_SQFILE *sqfp, ESL_SQ *sq, int len, uint64_t *nres)
{
  int     inx;
  int     off;
  int     size;

  char   *ptr;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  /* if we don't know the sequence length, figure it out */
  if (ncbi->seq_L == -1) ncbi->seq_L = sq->eoff - sq->doff - 1;

  /* check if we are at the end */
  if (sq->start + sq->n > ncbi->seq_L) {
    if (nres != NULL) *nres = 0;
    sq->L = ncbi->seq_L;
    return eslEOD;
  }

  /* figure out if the sequence is in digital mode or not */
  ptr = (sq->dsq != NULL) ? (char *)sq->dsq + 1 : sq->seq;
  ptr += sq->n;

  /* calculate where to start reading from */
  off   = sq->doff + sq->start + sq->n - 1;

  /* calculate the size to read */
  size = ncbi->seq_L - (sq->start + sq->n - 1);
  size = (size > len) ? len : size;

  /* seek to the windows location and read into the buffer */
  if (fseek(ncbi->fppsq, off, SEEK_SET) != 0) return eslESYS;
  if (fread(ptr, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;

  /* figure out if the sequence is in digital mode or not */
  for (inx = 0; inx < size; ++inx) {
    *ptr = sqfp->inmap[(int) *ptr];
    if (sq->dsq == NULL) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  *ptr = (sq->dsq == NULL) ? '\0' : eslDSQ_SENTINEL;

  sq->n = sq->n + size;

  if (nres != NULL) *nres = size;

  return eslOK;
}

/* Function:  correct_ambiguity()
 * Synopsis:  Read in the dna sequence
 * Incept:    MSF, Thu Feb 4, 2010 [Janelia]
 *
 * Purpose:   Correct any ambiguity characters.
 */
static int
correct_ambiguity(ESL_SQFILE *sqfp, ESL_SQ *sq, int len)
{
  int     alen;         /* ambiguity length      */
  int     soff;         /* starting offset       */
  int     eoff;         /* ending offset         */
  int     ainx;         /* ambiguity index       */
  int     size;         /* size of table read in */
  int     cnt;          /* repeat count          */
  int     off;
  int     n;

  char   *ptr;

  unsigned char c;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->seq_alen == 0) return eslOK;

  if (fseek(ncbi->fppsq, ncbi->amb_off, SEEK_SET) != 0) return eslESYS;

  ptr = (sq->dsq != NULL) ? (char *)sq->dsq + 1 : sq->seq;
  ptr += sq->n;

  /* calculate the starting and ending offsets */
  soff = sq->start + sq->n - 1;
  eoff = soff + len;

  off = 0;
  ainx = 0;
  size = 0;
  alen = ncbi->seq_alen - 4;
  while (off < eoff) {
    /* check if we need to read in more of the  abmiguity table */
    if (ainx == size) {
      size = alen;
      size = (size > INIT_HDR_BUFFER_SIZE) ? INIT_HDR_BUFFER_SIZE : size;
      if (fread(ncbi->hdr_buf, sizeof(char), size, ncbi->fppsq) != size) return eslEFORMAT;
      alen -= size;
      ainx = 0;
    }

    /* get the ambiguity character */
    n = ((ncbi->hdr_buf[ainx] >> 4) & 0x0f);
    c = sqfp->inmap[n];
    if (sq->dsq == NULL) c = ncbi->alphasym[(int) c];

    /* get the repeat count */
    cnt = (ncbi->hdr_buf[ainx] & 0x0f) + 1;

    /* get the offset */
    off = ncbi->hdr_buf[ainx+1];
    off = (off << 8) | ncbi->hdr_buf[ainx+2];
    off = (off << 8) | ncbi->hdr_buf[ainx+3];

    if (off + cnt >= soff && off < eoff) {
      int inx;
      int start = (off > soff) ? off - soff : 0;
      int end   = (off + cnt > eoff) ? eoff : off - soff + cnt;
      for (inx = start; inx < end; ++inx) ptr[inx] = c;
    }

    off += cnt;
    ainx += 4;
  }

  return eslOK;
}

/* Function:  read_nres_dna()
 * Synopsis:  Read in the dna sequence
 * Incept:    MSF, Wed Jan 27, 2010 [Janelia]
 *
 * Purpose:   Read and translate the dna sequence.
 */
static int
read_nres_dna(ESL_SQFILE *sqfp, ESL_SQ *sq, int len, uint64_t *nres)
{
  int     inx;
  int     off;
  int     cnt;
  int     ncnt;
  int     start;
  int     skip;
  int     text;
  int     status;

  int     remainder;
  int     length;
  int     ssize;
  int     n;

  char   *ptr;
  void   *t;

  unsigned char c;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if (ncbi->index >= ncbi->num_seq) return eslEOF;

  /* if we don't know the sequence length, figure it out */
  if (ncbi->seq_L == -1) {
    if (fseek(ncbi->fppsq, ncbi->seq_apos - 1, SEEK_SET) != 0) return eslESYS;
    if (fread(&c, sizeof(char), 1, ncbi->fppsq) != 1)          return eslEFORMAT;

    ssize       = ncbi->seq_apos - sq->doff - 1;
    remainder   = c & 0x03;
    length      = ssize * 4 + remainder;

    ncbi->seq_L = length;
  }

  /* check if we are at the end */
  if (sq->start + sq->n > ncbi->seq_L) {
    if (nres != NULL) *nres = 0;
    sq->L = ncbi->seq_L;
    return eslEOD;
  }

  /* calculate where to start reading from */
  start = sq->start + sq->n - 1;
  off   = sq->doff + start / 4;

  /* calculate bits to skip at the beginning and end */
  cnt   = ncbi->seq_L - (sq->start + sq->n - 1);
  cnt   = (cnt > len) ? len : cnt;

  skip      = start & 0x03;
  remainder = skip + cnt;
  remainder = (remainder & 0x03) ? (remainder & 0x03) : 4;

  /* calculate bytes need to read in the window */
  ssize = (cnt + skip + 3) / 4;

  /* seek to the windows location and read into the buffer */
  if (fseek(ncbi->fppsq, off, SEEK_SET) != 0) return eslESYS;
  if (ncbi->hdr_alloced < ssize) {
    while (ncbi->hdr_alloced < ssize) ncbi->hdr_alloced += ncbi->hdr_alloced;
    ESL_RALLOC(ncbi->hdr_buf, t, sizeof(char) * ncbi->hdr_alloced);
  }
  if (fread(ncbi->hdr_buf, sizeof(char), ssize, ncbi->fppsq) != ssize) return eslEFORMAT;

  /* figure out if the sequence is in digital mode or not */
  if (sq->dsq != NULL) {
    text = FALSE;
    ptr = (char *)sq->dsq + 1;
  } else {
    text = TRUE;
    ptr = sq->seq;
  }
  ptr += sq->n;

  inx = 0;
  c = ncbi->hdr_buf[inx];
  for (inx = skip; inx < 4; ++inx) {
    n = 1 << ((c >> (6 - inx * 2)) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  for (inx = 1; inx < ssize - 1; ++inx) {
    c = ncbi->hdr_buf[inx];
    n = 1 << ((c >> 6) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 4) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 2) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
    n = 1 << ((c >> 0) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  /* handle the remainder */
  c = ncbi->hdr_buf[inx];
  for (inx = 0; inx < remainder; ++inx) {
    n = 1 << ((c >> (6 - inx * 2)) & 0x03);
    *ptr = sqfp->inmap[n];
    if (text) *ptr = ncbi->alphasym[(int) *ptr];
    ++ptr;
  }

  *ptr = (text) ? '\0' : eslDSQ_SENTINEL;

  /* calculate the number of residues processed */
  ncnt = (ssize - 1) * 4 + remainder - skip;

  /* start processing the abmiguity table if there is one */
  if (ncbi->seq_alen > 0) {
    correct_ambiguity(sqfp, sq, ncnt);
  }

  sq->n = sq->n + ncnt;

  if (nres != NULL) *nres = ncnt;

  return eslOK;

 ERROR:
  return eslEMEM;
}


/* Function:  inmap_ncbi()
 * Synopsis:  Set up a translation map
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Initialize the translation map used to translate a ncbi
 *            sequences to the internal representation used in hmmer.
 *
 * Returns:   <eslOK> on success;
 * 
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
inmap_ncbi(ESL_SQFILE *sqfp)
{
  int status = eslOK;

  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  switch(ncbi->alphatype) {
  case eslDNA:
    status = inmap_ncbi_dna(sqfp);
    break;
  case eslAMINO:
    status = inmap_ncbi_amino(sqfp);
    break;
  default:
    ESL_EXCEPTION(eslEINVAL, "bad alphabet type: unrecognized");
  }

  return status;
}

/* Function:  inmap_ncbi_dna()
 * Synopsis:  Set up a translation map
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Initialize the translation map used to translate a ncbi
 *            protein sequence to the internal representation used in
 *            hmmer.
 *
 * Returns:   <eslOK> on success;
 * 
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
inmap_ncbi_dna(ESL_SQFILE *sqfp)
{
  int x, y;
  const char *ncbisym = "-ACMGRSVTWYHKDBN";

  ESL_ALPHABET    *abc  = NULL;
  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if ((abc = esl_alphabet_Create(eslDNA)) == NULL) return eslEMEM;

  for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;

  /* for each letter in the ncbi alphabet, find that letter in the
   * hmmer alphabet and map the translation.
   */
  for (x = 0; x < strlen(ncbisym); ++x) {
    for (y = 0; y < strlen(abc->sym); ++y) {
      if (ncbisym[x] == abc->sym[y]) {
	sqfp->inmap[x] = y;
	break;
      }
    }

    /* there is a problem if a translation does not exist */
    if (y >= strlen(abc->sym)) return eslEFORMAT;
  }

  if (ncbi->alphasym == NULL) esl_strdup(abc->sym, -1, &ncbi->alphasym);

  esl_alphabet_Destroy(abc);

  return eslOK;
}


/* Function:  inmap_ncbi_amino()
 * Synopsis:  Set up a translation map
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Initialize the translation map used to translate a ncbi
 *            protein sequence to the internal representation used in
 *            hmmer.
 *
 * Returns:   <eslOK> on success;
 * 
 * Throws:    <eslEMEM> on allocation failure;
 */
static int
inmap_ncbi_amino(ESL_SQFILE *sqfp)
{
  int x, y;
  const char *ncbisym = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

  ESL_ALPHABET    *abc  = NULL;
  ESL_SQNCBI_DATA *ncbi = &sqfp->data.ncbi;

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) return eslEMEM;

  for (x =  0;  x < 128;  x++) sqfp->inmap[x] = eslDSQ_ILLEGAL;

  /* for each letter in the ncbi alphabet, find that letter in the
   * hmmer alphabet and map the translation.
   */
  for (x = 0; x < strlen(ncbisym); ++x) {
    for (y = 0; y < strlen(abc->sym); ++y) {
      if (ncbisym[x] == abc->sym[y]) {
	sqfp->inmap[x] = y;
	break;
      }
    }

    /* there is a problem if a translation does not exist */
    if (y >= strlen(abc->sym)) return eslEFORMAT;
  }

  if (ncbi->alphasym == NULL) esl_strdup(abc->sym, -1, &ncbi->alphasym);

  esl_alphabet_Destroy(abc);

  return eslOK;
}



/*****************************************************************
 *# 5. Parsing routines
 *****************************************************************/ 

/* Function:  parse_expect()
 * Synopsis:  Expect the next bytes to parse match
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Match if the next <len> bytes to parse match the bytes
 *            in <str>.  If the bytes do not match, throw <eslEFORMAT>
 *            error.  Advance the parsers pointer.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header or if the data to parse does not match
 *            what is expected.
 */
static int
parse_expect(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %d : 0x%X(%d)\n",
	       len, (uint32_t) (ncbi->hdr_ptr - ncbi->hdr_buf), ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  /* check the buffer matches the token string */
  c = (unsigned char *) str;
  while (len--) {
    if (*ncbi->hdr_ptr != *c) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting 0x%X found 0x%X at %d : 0x%X(%d)\n",
	       *ncbi->hdr_ptr, *c, (uint32_t) (ncbi->hdr_ptr - ncbi->hdr_buf), ncbi->hdr_fpos, ncbi->hdr_size); 
    }
    ncbi->hdr_ptr++;
    c++;
  }

  return eslOK;
}

/* Function:  parse_accept()
 * Synopsis:  Check if the next bytes to parse match
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Check if the next <len> bytes to parse match the bytes
 *            in <str>.  If the bytes match, they are consumed and the
 *            parsers pointer is advanced.
 *
 * Returns:   <eslOK> on success
 *            <eslEFORMAT> if the bytes to not match.
 */
static int
parse_accept(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  int i;
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* check the buffer matches the token string */
  if (ncbi->hdr_ptr + len > limit)  return eslEFORMAT;

  /* verify the buffer matches the token string without advancing
   * the buffer pointers until we have a complete match.
   */
  c = (unsigned char *) str;
  for (i = 0; i < len; ++i) {
    if (ncbi->hdr_ptr[i] != c[i])   return eslEFORMAT;
  }
  ncbi->hdr_ptr += len;

  return eslOK;
}

/* Function:  parse_peek()
 * Synopsis:  Peek at the next byte to parse
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Return the next characer to be parsed without advancing the
 *            parsers pointer.
 *
 * Returns:   <eslOK> on success
 *            <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_peek(ESL_SQNCBI_DATA *ncbi, unsigned char *c)
{
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + 1 > limit)    return eslEFORMAT;

  *c = *ncbi->hdr_ptr;

  return eslOK;
}

/* Function:  parse_consume()
 * Synopsis:  Copies bytes from the header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Copies <len> bytes from the header to the buffer supplied by
 *            <str> if non-null.  Adcance the parser pointer.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_consume(ESL_SQNCBI_DATA *ncbi, void *str, int len)
{
  int i;
  unsigned char *c;
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %d : 0x%X(%d)\n",
	       len, (uint32_t) (ncbi->hdr_ptr - ncbi->hdr_buf), ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  /* copy the characters in the buffer to <str> */
  c = (unsigned char *) str;
  for (i = 0; i < len; ++i) {
    if (c != NULL) *c++ = *ncbi->hdr_ptr++;
  }

  return eslOK;
}

/* Function:  parse_advance()
 * Synopsis:  Advance the parser pointer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Advance the parser pointer <len> bytes.
 *
 * Returns:   <eslOK> on success
 * 
 * Throws:    <eslEFORMAT> if there are insufficient bytes remaining
 *            in the header.
 */
static int
parse_advance(ESL_SQNCBI_DATA *ncbi, int len)
{
  unsigned char *limit;

  limit = ncbi->hdr_buf + ncbi->hdr_size;

  /* verify the buffer has atleast len bytes remaining */
  if (ncbi->hdr_ptr + len > limit) {
      ESL_FAIL(eslEFORMAT, ncbi->errbuf, "Expecting %d bytes at %d : 0x%X(%d)\n",
	       len, (uint32_t) (ncbi->hdr_ptr - ncbi->hdr_buf), ncbi->hdr_fpos, ncbi->hdr_size); 
  }

  ncbi->hdr_ptr += len;

  return eslOK;
}

/* Function:  parse_header()
 * Synopsis:  Parse the ncbi db header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Parse a ncbi database header.  This routine implements
 *            a recursive descent parser for the ASN.1 definition of
 *            a blast database header filling in <sq>.
 *
 *            The blast db header can have multiple definitions defined
 *            within it.  Only the information from the first usable
 *            defition will be used.
 *
 * Returns:   <eslOK> on success; the new sequence is stored in <s>.
 * 
 *            Returns <eslEFORMAT> if there's a problem with the format,
 *            such as an illegal character.
 */
static int
parse_header(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;
  void *tmp;

  unsigned char c;

  /* read in the header data */
  if (ncbi->hdr_alloced < ncbi->hdr_size) {
    while (ncbi->hdr_alloced < ncbi->hdr_size) ncbi->hdr_alloced += ncbi->hdr_alloced;
    ESL_RALLOC(ncbi->hdr_buf, tmp, sizeof(char) * ncbi->hdr_alloced);
  }
  if (fread(ncbi->hdr_buf, sizeof(char), ncbi->hdr_size, ncbi->fpphr) != ncbi->hdr_size) return eslEFORMAT;
  ncbi->hdr_ptr = ncbi->hdr_buf;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  if (parse_peek(ncbi, &c) != eslOK)                          return eslEFORMAT;

  /* parse the different seq id structures */
  while (c != 0x00) {
    if ((status = parse_def_line(ncbi, sq)) != eslOK)         return status;

    if (parse_peek(ncbi, &c) != eslOK)                        return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
 ERROR:
  return eslEMEM;
}

/* Function:  parse_def_line()
 * Synopsis:  Parse the Blast-def-line definition
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Blast-def-line ::= SEQUENCE {
 * 	title       VisibleString       OPTIONAL,  -- simple title
 * 	seqid       SEQUENCE OF Seq-id,            -- Regular NCBI Seq-Id
 * 	taxid       INTEGER             OPTIONAL,  -- taxonomy id
 * 	memberships SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
 * 	links       SEQUENCE OF INTEGER OPTIONAL,  -- bit arrays
 * 	other-info  SEQUENCE OF INTEGER OPTIONAL   -- for future use (probably genomic sequences)
 * }
 */
static int
parse_def_line(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  char *buf;
  int   taxid;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional title */
  sq->desc[0] = 0;
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    free(sq->desc);
    sq->dalloc = strlen(buf) + 1;
    sq->desc   = buf;
  }

  /* look for sequence id structure */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_seq_id(ncbi, sq)) != eslOK)             return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional taxonomy id */
  sq->tax_id = -1;
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, &taxid)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    sq->tax_id = taxid;
  }

  /* look for an optional memberships */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional links */
  if (parse_accept(ncbi, "\xa4\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional other info */
  if (parse_accept(ncbi, "\xa5\x80", 2) == eslOK) {
    if ((status = ignore_sequence_of_integer(ncbi)) != eslOK) return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_seq_id()
 * Synopsis:  Parse the Blast-def-line definition
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Seq-id ::= CHOICE {
 *     local             Object-id ,       -- local use
 *     gibbsq            INTEGER ,         -- Geninfo backbone seqid
 *     gibbmt            INTEGER ,         -- Geninfo backbone moltype
 *     giim              Giimport-id ,     -- Geninfo import id
 *     genbank           Textseq-id ,
 *     embl              Textseq-id ,
 *     pir               Textseq-id ,
 *     swissprot         Textseq-id ,
 *     patent            Patent-seq-id ,
 *     other             Textseq-id ,      -- for historical reasons, 'other' = 'refseq'
 *     general           Dbtag ,           -- for other databases
 *     gi                INTEGER ,         -- GenInfo Integrated Database
 *     ddbj              Textseq-id ,      -- DDBJ
 *     prf               Textseq-id ,      -- PRF SEQDB
 *     pdb               PDB-seq-id ,      -- PDB sequence
 *     tpg               Textseq-id ,      -- Third Party Annot/Seq Genbank
 *     tpe               Textseq-id ,      -- Third Party Annot/Seq EMBL
 *     tpd               Textseq-id ,      -- Third Party Annot/Seq DDBJ
 *     gpipe             Textseq-id ,      -- Internal NCBI genome pipeline processing ID
 *     named-annot-track Textseq-id        -- Internal named annotation tracking ID
 * }
 */
static int
parse_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  unsigned char c;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)           return eslEFORMAT;

  if (parse_consume(ncbi, &c, 1) != eslOK)                  return eslEFORMAT;

  /* parse the different seq id structures */
  while (c != 0x00) {
    if (parse_expect(ncbi, "\x80", 1) != eslOK)             return eslEFORMAT;
    switch (c) {
    case 0xa0: /* LOCAL */
      status = parse_object_id(ncbi, sq);
      break;
    case 0xa1: /* GIBBSQ */
    case 0xa2: /* GIBBMT */
      status = parse_integer(ncbi, NULL);
      break;
    case 0xa3: /* GIIM */
      return eslEFORMAT;
      break;
    case 0xa4: /* GENBANK */
    case 0xa5: /* EMBL */
    case 0xa6: /* PIR */
    case 0xa7: /* SWISSPROT */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xa8: /* PATENT */
      status = parse_patent_seq_id(ncbi, sq);
      break;
    case 0xa9: /* OTHER */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xaa: /* GENERAL */
      status = parse_dbtag(ncbi, sq);
      break;
    case 0xab: /* GI */
      status = parse_integer(ncbi, NULL);
      break;
    case 0xac: /* DDBJ */
    case 0xad: /* PRF */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    case 0xae: /* PDB */
      status = parse_pdb_seq_id(ncbi, sq);
      break;
    case 0xaf: /* TPG */
    case 0xb0: /* TPE */
    case 0xb1: /* TPD */
    case 0xb2: /* GPIPE */
    case 0xb3: /* NAMED ANNOT TRACK */
      status = parse_textseq_id(ncbi, sq);
      sq = NULL;
      break;
    default:
      status = eslEFORMAT;
    }

    if (status != eslOK)                                    return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)         return eslEFORMAT;
    if (parse_consume(ncbi, &c, 1)        != eslOK)         return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (c != 0x00 || parse_expect(ncbi, "\x00", 1) != eslOK)  return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_textseq_id()
 * Synopsis:  Parse the general text header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Textseq-id ::= SEQUENCE {
 *     name      VisibleString OPTIONAL ,
 *     accession VisibleString OPTIONAL ,
 *     release   VisibleString OPTIONAL ,
 *     version   INTEGER       OPTIONAL
 * }
 */
static int
parse_textseq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  char *buf = NULL;
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional name */
  if (sq != NULL) sq->name[0] = 0;
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    if (sq != NULL) {
      free(sq->name);
      sq->nalloc = strlen(buf) + 1;
      sq->name   = buf;
    } else {
      free(buf);
    }
    buf = NULL;
  }

  /* look for an optional accession */
  if (sq != NULL) sq->acc[0] = 0;
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, -1, &buf)) != eslOK)     return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;

    if (sq != NULL) {
      free(sq->acc);
      sq->aalloc = strlen(buf) + 1;
      sq->acc    = buf;
    } else {
      free(buf);
    }
    buf = NULL;
  }

  /* look for an optional release */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional version */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_dbtag()
 * Synopsis:  Parse the a general db tag
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Dbtag ::= SEQUENCE {
 *     db  VisibleString ,     -- name of database or system
 *     tag Object-id           -- appropriate tag
 * }
 */
static int
parse_dbtag(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an db name */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;

  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for a tag object */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_object_id(ncbi, sq)) != eslOK)          return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_patent_seq_id()
 * Synopsis:  Parse the patent header
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Patent-seq-id ::= SEQUENCE {
 *     seqid INTEGER ,          -- number of sequence in patent
 *     cit   Id-pat             -- patent citation
 * }
 */
static int
parse_patent_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a seqid */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_integer(ncbi, NULL)) != eslOK)          return status;

  /* look for a patent citation object */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_id_pat(ncbi, sq)) != eslOK)             return status;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_id_pat()
 * Synopsis:  Parse the patent citation
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Id-pat ::= SEQUENCE {                         -- just to identify a patent
 *     country  VisibleString ,                  -- Patent Document Country
 *     id       CHOICE {
 *         number     VisibleString ,            -- Patent Document Number
 *         app-number VisibleString              -- Patent Doc Appl Number
 *     } ,
 *     doc-type VisibleString         OPTIONAL   -- Patent Doc Type
 * }
 */
static int
parse_id_pat(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a country */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;

  /* look for an id */
  if (parse_expect(ncbi, "\xa1\x80", 2) != eslOK)             return eslEFORMAT;

  /* the id is a choice of two strings */

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional taxonomy id */
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    status = parse_string(ncbi, 0, NULL);
  } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    status = parse_string(ncbi, 0, NULL);
  } else {
    status = eslEFORMAT;
  }
  if (status != eslOK)                                        return status;

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for a doc type */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_object_id()
 * Synopsis:  Parse a generic sequence id
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Object-id ::= CHOICE {
 *     id  INTEGER ,
 *     str VisibleString
 * }
 */
static int
parse_object_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  char  *buf = NULL;
  int    status;

  /* look for an optional taxonomy id */
  if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
    status = parse_integer(ncbi, NULL);
  } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    status = parse_string(ncbi, -1, &buf);

    if (sq != NULL) {
      free(sq->name);
      sq->nalloc = strlen(buf) + 1;
      sq->name   = buf;
    } else {
      free(buf);
    }
    buf = NULL;
  } else {
    status = eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (status == eslOK) {
    status = parse_expect(ncbi, "\x00\x00", 2);
  }

  return status;
}


/* Function:  parse_pdb_seq_id()
 * Synopsis:  Parse a PDB sequence
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * PDB-seq-id ::= SEQUENCE {
 *     mol   PDB-mol-id ,              -- the molecule name
 *     chain INTEGER ,                 -- a single ASCII character, chain id
 *     rel   Date         OPTIONAL }   -- release date, month and year
 *
 * Date ::= CHOICE {
 *     str   VisibleString ,           -- for those unparsed dates
 *     std   Date-std                  -- use this if you can
 * }
 */
static int
parse_pdb_seq_id(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for an pdb mol id */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_string(ncbi, 0, NULL)) != eslOK)        return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for chain */
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional date */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if (parse_accept(ncbi, "\xa0\x80", 2) == eslOK) {
      status = parse_string(ncbi, 0, NULL);
    } else if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
      status = parse_date_std(ncbi, sq);
    } else {
      status = eslEFORMAT;
    }
    if (status != eslOK)                                      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_date_std()
 * Synopsis:  Parse the data structure
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Date-std ::= SEQUENCE {              -- NOTE: this is NOT a unix tm struct
 *     year   INTEGER ,                 -- full year (including 1900)
 *     month  INTEGER       OPTIONAL ,  -- month (1-12)
 *     day    INTEGER       OPTIONAL ,  -- day of month (1-31)
 *     season VisibleString OPTIONAL ,  -- for "spring", "may-june", etc
 *     hour   INTEGER       OPTIONAL ,  -- hour of day (0-23)
 *     minute INTEGER       OPTIONAL ,  -- minute of hour (0-59)
 *     second INTEGER       OPTIONAL    -- second of minute (0-59)
 * }
 */
static int
parse_date_std(ESL_SQNCBI_DATA *ncbi, ESL_SQ *sq)
{
  int   status;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)             return eslEFORMAT;

  /* look for a year */
  if (parse_expect(ncbi, "\xa0\x80", 2) != eslOK)             return eslEFORMAT;
  if ((status = parse_integer(ncbi, NULL)) != eslOK)          return status;
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  /* look for an optional month */
  if (parse_accept(ncbi, "\xa1\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional day */
  if (parse_accept(ncbi, "\xa2\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional season */
  if (parse_accept(ncbi, "\xa3\x80", 2) == eslOK) {
    if ((status = parse_string(ncbi, 0, NULL)) != eslOK)      return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional hour */
  if (parse_accept(ncbi, "\xa4\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional minute */
  if (parse_accept(ncbi, "\xa5\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* look for an optional second */
  if (parse_accept(ncbi, "\xa6\x80", 2) == eslOK) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK)        return status;
    if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)           return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)             return eslEFORMAT;

  return eslOK;
}


/* Function:  parse_string()
 * Synopsis:  Parse a visible string
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads an string from the header stream.  The arguement <max>
 *            specified the maximum number of characters to save.  If <max>
 *            is -1, the entire string will be saved.
 *
 *            The string will always be zero terminated.
 *
 *            If <str> is non null, the parsed integer will be placed
 *            in the pointer.  The calling routine is responsible for
 *            freeing the allocated memory.
 *
 * Returns:   <eslOK> on success.
 *            <eslEMEM> if there's a memory allocation error.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
parse_string(ESL_SQNCBI_DATA *ncbi, int max, char **str)
{
  int n;
  int len;
  int status;

  char *v  = NULL;

  unsigned char  x;
  unsigned char  c;
  unsigned char *ptr;

  if (parse_expect(ncbi, "\x1a", 1) != eslOK)  return eslEFORMAT;

  /* the next byte is the length of the string.  if the length is
   * less than 128, then this is the true length; otherwise this
   * length describes the number of bytes necessary to hold the
   * true length of the string in the lower 7 bits.
   */
  if (parse_consume(ncbi, &c, 1) != eslOK)     return eslEFORMAT;
  if (c < 128) {
    n = c;
  } else {
    c = c & 0x7f;
    if (c > sizeof(n))                                 return eslEFORMAT;

    n = 0;
    while (c > 0) {
      if (parse_consume(ncbi, &x, 1) != eslOK) return eslEFORMAT;
      n = (n << 8) + (unsigned int) x;
      --c;
    }
  }

  /* validate the length of the string */
  ptr = ncbi->hdr_ptr;
  if (parse_advance(ncbi, n) != eslOK)         return eslEFORMAT;

  /* now that we have the length of the string, check how much
   * of it (if any) we need to save.
   */
  if (str != NULL && max != 0) {
    if (max == -1 || max > n)  len = n;
    else                       len = max - 1;

    ESL_ALLOC(v, sizeof(char) * (len + 1));
    memcpy(v, ptr, len);
    v[len] = 0;

    *str = v;
  }

  return eslOK;

 ERROR:
  if (v != NULL) free(v);
  return status;
}


/* Function:  parse_integer()
 * Synopsis:  Parse an integer
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Reads an integer from the header stream.  If the integer is
 *            more bytes than the native int format, the most significant
 *            bytes will be lost.
 *
 *            If <value> is non null, the parsed integer will be placed
 *            in the pointer.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
parse_integer(ESL_SQNCBI_DATA *ncbi, int *value)
{
  int n;

  unsigned char  c;
  unsigned char *ptr;

  if (parse_expect(ncbi, "\x02", 1) != eslOK) return eslEFORMAT;

  /* get the length of the integer */
  if (parse_peek(ncbi, &c) != eslOK)          return eslEFORMAT;
  ptr = ncbi->hdr_ptr + 1;

  /* advance past the integer to make sure the buffer holds all
   * of the integer.  the pointer <ptr> points the the start of
   * the integer.
   */
  if (parse_advance(ncbi, c + 1) != eslOK)    return eslEFORMAT;

  n = 0;
  while (c > 0) {
    n = (n << 8) + (unsigned int) *ptr++;
    --c;
  }

  if (value != NULL) *value = n;

  return eslOK;
}


/* Function:  ignore_sequence_of_integer()
 * Synopsis:  Skip over the sequence of integers
 * Incept:    MSF, Mon Dec 10, 2009 [Janelia]
 *
 * Purpose:   Skip over a sequence of integers.
 *
 * Returns:   <eslOK> on success.
 *            <eslEFORMAT> if there's a problem with the format.
 *
 */
static int
ignore_sequence_of_integer(ESL_SQNCBI_DATA *ncbi)
{
  int status;
  unsigned char c;

  /* verify we are at the beginning of a structure */
  if (parse_expect(ncbi, "\x30\x80", 2) != eslOK)      return eslEFORMAT;

  if (parse_peek(ncbi, &c) != eslOK)                   return eslEFORMAT;
  while (c == 0x02) {
    if ((status = parse_integer(ncbi, NULL)) != eslOK) return status;
    if (parse_peek(ncbi, &c) != eslOK)                 return eslEFORMAT;
  }

  /* verify we are at the end of the structure */
  if (parse_expect(ncbi, "\x00\x00", 2) != eslOK)      return eslEFORMAT;

  return eslOK;
}


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
