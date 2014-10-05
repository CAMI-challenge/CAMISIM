/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

#ifndef GSI64H_INCLUDED
#define GSI64H_INCLUDED
#ifdef USE_GSI64

/* gsi64.h
 * Database indexing (GSI64 format support)
 * CVS $Id: gsi64.h,v 1.2 2000/12/21 23:42:59 eddy Exp $
 *
 * A GSI64 (generic sequence index, 64 bit hack) file is composed of
 * recnum + nfiles + 1 records. Each record contains
 * three fields; key, file number, and disk offset.
 * Record 0 contains:
 *   [ "GSI64" ]  [ nfiles ]  [ recnum ]
 * Records 1..nfiles map file names to file numbers, and contain:
 *   [ filename ] [ file number, 1..nfiles ] [ 0 (unused) ]
 * Records nfiles+1 to recnum+nfiles+1 provide disk offset
 * and file number indices for every key:
 *   [ key ] [ file number ] [ offset]
 *
 * Because the file is binary, we take some (but not 
 * complete) care to improve portability amongst platforms.
 * This means using network order integers (see ntohl())
 * and defining types for 16 and 64 bit integers.
 * 
 * A short test program that verifies the sizes of these
 * data types would be a good idea...
 * 
 * Because we use 64-bit offsets, ftell64(), and fseek64(),
 * we rely on the OS actually providing these. This is
 * a temporary hack for human genome analysis.
 */
typedef unsigned long long  sqd_uint64; /* 64 bit integer. */

#define GSI64_KEYSIZE    32         /* keys are 32 bytes long */  
#define GSI64_RECSIZE    42	    /* 32 + 2 + 8 bytes       */
#define SQD_UINT16_MAX 65535	    /* 2^16-1 */
#define SQD_UINT64_MAX 18446744073709551615LU /* 2^64-1 */

struct gsi64_s {
  FILE        *gsifp;		/* open GSI index file            */
  sqd_uint16   nfiles;		/* number of files = 16 bit int   */
  sqd_uint64   recnum;		/* number of records = 64 bit int */
};
typedef struct gsi64_s GSI64FILE;

struct gsi64key_s {
  char       key[GSI64_KEYSIZE];
  sqd_uint16 filenum;
  sqd_uint64 offset;
};
struct gsi64index_s {
  char           **filenames;
  int             *fmt;
  sqd_uint16       nfiles;

  struct gsi64key_s *elems;
  sqd_uint64         nkeys;
};  



/* if ntohl() and friends are not available, you
 * can slip replacements in by providing sre_ntohl()
 * functions. (i.e., there is a possible portability problem here.)
 */
#if 0
#define sre_ntohl(x)  ntohl(x); 
#define sre_ntohs(x)  ntohs(x);
#define sre_htonl(x)  htonl(x); 
#define sre_htons(x)  htons(x);
#endif

/* from gsi64.c
 */
extern GSI64FILE *GSI64Open(char *gsifile);
extern int        GSI64GetRecord(GSI64FILE *gsi, char *f1, sqd_uint16 *f2, sqd_uint64 *f3);
extern int        GSI64GetOffset(GSI64FILE *gsi, char *key, char *sqfile, 
			         int *fmt, long long *ret_offset);
extern void     GSI64Close(GSI64FILE *gsi);
extern struct gsi64index_s *GSI64AllocIndex(void);
extern void     GSI64FreeIndex(struct gsi64index_s *g);
extern void     GSI64AddFileToIndex(struct gsi64index_s *g, char *filename, int fmt);
extern void     GSI64AddKeyToIndex(struct gsi64index_s *g, char *key, int filenum, long long offset);
extern void     GSI64SortIndex(struct gsi64index_s *g);
extern void     GSI64WriteIndex(FILE *fp, struct gsi64index_s *g);
extern void     GSI64WriteHeader(FILE *fp, int nfiles, long long nkeys);
extern int      GSI64WriteFileRecord(FILE *fp, char *fname, int idx, int fmt);
extern int      GSI64WriteKeyRecord(FILE *fp, char *key, int fileidx, long long offset);

#endif /* USE_GSI64 */
#endif /*GSIH_INCLUDED*/
