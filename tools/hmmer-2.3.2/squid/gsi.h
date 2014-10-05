/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

#ifndef GSIH_INCLUDED
#define GSIH_INCLUDED

/* gsi.h
 * Database indexing (GSI format support)
 * RCS $Id: gsi.h,v 1.3 2001/08/04 20:15:42 eddy Exp $
 *
 * A GSI (generic sequence index) file is composed of
 * recnum + nfiles + 1 records. Each record contains
 * three fields; key, file number, and disk offset.
 * Record 0 contains:
 *   [ "GSI" ]  [ nfiles ]  [ recnum ]
 * Records 1..nfiles map file names to file numbers, and contain:
 *   [ filename ] [ file number, 1..nfiles ] [ 0 (unused) ]
 * Records nfiles+1 to recnum+nfiles+1 provide disk offset
 * and file number indices for every key:
 *   [ key ] [ file number ] [ offset]
 *
 * Because the file is binary, we take some (but not 
 * complete) care to improve portability amongst platforms.
 * This means using network order integers (see ntohl())
 * and defining types for 16 and 32 bit integers.
 * 
 * Because we use 32-bit offsets, ftell(), and fseek(),
 * there is an implicit 2 Gb file size maximum.
 * AFAIK neither ANSI C nor POSIX provide a portable solution
 * to this problem. fsetpos(), fgetpos() use an
 * opaque fpos_t datatype that we can't write portably
 * to a disk file. Suggestions welcomed.
 */
#define GSI_KEYSIZE    32         /* keys are 32 bytes long */  
#define GSI_RECSIZE    38	  /* 32 + 2 + 4 bytes       */
#define SQD_UINT16_MAX 65535	  /* 2^16-1 */
#define SQD_UINT32_MAX 4294967295U/* 2^32-1 */

struct gsi_s {
  FILE        *gsifp;		/* open GSI index file            */
  sqd_uint16   nfiles;		/* number of files = 16 bit int   */
  sqd_uint32   recnum;		/* number of records = 32 bit int */
};
typedef struct gsi_s GSIFILE;

struct gsikey_s {
  char       key[GSI_KEYSIZE];
  sqd_uint16 filenum;
  sqd_uint32 offset;
};
struct gsiindex_s {
  char           **filenames;
  int             *fmt;
  sqd_uint16       nfiles;

  struct gsikey_s *elems;
  int              nkeys;
};  


/* from gsi.c
 */
extern GSIFILE *GSIOpen(char *gsifile);
extern int      GSIGetRecord(GSIFILE *gsi, char *f1, sqd_uint16 *f2, sqd_uint32 *f3);
extern int      GSIGetOffset(GSIFILE *gsi, char *key, char *sqfile, 
			  int *fmt, long *ret_offset);
extern void     GSIClose(GSIFILE *gsi);
extern struct gsiindex_s *GSIAllocIndex(void);
extern void     GSIFreeIndex(struct gsiindex_s *g);
extern void     GSIAddFileToIndex(struct gsiindex_s *g, char *filename, int fmt);
extern void     GSIAddKeyToIndex(struct gsiindex_s *g, char *key, int filenum, long offset);
extern void     GSISortIndex(struct gsiindex_s *g);
extern void     GSIWriteIndex(FILE *fp, struct gsiindex_s *g);
extern void     GSIWriteHeader(FILE *fp, int nfiles, long nkeys);
extern int      GSIWriteFileRecord(FILE *fp, char *fname, int idx, int fmt);
extern int      GSIWriteKeyRecord(FILE *fp, char *key, int fileidx, long offset);

#endif /*GSIH_INCLUDED*/
