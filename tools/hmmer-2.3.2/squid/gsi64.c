/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/ 


/* gsi64.c
 * Updated interfaces for GSI64 64-bit "generic sequence index" files.
 * See gsi.c for old interfaces.
 * This is a temporary hack! Needed for human genome project.
 */

/*    1 + <nfiles> + <nkeys> total records.
 *    Each record = 42 bytes.
 *
 *  one header record     :  <"GSI64"  (32)> <nfiles (2)> <nkeys (8)> 
 *  <nfiles> file records :  <filename (32)> <fileno (2)> <fmt   (8)> 
 *  <nkeys>  key records  :  <key      (32)> <fileno (2)> <offset(8)> 
 * 
 * CVS $Id: gsi64.c,v 1.3 2003/04/14 16:00:16 eddy Exp $
 */
#include "squidconf.h"

#ifdef USE_GSI64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SEEK_SET
#include <unistd.h>	/* needed for poor crippled SunOS */
#endif

#include "squid.h"
#include "gsi64.h"

/*****************************************************************
 * GSI64 index file access routines
 *****************************************************************/

/* Function: GSI64Open()
 * 
 * Purpose:  Open a GSI64 file. Returns the number of records in
 *           the file and a file pointer. Returns NULL on failure.
 *           The file pointer should be fclose()'d normally.
 */
GSI64FILE *
GSI64Open(char *gsifile)
{
  GSI64FILE  *gsi;
  char        magic[GSI64_KEYSIZE];

  gsi = (GSI64FILE *) MallocOrDie (sizeof(GSI64FILE));
  if ((gsi->gsifp = fopen(gsifile, "r")) == NULL)
    { free(gsi); squid_errno = SQERR_NOFILE; return NULL; }

  if (! fread(magic, sizeof(char), GSI64_KEYSIZE, gsi->gsifp))
    { free(gsi); squid_errno = SQERR_NODATA; return NULL; }
  if (strcmp(magic, "GSI64") != 0) 
    { free(gsi); squid_errno = SQERR_FORMAT; return NULL; }

  if (! fread(&(gsi->nfiles), sizeof(sqd_uint16), 1, gsi->gsifp))
    { free(gsi); squid_errno = SQERR_NODATA; return NULL; }
  if (! fread(&(gsi->recnum), sizeof(sqd_uint64), 1, gsi->gsifp))
    { free(gsi); squid_errno = SQERR_NODATA; return NULL; }

#if 0			/* HACK! we don't byteswap */
  gsi->nfiles = sre_ntohs(gsi->nfiles); /* convert from network short */
  gsi->recnum = sre_ntohl(gsi->recnum); /* convert from network long  */
#endif

  return gsi;
}

/* Function: GSI64GetRecord()
 * 
 * Purpose:  Each non-header record of a GSI64 index file consists
 *           of 42 bytes: 32 bytes of character string, a 2 byte
 *           short, and an 8 byte long long. This function returns the
 *           three values.
 *           
 * Args:     gsi  - open GSI64 index file, correctly positioned at a record
 *           f1   - char[32], allocated by caller (or NULL if unwanted)
 *           f2   - pointer to short (or NULL if unwanted)
 *           f3   - pointer to long long  (or NULL if unwanted)
 *                  
 * Return:   0 on failure and sets squid_errno.                  
 */
int
GSI64GetRecord(GSI64FILE *gsi, char *f1, sqd_uint16 *f2, sqd_uint64 *f3)
{
  if (f1 == NULL) fseek64(gsi->gsifp, GSI64_KEYSIZE, SEEK_CUR);
  else if (! fread(f1, GSI64_KEYSIZE, 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }

  if (f2 == NULL) fseek64(gsi->gsifp, sizeof(sqd_uint16), SEEK_CUR);
  else if (! fread(f2, sizeof(sqd_uint16), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }

  if (f3 == NULL) fseek64(gsi->gsifp, sizeof(sqd_uint64), SEEK_CUR);
  else if (! fread(f3, sizeof(sqd_uint64), 1, gsi->gsifp))
    { squid_errno = SQERR_NODATA; return 0; }

#if 0 /* no byteswap yet! HACK! */
  if (f2 != NULL) *f2 = sre_ntohs(*f2);
  if (f3 != NULL) *f3 = sre_ntohl(*f3);
#endif 

  return 1;
}


/* Function: GSI64GetOffset()
 * 
 * Purpose:  From a key (sequence name), find a disk offset
 *           in an open general sequence index file by binary
 *           search. Presumably GSI64 indexing could be even faster
 *           if we used hashing.
 *   
 * Args:     gsi         - GSI64 index file, opened by GSI64Open()
 *           key         - name of key to retrieve indices for
 *           ret_seqfile - pre-alloced char[32] array for seqfile name
 *           ret_fmt     - format of seqfile
 *           ret_offset  - return: disk offset in seqfile.         
 */
int
GSI64GetOffset(GSI64FILE *gsi, char *key, char *ret_seqfile, 
	       int *ret_format, long long *ret_offset)
{
  sqd_uint64  left, right, mid;
  int         cmp;
  char        name[GSI64_KEYSIZE + 1];
  sqd_uint64  offset;
  sqd_uint16  filenum;
  sqd_uint64  fmt;

  name[GSI64_KEYSIZE] = '\0';

  left  = gsi->nfiles + 1;
  right = gsi->nfiles + gsi->recnum;
  mid   = (left + right) / 2;
  fseek64(gsi->gsifp, mid * GSI64_RECSIZE, SEEK_SET);

  while (GSI64GetRecord(gsi, name, &filenum, &offset))
    {
      cmp = strcmp(name, key);
      if      (cmp == 0)      break;	       /* found it!              */
      else if (left >= right) return 0;        /* oops, missed it; fail. */
      else if (cmp < 0)       left = mid + 1;  /* it's right of mid      */
      else if (cmp > 0)	      right = mid - 1; /* it's left of mid       */ 
      mid = (left + right) / 2;
      fseek64(gsi->gsifp, mid * GSI64_RECSIZE, SEEK_SET);
    }

  /* Using file number, look up the sequence file and format.
   */
  fseek64(gsi->gsifp, filenum * GSI64_RECSIZE, SEEK_SET);
  GSI64GetRecord(gsi, ret_seqfile, NULL, &fmt);
  *ret_format =  (int) fmt;
  *ret_offset = (long long) offset;
  
  return 1;
}
    
/* Function: GSI64Close()
 * 
 * Purpose:  Close an open GSI64 sequence index file.
 */
void
GSI64Close(GSI64FILE *gsi)
{
  fclose(gsi->gsifp);
  free(gsi);
}


/*****************************************************************
 * GSI64 index construction routines
 * SRE, Wed Nov 10 11:49:14 1999 [St. Louis]
 * 
 * API:
 *       g = GSI64AllocIndex();
 *       
 *       [foreach filename, <32 char, no directory path]
 *          GSI64AddFileToIndex(g, filename);
 *          filenum++;
 *          [foreach key, <32 char, w/ filenum 1..nfiles, w/ 64bit offset]
 *             GSI64AddKeyToIndex(g, key, filenum, offset);
 *            
 *       GSI64SortIndex(g);
 *       GSI64WriteIndex(fp, g);
 *       GSI64FreeIndex(g);
 *****************************************************************/
struct gsi64index_s *
GSI64AllocIndex(void)
{
  struct gsi64index_s *g;
  
  g = MallocOrDie(sizeof(struct gsi64index_s));
  g->filenames = MallocOrDie(sizeof(char *) * 10);
  g->fmt       = MallocOrDie(sizeof(int) * 10); 
  g->elems     = MallocOrDie(sizeof(struct gsi64key_s) * 100);
  g->nfiles    = 0;
  g->nkeys     = 0;
  return g;
}
void
GSI64FreeIndex(struct gsi64index_s *g)
{
  int i;
  for (i = 0; i < g->nfiles; i++) free(g->filenames[i]);
  free(g->filenames);
  free(g->fmt);
  free(g->elems);
  free(g);
}
void
GSI64AddFileToIndex(struct gsi64index_s *g, char *filename, int fmt)
{
  int len;

  len = strlen(filename);
  if (len >= GSI64_KEYSIZE) Die("File name too long to be indexed.");
  g->filenames[g->nfiles] = sre_strdup(filename, len);
  g->fmt[g->nfiles]       = fmt;
  g->nfiles++;
  if (g->nfiles % 10 == 0) {
    g->filenames = ReallocOrDie(g->filenames, sizeof(char *) * (g->nfiles + 10)); 
    g->fmt       = ReallocOrDie(g->fmt,       sizeof(int)    * (g->nfiles + 10)); 
  }
}
void
GSI64AddKeyToIndex(struct gsi64index_s *g, char *key, int filenum, long long offset)
{
  if (strlen(key) >= GSI64_KEYSIZE) Die("key too long in GSI64 index");
  if (filenum > SQD_UINT16_MAX) Die("too many files in GSI64 index");
  if (offset  > SQD_UINT64_MAX) Die("offset too big in GSI64 index");
  
  strncpy(g->elems[g->nkeys].key, key, GSI64_KEYSIZE-1);
  g->elems[g->nkeys].key[GSI64_KEYSIZE-1] = '\0';
  g->elems[g->nkeys].filenum = (sqd_uint16) filenum;
  g->elems[g->nkeys].offset  = (sqd_uint64) offset;
  g->nkeys++;

  if (g->nkeys % 100 == 0)
    g->elems = ReallocOrDie(g->elems, sizeof(struct gsi64key_s) * (g->nkeys + 100));
}
static int 
gsi_keysorter(const void *k1, const void *k2)
{
  struct gsi64key_s *key1;
  struct gsi64key_s *key2;
  key1 = (struct gsi64key_s *) k1;
  key2 = (struct gsi64key_s *) k2;
  return strcmp(key1->key, key2->key);
}
void
GSI64SortIndex(struct gsi64index_s *g)
{
  qsort((void *) g->elems, g->nkeys, sizeof(struct gsi64key_s), gsi_keysorter); 
}
void
GSI64WriteIndex(FILE *fp, struct gsi64index_s *g)
{
  sqd_uint16 i;
  sqd_uint64 j;

  /* Range checking.
   */
  if (g->nfiles > SQD_UINT16_MAX) Die("Too many files in GSI64 index.");
  if (g->nkeys  > SQD_UINT64_MAX) Die("Too many keys in GSI64 index.");

  GSI64WriteHeader(fp, g->nfiles, g->nkeys);
  for (i = 0; i < g->nfiles; i++)
    GSI64WriteFileRecord(fp, g->filenames[i], i+1, g->fmt[i]);
  for (j = 0; j < g->nkeys; j++)
    GSI64WriteKeyRecord(fp, g->elems[j].key, g->elems[j].filenum, g->elems[j].offset);
}





/* Function: GSI64WriteHeader()
 * Date:     SRE, Wed Aug  5 10:36:02 1998 [St. Louis]
 *
 * Purpose:  Write the first record to an open GSI64 file:
 *           "GSI64" <nfiles> <nkeys>
 *
 * Args:     fp      - open file to write to.
 *           nfiles  - number of files indexed
 *           nkeys   - number of keys indexed          
 *
 * Returns:  void
 */
void
GSI64WriteHeader(FILE *fp, int nfiles, long long nkeys)
{
  char       key[GSI64_KEYSIZE];
  sqd_uint16 f1;
  sqd_uint64 f2;

  /* beware potential range errors!
   */
  if (nfiles > SQD_UINT16_MAX) Die("GSI64: nfiles out of range");
  if (nkeys > SQD_UINT64_MAX)  Die("GSI64: nkeys out of range");

  f1 = (sqd_uint16) nfiles;
  f2 = (sqd_uint64) nkeys;
#if 0 /* HACK no byteswap */
  f1 = sre_htons(f1);
  f2 = sre_htonl(f2);
#endif
  strcpy(key, "GSI64");

  if (fwrite(key,   1, GSI64_KEYSIZE, fp) < GSI64_KEYSIZE) PANIC;
  if (fwrite(&f1,   2,  1, fp) < 1)  PANIC;
  if (fwrite(&f2,   8,  1, fp) < 1)  PANIC;
}


/* Function: GSI64WriteFileRecord()
 * Date:     SRE, Wed Aug  5 10:45:51 1998 [St. Louis]
 *
 * Purpose:  Write a file record to an open GSI64 file.
 *
 * Args:     fp    - open GSI64 file
 *           fname - file name (max 31 characters)
 *           idx   - file number
 *           fmt   - file format (e.g. kPearson, etc.)
 *
 * Returns:  0 on failure. 1 on success.
 */
int
GSI64WriteFileRecord(FILE *fp, char *fname, int idx, int fmt)
{
  sqd_uint16 f1;
  sqd_uint64 f2;

  if (strlen(fname) >= GSI64_KEYSIZE) return 0;
  if (idx > SQD_UINT16_MAX) Die("GSI64: file index out of range");
  if (fmt > SQD_UINT64_MAX) Die("GSI64: format index out of range");

  f1 = (sqd_uint16) idx;
  f2 = (sqd_uint64) fmt;
#if 0 /* hack : no byteswap */
  f1 = sre_htons(f1);
  f2 = sre_htonl(f2);
#endif

  if (fwrite(fname, 1, GSI64_KEYSIZE, fp) < GSI64_KEYSIZE) PANIC;
  if (fwrite(&f1, 2, 1, fp) < 1)    PANIC;
  if (fwrite(&f2, 8, 1, fp) < 1)    PANIC;
  return 1;
}


/* Function: GSI64WriteKeyRecord()
 * Date:     SRE, Wed Aug  5 10:52:30 1998 [St. Louis]
 *
 * Purpose:  Write a key record to a GSI64 file.
 *
 * Args:     fp      - open GSI64 file for writing
 *           key     - key (max 31 char + \0)
 *           fileidx - which file number to find this key in
 *           offset  - offset for this key       
 * 
 * Returns:  1 on success, else 0.
 *           will fail if key >= 32 chars, for instance.
 */
int
GSI64WriteKeyRecord(FILE *fp, char *key, int fileidx, long long offset)
{
  sqd_uint16 f1;
  sqd_uint64 f2;

  if (strlen(key) >= GSI64_KEYSIZE) return 0;
  if (fileidx > SQD_UINT16_MAX) Die("GSI64: file index out of range");
  if (offset  > SQD_UINT64_MAX) Die("GSI64: offset out of range");

  f1 = (sqd_uint16) fileidx;
  f2 = (sqd_uint64) offset;
#if 0 /* HACK! */
  f1 = sre_htons(f1);
  f2 = sre_htonl(f2);
#endif 
  
  if (fwrite(key, 1, GSI64_KEYSIZE, fp) < GSI64_KEYSIZE) PANIC;
  if (fwrite(&f1, 2,  1, fp) < 1) PANIC;
  if (fwrite(&f2, 8,  1, fp) < 1) PANIC;
  return 1;
}

#endif /*USE_GSI64 */
