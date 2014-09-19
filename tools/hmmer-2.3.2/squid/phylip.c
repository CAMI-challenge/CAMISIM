/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* phylip.c
 * SRE, Mon Jun 14 14:08:33 1999 [St. Louis]
 * 
 * Import/export of PHYLIP interleaved multiple sequence alignment
 * format files.
 * 
 * CVS $Id: phylip.c,v 1.3 2003/04/14 16:00:16 eddy Exp $
 */

#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"

#ifdef TESTDRIVE_PHYLIP
/*****************************************************************
 * phylip.c test driver:
 * 
 */
int 
main(int argc, char **argv)
{
  MSAFILE *afp;
  MSA     *msa;
  char    *file;
  
  file = argv[1];

  if ((afp = MSAFileOpen(file, MSAFILE_UNKNOWN, NULL)) == NULL)
    Die("Couldn't open %s\n", file);

  printf("format %d\n", afp->format);

  while ((msa = ReadPhylip(afp)) != NULL)
    {
      WritePhylip(stdout, msa);
      MSAFree(msa); 
    }
  
  MSAFileClose(afp);
  exit(0);
}
/******************************************************************/
#endif /* testdrive_phylip */



/* Function: ReadPhylip()
 * Date:     SRE, Fri Jun 18 12:59:37 1999 [Sanger Centre]
 *
 * Purpose:  Parse an alignment from an open Phylip format
 *           alignment file. Phylip is a single-alignment format.
 *           Return the alignment, or NULL if we have no data.
 *
 * Args:     afp - open alignment file
 *
 * Returns:  MSA * - an alignment object
 *                   Caller responsible for an MSAFree()
 *           NULL if no more alignments        
 */
MSA *
ReadPhylip(MSAFILE *afp)
{
  MSA  *msa;
  char *s, *s1, *s2;
  char  name[11];		/* seq name max len = 10 char */
  int   nseq, alen;
  int   idx;			/* index of current sequence */
  int   slen;
  int   nblock;
  
  if (feof(afp->f)) return NULL;

  /* Skip until we see a nonblank line; it's the header,
   * containing nseq/alen
   */
  nseq = 0; alen = 0;
  while ((s = MSAFileGetLine(afp)) != NULL)
    {
      if ((s1 = sre_strtok(&s, WHITESPACE, NULL)) == NULL) continue;
      if ((s2 = sre_strtok(&s, WHITESPACE, NULL)) == NULL)
	Die("Failed to parse nseq/alen from first line of PHYLIP file %s\n", afp->fname);
      if (! IsInt(s1) || ! IsInt(s2))
	Die("nseq and/or alen not an integer in first line of PHYLIP file %s\n", afp->fname);
      nseq = atoi(s1);
      alen = atoi(s2);
      break;
    }

  msa = MSAAlloc(nseq, 0);
  idx    = 0;
  nblock = 0;
  while ((s = MSAFileGetLine(afp)) != NULL) 
    {
      /* ignore blank lines. nonblank lines start w/ nonblank char */
      if (isspace((int) *s)) continue;
				/* First block has seq names */
      if (nblock == 0) {
	strncpy(name, s, 10);
	name[10] = '\0';
	GKIStoreKey(msa->index, name);
	msa->sqname[idx] = sre_strdup(name, -1);
	s += 10;		
      }
				/* be careful of trailing whitespace on lines */
      if ((s1 = sre_strtok(&s, WHITESPACE, &slen)) == NULL)
	Die("Failed to parse sequence at line %d of PHYLIP file %s\n", 
	    afp->linenumber, afp->fname);
      msa->sqlen[idx] = sre_strcat(&(msa->aseq[idx]), msa->sqlen[idx], s1, slen);

      idx++;
      if (idx == nseq) { idx = 0; nblock++; }
    }
  msa->nseq = nseq;
  MSAVerifyParse(msa);		/* verifies; sets alen, wgt; frees sqlen[] */
  return msa;
}



/* Function: WritePhylip()
 * Date:     SRE, Fri Jun 18 12:07:41 1999 [Sanger Centre]
 *
 * Purpose:  Write an alignment in Phylip format to an open file.
 *
 * Args:     fp    - file that's open for writing.
 *           msa   - alignment to write. 
 *
 * Returns:  (void)
 */
void
WritePhylip(FILE *fp, MSA *msa)
{
  int    idx;			/* counter for sequences         */
  int    cpl = 50;		/* 50 seq char per line          */
  char   buf[51];		/* buffer for writing seq        */
  int    pos;

  /* First line has nseq, alen
   */
  fprintf(fp, " %d  %d\n", msa->nseq, msa->alen);

  /* Alignment section.
   * PHYLIP is a multiblock format, blocks (optionally) separated
   * by blanks; names only attached to first block. Names are
   * restricted to ten char; we achieve this by simple truncation (!).
   * (Do we need to convert gap characters from our ./- convention?)
   */
  for (pos = 0; pos < msa->alen; pos += cpl)
    {
      if (pos > 0) fprintf(fp, "\n");

      for (idx = 0; idx < msa->nseq; idx++)
	{
	  strncpy(buf, msa->aseq[idx] + pos, cpl);
	  buf[cpl] = '\0';	      
	  if (pos > 0) fprintf(fp, "%s\n", buf);
	  else         fprintf(fp, "%-10.10s%s\n", msa->sqname[idx], buf);
	}
    }
  return;
}
