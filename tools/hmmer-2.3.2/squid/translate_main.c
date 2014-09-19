/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* translate_main.c
 * 
 * translate - create a file of all possible protein ORFs, given
 *             an input nucleic acid sequence
 *
 * 
 * Not currently compliant w/ HMMER API.
 * 
 * 1.02 Thu Apr 20 16:12:41 1995
 *     + incorporated into squid
 *     +  -a, -s options added
 *
 * CVS $Id: translate_main.c,v 1.7 2003/04/14 16:00:16 eddy Exp $
 */

#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "squid.h"

#ifdef NEED_GETOPTH
#include <getopt.h>
#endif

#define OPTIONS "ahl:mo:qs:"

static char usage[] = "\
Usage: translate [-options] <seqfile>\n\
   Translate a nucleic acid sequence into protein ORFs.\n\
   Available options are:\n\
   -a            : translate in full, with stops; no individual ORFs\n\
   -h            : help; show brief usage and version info\n\
   -l <minlen>   : report only ORFs greater than minlen (default 20)\n\
   -m            : require ORFs to start with AUG/Met\n\
   -o <outfile>  : save results in output file\n\
   -q            : quiet; silence banner, for piping or redirection\n\
   -s <stopchar> : with -a, set stop character to <stopchar>\n";

int
main(int argc, char **argv)
{
  char        *seqfile;         /* name of seq file to read             */
  SQFILE      *seqfp;		/* ptr to opened seq file               */
  int          format;		/* format of sequence file              */
  char        *seq;             /* ptr to current sequence              */
  SQINFO       sqinfo;          /* sequence information                 */
  char        *revseq;		/* reverse complement of seq            */
  int          start, end;	/* coords of ORF in current seq         */
  int          orfnumber;	/* counter for ORFs in current seq      */
  char        *aaseq[6];        /* full translations in all 6 frames    */
  char        *orf;             /* ptr to translated ORF sequence       */
  char        *sptr;		/* ptr into orf                         */
  int          len;		/* length of an ORF                     */
  int          frame;		/* counter for frames (3..5 are reverse)*/

  int          minimum_len;	/* minimum length of ORFs to print out  */
  char        *outfile;         /* file to save output in               */
  FILE        *ofp;		/* where to direct output               */
  char         stopchar;	/* what to use as a stop character      */
  int          keepstops;	/* TRUE to do six big ORFs              */
  int          quiet;		/* TRUE to silence banner               */
  int          require_met;	/* TRUE to start orfs with M            */

  int          optchar;		/* option character */
  extern char *optarg;          /* for getopt() */
  extern int   optind;		/* for getopt() */

  /***********************************************
   * Parse the command line
   ***********************************************/

  format      = SQFILE_UNKNOWN;	/* autodetect by default */
  minimum_len = 20;
  outfile     = NULL;
  stopchar    = '*';
  keepstops   = FALSE;
  quiet       = FALSE;
  require_met = FALSE;

  while ((optchar = getopt(argc, argv, OPTIONS)) != -1)
    switch (optchar) {

    case 'a': keepstops   = TRUE;         break;
    case 'l': minimum_len = atoi(optarg); break;
    case 'm': require_met = TRUE;         break;
    case 'o': outfile     = optarg;       break;
    case 'q': quiet       = TRUE;         break;
    case 's': stopchar    = *optarg;      break;

    case 'h':
      printf("translate %s, %s\n%s\n", SQUID_VERSION, SQUID_DATE, usage);
      exit(EXIT_SUCCESS);
    default:
      Die("%s\n", usage);
    }

  if (argc - optind != 1)
    Die("Incorrect number of command line arguments\n%s\n", usage);

  seqfile = argv[optind];
  
  /***********************************************
   * Open sequence file and output file
   ***********************************************/

  seqfp = SeqfileOpen(seqfile, format, NULL);
  if (seqfp == NULL)
    Die("Failed to open sequence file %s\n%s\n", 
	seqfile, usage);

  if (outfile != NULL)
    {
      if ((ofp = fopen(outfile, "w")) == NULL)
	Die("Failed to open output file %s\n", outfile);
    }
  else
    ofp = stdout;
	

  /***********************************************
   * Main routine
   ***********************************************/

  if (! quiet) printf("translate %s, %s\n", SQUID_VERSION, SQUID_DATE);

  while (ReadSeq(seqfp, seqfp->format, &seq, &sqinfo))
    {
      s2upper(seq); 
      revseq = (char *) malloc (sqinfo.len + 1);
      revcomp(revseq, seq);
      orfnumber = 1;

				/* Translate seq in all six frames */
      aaseq[0] = Translate(seq, stdcode1);
      aaseq[1] = Translate(seq + 1, stdcode1);
      aaseq[2] = Translate(seq + 2, stdcode1);
      aaseq[3] = Translate(revseq, stdcode1);
      aaseq[4] = Translate(revseq + 1, stdcode1);
      aaseq[5] = Translate(revseq + 2, stdcode1);
      


      if (keepstops)
	{			/* full translation including stops */
	  for (frame = 0; frame < 6; frame++)
	    { 
	      fprintf(ofp, "> %s:%d", sqinfo.name, frame);
	      for (sptr = aaseq[frame]; *sptr; sptr++)
		{
		  if (*sptr == '*') *sptr = stopchar;
		  if (! ((sptr - aaseq[frame]) % 50)) putc('\n', ofp);
		  putc((int) *sptr, ofp);
		}
	      putc('\n', ofp);
	    }		  
	}
      else
	{			/* Print all decent ORF's in FASTA format */
	  for (frame = 0; frame < 6; frame++)
	    {
				/* initialize strtok on the first ORF;
				   termination codons are '*' symbols */
	      orf = strtok(aaseq[frame], "*");
	      while (orf != NULL && *orf != '\0')
		{
		  if (require_met) {
		    while (*orf != 'M' && *orf != '\0') orf++;
		  } 

		  if (*orf != '\0') {
		    len = strlen(orf);	      
		    if (len > minimum_len)
		      {
				/* calculate coords */
			start = (orf - aaseq[frame]) * 3 + 1;
			if (frame < 3) start += frame; /* frame corrections */
			else       start -= frame-3;
		      
			if (frame < 3) 
			  end = start + len * 3 - 1;
			else
			  {
			    start = -1 * (start - sqinfo.len - 1);
			    end = start - len * 3 + 1;
			  }
		  
			fprintf(ofp, "> %s.%d    length %d, nt %d..%d",
				sqinfo.name,
				orfnumber,
				len,
				start,
				end);

			for (sptr = orf; *sptr; sptr++)
			  {
			    if (! ((sptr - orf) % 50))
			      putc('\n', ofp);
			    putc((int) *sptr, ofp);
			  }
			putc('\n', ofp);
		  
			orfnumber++;
		      }
		  }
				/* pick off next orf */
		  orf = strtok(NULL, "*");
		  
		}
	    }
	}

      for (frame = 0; frame < 6; frame++)
	free(aaseq[frame]);
      FreeSequence(seq, &sqinfo);
      free(revseq);
    }

  SeqfileClose(seqfp);

  /**************************************************
   * Successful return to invocation environment
   **************************************************/
  return 0;
}

