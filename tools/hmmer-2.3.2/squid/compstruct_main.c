/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* compstruct_main.c
 * SRE, Tue Aug 30 10:35:31 1994
 * 
 * Compare RNA secondary structures. 
 * CVS $Id: compstruct_main.c,v 1.5 2003/04/14 16:00:16 eddy Exp $
 */

#include "squidconf.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "squid.h"
#include "msa.h"

static char banner[] = "compalign - compare test RNA secondary structure predictions to trusted set";

char usage[]  = "\
Usage: compstruct [-options] <trusted file> <test file>\n\
  Both files must contain secondary structure markup (e.g. Stockholm, SQUID,\n\
  SELEX formats), and sequences must occur in the same order in the two files.\n\
\n\
  Available options are:\n\
   -h : print short help and usage info\n\
";

static char experts[] = "\
   --informat <s> : specify that both alignments are in format <s> (SELEX, for instance)\n\
   --quiet        : suppress verbose header (used in regression testing)\n\
"; 

struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },     
  { "--informat", FALSE, sqdARG_STRING },
  { "--quiet",    FALSE, sqdARG_NONE   },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static int  KHS2ct(char *ss, int **ret_ct);
/* static void WriteCT(FILE *fp, char *seq, int *ct, int len); */

int
main(int argc, char **argv)
{
  char     *kfile, *tfile;      /* known, test structure file      */
  int       format;		/* expected format of kfile, tfile */
  SQFILE   *kfp, *tfp;          /* open kfile, tfile               */
  char     *kseq, *tseq;        /* known, test sequence            */
  SQINFO    kinfo, tinfo;       /* known, test info                */
  int      *kct, *tct;          /* known, test CT rep of structure */
  int       pos;
  int       nseq;

  int correct;			/* count of correct base pair predictions */
  int missedpair;		/* count of false negatives               */
  int falsepair;		/* count of false positives               */
  int tot_trusted;		/* total base pairs in trusted structure  */
  int tot_predicted;		/* total base pairs in predicted structure*/
  int tot_correct;		/* cumulative total correct pairs   */

  int dscorrect;		/* count of correct 2-state paired prediction   */
  int sscorrect;		/* count of correct 2-state unpaired prediction */
  int tot_dscorrect;
  int tot_sscorrect;
  int tot_positions;

  int quiet;			/* TRUE to silence verbose banner */

  char *optname; 
  char *optarg;
  int   optind;
  
  /***********************************************
   * Parse command line
   ***********************************************/

  format = MSAFILE_UNKNOWN;
  quiet  = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))
    {
      if      (strcmp(optname, "--quiet") == 0)  quiet  = TRUE; 
      else if (strcmp(optname, "--informat") == 0) {
	format = String2SeqfileFormat(optarg);
	if (format == MSAFILE_UNKNOWN) 
	  Die("unrecognized sequence file format \"%s\"", optarg);
	if (! IsAlignmentFormat(format))
	  Die("%s is an unaligned format, can't read as an alignment", optarg);
      }
      else if (strcmp(optname, "-h") == 0) {
	SqdBanner(stdout, banner);
	puts(usage);
	puts(experts);
        exit(EXIT_SUCCESS);
      }
    }

  if (argc - optind != 2) 
    Die("Incorrect number of command line arguments.\n%s\n", usage); 

  kfile = argv[optind++];
  tfile = argv[optind];
  
  if (! quiet) SqdBanner(stdout, banner);

  /***********************************************
   * Open the files
   ***********************************************/

  if ((kfp = SeqfileOpen(kfile, format, NULL)) == NULL)
    Die("Failed to open trusted structure file %s for reading", kfile);
  if ((tfp = SeqfileOpen(tfile, format, NULL)) == NULL)
    Die("Failed to open test structure file %s for reading", tfile);
  
  /***********************************************
   * Do structure comparisons, one seq at a time
   ***********************************************/

  tot_trusted = tot_predicted = tot_correct = 0;
  tot_dscorrect = tot_sscorrect = tot_positions = 0;
  nseq = 0;
  while (ReadSeq(kfp, kfp->format, &kseq, &kinfo) && ReadSeq(tfp, tfp->format, &tseq, &tinfo))
    {
      if (!quiet && strcmp(tinfo.name, kinfo.name) != 0)
	Warn("Trusted sequence %s, test sequence %s -- names not identical\n",
	     kinfo.name, tinfo.name);
      if (!quiet && strcmp(kseq, tseq) != 0)
	Warn("Trusted sequence %s, test sequence %s -- sequences not identical\n",
	     kinfo.name, tinfo.name);

      printf("%s %s\n", kinfo.name, (kinfo.flags & SQINFO_DESC) ? kinfo.desc : "");

      if (! (tinfo.flags & SQINFO_SS) && ! (kinfo.flags & SQINFO_SS))
	printf("[no test or trusted structure]\n\n");
      else if (! (tinfo.flags & SQINFO_SS))
	printf("[no test structure]\n\n");
      else if (! (kinfo.flags & SQINFO_SS))
	printf("[no trusted structure]\n\n");
      else
	{
	  if (! KHS2ct(kinfo.ss, &kct))
	    { printf("[bad trusted structure]\n"); goto CLEANUP;}
	  if (! KHS2ct(tinfo.ss, &tct))
	    { printf("[bad test structure]\n"); free(kct); goto CLEANUP; }
	  
/*	  WriteCT(stdout, tseq, tct, tinfo.len); */
/*	  WriteCT(stdout, tseq, kct, tinfo.len); */

	  correct = falsepair = missedpair = 0;
	  dscorrect = sscorrect = 0;
	  for (pos = 0; pos < kinfo.len; pos++)
	    {
				/* check if actual base pair is predicted */
	      if (kct[pos] >= 0 && kct[pos] == tct[pos])
		correct++;
	      else if (kct[pos] >= 0)
		missedpair++;

	      if (tct[pos] >= 0 && kct[pos] != tct[pos])
		falsepair++;

				/* 2 state prediction */
	      if (kct[pos] >= 0 && tct[pos] >= 0)
		dscorrect++;
	      else if (kct[pos] < 0 && tct[pos] < 0)
		sscorrect++;
	    }
	  nseq++;
	  tot_trusted   += correct + missedpair;
	  tot_predicted += correct + falsepair;
	  tot_correct   += correct;

	  tot_dscorrect += dscorrect;
	  tot_sscorrect += sscorrect;
	  tot_positions += kinfo.len;
	  
				/* print out per sequence info */
	  printf("   %d/%d trusted pairs predicted (%.2f%% sensitivity)\n", 
		 correct, correct+missedpair, 
		 100. * (float) correct/ (float) (correct + missedpair));
	  printf("   %d/%d predicted pairs correct (%.2f%% specificity)\n",
		 correct, correct + falsepair,
		 100. * (float) correct/ (float) (correct + falsepair));
	  
	  printf("   Two state: %d/%d positions correctly predicted (%.2f%% accuracy)\n", 
		 dscorrect + sscorrect,
		 kinfo.len,
		 100. * (float) (dscorrect + sscorrect) / (float) kinfo.len);
	  puts("");


	  free(kct);
	  free(tct);
	}

    CLEANUP:
      FreeSequence(kseq, &kinfo);
      FreeSequence(tseq, &tinfo);
    }

  /* And the final summary:
   */
  puts("");
  printf("Overall structure prediction accuracy (%d sequences, %d positions)\n",
	 nseq, tot_positions);
  printf("   %d/%d trusted pairs predicted (%.2f%% sensitivity)\n", 
	 tot_correct, tot_trusted, 
	 100. * (float) tot_correct/ (float) tot_trusted);
  printf("   %d/%d predicted pairs correct (%.2f%% specificity)\n",
	 tot_correct, tot_predicted, 
	 100. * (float) tot_correct/ (float) tot_predicted);
  printf("   Two state: %d/%d positions correctly predicted (%.2f%% accuracy)\n", 
	 tot_dscorrect + tot_sscorrect, tot_positions,
	 100. * (float) (tot_dscorrect + tot_sscorrect) / (float) tot_positions);
  puts("");

  SeqfileClose(tfp);
  SeqfileClose(kfp);
  return 0;
}


/* Function: KHS2ct()
 * 
 * Purpose:  Convert a secondary structure string to an array of integers
 *           representing what position each position is base-paired 
 *           to (0..len-1), or -1 if none. This is off-by-one from a
 *           Zuker .ct file representation.
 *           
 *           The .ct representation can accomodate pseudoknots but the 
 *           secondary structure string cannot easily; the string contains
 *           "Aa", "Bb", etc. pairs as a limited representation of
 *           pseudoknots. The string contains "><" for base pairs.
 *           Other symbols are ignored.
 *           
 * Return:   ret_ct is allocated here and must be free'd by caller.
 *           Returns 1 on success, 0 if ss is somehow inconsistent.
 */
static int 
KHS2ct(char *ss, int **ret_ct)
{
  struct intstack_s *dolist[27];
  int *ct;
  int  i;
  int  pos, pair;
  int  status = 1;		/* success or failure return status */
  int  len;

  for (i = 0; i < 27; i++)
    dolist[i] = InitIntStack();
  len = strlen(ss);

  if ((ct = (int *) malloc (len * sizeof(int))) == NULL)
    Die("malloc failed");
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; ss[pos] != '\0'; pos++)
    {
      if (ss[pos] == '>')	/* left side of a pair: push onto stack 0 */
	PushIntStack(dolist[0], pos);
      else if (ss[pos] == '<')	/* right side of a pair; resolve pair */
	{
	  if (! PopIntStack(dolist[0], &pair))
	    { status = 0; }
	  else
	    {
	      ct[pos]  = pair;
	      ct[pair] = pos;
	    }
	}
				/* same stuff for pseudoknots */
      else if (isupper((int) ss[pos]))
	PushIntStack(dolist[ss[pos] - 'A' + 1], pos);
      else if (islower((int) ss[pos]))
	{
	  if (! PopIntStack(dolist[ss[pos] - 'a' + 1], &pair))
	    { status = 0; }
	  else
	    {
	      ct[pos]  = pair;
	      ct[pair] = pos;
	    }
	}
      else if (!isgap(ss[pos])) status = 0; /* bad character */
    }

  for (i = 0; i < 27; i++)
    if ( FreeIntStack(dolist[i]) > 0)
      status = 0;

  *ret_ct = ct;
  return status;
}


#ifdef SRE_REMOVED
/* Function: WriteCT()
 * 
 * Purpose:  Write a CT representation of a structure.
 *           Written in 1..len sense, with 0 for unpaired
 *           positions.
 */
static void
WriteCT(FILE *fp, char *seq, int *ct, int len)
{
  int pos;
  for (pos = 0; pos < len; pos++)
    fprintf(fp, "%d %c %d\n", pos+1, seq[pos], ct[pos]+1);
}
#endif
