/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* sreformat_main.c
 * Mon Sep 13 13:06:51 1993
 * 
 * sreformat - reformat sequence files.
 * renamed sreformat from reformat, Tue Jun 30 10:53:38 1998
 *
 * CVS $Id: sreformat_main.c,v 1.19 2003/04/14 16:00:16 eddy Exp $
 */


#include "squidconf.h"

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "squid.h"
#include "msa.h"

static char banner[] = "sreformat - convert between sequence formats";

static char usage[] = "\
Usage: sreformat [-options] <format> <seqfile>\n\
  Output format choices: Unaligned      Aligned\n\
                         -----------    -------\n\
                         fasta          stockholm\n\
                         embl           msf\n\
                         genbank        a2m\n\
                         gcg            phylip\n\
                         gcgdata        clustal\n\
                         pir            selex\n\
                         raw            eps\n\n\
  Available options are:\n\
    -h : help; print brief help on version and usage\n\
    -d : force DNA alphabet for nucleic acid sequence\n\
    -r : force RNA alphabet for nucleic acid sequence\n\
    -l : force lower case\n\
    -u : force upper case\n\
    -x : convert non-IUPAC chars in DNA to N's for IUPAC/BLAST compatibility\n\
";

static char experts[] = "\
  Expert options:\n\
    --informat <s>: input sequence file is in format <s>\n\
    --mingap      : remove columns containing all gaps (seqfile=alignment)\n\
    --nogap       : remove columns containing any gaps (seqfile=alignment)\n\
    --pfam        : modify Stockholm format output to be in PFAM style (1 line/seq)\n\
    --sam         : try to convert gaps to SAM style (seqfile=alignment)\n\
    --samfrac <x> : convert to SAM convention; cols w/ gapfrac > x are inserts\n\
    --gapsym <c>  : convert all gaps to character '<c>'\n\
";

static struct opt_s OPTIONS[] = {
  { "-d", TRUE, sqdARG_NONE },
  { "-h", TRUE, sqdARG_NONE },
  { "-l", TRUE, sqdARG_NONE },
  { "-r", TRUE, sqdARG_NONE },
  { "-u", TRUE, sqdARG_NONE },
  { "-x", TRUE, sqdARG_NONE },
  { "--gapsym",  FALSE, sqdARG_CHAR },
  { "--informat",FALSE, sqdARG_STRING }, 
  { "--mingap",  FALSE, sqdARG_NONE },
  { "--nogap",   FALSE, sqdARG_NONE },
  { "--pfam",    FALSE, sqdARG_NONE },
  { "--sam",     FALSE, sqdARG_NONE },
  { "--samfrac", FALSE, sqdARG_FLOAT },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  char     *seqfile;            /* name of sequence file */
  char     *format;
  SQFILE   *dbfp;		/* open sequence file */
  int       fmt;		/* format of seqfile  */
  int       outfmt;		/* output format */
  char     *seq;		/* sequence */
  SQINFO    sqinfo;
  int       i;

  int    force_rna;		/* TRUE to force RNA alphabet */
  int    force_dna;		/* TRUE to force DNA alphabet */
  int    force_lower;		/* TRUE to force lower case   */
  int    force_upper;		/* TRUE to force upper case   */
  int    x_is_bad;		/* TRUE to convert X to N     */
  int    do_mingap;		/* TRUE to remove columns containing all gaps */
  int    do_nogap;		/* TRUE to remove columns containing any gaps */
  int    do_pfam;		/* TRUE to make SELEX -> PFAM */
  int    samize;		/* TRUE to SAMize an A2M conversion */
  float  samfrac;		/* -1, or gap fraction for a SAM conversion */
  int    expect_alignment;	/* TRUE to expect an input alignment to convert */
  char   gapsym;		/* 0 if unset; else = character to use for gaps */

  char *optname;                /* name of option found by Getopt()      */
  char *optarg;                 /* argument found by Getopt()            */
  int   optind;                 /* index in argv[]                       */

  /***********************************************
   * Parse command line
   ***********************************************/

  force_rna        = FALSE;
  force_dna        = FALSE;
  force_upper      = FALSE;
  force_lower      = FALSE;
  x_is_bad         = FALSE;
  do_mingap        = FALSE;
  do_nogap         = FALSE;
  do_pfam          = FALSE;   
  samize           = FALSE;
  samfrac          = -1.0;
  fmt              = SQFILE_UNKNOWN;
  expect_alignment = FALSE;
  gapsym           = 0;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-a")        == 0) expect_alignment= TRUE;
    else if (strcmp(optname, "-d")        == 0) force_dna   = TRUE;
    else if (strcmp(optname, "-l")        == 0) force_lower = TRUE;
    else if (strcmp(optname, "-r")        == 0) force_rna   = TRUE;
    else if (strcmp(optname, "-u")        == 0) force_upper = TRUE;
    else if (strcmp(optname, "-x")        == 0) x_is_bad    = TRUE;
    else if (strcmp(optname, "--gapsym")  == 0) gapsym      = *optarg;
    else if (strcmp(optname, "--mingap")  == 0) do_mingap   = TRUE;
    else if (strcmp(optname, "--nogap")   == 0) do_nogap    = TRUE;
    else if (strcmp(optname, "--pfam")    == 0) do_pfam     = TRUE;
    else if (strcmp(optname, "--sam")     == 0) samize      = TRUE;
    else if (strcmp(optname, "--samfrac") == 0) samfrac     = atof(optarg);
    else if (strcmp(optname, "--informat") == 0) {
      fmt = String2SeqfileFormat(optarg);
      if (fmt == SQFILE_UNKNOWN) 
	Die("unrecognized sequence file format \"%s\"", optarg);
    }
    else if (strcmp(optname, "-h") == 0) {
      SqdBanner(stdout, banner);
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind != 2)
    Die("%s\n", usage); 
  if (force_lower && force_upper)
    Die("Can't force both upper case and lower case. Stop trying to confuse me.\n%s", 
	usage);
  if (force_rna && force_dna)
    Die("Can't force both RNA and DNA. Stop trying to find bugs. You'll be sorry.\n%s", 
	usage);

  format  = argv[optind]; optind++;
  seqfile = argv[optind]; optind++;
  
  /* Try to work around inability to autodetect from a pipe or .gz:
   * assume FASTA format
   */
  if (fmt == SQFILE_UNKNOWN &&
      (Strparse("^.*\\.gz$", seqfile, 0) || strcmp(seqfile, "-") == 0))
    fmt = SQFILE_FASTA;

  /***********************************************
   * Figure out what format we're supposed to write
   ***********************************************/

  if ((outfmt = String2SeqfileFormat(format)) == SQFILE_UNKNOWN)
    Die("Unknown output format %s\n%s", format, usage);

  /***********************************************
   * Reformat the file, printing to stdout.
   ***********************************************/

  /* If the output format is an alignment, then the input format
   * has to be an alignment.
   */
  if (IsAlignmentFormat(outfmt))
    {
      MSAFILE *afp;
      MSA     *msa;

      if ((afp = MSAFileOpen(seqfile, fmt, NULL)) == NULL)
	Die("Alignment file %s could not be opened for reading", seqfile);

      while ((msa = MSAFileRead(afp)) != NULL)
	{
	  /* If asked, convert upper/lower convention and
	   * gap character conventions now
	   */
	  if (do_mingap)    MSAMingap(msa);
	  if (do_nogap)     MSANogap(msa);
	  if (gapsym)       AlignmentHomogenousGapsym(msa->aseq, msa->nseq, msa->alen, gapsym);
	  if (samize)       SAMizeAlignment(msa->aseq, msa->nseq, msa->alen);
	  if (samfrac >= 0) SAMizeAlignmentByGapFrac(msa->aseq, msa->nseq, msa->alen, samfrac);

	  for (i = 0; i < msa->nseq; i++)
	    {
	      if (force_dna)   ToDNA(msa->aseq[i]);
	      if (force_rna)   ToRNA(msa->aseq[i]);
	      if (x_is_bad)    ToIUPAC(msa->aseq[i], TRUE);
	      if (force_lower) s2lower(msa->aseq[i]);
	      if (force_upper) s2upper(msa->aseq[i]);
	    }
      
	  /* This code block can be replaced with a 
	   * MSAFileWrite() call someday... SRE Sun Apr 22 19:17:19 2001
	   */
	  switch (outfmt) {
	  case MSAFILE_A2M:       WriteA2M(stdout, msa);         break;
	  case MSAFILE_CLUSTAL:   WriteClustal(stdout, msa);     break;
	  case MSAFILE_MSF:       WriteMSF(stdout, msa);         break;
	  case MSAFILE_PHYLIP:    WritePhylip(stdout, msa);      break;
	  case MSAFILE_SELEX:     
	    if (do_pfam) WriteSELEXOneBlock(stdout, msa);
	    else         WriteSELEX(stdout, msa);       
	    break;
	  case MSAFILE_EPS:       EPSWriteSmallMSA(stdout, msa); break;
	  case MSAFILE_STOCKHOLM:
	    if (do_pfam) WriteStockholmOneBlock(stdout, msa);
	    else         WriteStockholm(stdout, msa);
	    break;
	  default:
	    Die("can't write. no such alignment format %d\n", outfmt);
	  }
	  
	  MSAFree(msa);
	}
      MSAFileClose(afp);
    }
  else
    {
      if ((dbfp = SeqfileOpen(seqfile, fmt, NULL)) == NULL)
	Die("Failed to open sequence file %s for reading", seqfile);
  
      while (ReadSeq(dbfp, fmt, &seq, &sqinfo))
	{
	  if (force_dna)   ToDNA(seq);
	  if (force_rna)   ToRNA(seq);
	  if (x_is_bad)    ToIUPAC(seq, FALSE);
	  if (force_lower) s2lower(seq);
	  if (force_upper) s2upper(seq);

	  WriteSeq(stdout, outfmt, seq, &sqinfo);
	  FreeSequence(seq, &sqinfo);
	}
      SeqfileClose(dbfp);
    }

  return 0;
}

