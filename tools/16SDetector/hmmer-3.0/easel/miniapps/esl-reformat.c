/* Convert sequence file formats
 *             
 * SRE, Sun Feb 27 08:24:33 2005
 * from squid's sreformat (1993).
 * SVN $Id: esl-reformat.c 545 2010-03-02 17:06:54Z eddys $            
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_sqio.h"
#include "esl_sq.h"
#include "esl_msa.h"
#include "esl_wuss.h"

static char banner[] = "convert sequence file formats";

static char usage[] = "[-options] <format> <seqfile>\n\
  Output format choices: Unaligned      Aligned\n\
                         -----------    -------\n\
                         fasta          stockholm\n\
                                        pfam\n\
                                        a2m\n\
                                        psiblast\n\
                                        afa\n\
\n";

static ESL_OPTIONS options[] = {
   /* name          type        default env   range togs  reqs  incompat                     help                                      docgroup */
  { "-d",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-r",                  "convert to DNA alphabet (U->T)",                     0 },
  { "-h",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL,                  "help; print brief info on version and usage",        0 },
  { "-l",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-u",                  "convert to lower case",                              0 },
  { "-n",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-x",                  "remove DNA IUPAC codes; convert ambig chars to N",   0 },
  { "-o",         eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, NULL,                  "send output to file <f>, not stdout",                0 },
  { "-r",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-d",                  "convert to RNA alphabet (T->U)",                     0 }, 
  { "-u",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-l",                  "convert to upper case",                              0 },
  { "-x",         eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "-n",                  "convert non-IUPAC chars (e.g. X) in DNA to N",       0 },
  { "--gapsym",   eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, "--mingap,--nogap",    "convert all gaps to character <c>",                  0 },
  { "--informat", eslARG_STRING,  NULL, NULL, NULL, NULL, NULL, NULL,                  "input sequence file is in format <s>",               0 },
  { "--mingap",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--nogap",             "remove columns containing all gaps (seqfile=MSA)",   0 },
  { "--nogap",    eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--mingap,--gapsym",   "remove columns containing any gaps (seqfile=MSA)",   0 },
  { "--wussify",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--dewuss,--fullwuss", "convert old RNA structure markup lines to WUSS",     0 },
  { "--dewuss",   eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--wussify,--fullwuss","convert WUSS RNA structure markup to old format",    0 },
  { "--fullwuss", eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, "--wussify,--dewuss",  "convert simple WUSS notation to full (output) WUSS", 0 },
  { "--ignore",   eslARG_STRING, FALSE, NULL, NULL, NULL, NULL, NULL,                  "ignore input seq characters listed in string <s>",   0 },
  { "--acceptx",  eslARG_STRING, FALSE, NULL, NULL, NULL, NULL, NULL,                  "accept input seq chars in string <s> as X",          0 },
  { "--rename",   eslARG_STRING, FALSE, NULL, NULL, NULL, NULL, NULL,                  "rename and number each sequence <s>.<n>",            0 },
  { 0,0,0,0,0,0,0,0 },
};

static void symconvert(char *s, char *oldsyms, char *newsyms);

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;	        /* application configuration               */
  char        *outformat;	/* output format as a string               */
  char        *infile;          /* name of input sequence file             */
  int          infmt;		/* input format as a code; eslSQFILE_FASTA */
  int          outfmt;		/* output format as a code                 */
  int          status;		/* return code from an Easel call          */
  FILE        *ofp;		/* output stream                           */

  char  *informat;		/* input format as string; "fasta"           */
  char  *outfile;		/* output file, or NULL                      */
  int    force_rna;		/* TRUE to force RNA alphabet                */
  int    force_dna;		/* TRUE to force DNA alphabet                */
  int    force_lower;		/* TRUE to force lower case                  */
  int    force_upper;		/* TRUE to force upper case                  */
  int    iupac_to_n;            /* TRUE to convert ambiguities all to N's    */
  int    x_is_bad;		/* TRUE to convert X to N                    */
  int    do_mingap;		/* TRUE to remove cols containing all gaps   */
  int    do_nogap;		/* TRUE to remove cols containing any gaps   */
  char  *gapsym;		/* NULL if unset; else, char for gaps        */
  int    wussify;		/* TRUE to convert old KH SS markup to WUSS  */
  int    dewuss;		/* TRUE to convert WUSS back to old KH       */
  int    fullwuss;		/* TRUE to convert simple WUSS to full WUSS  */
  char  *rename; 		/* if non-NULL rename seqs to <s>.<n>        */
  int    idx;
  char   errbuf[eslERRBUFSIZE]; /* for error messages                        */

  /*****************************************************************
   * Parse the command line
   *****************************************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "-h"))
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("  where options are:\n");
      esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0= group; 2 = indentation; 80=textwidth*/
      exit(EXIT_SUCCESS);
    }

  force_dna   = esl_opt_GetBoolean(go, "-d");
  force_lower = esl_opt_GetBoolean(go, "-l");
  iupac_to_n  = esl_opt_GetBoolean(go, "-n");
  outfile     = esl_opt_GetString (go, "-o");
  force_rna   = esl_opt_GetBoolean(go, "-r");
  force_upper = esl_opt_GetBoolean(go, "-u");
  x_is_bad    = esl_opt_GetBoolean(go, "-x");
  gapsym      = esl_opt_GetString( go, "--gapsym");
  informat    = esl_opt_GetString( go, "--informat");
  do_mingap   = esl_opt_GetBoolean(go, "--mingap");
  do_nogap    = esl_opt_GetBoolean(go, "--nogap");
  wussify     = esl_opt_GetBoolean(go, "--wussify");
  dewuss      = esl_opt_GetBoolean(go, "--dewuss");
  fullwuss    = esl_opt_GetBoolean(go, "--fullwuss");
  rename      = esl_opt_GetString (go, "--rename");

  if (esl_opt_ArgNumber(go) != 2) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  outformat = esl_opt_GetArg(go, 1);
  infile    = esl_opt_GetArg(go, 2);

  infmt = eslSQFILE_UNKNOWN;
  if (informat != NULL)
    {
      infmt = esl_sqio_EncodeFormat(informat);
      if (infmt == eslSQFILE_UNKNOWN)
	esl_fatal("%s is not a recognized input seqfile format\n");
    }
    
  outfmt = esl_sqio_EncodeFormat(outformat);
  if (outfmt == eslSQFILE_UNKNOWN)
    esl_fatal("%s is not a recognized output seqfile format\n");

  if (gapsym != NULL && strlen(gapsym) != 1)
    esl_fatal("Argument to --gapsym must be a single character.");
  
  if (outfile == NULL) ofp = stdout;
  else if ((ofp = fopen(outfile, "w")) == NULL)
    esl_fatal("Failed to open output file %s\n", outfile);


  /***********************************************
   * Reformat the file, printing to stdout.
   ***********************************************/

  /* If the output format is an alignment, then the input format
   * has to be an alignment.
   */
  if (esl_sqio_IsAlignment(outfmt))
    {
      ESL_MSAFILE *afp;
      ESL_MSA     *msa;

      status = esl_msafile_Open(infile, infmt, NULL, &afp);
      if (status == eslENOTFOUND)
	esl_fatal("Alignment file %s not readable\n", infile);
      else if (status == eslEFORMAT) 
	esl_fatal("Couldn't determine format of alignment %s\n", infile);
      else if (status == eslEINVAL)
	esl_fatal("Can't autodetect format of stdin or .gz; use --informat\n");
      else if (status != eslOK) 
	esl_fatal("Alignment file open failed with error %d\n", status);

      if ( esl_opt_IsOn(go, "--ignore"))  esl_fatal("The --ignore option is unimplemented for alignment reformatting.");
      if ( esl_opt_IsOn(go, "--acceptx")) esl_fatal("The --acceptx option is unimplemented for alignment reformatting.");

      while ((status = esl_msa_Read(afp, &msa)) != eslEOF)
	{
	  if      (status == eslEFORMAT) esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
	  else if (status == eslEINVAL)  esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
	  else if (status != eslOK)      esl_fatal("Alignment file read failed with error code %d\n", status);

	  if (do_mingap)    if((status = esl_msa_MinimGaps(msa, errbuf, "-_.")) != eslOK) esl_fatal(errbuf);
	  if (do_nogap)     if((status = esl_msa_NoGaps   (msa, errbuf, "-_.")) != eslOK) esl_fatal(errbuf);
	  if (gapsym!=NULL) esl_msa_SymConvert(msa, "-_.", gapsym);
	  if (force_lower)  esl_msa_SymConvert(msa,
					       "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
					       "abcdefghijklmnopqrstuvwxyz");
	  if (force_upper)  esl_msa_SymConvert(msa,
					       "abcdefghijklmnopqrstuvwxyz",
					       "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	  if (force_rna)    esl_msa_SymConvert(msa, "Tt", "Uu");
	  if (force_dna)    esl_msa_SymConvert(msa, "Uu", "Tt");
	  if (iupac_to_n)   esl_msa_SymConvert(msa, 
					       "RYMKSWHBVDrymkswhbvd",
					       "NNNNNNNNNNnnnnnnnnnn");
	  if (x_is_bad)     esl_msa_SymConvert(msa, "Xx", "Nn");
	  
	  if (rename)
	    {
	      for (idx = 0; idx < msa->nseq; idx++)
		esl_msa_FormatSeqName(msa, idx, "%s.%d", rename, idx+1);
	    }

	  if (wussify)
	    {
	      if (msa->ss_cons != NULL) 
		esl_kh2wuss(msa->ss_cons, msa->ss_cons);
	      if (msa->ss != NULL)
		for (idx = 0; idx < msa->nseq; idx++)
		  if (msa->ss[idx] != NULL)
		    esl_kh2wuss(msa->ss[idx], msa->ss[idx]);
	    }

	  if (dewuss)
	    {
	      if (msa->ss_cons != NULL)
		esl_wuss2kh(msa->ss_cons, msa->ss_cons);
	      if (msa->ss != NULL)
		for (idx = 0; idx < msa->nseq; idx++)
		  if (msa->ss[idx] != NULL)
		    esl_wuss2kh(msa->ss[idx], msa->ss[idx]);
	    }

	  if (fullwuss)
	    {
	      if (msa->ss_cons != NULL)
		{
		  status = esl_wuss_full(msa->ss_cons, msa->ss_cons);
		  if (status == eslESYNTAX) 
		    esl_fatal("Bad consensus SS: not in WUSS format\n");
		  else if (status != eslOK)
		    esl_fatal("Conversion of SS_cons failed, code %d\n", status);
		}
	      if (msa->ss != NULL)
		for (idx = 0; idx < msa->nseq; idx++)
		  if (msa->ss[idx] != NULL)
		    {
		      status = esl_wuss_full(msa->ss[idx], msa->ss[idx]);
		      if (status == eslESYNTAX) 
			esl_fatal("Bad SS for %s: not in WUSS format\n",
				  msa->sqname[idx]);
		      else if (status != eslOK)
			esl_fatal("Conversion of SS for %s failed, code %d\n", 
				  msa->sqname[idx], status);
		    }
	    }

	  esl_msa_Write(ofp, msa, outfmt);
	  esl_msa_Destroy(msa);
	}
      esl_msafile_Close(afp);
    } /* end of alignment->alignment conversion */
  else
    { /* else: conversion to unaligned file formats */
      ESL_SQFILE  *sqfp;	/* open input sequence file                */
      ESL_SQ      *sq;		/* an input sequence                       */

      status = esl_sqfile_Open(infile, infmt, NULL, &sqfp);
      if (status == eslENOTFOUND)
	esl_fatal("Couldn't open seqfile %s\n", infile);
      else if (status == eslEFORMAT)
	esl_fatal("Couldn't determine format of seqfile %s\n", infile);      
      else if (status == eslEINVAL)
	esl_fatal("Can't autodetect format of stdin or .gz; use --informat\n");
      else if (status != eslOK)
	esl_fatal("Open of seqfile %s failed, code %d\n", infile, status);
      
      if ( esl_opt_IsOn(go, "--ignore"))  esl_sqio_Ignore  (sqfp, esl_opt_GetString(go, "--ignore"));
      if ( esl_opt_IsOn(go, "--acceptx")) esl_sqio_AcceptAs(sqfp, esl_opt_GetString(go, "--acceptx"), 'X');

      sq  = esl_sq_Create();
      idx = 0;
      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  if (force_lower) symconvert(sq->seq, 
				      "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
				      "abcdefghijklmnopqrstuvwxyz");
	  if (force_upper) symconvert(sq->seq, 
				      "abcdefghijklmnopqrstuvwxyz",
				      "ABCDEFGHIJKLMNOPQRSTUVWXYZ");
	  if (force_rna)   symconvert(sq->seq, "Tt", "Uu");
	  if (force_dna)   symconvert(sq->seq, "Uu", "Tt");
	  if (iupac_to_n)  symconvert(sq->seq, 
				      "RYMKSWHBVDrymkswhbvd",
				      "NNNNNNNNNNnnnnnnnnnn");
	  if (x_is_bad)    symconvert(sq->seq, "Xx", "Nn");
	  
	  if (wussify && sq->ss != NULL) esl_kh2wuss(sq->ss, sq->ss);	    
	  if (dewuss  && sq->ss != NULL) esl_wuss2kh(sq->ss, sq->ss);	    

	  if (fullwuss && sq->ss != NULL)
	    {
	      status = esl_wuss_full(sq->ss, sq->ss);
	      if (status == eslESYNTAX) 
		esl_fatal("Bad SS for %s: not in WUSS format\n", sq->name);
	      else if (status != eslOK)
		esl_fatal("Conversion of SS for %s failed, code %d\n", 
			  sq->name, status);
	    }

	  if (rename) esl_sq_FormatName(sq, "%s.%d", rename, idx+1);

	  esl_sqio_Write(ofp, sq, outfmt, FALSE);
	  esl_sq_Reuse(sq);
	  idx++;
	}
      /* status should be eslEOF on normal end; if it isn't, deal w/ error */
      if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					       sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					       status, sqfp->filename);
      
      esl_sq_Destroy(sq);
      esl_sqfile_Close(sqfp);
    } /* end of unaligned seq conversion */

  esl_getopts_Destroy(go);
  exit(0);
}

/* symconvert()
 * 
 * single seq version of esl_msa_SymConvert(); see
 * documentation there.
 * 
 * no reason yet to include in sqio API, but that may change.
 * 
 * inefficient to use this for upper/lower case conversion,
 * prob by an order of magnitude (because of the strchr() call,
 * which could be replaced by a range test), but I bet it's
 * unnoticeable.
 */
static void
symconvert(char *s, char *oldsyms, char *newsyms)
{
  int   pos;
  char *sptr;
  int   special;

  special = (strlen(newsyms) == 1 ? TRUE : FALSE);

  for (pos = 0; s[pos] != '\0'; pos++)
    if ((sptr = strchr(oldsyms, s[pos])) != NULL)
      s[pos] = (special ? *newsyms : newsyms[sptr-oldsyms]);
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
