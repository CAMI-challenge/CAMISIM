/* Fetch a sequence (or part of one) from a sequence flatfile.
 * 
 * From squid's sfetch and ffetch
 * SRE, Mon Mar 31 16:12:50 2008 [Janelia] 
 * SVN $Id: esl-sfetch.c 509 2010-02-07 22:56:55Z eddys $
 */
#include "esl_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_fileparser.h"
#include "esl_keyhash.h"
#include "esl_regexp.h"
#include "esl_ssi.h"
#include "esl_sq.h"
#include "esl_sqio.h"

static char banner[] = "retrieve sequence(s) from a file";
static char usage1[] = "[options] <sqfile> <name>        (one seq named <name>)";
static char usage2[] = "[options] -f <sqfile> <namefile> (all seqs in <namefile>)";
static char usage3[] = "[options] --index <sqfile>       (index <sqfile>)";

static void
cmdline_failure(char *argv0, char *format, ...) 
{
  va_list argp;
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage1);
  esl_usage(stdout, argv0, usage2);
  esl_usage(stdout, argv0, usage3);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go) 
{
  esl_banner(stdout, argv0, banner);
  esl_usage (stdout, argv0, usage1);
  esl_usage (stdout, argv0, usage2);
  esl_usage (stdout, argv0, usage3);
  puts("\n where general options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  puts("\n Options for retrieving subsequences:");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 80);
  puts("\n  On command line, subseq coords are separated by any nonnumeric, nonspace character(s).");
  puts("  for example, -c 23..100 or -c 23/100 or -c 23-100 all work.\n");
  puts("  Additionally, to retrieve a suffix to the end, omit the end coord; -c 23: will work.");
  puts("  By default, the subseq will be named <source name>/<from>-<to>. To assign a name of");
  puts("  your choice, use -n <newname>.\n");
  puts("  In retrieving subsequences listed in a file (-C -f, or just -Cf), each line of the file");
  puts("  is in GDF format: <newname> <from> <to> <source seqname>, space/tab delimited.\n");
  puts("  When <start> coordinate is greater than <end>, for DNA or RNA, the reverse complement is");
  puts("  retrieved; in protein sequence, this is an error. The -r option is another way to revcomp.");
  puts("\n other options:");
  esl_opt_DisplayHelp(stdout, go, 3, 2, 80);
  puts("\n options for retreiving subsequences from cmsearch tab file (require -C and -f):");
  esl_opt_DisplayHelp(stdout, go, 4, 2, 80);
  exit(0);
}

static ESL_OPTIONS options[] = {
  /* name          type           default env   range togs  reqs               incomp                help                                                 docgroup */
  { "-h",          eslARG_NONE,   FALSE,  NULL, NULL, NULL, NULL,              NULL,                 "help; show brief info on version and usage",        1 },
  { "-o",          eslARG_OUTFILE,FALSE,  NULL, NULL, NULL, NULL,              "-O,--index",         "output sequences to file <f> instead of stdout",    1 },
  { "-O",          eslARG_NONE,   FALSE,  NULL, NULL, NULL, NULL,              "-o,-f,--index",      "output sequence to file named <key>",               1 },
  { "-n",          eslARG_STRING, FALSE,  NULL, NULL, NULL, NULL,              "-f,--index",         "rename the sequence <s>",                           1 },
  { "-r",          eslARG_NONE,   FALSE,  NULL, NULL, NULL, NULL,              "--index",            "reverse complement the seq(s)",                     1 },

  { "-c",          eslARG_STRING, FALSE,  NULL, NULL, NULL, NULL,              "-f,--index",         "retrieve subsequence coords <from>..<to>",          2 },
  { "-C",          eslARG_NONE,   FALSE,  NULL, NULL, NULL, "-f",              "--index",            "<namefile> in <f> contains subseq coords too",      2 },

  { "--informat",  eslARG_STRING, FALSE,  NULL, NULL, NULL, NULL,              NULL,                 "specify that input file is in format <s>",          3 },
  { "--tabfile",   eslARG_NONE,   FALSE,  NULL, NULL, NULL, "-C,-f",           "--index",            "<namefile> in <f> is Infernal cmsearch tab file",   4 },
  { "--shortname", eslARG_NONE,   FALSE,  NULL, NULL, NULL, "-C,-f,--tabfile", "--index",            "w/--tabfile, do not add bit score, E value, GC to name", 4 },
  { "--Tmin",      eslARG_REAL,   NULL,   NULL, NULL, NULL, "-C,-f,--tabfile", "--index",            "w/--tabfile, only fetch sequences with bit scores above <x>", 4},
  { "--Emax",      eslARG_REAL,   NULL,   NULL, "x>0.",NULL,"-C,-f,--tabfile", "--index",            "w/--tabfile, only fetch sequences with E-values below <x>", 4},

  /* undocumented as options, because they're documented as alternative invocations: */
  { "-f",          eslARG_NONE,  FALSE,   NULL, NULL, NULL, NULL,              "--index",           "second cmdline arg is a file of names to retrieve", 99 },
  { "--index",     eslARG_NONE,  FALSE,   NULL, NULL, NULL, NULL,               NULL,               "index <sqfile>, creating <sqfile>.ssi",             99 },

 { 0,0,0,0,0,0,0,0,0,0 },
};

static void create_ssi_index(ESL_GETOPTS *go, ESL_SQFILE *sqfp);
static void multifetch(ESL_GETOPTS *go, FILE *ofp, char *keyfile, ESL_SQFILE *sqfp);
static void onefetch(ESL_GETOPTS *go, FILE *ofp, char *key, ESL_SQFILE *sqfp);
static void multifetch_subseq(ESL_GETOPTS *go, FILE *ofp, char *keyfile, ESL_SQFILE *sqfp);
static void multifetch_subseq_infernal(ESL_GETOPTS *go, FILE *ofp, char *tabfile, ESL_SQFILE *sqfp);
static void onefetch_subseq(ESL_GETOPTS *go, FILE *ofp, ESL_SQFILE *sqfp, char *newname, 
			    char *key, uint32_t given_start, uint32_t given_end);
static int  parse_coord_string(const char *cstring, uint32_t *ret_start, uint32_t *ret_end);
static void infernal_name_subseq(char **ret_name, const char *name, ...);

int
main(int argc, char **argv)
{
  ESL_GETOPTS  *go      = NULL;	                        /* application configuration       */
  char         *seqfile = NULL;	                        /* sequence file name              */
  int           infmt   = eslSQFILE_UNKNOWN;		/* format code for seqfile         */
  ESL_SQFILE   *sqfp    = NULL;                         /* open sequence file              */
  FILE         *ofp     = NULL;	                        /* output stream for sequences     */
  int           status;		                        /* easel return code               */

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h") )                   cmdline_help   (argv[0], go);
  if (esl_opt_ArgNumber(go) < 1)                       cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        

  /* Open the sequence file */
  seqfile = esl_opt_GetArg(go, 1);
  if (esl_opt_GetString(go, "--informat") != NULL) {
    infmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (infmt == eslSQFILE_UNKNOWN) esl_fatal("%s is not a valid input sequence file format for --informat"); 
  }
  status = esl_sqfile_Open(seqfile, infmt, NULL, &sqfp);
  if      (status == eslENOTFOUND) cmdline_failure(argv[0], "Sequence file %s not found.\n",     seqfile);
  else if (status == eslEFORMAT)   cmdline_failure(argv[0], "Format of file %s unrecognized.\n", seqfile);
  else if (status == eslEINVAL)    cmdline_failure(argv[0], "Can't autodetect stdin or .gz.\n");
  else if (status != eslOK)        cmdline_failure(argv[0], "Open failed, code %d.\n", status);

  /* Open the output file, if any */
  if (esl_opt_GetBoolean(go, "-O")) 
    {
      if ((ofp = fopen(esl_opt_GetArg(go, 2), "w")) == NULL)
	cmdline_failure(argv[0], "Failed to open output file %s\n", esl_opt_GetArg(go, 2));
    }
  else if (esl_opt_GetString(go, "-o") != NULL)
    {
      if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL)
	cmdline_failure(argv[0], "Failed to open output file %s\n", esl_opt_GetString(go, "-o"));
    }
  else ofp = stdout;

  /* Indexing  mode */
  if (esl_opt_GetBoolean(go, "--index")) 
    {
      if (esl_opt_ArgNumber(go) != 1) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      if (sqfp->data.ascii.do_gzip)  cmdline_failure(argv[0], "Can't index a .gz compressed file");
      if (sqfp->data.ascii.do_stdin) cmdline_failure(argv[0], "Can't index a standard input pipe");

      create_ssi_index(go, sqfp);
    }

  /* List retrieval mode */
  else if (esl_opt_GetBoolean(go, "-f"))
    {
      if (esl_opt_ArgNumber(go) != 2) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        

      /* Open the SSI index for retrieval */
      if (! sqfp->data.ascii.do_gzip && ! sqfp->data.ascii.do_stdin &&  ! esl_sqio_IsAlignment(sqfp->format)) 
	{
	  status = esl_sqfile_OpenSSI(sqfp, NULL);
	  if      (status == eslEFORMAT)   cmdline_failure(argv[0], "SSI index is in incorrect format\n");
	  else if (status == eslERANGE)    cmdline_failure(argv[0], "SSI index is in 64-bit format and we can't read it\n");
	  else if (status != eslOK)        cmdline_failure(argv[0], "Failed to open SSI index\n");
	}

      if (esl_opt_GetBoolean(go, "-C")) { 
	if(esl_opt_GetBoolean(go, "--tabfile")) multifetch_subseq_infernal(go, ofp, esl_opt_GetArg(go, 2), sqfp);
	else                            	multifetch_subseq         (go, ofp, esl_opt_GetArg(go, 2), sqfp);
      }
      else                        	        multifetch       (go, ofp, esl_opt_GetArg(go, 2), sqfp);
    }

  /* Single sequence retrieval mode */
  else 
    {
      if (esl_opt_ArgNumber(go) != 2) cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");        
      char *key     = esl_opt_GetArg(go, 2);
      char *cstring = esl_opt_GetString(go, "-c");
      char *newname = esl_opt_GetString(go, "-n");

      /* Open the SSI index for retrieval */
      if (! sqfp->data.ascii.do_gzip && ! sqfp->data.ascii.do_stdin &&  ! esl_sqio_IsAlignment(sqfp->format)) 
	{
	  status = esl_sqfile_OpenSSI(sqfp, NULL);
	  if      (status == eslEFORMAT)   cmdline_failure(argv[0], "SSI index is in incorrect format\n");
	  else if (status == eslERANGE)    cmdline_failure(argv[0], "SSI index is in 64-bit format and we can't read it\n");
	  else if (status != eslOK)        cmdline_failure(argv[0], "Failed to open SSI index\n");
	}

      /* -c: subsequence retrieval; else full sequence retrieval */
      if (cstring != NULL)
	{
	  uint32_t start, end;

	  parse_coord_string(cstring, &start, &end);
	  onefetch_subseq(go, ofp, sqfp, newname, key, start, end);
	  if (ofp != stdout) printf("\n\nRetrieved subsequence %s/%d-%d.\n",  key, start, end);
	}
      else 
	{
	  onefetch(go, ofp, esl_opt_GetArg(go, 2), sqfp);
	  if (ofp != stdout) printf("\n\nRetrieved sequence %s.\n",  esl_opt_GetArg(go, 2));
	}
    }

  esl_sqfile_Close(sqfp);
  esl_getopts_Destroy(go);
  return 0;
}


/* Create an SSI index file for open sequence file <sqfp>.
 * Both name and accession of sequences are stored as keys.
 */
static void
create_ssi_index(ESL_GETOPTS *go, ESL_SQFILE *sqfp)
{
  ESL_NEWSSI *ns      = NULL;
  ESL_SQ     *sq      = esl_sq_Create();
  int         nseq    = 0;
  char       *ssifile = NULL;
  uint16_t    fh;
  int         status;

  esl_strdup(sqfp->filename, -1, &ssifile);
  esl_strcat(&ssifile, -1, ".ssi", 4);
  status = esl_newssi_Open(ssifile, FALSE, &ns);
  if      (status == eslENOTFOUND)   esl_fatal("failed to open SSI index %s", ssifile);
  else if (status == eslEOVERWRITE)  esl_fatal("SSI index %s already exists; delete or rename it", ssifile);
  else if (status != eslOK)          esl_fatal("failed to create a new SSI index");

  if (esl_newssi_AddFile(ns, sqfp->filename, sqfp->format, &fh) != eslOK)
    esl_fatal("Failed to add sequence file %s to new SSI index\n", sqfp->filename);

  printf("Creating SSI index for %s...    ", sqfp->filename); 
  fflush(stdout);
  
  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      nseq++;
      if (sq->name == NULL) esl_fatal("Every sequence must have a name to be indexed. Failed to find name of seq #%d\n", nseq);

      if (esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->L) != eslOK)
	esl_fatal("Failed to add key %s to SSI index", sq->name);

      if (sq->acc[0] != '\0') {
	if (esl_newssi_AddAlias(ns, sq->acc, sq->name) != eslOK)
	  esl_fatal("Failed to add secondary key %s to SSI index", sq->acc);
      }
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					   sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					    status, sqfp->filename);

  /* Determine if the file was suitable for fast subseq lookup. */
  if (sqfp->data.ascii.bpl > 0 && sqfp->data.ascii.rpl > 0) {
    if ((status = esl_newssi_SetSubseq(ns, fh, sqfp->data.ascii.bpl, sqfp->data.ascii.rpl)) != eslOK) 
      esl_fatal("Failed to set %s for fast subseq lookup.");
  }

  /* Save the SSI file to disk */
  if (esl_newssi_Write(ns) != eslOK)  esl_fatal("Failed to write keys to ssi file %s\n", ssifile);

  /* Done - output and exit. */
  printf("done.\n");
  if (ns->nsecondary > 0) 
    printf("Indexed %d sequences (%ld names and %ld accessions).\n", nseq, (long) ns->nprimary, (long) ns->nsecondary);
  else 
    printf("Indexed %d sequences (%ld names).\n", nseq, (long) ns->nprimary);
  printf("SSI index written to file %s\n", ssifile);

  free(ssifile);
  esl_sq_Destroy(sq);
  esl_newssi_Close(ns);
  return;
}

/* multifetch:
 * given a file containing lines with one name or key per line;
 * parse the file line-by-line;
 * if we have an SSI index available, retrieve the seqs by key
 * as we see each line;
 * else, without an SSI index, store the keys in a hash, then
 * read the entire seq file in a single pass, outputting seqs
 * that are in our keylist. 
 * 
 * Note that with an SSI index, you get the seqs in the order they
 * appear in the <keyfile>, but without an SSI index, you get seqs in
 * the order they occur in the seq file.
 */
static void
multifetch(ESL_GETOPTS *go, FILE *ofp, char *keyfile, ESL_SQFILE *sqfp)
{
  ESL_KEYHASH    *keys   = esl_keyhash_Create();
  ESL_FILEPARSER *efp    = NULL;
  int             nseq   = 0;
  int             nkeys  = 0;
  char           *key;
  int             keylen;
  int             keyidx;
  int             status;

  
  if (esl_fileparser_Open(keyfile, NULL, &efp) != eslOK)  esl_fatal("Failed to open key file %s\n", keyfile);
  esl_fileparser_SetCommentChar(efp, '#');

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &key, &keylen) != eslOK)
	esl_fatal("Failed to read seq name on line %d of file %s\n", efp->linenumber, keyfile);
      
      status = esl_key_Store(keys, key, &keyidx);
      if (status == eslEDUP) esl_fatal("seq key %s occurs more than once in file %s\n", key, keyfile);
	
      /* if we have an SSI index, just fetch them as we go. */
      if (sqfp->data.ascii.ssi != NULL) { onefetch(go, ofp, key, sqfp);  nseq++; }
      nkeys++;
    }

  /* If we don't have an SSI index, we haven't fetched anything yet; do it now. */
  if (sqfp->data.ascii.ssi == NULL) 
    {
      ESL_SQ *sq     = esl_sq_Create();

      while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
	{
	  if ( (sq->name[0] != '\0' && esl_key_Lookup(keys, sq->name, NULL) == eslOK) ||
	       (sq->acc[0]  != '\0' && esl_key_Lookup(keys, sq->acc,  NULL) == eslOK))
	    {
	      if (esl_opt_GetBoolean(go, "-r") )
		if (esl_sq_ReverseComplement(sq) != eslOK) 
		  esl_fatal("Failed to reverse complement %s\n", sq->name);
	      esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE);
	      nseq++;
	    }
	  esl_sq_Reuse(sq);
	}
      if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					       sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s",
					       status, sqfp->filename);
      esl_sq_Destroy(sq);
    }
  
  if (nkeys != nseq) esl_fatal("Tried to retrieve %d keys, but only retrieved %d sequences\n", nkeys, nseq);

  if (ofp != stdout) printf("\nRetrieved %d sequences.\n", nseq);

  esl_keyhash_Destroy(keys);
  esl_fileparser_Close(efp);
  return;
}
  


/* onefetch():
 * Given one <key> (a seq name or accession), retrieve the corresponding sequence.
 * In SSI mode, we can do this quickly by positioning the file, then regurgitating
 * every line until the end-of-record marker; we don't even have to parse.
 * Without an SSI index, we have to parse the file sequentially 'til we find
 * the one we're after.
 */
static void
onefetch(ESL_GETOPTS *go, FILE *ofp, char *key, ESL_SQFILE *sqfp)
{
  ESL_SQ  *sq            = esl_sq_Create();
  int      do_revcomp    = esl_opt_GetBoolean(go, "-r");
  char    *newname       = esl_opt_GetString(go, "-n");
  int      status;

  /* Try to position the file at the desired sequence with SSI. */
  if (sqfp->data.ascii.ssi != NULL)	
    {
      status = esl_sqfile_PositionByKey(sqfp, key);
      if      (status == eslENOTFOUND) esl_fatal("seq %s not found in SSI index for file %s\n", key, sqfp->filename);
      else if (status == eslEFORMAT)   esl_fatal("Failed to parse SSI index for %s\n", sqfp->filename);
      else if (status != eslOK)        esl_fatal("Failed to look up location of seq %s in SSI index of file %s\n", key, sqfp->filename);

      status = esl_sqio_Read(sqfp, sq);
      if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
					       sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
      else if (status == eslEOF)     esl_fatal("Unexpected EOF reading sequence file %s",
					       status, sqfp->filename);
      else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s",
					       status, sqfp->filename);

      if (strcmp(key, sq->name) != 0 && strcmp(key, sq->acc) != 0) 
	esl_fatal("whoa, internal error; found the wrong sequence %s, not %s", sq->name, key);
    }  
  else 
    { /* Else, we have to read the whole damn file sequentially until we find the seq */
      while ((status = esl_sqio_Read(sqfp, sq)) != eslEOF) {
	if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n",
						 sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
	else if (status != eslOK)      esl_fatal("Unexpected error %d reading sequence file %s",
						 status, sqfp->filename);

	if (strcmp(key, sq->name) == 0 || strcmp(key, sq->acc) == 0) break;
	esl_sq_Reuse(sq);
      }
      if (status == eslEOF) esl_fatal("Failed to find sequence %s in file %s\n", key, sqfp->filename);

    }

  if (do_revcomp == FALSE && newname == NULL && ! esl_sqio_IsAlignment(sqfp->format)) 
    { /* If we're not manipulating the sequence in any way, and it's not from an alignment file, we can Echo() it. */
      if (esl_sqio_Echo(sqfp, sq, ofp) != eslOK) esl_fatal("Echo failed: %s\n", esl_sqfile_GetErrorBuf(sqfp));
    }
  else
    { /* Otherwise we Write() the parsed version. */
      if (do_revcomp && esl_sq_ReverseComplement(sq) != eslOK) esl_fatal("Failed to reverse complement %s; is it a protein?\n", sq->name);
      if (newname != NULL) esl_sq_SetName(sq, newname);
      esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE);
    }

  esl_sq_Destroy(sq);
}

static void
multifetch_subseq(ESL_GETOPTS *go, FILE *ofp, char *gdffile, ESL_SQFILE *sqfp)
{
  ESL_FILEPARSER *efp    = NULL;
  char           *newname;
  char           *s;
  int             n1, n2;
  int             start, end;
  char           *source;
 
  if (esl_fileparser_Open(gdffile, NULL, &efp) != eslOK)  esl_fatal("Failed to open key file %s\n", gdffile);
  esl_fileparser_SetCommentChar(efp, '#');

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      if (esl_fileparser_GetTokenOnLine(efp, &newname, &n1) != eslOK)
	esl_fatal("Failed to read subseq name on line %d of file %s\n", efp->linenumber, gdffile);
      if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	esl_fatal("Failed to read start coord on line %d of file %s\n", efp->linenumber, gdffile);
      start = atoi(s);
      if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	esl_fatal("Failed to read end coord on line %d of file %s\n", efp->linenumber, gdffile);
      end   = atoi(s);
      if (esl_fileparser_GetTokenOnLine(efp, &source, &n2) != eslOK)
	esl_fatal("Failed to read source seq name on line %d of file %s\n", efp->linenumber, gdffile);

      onefetch_subseq(go, ofp, sqfp, newname, source, start, end);
    }
  esl_fileparser_Close(efp);
}

static void
onefetch_subseq(ESL_GETOPTS *go, FILE *ofp, ESL_SQFILE *sqfp, char *newname, char *key, uint32_t given_start, uint32_t given_end)
{
  int    start, end;
  int    do_revcomp;
  ESL_SQ *sq = esl_sq_Create();

  if (sqfp->data.ascii.ssi == NULL) esl_fatal("no ssi index");

  /* reverse complement indicated by coords. */
  /* -c 52: would be 52,0, so watch out for given_end = 0 case */
  if (given_end != 0 && given_start > given_end)
    { start = given_end;   end = given_start; do_revcomp = TRUE;  }
  else
    { start = given_start; end = given_end;   do_revcomp = FALSE; }

  if (esl_sqio_FetchSubseq(sqfp, key, start, end, sq) != eslOK) esl_fatal(esl_sqfile_GetErrorBuf(sqfp));

  if      (newname != NULL) esl_sq_SetName(sq, newname);
  else                      esl_sq_FormatName(sq, "%s/%d-%d", key, given_start, (given_end == 0) ? sq->L : given_end);

  /* Two ways we might have been asked to revcomp: by coord, or by -r option */
  /* (If both happen, they'll cancel each other out) */
  if (do_revcomp) 
    if (esl_sq_ReverseComplement(sq) != eslOK) esl_fatal("Failed to reverse complement %s; is it a protein?\n", sq->name);
  if (esl_opt_GetBoolean(go, "-r"))
    if (esl_sq_ReverseComplement(sq) != eslOK) esl_fatal("Failed to reverse complement %s; is it a protein?\n", sq->name);

  esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE);
  esl_sq_Destroy(sq);
}


static int
parse_coord_string(const char *cstring, uint32_t *ret_start, uint32_t *ret_end)
{
  ESL_REGEXP *re = esl_regexp_Create();
  char        tok1[32];
  char        tok2[32];

  if (esl_regexp_Match(re, "^(\\d+)\\D+(\\d*)$", cstring) != eslOK) esl_fatal("-c takes arg of subseq coords <from>..<to>; %s not recognized", cstring);
  if (esl_regexp_SubmatchCopy(re, 1, tok1, 32)            != eslOK) esl_fatal("Failed to find <from> coord in %s", cstring);
  if (esl_regexp_SubmatchCopy(re, 2, tok2, 32)            != eslOK) esl_fatal("Failed to find <to> coord in %s",   cstring);
  
  *ret_start = atol(tok1);
  *ret_end   = (tok2[0] == '\0') ? 0 : atol(tok2);
  
  esl_regexp_Destroy(re);
  return eslOK;
}

static void
multifetch_subseq_infernal(ESL_GETOPTS *go, FILE *ofp, char *tabfile, ESL_SQFILE *sqfp)
{
  int status;
  ESL_FILEPARSER *efp    = NULL;
  char           *s;
  int             start, end, n;
  char           *tname = NULL;
  char           *mname = NULL;
  char           *newname = NULL;
  double          bit, E;
  int             gc;
  int             has_E = FALSE;
  char           *tok1 = NULL;
  char           *tok2 = NULL;
  char           *tok3 = NULL;
  char           *tok4 = NULL;
  char           *tok5 = NULL;
  char           *tok6 = NULL;
  char           *tok7 = NULL;
  char           *tok8 = NULL;
  char           *tok9 = NULL;
  float           Tmin = -eslINFINITY;
  float           Emax = eslINFINITY;
  int             Emax_is_on = FALSE;
  int             has_model_name = FALSE;

  if (esl_fileparser_Open(tabfile, NULL, &efp) != eslOK)  esl_fatal("Failed to open Infernal tab file %s\n", tabfile);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if (esl_opt_IsOn(go, "--Tmin")) { Tmin = esl_opt_GetReal(go, "--Tmin"); }
  if (esl_opt_IsOn(go, "--Emax")) { Emax = esl_opt_GetReal(go, "--Emax"); Emax_is_on = TRUE; }

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      /* infernal 1.01 introduced a new tab file format.
       * This function recognizes both formats (version 1.0 and 1.01+)
       *
       * Example line from infernal 1.0:
       *se-1-2                1          24      1     24     18.93      0.234   54
       *<target (t) name>     <t start>  <t end> <q s> <q e>  <bit sc>   <E-val> <GC content>
       *
       * Example line from infernal 1.01+:
       *SEED          se-1-2                1          24      1     24     18.93      0.234   54
       *<model name>  <target (t) name>     <t start>  <t end> <q s> <q e>  <bit sc>   <E-val> <GC content>
       */
      has_E = FALSE;
      has_model_name = FALSE;

      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read first token on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok1) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 2 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok2) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 3 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok3) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 4 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok4) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 5 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok5) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 6 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok6) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 7 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok7) != eslOK) goto ERROR;
      if (esl_fileparser_GetTokenOnLine(efp, &s, &n)  != eslOK) esl_fatal("Failed to read token 8 on line %d of Infernal tab file %s\n", efp->linenumber, tabfile);
      if (esl_strdup(s, n, &tok8) != eslOK) goto ERROR;
      /* if we're in an infernal 1.0 file, we should be at the end of the line, we can tell if this is the case if the next esl_fileparser_GetTokenOnLine() call returns eslEOL. */
      if ((status = esl_fileparser_GetTokenOnLine(efp, &s, &n)) == eslEOL) { 
	/* we're in an infernal 1.0 cmsearch tabfile, parse the line */
	tname = tok1;
	start = atoi(tok2);
	end   = atoi(tok3);
	/* tok4, query start is not used */
	/* tok5, query end   is not used */
	bit   = atof(tok6);
	if(strcmp(tok7, "-") != 0) { has_E = TRUE; E = atof(tok7); } 
	else                       { has_E = FALSE; }
	gc = atoi(tok8);
      }
      else if (status == eslOK) { /* 9th token was read */
	/* we're in an infernal 1.01+ cmsearch tabfile, we have the model name */
	if (esl_strdup(s, n, &tok9) != eslOK) goto ERROR;
	mname = tok1;
	has_model_name = TRUE;
	tname = tok2;
	start = atoi(tok3);
	end   = atoi(tok4);
	/* tok5, query start is not used */
	/* tok6, query end   is not used */
	bit   = atof(tok7);
	if(strcmp(tok8, "-") != 0) { has_E = TRUE; E = atof(tok8); } 
	else                       { has_E = FALSE; }
	gc = atoi(tok9);
	/* ensure next token is end of line */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) == eslOK) { esl_fatal("Read more than 9 tokens on line %d of Infernal tab file %s\n", efp->linenumber, tabfile); }
      }
      else { /* this should never happen, esl_fileparser_GetTokenOnLine returns only either eslOK or eslEOL */
	esl_fatal("Unexpected error status %d, reading line %d of Infernal tab file %s\n", status, efp->linenumber, tabfile); 
      }

      if(! esl_opt_GetBoolean(go, "--shortname")) { 
	if(has_model_name) { 
	  if(has_E) infernal_name_subseq(&newname, "%s/%d-%d/%s/B%.2f/E%.2g/GC%d", tname, start, end, mname, bit, E, gc);
	  else      infernal_name_subseq(&newname, "%s/%d-%d/%s/B%.2f/GC%d",       tname, start, end, mname, bit, gc);
	}
	else { /* we don't have model name (infernal 1.0 file) */
	  if(has_E) infernal_name_subseq(&newname, "%s/%d-%d/B%.2f/E%.2g/GC%d", tname, start, end, bit, E, gc);
	  else      infernal_name_subseq(&newname, "%s/%d-%d/B%.2f/GC%d",       tname, start, end, bit, gc);
	}
      }
      /* else use default esl-sfetch name style, pass newname as NULL to onefetch_subseq (note: --modelname and --shortname are incompatible)*/

      /* make sure that if Emax is enabled, we've read an E-value */
      if(Emax_is_on && (!has_E)) { esl_fatal("--Emax enabled, but no E-value read from line %d of Infernal tab file %s\n", efp->linenumber, tabfile); }

      /* check if the hit meets our minimum score criteria (these are -inf bit score and inf E-value unless --Tmin and/or --Emax enabled) */
      if(((has_E && E <= Emax) || (! has_E)) && (bit >= Tmin)) { 
	onefetch_subseq(go, ofp, sqfp, newname, tname, start, end);
      }

      if(newname != NULL) { free(newname); newname = NULL; }
      if(tok1    != NULL) { free (tok1);   tok1    = NULL; }
      if(tok2    != NULL) { free (tok2);   tok2    = NULL; }
      if(tok3    != NULL) { free (tok3);   tok3    = NULL; }
      if(tok4    != NULL) { free (tok4);   tok4    = NULL; }
      if(tok5    != NULL) { free (tok5);   tok5    = NULL; }
      if(tok6    != NULL) { free (tok6);   tok6    = NULL; }
      if(tok7    != NULL) { free (tok7);   tok7    = NULL; }
      if(tok8    != NULL) { free (tok8);   tok8    = NULL; }
      if(tok9    != NULL) { free (tok9);   tok9    = NULL; }
    }
  esl_fileparser_Close(efp);
  return;
  
 ERROR:
  if(newname != NULL) { free(newname); newname = NULL; }
  if(tok1    != NULL) { free (tok1);   tok1    = NULL; }
  if(tok2    != NULL) { free (tok2);   tok2    = NULL; }
  if(tok3    != NULL) { free (tok3);   tok3    = NULL; }
  if(tok4    != NULL) { free (tok4);   tok4    = NULL; }
  if(tok5    != NULL) { free (tok5);   tok5    = NULL; }
  if(tok6    != NULL) { free (tok6);   tok6    = NULL; }
  if(tok7    != NULL) { free (tok7);   tok7    = NULL; }
  if(tok8    != NULL) { free (tok8);   tok8    = NULL; }
  if(tok9    != NULL) { free (tok9);   tok9    = NULL; }
  if(efp != NULL) esl_fileparser_Close(efp);
  esl_fatal("Error, status code %d, probably out of memory.", status);
  return;
}


static void
infernal_name_subseq(char **ret_name, const char *name, ...)
{
  va_list argp;
  va_list argp2;
  int   n;
  void *tmp;
  int   status;
  char *newname;
  int nalloc = 1;
  ESL_ALLOC(newname, sizeof(char) * nalloc);
  if (name == NULL) { newname[0] = '\0'; *ret_name = newname; return; }

  va_start(argp, name);
  va_copy(argp2, argp);
  if ((n = vsnprintf(newname, nalloc, name, argp)) >= nalloc)
    {
      ESL_RALLOC(newname, tmp, sizeof(char) * (n+1)); 
      nalloc = n+1;
      vsnprintf(newname, nalloc, name, argp2);
    }
  va_end(argp);
  va_end(argp2);
  *ret_name = newname;
  return;

 ERROR:
  esl_fatal("Memory allocation error in infernal_name_subseq().\n");
  return;
}
