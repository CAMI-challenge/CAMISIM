/* A simple token-based file parsing system.
 * 
 * SRE, Tue Jul 13 14:40:35 2004 [St. Louis]
 * SVN $Id: esl_fileparser.h 309 2008-12-17 13:49:36Z eddys $
 */

#ifndef ESL_FILEPARSER_INCLUDED
#define ESL_FILEPARSER_INCLUDED

#include <stdio.h>
#include "easel.h"

typedef struct {
  FILE *fp;			/* open file pointer, for reading                  */
  char *buf;			/* current line; will be modified by esl_strtok(). */
  int   buflen;			/* current allocated length of buf                 */
  char *s;			/* used by esl_strtok(); current position in buf.  */
  char  commentchar;		/* often '#'                                       */

  char *tok;			/* _NextLine() may remember a token...             */
  int   toklen;			/* ... and its length                              */
  char  tokchar;		/* ... and char that got overwritten by \0, if any */

  char *filename;		/* name of opened file; or NULL (if just a stream) */
  int   linenumber;		/* what line is loaded into buf; 1..nlines         */
  char  errbuf[eslERRBUFSIZE];  /* for holding error diagnostics                   */
} ESL_FILEPARSER;

extern int  esl_fileparser_Open(const char *filename, const char *envvar, ESL_FILEPARSER **ret_efp);
extern ESL_FILEPARSER *esl_fileparser_Create(FILE *fp);
extern int  esl_fileparser_SetCommentChar  (ESL_FILEPARSER *efp, char c);
extern int  esl_fileparser_GetToken        (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
extern int  esl_fileparser_NextLine        (ESL_FILEPARSER *efp);
extern int  esl_fileparser_NextLinePeeked  (ESL_FILEPARSER *efp, char *prefix, int plen);
extern int  esl_fileparser_GetTokenOnLine  (ESL_FILEPARSER *efp, char **opt_tok, int *opt_toklen);
extern int  esl_fileparser_GetRemainingLine(ESL_FILEPARSER *efp, char **ret_s);
extern void esl_fileparser_Destroy         (ESL_FILEPARSER *efp);
extern void esl_fileparser_Close           (ESL_FILEPARSER *efp);

#endif /*ESL_FILEPARSER_INCLUDED */
