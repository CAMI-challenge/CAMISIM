/* parse-t1
 * Test of easel's file parsing
 * 
 * 
 * CVS $Id$
 * SRE, Sat Jul 10 08:59:48 2004 [St. Louis]
 */

#include <stdio.h>
#include <stdlib.h>
#include <easel/esl_core.h>
#include <easel/esl_parse.h>


int
main(int argc, char **argv)
{
  FILE           *fp;
  ESL_FILEPARSER *efp;
  char           *testfile = "parse-t1.dat";
  int             ntok_expected;
  int             lastval_expected;
  int             val;
  int             ntok;
  char           *tok;

  if ((fp = fopen(testfile, "r")) == NULL)
    esl_Die("couldn't open %s", testfile);

  esl_fileparse_create(fp, &efp);
  esl_fileparse_set_commentchar(efp, '#');

  esl_fileparse_token(efp, &tok, NULL);
  ntok_expected = atoi(tok);

  esl_fileparse_token(efp, &tok, NULL);
  lastval_expected = atoi(tok);
  
  ntok = 0;
  val  = -1;
  while (esl_fileparse_token(efp, &tok, NULL) == ESL_OK)
    {
      val = atoi(tok);
      ntok++;
    }

  if (val != lastval_expected)
    esl_Die("expected to see %d, but read %d", lastval_expected, val);
  if (ntok != ntok_expected)
    esl_Die("expected %d tokens, but saw %d", ntok_expected, ntok);
  
  esl_fileparse_free(efp);
  return ESL_OK;
}
