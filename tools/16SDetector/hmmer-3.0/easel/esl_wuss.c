/* wuss.c
 * RNA secondary structure markup in WUSS notation.
 * 
 * SRE, Tue Feb 15 08:43:23 2005
 * SVN $Id: esl_wuss.c 241 2008-04-01 19:01:52Z eddys $
 * xref squid wuss.c.
 */

#include <esl_config.h>

#include <string.h>
#include <ctype.h>

#include "easel.h"
#include "esl_stack.h"
#include "esl_wuss.h"

/* Function:  esl_wuss2ct()
 * Incept:    SRE, Tue Feb 15 08:44:54 2005 [St. Louis]
 *
 * Purpose:   Given a secondary structure string <ss>, <0..len-1>,
 *            in WUSS notation, convert it to a CT array, <1..len>,
 *            in <ct>. Caller provides a <ct> allocated for at least 
 *            <len+1> ints. <ct[i]> is the position that residue i
 *            base pairs to, or 0 if i is unpaired. <ct[0]> is undefined
 *            (but if you care: it is set to 0).
 *            
 *            WUSS notation is interpreted loosely here, as input
 *            WUSS.  Any matching bracket pair or upper/lower case
 *            alphabetic pair is interpreted as a base pair; any other
 *            WUSS annotation is interpreted as unpaired.
 *            
 * Returns:   <eslOK> on success. Returns <eslESYNTAX> if the WUSS
 *            string isn't valid.
 *            
 * Throws:    <eslEMEM> on allocation failure.           
 */
int 
esl_wuss2ct(char *ss, int len, int *ct)
{
  ESL_STACK *pda[27];     /* 1 secondary structure + up to 26 levels of pk's */
  int        i;
  int        pos, pair;
  int        status;      /* success or failure return status */

 /* Initialization: always initialize the main pda (0);
  * we'll init the pk pda's on demand.
  */
  if ((pda[0] = esl_stack_ICreate()) == NULL) goto FINISH;
  for (i = 1; i <= 26; i++) pda[i] = NULL;

  for (pos = 0; pos <= len; pos++) ct[pos] = 0;

  for (pos = 1; pos <= len; pos++)
    {
      if (!isprint((int) ss[pos-1]))  /* armor against garbage */
	{ status = eslESYNTAX; goto FINISH; }

      /* left side of a pair: push position onto stack 0 (pos = 1..L) */
      else if (ss[pos-1] == '<' ||
	       ss[pos-1] == '(' ||
	       ss[pos-1] == '[' ||
	       ss[pos-1] == '{')
	{
	  if ((status = esl_stack_IPush(pda[0], pos)) != eslOK) goto FINISH;
	}
      
      /* right side of a pair; resolve pair; check for agreement */
      else if (ss[pos-1] == '>' || 
	       ss[pos-1] == ')' ||
	       ss[pos-1] == ']' ||
	       ss[pos-1] == '}')
        {
          if (esl_stack_IPop(pda[0], &pair) == eslEOD)
            { status = eslESYNTAX; goto FINISH;} /* no closing bracket */
          else if ((ss[pair-1] == '<' && ss[pos-1] != '>') ||
		   (ss[pair-1] == '(' && ss[pos-1] != ')') ||
		   (ss[pair-1] == '[' && ss[pos-1] != ']') ||
		   (ss[pair-1] == '{' && ss[pos-1] != '}'))
	    { status = eslESYNTAX; goto FINISH; }  /* brackets don't match */
	  else
	    {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
                                /* same stuff for pseudoknots */
      else if (isupper((int) ss[pos-1])) 
	{
	  /* Create the PK stacks on demand.
	   */
	  i = ss[pos-1] - 'A' + 1;
	  if (pda[i] == NULL) 
	    if ((pda[i] = esl_stack_ICreate()) == NULL) 
	      { status = eslEMEM; goto FINISH; }

	  if ((status = esl_stack_IPush(pda[i], pos)) != eslOK) goto FINISH;
	}
      else if (islower((int) ss[pos-1])) 
	{
	  i = ss[pos-1] - 'a' + 1;
	  if (pda[i] == NULL || 
	      esl_stack_IPop(pda[i], &pair) == eslEOD)
            { status = eslESYNTAX; goto FINISH;}
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
	}
      else if (strchr(":,_-.~", ss[pos-1]) == NULL)
	{ status = eslESYNTAX; goto FINISH; } /* bogus character */
    }
  status = eslOK;

 FINISH:
  for (i = 0; i <= 26; i++)
    if (pda[i] != NULL) 
      { /* nothing should be left on stacks */
	if (esl_stack_ObjectCount(pda[i]) != 0)
	  status = eslESYNTAX;
	esl_stack_Destroy(pda[i]);
      }
  return status;
}


/* Function:  esl_ct2wuss()
 * Incept:    SRE, Wed Feb 16 11:22:53 2005 [St. Louis]
 *
 * Purpose:   Convert a CT array <ct> for <n> residues (1..n) to a WUSS
 *            format string <ss>. <ss> must be allocated for at least
 *            n+1 chars (+1 for the terminal NUL). 
 *
 *            Currently limited to nonpseudoknotted structures. Attempting
 *            to convert a pseudoknot-containing <ct> will return an
 *            <eslEINVAL> error.
 *
 * Returns:   <eslOK> on success.
 *            <eslEINVAL> if <ct> contains a pseudoknot.
 * 
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal failure.
 */
int
esl_ct2wuss(int *ct, int n, char *ss)
{
  ESL_STACK *pda    = NULL;	/* main stack  */
  ESL_STACK *aux    = NULL;	/* aux storage */
  int        status = eslEMEM;	/* exit status 'til proven otherwise */
  int        i,j;
  int        nfaces;
  int        minface;

  ss[0] = '\0';	/* in case we abort, and caller does something dumb w/ ss */  

  if ((pda = esl_stack_ICreate()) == NULL) goto FINISH;
  if ((aux = esl_stack_ICreate()) == NULL) goto FINISH;

  for (j = 1; j <= n; j++)
    {
      if (ct[j] == 0)	/* unpaired: push j. */
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else if (ct[j] > j) /* left side of a bp: push j. */
	{
	  if (esl_stack_IPush(pda, j) != eslOK) goto FINISH;
	}
      else   /* right side of a bp; main routine: resolve a subseq */
	{
	  /* Pop back until we find the left partner of i;
           * store SS residues in aux;
           * keep track of #faces and the maximum face depth.
	   */
	  nfaces  = 0;
	  minface = -1;
	  while (1) 
	    {
	      if (esl_stack_IPop(pda, &i) != eslOK) goto FINISH;

	      if (i < 0) 		/* a face counter */
		{
		  nfaces++;
		  if (i < minface) minface = i;
		}
	      else if (ct[i] == j) 
		break;		/* we found the i,j pair. */
	      else if (ct[i] == 0) 
		{
		  if (esl_stack_IPush(aux, i) != eslOK) goto FINISH;
		}
	      else /* ct[i]>0, != j: i is paired, but not to j: pseudoknot! */
		{
		  esl_stack_Destroy(pda); esl_stack_Destroy(aux);	 
		  ESL_EXCEPTION(eslEINVAL, "pseudoknots not permitted yet");
		}
	    }
	  
	  /* Now we know i,j pair; and we know how many faces are
	   * above them; and we know the max depth of those faces.
	   * That's enough to label the pair in WUSS notation.
	   * if nfaces == 0, minface is -1; <> a closing bp of a hairpin.
	   * if nfaces == 1, inherit minface, we're continuing a stem.
	   * if nfaces > 1, bump minface in depth; we're closing a bifurc.
	   */
	  if (nfaces > 1 && minface > -4) minface--;
	  switch (minface) {
	  case -1: ss[i-1] = '<'; ss[j-1] = '>'; break;
	  case -2: ss[i-1] = '('; ss[j-1] = ')'; break;
	  case -3: ss[i-1] = '['; ss[j-1] = ']'; break;
	  case -4: ss[i-1] = '{'; ss[j-1] = '}'; break;
	  default:
	    esl_stack_Destroy(pda); esl_stack_Destroy(aux);
	    ESL_EXCEPTION(eslEINCONCEIVABLE, "no such face code");
	  }
	  if (esl_stack_IPush(pda, minface) != eslOK) goto FINISH;

	  /* Now, aux contains all the unpaired residues we need to label,
	   * according to the # of faces "above" them:
	   *  nfaces = 0: hairpin loop
	   *  nfaces = 1: bulge or interior loop
	   *  nfaces > 1: multifurc
	   */
	  while (esl_stack_IPop(aux, &i) == eslOK)
	    {
	      switch (nfaces) {
	      case 0:  ss[i-1] = '_'; break;
	      case 1:  ss[i-1] = '-'; break;
	      default: ss[i-1] = ','; break; /* nfaces > 1 */
	      }
	    }
	  
	} /* finished processing a subseq enclosed by a bp */
    } /* finished loop over j: end position on seq, 1..n*/

  /* Anything that's left in the pda is either a face counter
   * or external single-strand. Face counters are negative; 
   * position indices are positive.
   */
  while (esl_stack_IPop(pda, &i) == eslOK) 
    if (i > 0) ss[i-1] = ':';
  ss[n] = '\0';
  status = eslOK;

 FINISH:
  if (pda != NULL) esl_stack_Destroy(pda);
  if (aux != NULL) esl_stack_Destroy(aux);
  return status;
}



/* Function:  esl_wuss2kh()
 * Incept:    SRE, Tue Feb 15 10:05:35 2005 [St. Louis]
 *
 * Purpose:   Converts a secondary structure string <ss> in 
 *            WUSS notation back to old KHS format in <kh>.
 *            <kh> must be allocated for at least as much
 *            space as <ss>. <kh> may be the same as <ss>,
 *            in which case the conversion is done in-place.
 *
 * Note:      Left bp chars  are converted to >   (left base of base pairs)
 *            Right bp chars are converted to <   (right base of base pairs)
 *            Characters _-,:~ are converted to . (unpaired bases)
 *            Character  .     is untouched       (unpaired)
 *            Everything else is untouched, including any pseudoknot notation.
 * 
 * Returns:   <eslOK> on success.
 */
int
esl_wuss2kh(char *ss, char *kh)
{
  while (*ss != '\0')
    {
      if       (*ss == '<') *kh = '>';
      else if  (*ss == '(') *kh = '>';
      else if  (*ss == '[') *kh = '>';
      else if  (*ss == '{') *kh = '>';
      else if  (*ss == '>') *kh = '<';
      else if  (*ss == ')') *kh = '<';
      else if  (*ss == ']') *kh = '<';
      else if  (*ss == '}') *kh = '<';
      else if  (*ss == '_') *kh = '.';
      else if  (*ss == '-') *kh = '.';
      else if  (*ss == ',') *kh = '.';
      else if  (*ss == ':') *kh = '.';
      else if  (*ss == '~') *kh = '.';
      else *kh = *ss;
      ss++;
      kh++;
    }
  *kh = '\0';
  return eslOK;
}


/* Function:  esl_kh2wuss()
 * Incept:    SRE, Tue Feb 15 10:10:40 2005 [St. Louis]
 *
 * Purpose:   Converts an old format secondary structure string <kh>
 *            to shorthand WUSS format <ss>. <ss> must be allocated at least
 *            as large as <kh>. <ss> can be identical to <kh>, in which
 *            case the conversion is done in-place.
 *
 * Note:      Character > is converted to <  (left base of base pairs)
 *            Character < is converted to >  (right base of base pairs)
 *            A space is converted to .      (just in case)      
 *
 * Returns:   <eslOK> on success.
 */
int
esl_kh2wuss(char *kh, char *ss)
{
  while (*kh != '\0')
    {
      if      (*kh == '>') *ss = '<';
      else if (*kh == '<') *ss = '>';
      else if (*kh == ' ') *ss = '.';
      else *ss = *kh;
      kh++;
      ss++;
    }
  *ss = '\0';
  return eslOK;
}


/* Function:  esl_wuss_full()
 * Incept:    SRE, Mon Feb 28 09:44:40 2005 [St. Louis]
 *
 * Purpose:   Given a simple ("input") WUSS format annotation string <oldss>,
 *            convert it to full ("output") WUSS format in <newss>.
 *            <newss> must be allocated by the caller to be at least as 
 *            long as <oldss>. <oldss> and <newss> can be the same,
 *            to convert a secondary structure string in place.
 *            
 *            Pseudoknot annotation is preserved, if <oldss> had it.
 *
 * Returns:   <eslSYNTAX> if <oldss> isn't in valid WUSS format.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslEINCONCEIVABLE> on internal error that can't happen.
 */
int
esl_wuss_full(char *oldss, char *newss)
{
  char *tmp = NULL;
  int  *ct  = NULL;
  int   n;
  int   i;
  int   status;

  /* We can use the ct2wuss algorithm to generate a full WUSS string -
   * convert to ct, then back to WUSS.  ct2wuss doesn't deal with pk's
   * though, and we want to propagate pk annotation if it's there.  So
   * we need two workspaces: ct array, and a temporary ss string that
   * we use to hold non-pk annotation.  As a final step, we overlay
   * the pk annotation from the original oldss annotation.
   */
  n = strlen(oldss);
  ESL_ALLOC(ct,  sizeof(int)  * (n+1));
  ESL_ALLOC(tmp, sizeof(char) * (n+1));
  
  esl_wuss_nopseudo(oldss, tmp);/* tmp = nonpseudoknotted oldss */

  status = esl_wuss2ct(tmp, n, ct);   /* ct  = oldss in ct format, no pks */
  if (status != eslOK) goto ERROR;

  status = esl_ct2wuss(ct, n, tmp);   /* now tmp is a full WUSS string */
  if (status == eslEINVAL) { status = eslEINCONCEIVABLE; goto ERROR; }/* we're sure, no pk's */
  else if (status != eslOK) goto ERROR; /* EMEM, EINCONCEIVABLE  */
  
  for (i = 0; i < n; i++)
    if (isalpha(oldss[i])) newss[i] = oldss[i];	/* transfer pk annotation */
    else newss[i] = tmp[i];                     /* transfer new WUSS      */

  free(ct);
  free(tmp);
  return eslOK;

 ERROR:
  free(ct);
  free(tmp);
  return status;
}



/* Function:  esl_wuss_nopseudo()
 * Incept:    SRE, Tue Feb 15 11:02:43 2005 [St. Louis]
 *
 * Purpose:   Given a WUSS format annotation string <ss1>,
 *            removes all pseudoknot annotation to create a new 
 *            WUSS string <ss2> that contains only a "canonical"
 *            (nonpseudoknotted) structure. <ss2> must be allocated to
 *            be at least as large as <ss1>. <ss1> and <ss2>
 *            may be the same, in which case the conversion is
 *            done in place. Pseudoknot annotation in <ss1> is
 *            simply replaced by <.> in <ss2>; the resulting
 *            <ss2> WUSS string is therefore in valid simplified format,
 *            but may not be valid full format WUSS.
 *
 * Returns:   <eslOK>.
 */
int
esl_wuss_nopseudo(char *ss1, char *ss2)
{
  while (*ss1 != '\0') 
    {
      if (isalpha(*ss1)) *ss2 = '.';
      else *ss2 = *ss1;
      ss1++;
      ss2++;
    }
  *ss2 = '\0';
  return eslOK;
}


#ifdef eslWUSS_TESTDRIVE
/* gcc -g -Wall -o test -I. -DeslWUSS_TESTDRIVE wuss.c stack.c easel.c
 * ./test
 */
#include <stdlib.h>

#include "easel.h"
#include "esl_stack.h"
#include "esl_wuss.h"

int
main(int argc, char **argv)
{
  /* The example is E. coli RNase P, w/ and w/o pks. 
   * J Brown figure 10.3.00 shows 1 too many bp for pk stem A. 
   */
  char ss[] = "\
{{{{{{{{{{{{{{{{{{,<<<<<<<<<<<<<-<<<<<____>>>>>>>>>->>>>>>>>\
>,,,,AAA-AAAAA[[[[---BBBB-[[[[[<<<<<_____>>>>><<<<____>>>->(\
(---(((((,,,,,,,,,,,,<<<<<--<<<<<<<<____>>>>>->>>>>>-->>,,,,\
,,,<<<<<<_______>>>>>><<<<<<<<<____>>>->>>>>->,,)))--))))]]]\
]]]]]],,,<<<<------<<<<<<----<<<<<_bbbb>>>>>>>>>>>----->>>>,\
,,,,,<<<<<<<<____>>>>>>>>,,,,,,,,,,}}}}}}}----------aaaaaaaa\
-}-}}}}}}}}}}::::";
  char ss_nopk[] = "\
{{{{{{{{{{{{{{{{{{,<<<<<<<<<<<<<-<<<<<____>>>>>>>>>->>>>>>>>\
>,,,,,,,,,,,,,[[[[--------[[[[[<<<<<_____>>>>><<<<____>>>->(\
(---(((((,,,,,,,,,,,,<<<<<--<<<<<<<<____>>>>>->>>>>>-->>,,,,\
,,,<<<<<<_______>>>>>><<<<<<<<<____>>>->>>>>->,,)))--))))]]]\
]]]]]],,,<<<<------<<<<<<----<<<<<_____>>>>>>>>>>>----->>>>,\
,,,,,<<<<<<<<____>>>>>>>>,,,,,,,,,,}}}}}}}------------------\
-}-}}}}}}}}}}::::";
  int  len;
  int  *ct1 = NULL;
  int  *ct2 = NULL;
  char *ss2 = NULL;
  int  i;
  int  nbp, nbp_true, npk;
  int  status;

  len = strlen(ss);
  ESL_ALLOC(ct1, sizeof(int)  * (len+1));
  ESL_ALLOC(ct2, sizeof(int)  * (len+1));
  ESL_ALLOC(ss2, sizeof(char) * (len+1));
  nbp_true = npk = 0;
  for (i = 0; i < len; i++)
    {
      if (strchr("{[(<", ss[i]) != NULL)
	nbp_true++;
      if (isupper(ss[i]))
	npk++;
    }
	
  if (esl_wuss2ct(ss, len, ct1) != eslOK) abort();
  nbp = 0;
  for (i = 1; i <= len; i++)
    if (ct1[i] > i) nbp++;
  if (nbp != nbp_true + npk) abort();

  if (esl_wuss2kh(ss, ss)       != eslOK) abort();
  if (esl_kh2wuss(ss, ss)      != eslOK) abort();
  if (esl_wuss2ct(ss, len, ct2) != eslOK) abort();
  
  for (i = 1; i <= len; i++)
    if (ct1[i] != ct2[i]) abort();

  if (esl_wuss_nopseudo(ss, ss)      != eslOK) abort();
  if (esl_wuss2ct(ss, len, ct1)      != eslOK) abort();
  if (esl_wuss2ct(ss_nopk, len, ct2) != eslOK) abort();
  for (i = 1; i <= len; i++)
    if (ct1[i] != ct2[i]) abort();

  if (esl_wuss2ct(ss_nopk, len, ct1) != eslOK) abort();
  if (esl_ct2wuss(ct1, len, ss2)     != eslOK) abort();
  if (strcmp(ss_nopk, ss2) != 0) abort();

  free(ct1);
  free(ct2);
  free(ss2);
  return 0;

 ERROR:
  free(ct1);
  free(ct2);
  free(ss2);
  return status;
}
#endif /*eslWUSS_TESTDRIVE*/



/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
