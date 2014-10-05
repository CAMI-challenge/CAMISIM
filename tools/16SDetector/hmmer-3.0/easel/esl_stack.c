/* Pushdown stacks for integers, pointers, and characters.
 *
 * Contents:
 *   1. The <ESL_STACK> object.
 *   2. Other functions in the API.
 *   3. Shuffling stacks.      [eslAUGMENT_RANDOM]
 *   4. Unit tests.
 *   5. Test driver.
 *   6. Example.
 *   7. Copyright and license.
 *
 * Augmentations:
 *   eslAUGMENT_RANDOM  : adds function for shuffling a stack. 
 * 
 * SRE 1 March 2000 [Seattle]
 * Incorp into Easel SRE, Sun Dec 26 07:42:12 2004 [Zaragoza]
 * SVN $Id: esl_stack.c 249 2008-04-24 19:19:50Z eddys $
 */ 
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_stack.h"
#ifdef eslAUGMENT_RANDOM
#include "esl_random.h"
#endif

/*****************************************************************
 *# 1. The <ESL_STACK> object.
 *****************************************************************/

/* Function:  esl_stack_ICreate()
 * Synopsis:  Create an integer stack.
 * Incept:    SRE, Sun Dec 26 09:11:50 2004 [Zaragoza]
 *
 * Purpose:   Creates an integer stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_ICreate(void)
{
  int status;
  ESL_STACK *ns = NULL;
  
  ESL_ALLOC(ns, sizeof(ESL_STACK));
  ns->nalloc   = ESL_STACK_INITALLOC;
  ns->pdata    = NULL;
  ns->cdata    = NULL;
  ESL_ALLOC(ns->idata, sizeof(int) * ns->nalloc);
  ns->n        = 0;
  return ns;

 ERROR:
  esl_stack_Destroy(ns);
  return NULL;
}

/* Function:  esl_stack_CCreate()
 * Synopsis:  Create a character stack.
 * Incept:    SRE, Sun Dec 26 09:15:35 2004 [Zaragoza]
 *
 * Purpose:   Creates a character stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_CCreate(void)
{
  int status;
  ESL_STACK *cs = NULL;
  
  ESL_ALLOC(cs, sizeof(ESL_STACK));
  cs->nalloc   = ESL_STACK_INITALLOC;
  cs->idata    = NULL;
  cs->pdata    = NULL;
  ESL_ALLOC(cs->cdata, sizeof(char) * cs->nalloc);
  cs->n        = 0;
  return cs;

 ERROR:
  esl_stack_Destroy(cs);
  return NULL;
}

/* Function:  esl_stack_PCreate()
 * Synopsis:  Create a pointer stack.
 * Incept:    SRE, Sun Dec 26 09:16:07 2004 [Zaragoza]
 *
 * Purpose:   Creates a pointer stack.
 *
 * Returns:   a pointer to the new stack.
 *
 * Throws:    <NULL> on an allocation failure.
 */
ESL_STACK *
esl_stack_PCreate(void)
{
  int status;
  ESL_STACK *ps = NULL;
  
  ESL_ALLOC(ps, sizeof(ESL_STACK));
  ps->nalloc   = ESL_STACK_INITALLOC;
  ps->idata    = NULL;
  ps->cdata    = NULL;
  ESL_ALLOC(ps->pdata, sizeof(void *) * ps->nalloc);
  ps->n        = 0;
  return ps;

 ERROR:
  esl_stack_Destroy(ps);
  return NULL;
}

/* Function:  esl_stack_Reuse()
 * Synopsis:  Reuse a stack.
 * Incept:    SRE, Tue Dec 28 04:21:36 2004 [Zaragoza]
 *
 * Purpose:   Empties stack <s> so it can be reused without
 *            creating a new one. The stack <s>
 *            can be of any data type; it retains its original
 *            type.
 *
 * Returns:   <eslOK>
 */
int
esl_stack_Reuse(ESL_STACK *s)
{
  s->n = 0;	/* it's that simple in this implementation */
  return eslOK;
}

/* Function:  esl_stack_Destroy()
 * Synopsis:  Free a stack.
 * Incept:    SRE, Sun Dec 26 09:16:24 2004 [Zaragoza]
 *
 * Purpose:   Destroys a created stack <s>, of any data type.
 */
void
esl_stack_Destroy(ESL_STACK *s)
{
  if (s->idata != NULL) free(s->idata);
  if (s->cdata != NULL) free(s->cdata);
  if (s->pdata != NULL) free(s->pdata);
  free(s);
}


/*****************************************************************
 *# 2. Other functions in the API.
 *****************************************************************/

/* Function:  esl_stack_IPush()
 * Synopsis:  Push an integer onto a stack.
 * Incept:    SRE, Sun Dec 26 09:17:17 2004 [Zaragoza]
 *
 * Purpose:   Push an integer <x> onto an integer stack <ns>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
esl_stack_IPush(ESL_STACK *ns, int x)
{
  int  status;
  int *ptr;

  if (ns->n == ns->nalloc) {
    ESL_RALLOC(ns->idata, ptr, sizeof(int) * ns->nalloc * 2);
    ns->nalloc += ns->nalloc;	/* reallocate by doubling */
  }
  ns->idata[ns->n] = x;
  ns->n++;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_stack_CPush()
 * Synopsis:  Push a char onto a stack.
 * Incept:    SRE, Sun Dec 26 09:18:24 2004 [Zaragoza]
 *
 * Purpose:   Push a character <c> onto a character stack <cs>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
esl_stack_CPush(ESL_STACK *cs, char c)
{
  int  status;
  char *ptr;

  if (cs->n == cs->nalloc) {
    ESL_RALLOC(cs->cdata, ptr, sizeof(char) * cs->nalloc * 2);
    cs->nalloc += cs->nalloc;	/* reallocate by doubling */
  }
  cs->cdata[cs->n] = c;
  cs->n++;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_stack_PPush()
 * Synopsis:  Push a pointer onto a stack.
 * Incept:    SRE, Sun Dec 26 09:18:49 2004 [Zaragoza]
 *
 * Purpose:   Push a pointer <p> onto a pointer stack <ps>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure.
 */
int
esl_stack_PPush(ESL_STACK *ps, void *p)
{
  int status;
  void *ptr;

  if (ps->n == ps->nalloc) {
    ESL_RALLOC(ps->pdata, ptr, sizeof(void *) * ps->nalloc * 2);
    ps->nalloc += ps->nalloc;	/* reallocate by doubling */
  }
  ps->pdata[ps->n] = p;
  ps->n++;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  esl_stack_IPop()
 * Synopsis:  Pop an integer off a stack.
 * Incept:    SRE, Sun Dec 26 09:19:12 2004 [Zaragoza]
 *
 * Purpose:   Pops an integer off the integer stack <ns>, and returns
 *            it through <ret_x>.
 *
 * Returns:   <eslOK> on success. <eslEOD> if stack is empty.
 */
int
esl_stack_IPop(ESL_STACK *ns, int *ret_x)
{
  if (ns->n == 0) {*ret_x = 0; return eslEOD;}
  ns->n--;
  *ret_x = ns->idata[ns->n];
  return eslOK;
}

/* Function:  esl_stack_CPop()
 * Synopsis:  Pop a char off a stack.
 * Incept:    SRE, Sun Dec 26 09:21:27 2004 [Zaragoza]
 *
 * Purpose:   Pops a character off the character stack <cs>, and returns
 *            it through <ret_c>.
 *
 * Returns:   <eslOK> on success. <eslEOD> if stack is empty.
 */
int
esl_stack_CPop(ESL_STACK *cs, char *ret_c)
{
  if (cs->n == 0) {*ret_c = 0; return eslEOD;}
  cs->n--;
  *ret_c = cs->cdata[cs->n];
  return eslOK;
}

/* Function:  esl_stack_PPop()
 * Synopsis:  Pop a pointer off a stack.
 * Incept:    SRE, Sun Dec 26 09:21:56 2004 [Zaragoza]
 *
 * Purpose:   Pops a pointer off the pointer stack <ps>, and returns
 *            it through <ret_p>.
 *
 * Returns:   <eslOK> on success. <eslEOD> if stack is empty.
 */
int
esl_stack_PPop(ESL_STACK *ps, void **ret_p)
{
  if (ps->n == 0) {*ret_p = 0; return eslEOD;}
  ps->n--;
  *ret_p = ps->pdata[ps->n];
  return eslOK;
}

/* Function:  esl_stack_ObjectCount()
 * Synopsis:  Return the number of objects in a stack.
 * Incept:    SRE, Sun Dec 26 09:22:41 2004 [Zaragoza]
 *
 * Purpose:   Returns the number of data objects stored in the
 *            stack <s>. The stack may be of any datatype.
 */
int 
esl_stack_ObjectCount(ESL_STACK *s)
{
  return s->n;
}

/* Function:  esl_stack_Convert2String()
 * Synopsis:  Convert a char stack to a string.
 * Incept:    SRE, Sun Dec 26 09:23:36 2004 [Zaragoza]
 *
 * Purpose:   Converts a character stack <cs> to a NUL-terminated
 *            string, and returns a pointer to the string. The
 *            characters in the string are in the same order they
 *            were pushed onto the stack.  The stack is destroyed by
 *            this operation, as if <esl_stack_Destroy()> had been
 *            called on it. The caller becomes responsible for
 *            free'ing the returned string.
 *
 * Returns:   Pointer to the string; caller must <free()> this.
 *
 * Throws:    NULL if a reallocation fails.
 */
char *
esl_stack_Convert2String(ESL_STACK *cs)
{
  char *s;

  if (esl_stack_CPush(cs, '\0') != eslOK)
    { free(cs->cdata); free(cs); return NULL; } /* nul-terminate the data or self-destruct */
  s = cs->cdata;		           /* data is already just a string - just return ptr to it */
  free(cs);			           /* free the stack around it. */
  return s;
}

/* Function:  esl_stack_DiscardTopN()
 * Synopsis:  Discard the top elements on a stack.
 * Incept:    SRE, Tue Dec 28 04:33:06 2004 [St. Louis]
 *
 * Purpose:   Throw away the top <n> elements on stack <s>.
 *            Equivalent to <n> calls to a <Pop()> function.
 *            If <n> equals or exceeds the number of elements 
 *            currently in the stack, the stack is emptied
 *            as if <esl_stack_Reuse()> had been called.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stack_DiscardTopN(ESL_STACK *s, int n)
{
  if (n <= s->n) s->n -= n;
  else           s->n = 0;
  return eslOK;
}

/*****************************************************************
 *# 3. Shuffling stacks [with <eslAUGMENT_RANDOM>]
 *****************************************************************/
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_stack_Shuffle()
 * Synopsis:  Randomly shuffle the elements in a stack.
 * Incept:    SRE, Mon Mar 31 11:01:06 2008 [Janelia]
 *
 * Purpose:   Randomly shuffle the elements in stack <s>, using
 *            random numbers from generator <r>.
 *
 * Returns:   <eslOK> on success, and the stack is randomly 
 *            shuffled.
 */
int
esl_stack_Shuffle(ESL_RANDOMNESS *r, ESL_STACK *s)
{
  int   n = s->n;
  int   w;

  while (n > 1) {
    w = esl_rnd_Roll(r, n);	/* shuffling algorithm: swap last elem with w, decrement n. */
    if      (s->idata != NULL)  ESL_SWAP(s->idata[w], s->idata[n-1], int);
    else if (s->cdata != NULL)  ESL_SWAP(s->cdata[w], s->cdata[n-1], char);
    else if (s->pdata != NULL)  ESL_SWAP(s->pdata[w], s->pdata[n-1], void *);
    n--;
  }
  return eslOK;
}
#endif /*eslAUGMENT_RANDOM*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef eslSTACK_TESTDRIVE

static void
utest_integer(void)
{
  char      *msg = "integer stack basic unit test failed";
  ESL_STACK *s   = NULL;
  int        n1  = ESL_STACK_INITALLOC*2+1;		/* force two reallocations */
  int        n2  = 0;
  int        i;
  int        val;

  if ((s = esl_stack_ICreate())                      == NULL)   esl_fatal(msg);
  for (i = 0; i < n1; i++) if (esl_stack_IPush(s, i) != eslOK)  esl_fatal(msg);
  while (esl_stack_IPop(s, &val) != eslEOD) n2++;
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Reuse(s);

  /* same again, with ObjectCount instead of EOD for popping */
  for (i = 0; i < n1; i++) if (esl_stack_IPush(s, i) != eslOK) esl_fatal(msg);
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_IPop(s, &val) != eslOK) esl_fatal(msg);
    n2++; 
  }
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Destroy(s);
}

static void
utest_char(void)  
{
  char      *msg = "char stack basic unit test failed";
  ESL_STACK *s   = NULL;
  int        n1  = ESL_STACK_INITALLOC*2+1;		/* force two reallocations */
  int        n2  = 0;
  int        i;
  char       c;

  if ((s = esl_stack_CCreate())                        == NULL)   esl_fatal(msg);
  for (i = 0; i < n1; i++) if (esl_stack_CPush(s, 'X') != eslOK)  esl_fatal(msg);
  while (esl_stack_CPop(s, &c) != eslEOD) {
    if (c != 'X') esl_fatal(msg);
    n2++; 
  }
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Reuse(s);

  /* same again, with ObjectCount instead of EOD for popping */
  for (i = 0; i < n1; i++) if (esl_stack_CPush(s, 'X') != eslOK) esl_fatal(msg);
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_CPop(s, &c) != eslOK) esl_fatal(msg);
    n2++; 
  }
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Destroy(s);
}
  
static void
utest_pointer(void)
{
  char      *msg = "pointer stack basic unit test failed";
  ESL_STACK *s   = NULL;
  int        n1  = ESL_STACK_INITALLOC*2+1;		/* force two reallocations */
  int        n2  = 0;
  int        i;
  void      *p;

  if ((s = esl_stack_PCreate())                        == NULL)   esl_fatal(msg);
  for (i = 0; i < n1; i++) {
    p = malloc(sizeof(int) * 64);
    if (esl_stack_PPush(s, p) != eslOK)  esl_fatal(msg);
  }
  while (esl_stack_PPop(s, &p) != eslEOD) { free(p); n2++; }
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Reuse(s);

  /* same again, with ObjectCount instead of EOD for popping */
  for (i = 0; i < n1; i++) {
    p = malloc(sizeof(int) * 64);
    if (esl_stack_PPush(s, p) != eslOK) esl_fatal(msg);
  }
  n2 = 0;
  while (esl_stack_ObjectCount(s)) {
    if (esl_stack_PPop(s, &p) != eslOK) esl_fatal(msg);
    free(p);
    n2++; 
  }
  if (n1 != n2) esl_fatal(msg);
  esl_stack_Destroy(s);
}  

static void
utest_convert2string(void)
{
  char      *msg = "stack::convert2string unit test failed";
  char      *str = "ABCDEF";
  ESL_STACK *s   = NULL;
  int        n   = strlen(str);
  int        i;
  char      *result = NULL;

  if ((s = esl_stack_CCreate())                          == NULL)   esl_fatal(msg);
  for (i = 0; i < n; i++) if (esl_stack_CPush(s, str[i]) != eslOK)  esl_fatal(msg);
  result = esl_stack_Convert2String(s);
  if (strcmp(result, str) != 0) esl_fatal(msg);
  free(result);	/* after Convert2String, only the string itself remains to be free'd */
}


#ifdef eslAUGMENT_RANDOM
static void
utest_shuffle(void)
{
  char           *msg  = "stack shuffle unit test failed";
  ESL_RANDOMNESS *r    = esl_randomness_Create(42);
  ESL_STACK      *s    = esl_stack_ICreate();
  int             n    = ESL_STACK_INITALLOC*2+1;      /* exercises reallocation */
  int            *seen = malloc(sizeof(int) * n);
  int             i;
  int             val;
  int             appears_shuffled = FALSE;

  for (i = 0; i < n; i++) esl_stack_IPush(s, i);
  esl_stack_Shuffle(r, s);
  
  for (i = 0; i < n; i++) seen[i] = 0;
  i = n-1;
  while (esl_stack_IPop(s, &val) != eslEOD) {
    seen[val]++;
    if (val != i--) appears_shuffled = TRUE;
  }
  for (i = 0; i < n; i++) if (seen[i] != 1) esl_fatal(msg);
  
  free(seen);
  esl_stack_Destroy(s);
  esl_randomness_Destroy(r);
}
#endif /*eslAUGMENT_RANDOM*/


#endif /*eslSTACK_TESTDRIVE*/
/*---------------- end of unit tests ----------------------------*/




/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

/*****************************************************************
 * Test driver and API example for the pushdown stack module.
 * To compile:
 *    gcc -g -Wall -I. -L. -DeslSTACK_TESTDRIVE -o testdrive esl_stack.c -leasel -lm
 * To run:
 *    ./testdrive
 * Returns 0 (success) w/ no output, or returns nonzero and says why.
 *****************************************************************/

/* why Pop() into a void *obj_p, instead of directly into int *obj, in
 * the test of the pointer stack? On PowerPC/Linux, using gcc -O3,
 * trying to Pop() into *obj causes a "dereferencing type-punned
 * pointer will break strict-aliasing rules" warning, and the test
 * driver crashes with a double free/corruption error in glibc.  Lower
 * optimization levels don't cause the problem; adding
 * -fno-strict-aliasing to the CFLAGS also avoids the problem. I'm
 * suspicious that it's a gcc optimizer bug. Pop()'ing into a void *
 * avoids the issue altogether. (SRE, Feb 22 2008 J2/119)
 */
#ifdef eslSTACK_TESTDRIVE
int 
main(void)
{
  utest_integer();
  utest_char();
  utest_pointer();
  utest_convert2string();

#ifdef eslAUGMENT_RANDOM
  utest_shuffle();
#endif

  return eslOK;
}
#endif /*eslSTACK_TESTDRIVE*/
/*-------------------- end of test driver -----------------------*/




/*****************************************************************
 * 6. Example.
 *****************************************************************/
#ifdef eslSTACK_EXAMPLE
/*::cexcerpt::stack_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslSTACK_EXAMPLE esl_stack.c easel.c -lm
 * run:     ./example
 */
#include "easel.h"
#include "esl_stack.h"

int
main(void)
{
  ESL_STACK *ns;
  int        x;

  ns = esl_stack_ICreate();
  esl_stack_IPush(ns, 42);
  esl_stack_IPush(ns, 7);
  esl_stack_IPush(ns, 3);
  while (esl_stack_IPop(ns, &x) != eslEOD) 
    printf("%d\n", x);
  esl_stack_Destroy(ns);   
}
/*::cexcerpt::stack_example::end::*/
#endif /*eslSTACK_EXAMPLE*/
/*------------------------ end of example -----------------------*/




/*****************************************************************  
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
