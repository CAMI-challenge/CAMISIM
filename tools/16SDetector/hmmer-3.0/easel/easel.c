/* Easel's foundation.
 * 
 * Contents:
 *    1. Error handling conventions.
 *    2. Memory allocation/deallocation conventions.
 *    3. Standard banner for Easel miniapplications.
 *    4. Replacements for some C library functions.
 *    5. File path/name manipulation, including tmpfiles.
 *    6. Typed comparison functions.
 *    7. Commonly used background composition (iid) frequencies.
 *    8. Unit tests [need to be written]
 *    9. Test driver [needs to be written]
 *   10. Examples. [need to be written]
 *   11. Copyright and license. 
 * 
 * SRE, Tue Oct 28 08:29:17 2003 [St. Louis]
 * SVN $Id: easel.c 509 2010-02-07 22:56:55Z eddys $
 */
#include "esl_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "easel.h"


/*****************************************************************
 * 1. Error handling.
 *****************************************************************/
static esl_exception_handler_f esl_exception_handler = NULL;

void
esl_exception_SetHandler(void (*handler)(int code, char *file, int line, 
					 char *format, va_list argp))
{
  esl_exception_handler = handler;
}

void
esl_exception_ResetDefaultHandler(void)
{
  esl_exception_handler = NULL;
}

void
esl_exception(int code, char *file, int line, char *format, ...)
{
  va_list argp;

  if (esl_exception_handler != NULL) {
    va_start(argp, format);
    (*esl_exception_handler)(code, file, line, format, argp);
    va_end(argp);
    return;
  } else {
    fprintf(stderr, "Fatal exception (source file %s, line %d):\n", file, line);
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);
    abort();
  }
}

void
esl_fatal(const char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
}

/* esl_nonfatal_handler()
 * SRE, Fri Sep  8 10:59:14 2006 [Janelia]
 * 
 * This stub is here to support the test harnesses, when they 
 * have eslTEST_THROWING turned on to test that thrown errors
 * are handled properly when a nonfatal error handler is
 * registered by the application.
 */
void
esl_nonfatal_handler(int code, char *file, int line, char *format, va_list argp)
{
  return;
}
/*---------------- end, error handling conventions --------------*/




/*****************************************************************
 * 2. Memory allocation/deallocation conventions.
 *****************************************************************/

/* Function:  esl_Free2D()
 * Incept:    squid's Free2DArray(), 1999.
 *
 * Purpose:   Free a 2D pointer array <p>, where first dimension is
 *            <dim1>. (That is, the array is <p[0..dim1-1][]>.)
 *            Tolerates any of the pointers being NULL, to allow
 *            sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free2D(void **p, int dim1)
{
  int i;
  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) free(p[i]);
    free(p);
  }
  return;
}

/* Function:  esl_Free3D()
 * Incept:    squid's Free3DArray(), 1999.
 *
 * Purpose:   Free a 3D pointer array <p>, where first and second
 *            dimensions are <dim1>,<dim2>. (That is, the array is
 *            <p[0..dim1-1][0..dim2-1][]>.) Tolerates any of the
 *            pointers being NULL, to allow sparse arrays.
 *
 * Returns:   void.
 */
void
esl_Free3D(void ***p, int dim1, int dim2)
{
  int i, j;

  if (p != NULL) {
    for (i = 0; i < dim1; i++)
      if (p[i] != NULL) {
        for (j = 0; j < dim2; j++)
          if (p[i][j] != NULL) free(p[i][j]);
        free(p[i]);
      }
    free(p);
  }
}
/*------------- end, memory allocation conventions --------------*/



/*****************************************************************
 * 3. Standard banner for Easel miniapplications.
 *****************************************************************/

/* Function:  esl_banner()
 * Synopsis:  print standard Easel application output header
 * Incept:    SRE, Mon Feb 14 11:26:56 2005 [St. Louis]
 *
 * Purpose:   Print the standard Easel command line application banner
 *            to <fp>, constructing it from <progname> (the name of the
 *            program) and a short one-line description <banner>.
 *            For example, 
 *            <esl_banner(stdout, "compstruct", "compare RNA structures");>
 *            might result in:
 *            
 *            \begin{cchunk}
 *            # compstruct :: compare RNA structures
 *            # Easel 0.1 (February 2005)
 *            # Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus
 *            # Freely licensed under the Janelia Software License.
 *            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *            \end{cchunk}
 *              
 *            <progname> would typically be an application's
 *            <argv[0]>, rather than a fixed string. This allows the
 *            program to be renamed, or called under different names
 *            via symlinks. Any path in the <progname> is discarded;
 *            for instance, if <progname> is "/usr/local/bin/esl-compstruct",
 *            "esl-compstruct" is used as the program name.
 *            
 * Note:    
 *    Needs to pick up preprocessor #define's from easel.h,
 *    as set by ./configure:
 *            
 *    symbol          example
 *    ------          ----------------
 *    EASEL_VERSION   "0.1"
 *    EASEL_DATE      "May 2007"
 *    EASEL_COPYRIGHT "Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
 *    EASEL_LICENSE   "Freely licensed under the Janelia Software License."
 *
 * Returns:   (void)
 */
void
esl_banner(FILE *fp, char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) appname = progname;

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# Easel %s (%s)\n", EASEL_VERSION, EASEL_DATE);
  fprintf(fp, "# %s\n", EASEL_COPYRIGHT);
  fprintf(fp, "# %s\n", EASEL_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname != NULL) free(appname);
  return;
}


/* Function:  esl_usage()
 * Synopsis:  print standard Easel application usage help line
 * Incept:    SRE, Wed May 16 09:04:42 2007 [Janelia]
 *
 * Purpose:   Given a usage string <usage> and the name of the program
 *            <progname>, output a standardized usage/help
 *            message. <usage> is minimally a one line synopsis like
 *            "[options] <filename>", but it may extend to multiple
 *            lines to explain the command line arguments in more 
 *            detail. It should not describe the options; that's the
 *            job of the getopts module, and its <esl_opt_DisplayHelp()> 
 *            function.
 *            
 *            This is used by the Easel miniapps, and may be useful in
 *            other applications as well.
 *
 *            As in <esl_banner()>, the <progname> is typically passed
 *            as <argv[0]>, and any path prefix is ignored.
 *            
 *            For example, if <argv[0]> is </usr/local/bin/esl-compstruct>,
 *            then 
 *            
 *            \begin{cchunk}
 *              esl_usage(stdout, argv[0], "[options] <trusted file> <test file>">
 *            \end{cchunk}
 *            
 *            produces
 *            
 *            \begin{cchunk}
 *              Usage: esl-compstruct [options] <trusted file> <test file>
 *            \end{cchunk}  
 *              
 * Returns:   (void).
 */
void
esl_usage(FILE *fp, char *progname, char *usage)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) appname = progname;
  fprintf(fp, "Usage: %s %s\n", appname, usage);
  if (appname != NULL) free(appname);
  return;
}


/*-------------------- end, standard miniapp banner --------------------------*/




/******************************************************************************
 * 4. Replacements for C library functions
 *  fgets()   ->  esl_fgets()     fgets() with dynamic allocation
 *  strdup()  ->  esl_strdup()    strdup() is not ANSI
 *  strcat()  ->  esl_strcat()    strcat() with dynamic allocation
 *  strtok()  ->  esl_strtok()    threadsafe strtok()
 *  sprintf() ->  esl_sprintf()   sprintf() with dynamic allocation
 *  strcmp()  ->  esl_strcmp()    strcmp() tolerant of NULL strings
 *  free()    ->  esl_free()      free() tolerant of NULL pointers
 *****************************************************************************/

/* Function: esl_fgets()
 * Date:     SRE, Thu May 13 10:56:28 1999 [St. Louis]
 *
 * Purpose:  Dynamic allocation version of fgets(),
 *           capable of reading almost unlimited line lengths.
 *
 * Args:     buf - ptr to a string (may be reallocated)
 *           n   - ptr to current allocated length of buf,
 *                 (may be changed)
 *           fp  - open file ptr for reading
 *           
 *           Before the first call to esl_fgets(), 
 *           initialize buf to NULL and n to 0.
 *           They're a linked pair, so don't muck with the
 *           allocation of buf or the value of n while
 *           you're still doing esl_fgets() calls with them.
 *
 * Returns:  <eslOK> on success. 
 *           Returns <eslEOF> on normal end-of-file.
 *
 *           When <eslOK>:
 *           <*buf> points to a <NUL>-terminated line from the file.
 *           <*n> contains the current allocated length for <*buf>.
 * 
 *           Caller must free <*buf> eventually. 
 *
 * Throws:   <eslEMEM> on an allocation failure.
 *
 * Example:  char *buf;
 *           int   n;
 *           FILE *fp;
 *           
 *           fp  = fopen("my_file", "r");
 *           buf = NULL;
 *           n   = 0;
 *           while (esl_fgets(&buf, &n, fp) == eslOK) 
 *           {
 *             do stuff with buf;
 *           }
 *           if (buf != NULL) free(buf);
 */
int
esl_fgets(char **buf, int *n, FILE *fp)
{
  int   status;
  void *p;
  char *s;
  int   len;
  int   pos;

  if (*n == 0) 
    {
      ESL_ALLOC(*buf, sizeof(char) * 128);
      *n   = 128;
    }

  /* Simple case 1. We're sitting at EOF, or there's an error.
   *                fgets() returns NULL, so we return EOF.
   */
  if (fgets(*buf, *n, fp) == NULL) return eslEOF;

  /* Simple case 2. fgets() got a string, and it reached EOF doing it.
   *                return success status, so caller can use
   *                the last line; on the next call we'll
   *                return the 0 for the EOF.
   */
  if (feof(fp)) return eslOK;

  /* Simple case 3. We got a complete string, with \n,
   *                and don't need to extend the buffer.
   */
  len = strlen(*buf);
  if ((*buf)[len-1] == '\n') return eslOK;

  /* The case we're waiting for. We have an incomplete string,
   * and we have to extend the buffer one or more times. Make
   * sure we overwrite the previous fgets's \0 (hence +(n-1)
   * in first step, rather than 128, and reads of 129, not 128).
   */
  pos = (*n)-1;
  while (1) {
    ESL_RALLOC(*buf, p, sizeof(char) * (*n+128));
    *n  += 128;
    s = *buf + pos;
    if (fgets(s, 129, fp) == NULL) return eslOK;
    len = strlen(s);
    if (s[len-1] == '\n') return eslOK;
    pos += 128;
  } 
  /*NOTREACHED*/
  return eslOK;

 ERROR:
  if (*buf != NULL) free(*buf);
  *buf = NULL;
  *n   = 0;
  return status;
}

/* Function: esl_strdup()
 * Date:     SRE, Wed May 19 17:57:28 1999 [St. Louis]
 *
 * Purpose: Makes a duplicate of string <s>, puts it in <ret_dup>.
 *          Caller can pass string length <n>, if it's known,
 *          to save a strlen() call; else pass -1 to have the string length
 *          determined.
 *          
 *          Tolerates <s> being <NULL>; in which case,
 *          returns <eslOK> with <*ret_dup> set to <NULL>.
 *
 * Args:     s       - string to duplicate (NUL-terminated)
 *           n       - length of string, if known; -1 if unknown.
 *           ret_dup - RETURN: duplicate of <s>.
 *                
 * Returns:  <eslOK> on success, and <ret_dup> is valid.
 *
 * Throws:   <eslEMEM> on allocation failure.
 */
int
esl_strdup(const char *s, int64_t n, char **ret_dup)
{
  int   status;
  char *new = NULL;

  if (ret_dup != NULL) *ret_dup = NULL;
  if (s == NULL) return eslOK;
  if (n < 0) n = strlen(s);

  ESL_ALLOC(new, sizeof(char) * (n+1));
  strcpy(new, s);

  if (ret_dup != NULL) *ret_dup = new; else free(new);
  return eslOK;

 ERROR:
  if (new     != NULL) free(new);
  if (ret_dup != NULL) *ret_dup = NULL;
  return status;
}


/* Function: esl_strcat()
 * Date:     SRE, Thu May 13 09:36:32 1999 [St. Louis]
 *
 * Purpose:  Dynamic memory version of strcat().
 *           Appends <src> to the string that <dest> points to,
 *           extending allocation for dest if necessary. Caller
 *           can optionally provide the length of <*dest> in
 *           <ldest>, and the length of <src> in <lsrc>; if 
 *           either of these is -1, <esl_strcat()> calls <strlen()>
 *           to determine the length. Providing length information,
 *           if known, accelerates the routine.
 *           
 *           <*dest> may be <NULL>, in which case this is equivalent
 *           to a <strdup()> of <src> (that is, <*dest> is allocated
 *           rather than reallocated). 
 *           
 *           <src> may be <NULL>, in which case <dest> is unmodified.
 *           
 * Note:     One timing experiment (100 successive appends of 
 *           1-255 char) shows esl_strcat() has about a 20%
 *           overhead relative to strcat(). If optional
 *           length info is passed, esl_strcat() is about 30%
 *           faster than strcat().
 *           
 * Args:     dest  - ptr to string (char **), '\0' terminated
 *           ldest - length of dest, if known; or -1 if length unknown.
 *           src   - string to append to dest, '\0' terminated       
 *           lsrc  - length of src, if known; or -1 if length unknown.
 *
 * Returns:  <eslOK> on success; <*dest> is (probably) reallocated, 
 *           modified, and nul-terminated.
 *           
 * Throws:   <eslEMEM> on allocation failure; initial state of <dest> 
 *           is unaffected.
 */
int
esl_strcat(char **dest, int64_t ldest, const char *src, int64_t lsrc)
{
  void     *p;
  int       status;
  int64_t   len1, len2;

  if (ldest < 0) len1 = ((*dest == NULL) ? 0 : strlen(*dest));
  else           len1 = ldest;

  if (lsrc < 0)  len2 = ((  src == NULL) ? 0 : strlen(src)); 
  else           len2 = lsrc;

  if (len2 == 0) return eslOK;

  if (*dest == NULL) ESL_ALLOC(*dest, sizeof(char) * (len2+1));
  else               ESL_RALLOC(*dest, p, sizeof(char) * (len1+len2+1));

  memcpy((*dest)+len1, src, len2+1);
  return eslOK;

 ERROR:
  return status;
}


/* Function: esl_strtok()
 * Synopsis: Threadsafe version of C's <strtok()>
 * Date:     SRE, Wed May 19 16:30:20 1999 [St. Louis]
 *
 * Purpose:  Thread-safe version of <strtok()> for parsing next token in
 *           a string.
 *          
 *           Increments <*s> while <**s> is a character in <delim>,
 *           then stops; the first non-<delim> character defines the
 *           beginning of a token. Increments <*s> until it reaches
 *           the next delim character (or \verb+\0+); this defines the end
 *           of the token, and this character is replaced with
 *           \verb+\0+. <*s> is then reset to point to the next character
 *           after the \verb+\0+ that was written, so successive calls can
 *           extract tokens in succession. Sets <*ret_tok> to point at
 *           the beginning of the token, and returns <eslOK>.
 *
 *           If a token is not found -- if <*s> already points to
 *           \verb+\0+, or to a string composed entirely of characters in
 *           <delim> -- then returns <eslEOL>, with <*ret_tok> set to
 *           <NULL>.
 *           
 *           <*s> cannot be a constant string, since we write \verb+\0+'s
 *           to it; caller must be willing to have this string
 *           modified. And since we walk <*s> through the string as we
 *           parse, the caller wants to use a tmp pointer <*s>, not
 *           the original string itself.
 *                      
 * Example:  
 *           char *tok;
 *           char *s;             
 *           char  buf[50] = "This is  a sentence.";
 *           
 *           s = buf;  
 *           esl_strtok(&s, " ", &tok);
 *                tok is "This"; s is "is  a sentence."
 *           esl_strtok(&s, " ", &tok);
 *                tok is "is"; s is " a sentence.".
 *           esl_strtok(&s, " ", &tok);
 *                tok is "a"; s is "sentence.".
 *           esl_strtok(&s, " ", &tok, &len);
 *                tok is "sentence."; s is "\0".
 *           esl_strtok(&s, " ", &tok, &len);
 *                returned eslEOL; tok is NULL; s is "\0".
 *       
 * Args:     s        - a tmp, modifiable ptr to a string
 *           delim    - characters that delimits tokens
 *           ret_tok  - RETURN: ptr to \0-terminated token 
 *
 * Returns:  <eslOK> on success, <*ret_tok> points to next token, and
 *           <*s> points to next character following the token. 
 *
 *           Returns <eslEOL> on end of line; in which case <*s>
 *           points to the terminal \verb+\0+ on the line, and <*ret_tok>
 *           is <NULL>.
 */
int
esl_strtok(char **s, char *delim, char **ret_tok)
{
  return esl_strtok_adv(s, delim, ret_tok, NULL, NULL);
}


/* Function: esl_strtok_adv()
 * Synopsis: More advanced interface to <esl_strtok()>
 * Date:     SRE, Mon Oct 13 10:16:26 2008
 *
 * Purpose:  Same as <esl_strtok()>, except the caller may also 
 *           optionally retrieve the length of the token in <*opt_toklen>,
 *           and the token-ending character that was replaced by \verb+\0+ 
 *           in <*opt_endchar>. 
 *           
 * Args:     s           - a tmp, modifiable ptr to string
 *           delim       - characters that delimits tokens
 *           ret_tok     - RETURN: ptr to \0-terminated token string
 *           opt_toklen  - optRETURN: length of token; pass NULL if not wanted
 *           opt_endchar - optRETURN: character that was replaced by <\0>.
 *
 * Returns:  <eslOK> on success, <*ret_tok> points to next token, <*s>
 *           points to next character following the token,
 *           <*opt_toklen> is the <strlen()> length of the token in
 *           characters (excluding its terminal \verb+\0+), and <*opt_endchar>
 *           is the character that got replaced by \verb+\0+ to form the token.
 *           
 *           Returns <eslEOL> if no token is found (end of line); in
 *           which case <*s> points to the terminal \verb+\0+ on the line,
 *           <*ret_tok> is <NULL>, <*opt_toklen> is 0 and
 *           <*opt_endchar> is \verb+\0+.
 */
int
esl_strtok_adv(char **s, char *delim, char **ret_tok, int *opt_toklen, char *opt_endchar)
{
  char *end;
  char *tok    = *s;
  char  c      = '\0';
  int   n      = 0;
  int   status = eslEOL;  /* unless proven otherwise */

  tok += strspn(tok, delim);
  if (! *tok) tok = NULL;         /* if *tok = 0, EOL, no token left */
  else
    {
      n    = strcspn(tok, delim);
      end  = tok + n;
      if (*end == '\0') *s = end; /* a final token that extends to end of string */
      else 
	{
	  c     = *end;		 /* internal token: terminate with \0 */
	  *end  = '\0';
	  *s    = end+1;
	}
      status = eslOK;
    }

  *ret_tok = tok;
  if (opt_toklen  != NULL) *opt_toklen  = n;
  if (opt_endchar != NULL) *opt_endchar = c;
  return status;
}


/* Function:  esl_sprintf()
 * Synopsis:  Dynamic allocation version of sprintf().
 * Incept:    SRE, Mon Oct 20 09:35:57 2008 [Janelia]
 *
 * Purpose:   Like ANSI C's <sprintf()>, except the string
 *            result is dynamically allocated, and returned
 *            through <*ret_s>. 
 *
 *            Caller is responsible for free'ing <*ret_s>.
 *            
 *            As a special case to facilitate some optional string
 *            initializations, if <format> is <NULL>, <*ret_s> is set
 *            to <NULL>.
 *
 * Returns:   <eslOK> on success, and <*ret_s> is the resulting
 *            string.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_sprintf(char **ret_s, const char *format, ...)
{
  va_list ap;
  int     status;

  va_start(ap, format);
  status = esl_vsprintf(ret_s, format, &ap);
  va_end(ap);
  return status;
}

/* Function:  esl_vsprintf()
 * Synopsis:  Dynamic allocation version of vsprintf()
 * Incept:    SRE, Wed Oct 22 14:48:44 2008 [Janelia]
 *
 * Purpose:   Like ANSI C's <vsprintf>, except the string
 *            result is dynamically allocated, and returned
 *            through <*ret_s>.
 *
 *            Caller is responsible for free'ing <*ret_s>.
 *            
 *            As a special case to facilitate some optional string
 *            initializations, if <format> is <NULL>, <*ret_s> is set
 *            to <NULL>.
 *            
 * Returns:   <eslOK> on success, and <*ret_s> is the resulting
 *            string.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *            <eslESYS> if a <*printf()> library call fails.
 */
int
esl_vsprintf(char **ret_s, const char *format, va_list *ap)
{
  char   *s = NULL;
  void   *p = NULL;
  va_list ap2;
  int     n1, n2;
  int     status;

  if (format == NULL) { *ret_s = NULL; return eslOK; }

  va_copy(ap2, *ap);
  n1 = strlen(format) * 2;	/* initial guess at string size needed */
  ESL_ALLOC(s, sizeof(char) * (n1+1));
  if ((n2 = vsnprintf(s, n1+1, format, *ap)) >= n1) 
    {
      ESL_RALLOC(s, p, sizeof(char) * (n2+1));
      if (vsnprintf(s, n2+1, format, ap2) == -1) ESL_EXCEPTION(eslESYS, "vsnprintf() failed");
    }
  else if (n2 == -1) ESL_EXCEPTION(eslESYS, "vsnprintf() failed");

  va_end(ap2);
  *ret_s = s;
  return eslOK;

 ERROR:
  if (s != NULL) free(s);
  va_end(ap2);
  *ret_s = NULL;
  return status;
}
   

/* Function:  esl_strcmp()
 * Synopsis:  a strcmp() that treats NULL as empty string.
 * Incept:    SRE, Wed Jan 21 14:11:48 2009 [Janelia]
 *
 * Purpose:   A version of <strcmp()> that accepts <NULL>
 *            strings. If both <s1> and <s2> are non-<NULL>
 *            they are compared by <strcmp()>. If both are
 *            <NULL>, return 0 (as if they are identical
 *            strings). If only <s1> (or <s2>) is non-<NULL>,
 *            return 1 (or -1), corresponding to ordering
 *            any non-<NULL> string as greater than a <NULL>
 *            string.
 *
 *            (Easel routinely uses NULL to mean an unset optional
 *            string, and often needs to compare two strings for
 *            equality.)
 *
 * Returns:   0 if <s1 == s2>; 1 if <s1 > s2>; -1 if <s1 < s2>.
 */
int 
esl_strcmp(const char *s1, const char *s2)
{
  if      (s1 && s2) return strcmp(s1, s2);
  else if (s1)       return 1;
  else if (s2)       return -1;
  else               return 0;
}



/*****************************************************************
 * Easel's optional replacements for common but non-ANSI C functions.
 * These alternatives are only compiled in when we need them,
 * and their inclusion is controlled by #define's in easel.h.
 *     strcasecmp() -> may be define'd to be esl_strcasecmp()
 */

#ifndef HAVE_STRCASECMP
/* Function:  esl_strcasecmp()
 * Incept:    SRE, Sat Dec 10 09:44:13 2005 [St. Louis]
 *
 * Purpose:   Compare strings <s1> and <s2>. Return -1 if 
 *            <s1> is alphabetically less than <s2>, 0 if they
 *            match, and 1 if <s1> is alphabetically greater
 *            than <s2>. All matching is case-insensitive.
 *
 * Args:      s1  - string 1, \0 terminated
 *            s2  - string 2, \0 terminated      
 *
 * Returns:   -1, 0, or 1, if <s1> is less than, equal, or 
 *            greater than <s2>, case-insensitively.
 *
 * Throws:    (no abnormal error conditions)
 */
int
esl_strcasecmp(const char *s1, const char *s2)
{
  int i, c1, c2;

  for (i = 0; s1[i] != '\0' && s2[i] != '\0'; i++)
    {
      c1 = s1[i];	/* total paranoia. don't trust toupper() to    */
      c2 = s2[i];       /* leave the original unmodified; make a copy. */
  
      if (islower(c1)) c1 = toupper(c1);
      if (islower(c2)) c2 = toupper(c2);
      
      if      (c1 < c2) return -1;
      else if (c1 > c2) return 1;
    }

  if      (s1[i] != '\0') return 1;   /* prefixes match, but s1 is longer */
  else if (s2[i] != '\0') return -1;  /* prefixes match, s2 is longer */

  return 0;  /* else, a case-insensitive match. */
}
#endif /* ! HAVE_STRCASECMP */


/*****************************************************************
 * and some extra str*() functions...
 *****************************************************************/ 

/* Function:  esl_strchop()
 * Incept:    SRE, Mon Apr  3 10:24:14 2006 [St. Louis]
 *
 * Purpose:   Chops trailing whitespace off of a string <s> (or if <s>
 *            is NULL, do nothing).
 *            <n> is the length of the input string, if known; or pass <n=-1>
 *            if length is unknown. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      from squid's StringChop().
 */
int
esl_strchop(char *s, int64_t n)
{
  int i;
  if (s == NULL) return eslOK;
  if (n < 0) n = strlen(s);
  for (i = n-1; i>=0 && isspace((int) s[i]); i--); 
  s[i+1] = '\0';
  return eslOK;
}


/* Function:  esl_strdealign()
 * Synopsis:  Dealign a string according to gaps in a reference aseq.
 * Incept:    SRE, Thu Feb 17 15:12:26 2005 [St. Louis]
 *
 * Purpose:   Dealign string <s> in place, by removing any characters 
 *            aligned to gaps in <aseq>. Gap characters are defined in the 
 *            string <gapstring>; for example, <-_.>. Optionally return the
 *            unaligned length of <s> in characters in <*opt_rlen>.
 *            
 *            By providing a reference <aseq> to dealign against, this
 *            function can dealign aligned annotation strings, such as
 *            secondary structure or surface accessibility strings.
 *            If <s> is the same as <aseq>, then the aligned sequence
 *            itself is dealigned in place. 
 *            
 *            To dealign both annotations and sequence, do the
 *            sequence last, since you need it as the reference <aseq>
 *            when doing the annotations.
 *           
 *            It is safe to pass a <NULL> <s> (an unset optional
 *            annotation), in which case the function no-ops and
 *            returns <eslOK>.
 *            
 * Args:      s        - string to dealign
 *            aseq     - reference aligned sequence seq
 *            gapchars - definition of gap characters ("-_." for example)
 *            opt_rlen - optRETURN: number of residues in <s> after dealign
 *
 * Returns:   <eslOK> on success.
 */
int
esl_strdealign(char *s, const char *aseq, const char *gapchars, int64_t *opt_rlen)
{
  int64_t n = 0;
  int64_t apos;

  if (s == NULL) return eslOK;

  for (apos = 0; aseq[apos] != '\0'; apos++)
    if (strchr(gapchars, aseq[apos]) == NULL)
      s[n++] = s[apos];
  s[n] = '\0';
  
  if (opt_rlen != NULL) *opt_rlen = n;
  return eslOK;
}





/*----------------- end, C library replacements  -------------------------*/




/******************************************************************************
 * 5. File path/name manipulation functions, including tmpfiles                             
 *                                                                      
 * Sufficiently widespread in the modules that we make them core.       
 * (Should be moved to their own module eventually)                     
 *****************************************************************************/

/* Function:  esl_FileExists()
 * Incept:    SRE, Sat Jan 22 09:07:24 2005 [St. Louis]
 *
 * Purpose:   Returns TRUE if <filename> exists, else FALSE.
 *     
 * Note:      Testing a read-only fopen() is the only portable ANSI C     
 *            I'm aware of. We could also use a POSIX func here, since
 *            we have a ESL_POSIX_AUGMENTATION flag in the code.
 *            
 * Xref:      squid's FileExists().
 */
int
esl_FileExists(const char *filename)
{
  FILE *fp;
  if ((fp = fopen(filename, "r"))) { fclose(fp); return TRUE; }
  return FALSE;
}


/* Function:  esl_FileTail()
 * Incept:    SRE, Tue Mar  7 08:30:00 2006 [St. Louis]
 *
 * Purpose:   Given a full pathname <path>, extract the filename
 *            without the directory path; return it via  
 *            <ret_filename>. <ret_filename> space is allocated
 *            here, and must be free'd by the caller.
 *            For example: 
 *               </foo/bar/baz.1> becomes <baz.1>;
 *               <foo/bar>        becomes <bar>; 
 *               <foo>            becomes <foo>; and
 *               </>              becomes the empty string.
 *            
 *            If <nosuffix> is <TRUE>, the rightmost trailing ".foo" extension
 *            is removed too. The suffix is defined as everything following
 *            the rightmost period in the filename in <path>:
 *            with <nosuffix> <TRUE>, 
 *                <foo.2/bar.idx> becomes <bar>,
 *                <foo.2/bar>     becomes <bar>, and
 *                <foo.2/bar.1.3> becomes <bar.1>.  
 *            
 * Args:      path     - full pathname to process, "/foo/bar/baz.1"
 *            nosuffix - TRUE to remove rightmost suffix from the filename
 *            ret_file - RETURN: filename portion of the path.
 *                     
 * Returns:   <eslOK> on success, and <ret_file> points to a newly
 *            allocated string containing the filename.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
esl_FileTail(const char *path, int nosuffix, char **ret_file)
{
  int   status;
  char *tail = NULL;
  char *lastslash;
  char *lastdot;
				/* remove directory prefix */
  lastslash = strrchr(path, eslDIRSLASH);
  ESL_ALLOC(tail, sizeof(char) * (strlen(path)+1)); /* a little overkill */
  if (lastslash == NULL) strcpy(tail, path);
  else                   strcpy(tail, lastslash+1);
				/* remove trailing suffix */
  if (nosuffix) {
    if ((lastdot = strrchr(tail, '.')) != NULL)
      *lastdot = '\0';
  }
  *ret_file = tail;
  return eslOK;

 ERROR:
  if (tail != NULL) free(tail);
  *ret_file = NULL;
  return status;
}

/* Function:  esl_FileConcat()
 * Incept:    SRE, Sat Jan 22 07:28:46 2005 [St. Louis]
 *
 * Purpose:   Concatenates directory path prefix <dir> and a filename
 *            <file>, and returns the new full pathname through
 *            <ret_path>. If <dir> does not already end in the
 *            appropriate delimiter (e.g. / for UNIX), one is added.
 *            
 *            If <dir> is NULL, then <ret_path> is just the same as
 *            <file>. Similarly, if <file> already appears to be a
 *            full path (because its first character is a /), then
 *            <dir> is ignored and <ret_path> is the same as
 *            <file>. It wouldn't normally make sense for a caller to
 *            call this function with such arguments.
 *            
 *            <file> may be a relative path. For example, 
 *            if <dir> is "/usr/local" and <file> is "lib/myapp/data",
 *            <ret_path> will be "/usr/local/lib/myapp/data".
 *
 * Returns:   <eslOK> on success, and puts the path
 *            in <ret_path>; this string is allocated here, 
 *            and must be free'd by caller with <free()>.
 *
 * Throws:    <eslEMEM>   on allocation failure.
 *            <eslEINVAL> on bad argument.
 *            In either case, <ret_path> is returned NULL.
 *
 * Xref:      squid's FileConcat().
 */
int
esl_FileConcat(const char *dir, const char *file, char **ret_path)
{
  char *path = NULL;
  int   nd, nf;
  int   status;

  if (ret_path != NULL) *ret_path = NULL;
  if (file == NULL) ESL_EXCEPTION(eslEINVAL, "null file");

  nd   = (dir  != NULL)? strlen(dir)  : 0;
  nf   = strlen(file);
  ESL_ALLOC(path, sizeof(char) * (nd+nf+2));

  if (dir == NULL)		     /* 1. silly caller didn't give a path */
    strcpy(path, file);
  else if (*file == eslDIRSLASH)     /* 2. <file> is already a path?   */
    strcpy(path, file); 
  else if (dir[nd-1] == eslDIRSLASH) /* 3. <dir><file> (dir is / terminated) */
    sprintf(path, "%s%s", dir, file);
  else				     /* 4. <dir>/<file> (usual case)   */
    sprintf(path, "%s%c%s", dir, eslDIRSLASH, file);	

  *ret_path = path;
  return eslOK;

 ERROR:
  if (path     != NULL) free(path);
  if (ret_path != NULL) *ret_path = NULL;
  return status;
}


/* Function:  esl_FileNewSuffix()
 * Incept:    SRE, Sat Jan 22 10:04:08 2005 [St. Louis]
 *
 * Purpose:   Add a file suffix <sfx> to <filename>; or if <filename>
 *            already has a suffix, replace it with <sfx>. A suffix is
 *            usually 2-4 letters following a '.' character. Returns
 *            an allocated string containing the result in <ret_newpath>.
 *            
 *            For example, if <filename> is "foo" and <sfx> is "ssi",
 *            returns "foo.ssi". If <filename> is "foo.db" and <sfx>
 *            is "idx", returns "foo.idx". 
 *
 * Returns:   <eslOK> on success, and <ret_newpath> is set
 *            string "<base_filename>.<sfx>". Caller must <free()>
 *            this string.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      squid's FileAddSuffix().
 */
int 
esl_FileNewSuffix(const char *filename, const char *sfx, char **ret_newpath)
{
  char *new = NULL;
  char *lastdot;
  int   nf;
  int   status;

  if (ret_newpath != NULL) *ret_newpath = NULL;

  lastdot   = strrchr(filename, '.'); /* check for suffix to replace */
  if (lastdot != NULL && 
      strchr(lastdot, eslDIRSLASH) != NULL) 
    lastdot = NULL; /*foo.1/filename case - don't be fooled.*/
  nf = (lastdot == NULL)? strlen(filename) : lastdot-filename;
  
  ESL_ALLOC(new, sizeof(char) * (nf+strlen(sfx)+2)); /* '.' too */
  strncpy(new, filename, nf);
  *(new+nf) = '.';
  strcpy(new+nf+1, sfx);

  if (ret_newpath != NULL) *ret_newpath = new; else free(new);
  return eslOK;

 ERROR:
  if (new         != NULL) free(new);
  if (ret_newpath != NULL) *ret_newpath = NULL;
  return status;
}



/* Function:  esl_FileEnvOpen()
 * Incept:    SRE, Sat Jan 22 08:41:48 2005 [St. Louis]
 *
 * Purpose:   Looks for a file <fname> in a colon-separated list of
 *            directories that is configured in an environment variable
 *            <env>. The first occurrence of file <fname> in this directory 
 *            list is opened read-only. The open file ptr is returned
 *            through <opt_fp>, and the full path name to the file
 *            that was opened is returned through <opt_path>. 
 *            Caller can pass NULL in place of <opt_fp> or <opt_path>
 *            if it is not interested in one or both of these. 
 *            
 *            Does not look in the current directory unless "." is
 *            explicitly in the directory list provided by <env>.
 *            
 * Note:      One reason to pass <opt_path> back to the caller is that
 *            sometimes we're opening the first in a group of files
 *            (for example, a database and its SSI index), and we want
 *            to make sure that after we find the main file, the
 *            caller can look for the auxiliary file(s) in exactly the
 *            same directory.
 *            
 * Examples:  % setenv BLASTDB /nfs/databases/blast-db:/nfs/databases/nr/
 *           
 *            FILE *fp;
 *            char *path;
 *            int   status;
 *            status = esl_FileEnvOpen("swiss42", "BLASTDB", &fp, &path);
 * 
 * Returns:   <eslOK> on success, and provides <opt_fp> and <opt_path>;
 *            <opt_fp> is opened here, and must be <fclose()>'d by caller;
 *            <opt_path> is allocated here, and must be <free()>'d by caller.
 *
 *            Returns <eslENOTFOUND> if the file not found in any directory,
 *            or if <env> does not contain any directories to look in.
 *
 * Throws:    <eslEMEM> on allocation error.
 *
 * Xref:      squid's EnvFileOpen().
 */
int
esl_FileEnvOpen(const char *fname, const char *env, FILE **opt_fp, char **opt_path)
{
  FILE *fp;
  char *dirlist;		/* :-separated list of directories */
  char *s, *s2;                 /* ptrs into elems in env list */
  char *path = NULL;
  int   np;
  int   status;

  fp = NULL;
  if (opt_fp   != NULL) *opt_fp   = NULL;
  if (opt_path != NULL) *opt_path = NULL;

  if (env == NULL)               return eslENOTFOUND;
  if ((s = getenv(env)) == NULL) return eslENOTFOUND;
  if (esl_strdup(s, -1, &dirlist) != eslOK) return eslEMEM;

  np   = strlen(fname) + strlen(s) + 2; /* upper bound on full path len */
  ESL_ALLOC(path, sizeof(char) * np);

  s  = dirlist;
  while (s != NULL) 
    {
      if ((s2 = strchr(s, ':')) != NULL) { *s2 = '\0'; s2++;} /* ~=strtok() */
      sprintf(path, "%s%c%s", s, eslDIRSLASH, fname); /* // won't hurt */
      if ((fp = fopen(path, "r")) != NULL) break;      
      s = s2;
    }
  if (fp == NULL) { free(path); free(dirlist); return eslENOTFOUND; }

  if (opt_path != NULL) { *opt_path = path; } else free(path);
  if (opt_fp   != NULL) { *opt_fp   = fp; }   else fclose(fp);
  free(dirlist);
  return eslOK;

 ERROR:
  if (path     != NULL) free(path);
  if (fp       != NULL) fclose(fp);
  if (dirlist  != NULL) free(dirlist);
  if (opt_path != NULL) *opt_path = NULL;
  if (opt_fp   != NULL) *opt_fp   = NULL;
  return status;
}

/* Function:  esl_tmpfile()
 * Incept:    SRE, Wed Sep  6 08:15:15 2006 [Janelia]
 *
 * Purpose:   Open a secure temporary <FILE *> handle and return it in
 *            <ret_fp>. The file is opened in read-write mode (<w+b>)
 *            with permissions 0600, as an atomic operation using the
 *            POSIX <mkstemp()> function.
 * 
 *            The <basename6X> argument is a modifiable string that must
 *            end in "XXXXXX" (for example, "esltmpXXXXXX"). The
 *            <basename6X> is used to construct a unique tmpfile name.
 *            
 *            Note that this string must be modifiable; do not declare
 *            it <char *tmpfile = "esltmpXXXXXX";> nor <char tmpfile[]
 *            = "esltmpXXXXXX";> because these will not work on some
 *            compilers. Something like <char tmpfile[16] =
 *            "esltmpXXXXXX";> that explicitly allocates storage will
 *            suffice.
 *            
 *            The file is opened in a standard temporary file
 *            directory. The path is obtained from the environment
 *            variable <TMPDIR>; failing that, from the environment
 *            variable <TMP>; and failing that, </tmp> is used. If the
 *            process is running <setuid> or <setgid>, then the
 *            environment variables are ignored, and the temp file is
 *            always created in </tmp>.
 *            
 *            The created tmpfile is not persistent and is not visible
 *            to a directory listing. The caller may <rewind()> the
 *            <ret_fp> and do cycles of reading and/or writing, but
 *            once the <ret_fp> is closed, the file disappears.  The
 *            caller does not need to <remove()> or <unlink()> it (and
 *            in fact, cannot do so, because it does not know the
 *            tmpfile's name).
 *            
 *            This function is a secure replacement for ANSI C
 *            <tmpfile()>, which is said to be insecurely implemented on
 *            some platforms.
 *
 * Returns:   <eslOK> on success, and now <ret_fp> points to a new <FILE *>
 *            stream for the opened tempfile. 
 *
 * Throws:    <eslESYS> if a system call (including the <mkstemp()> call)
 *            fails, and and <ret_fp> is returned NULL. One possible
 *            problem is if the temporary directory doesn't exist or
 *            is not writable. This is considered to be a system
 *            error, not a user error, so Easel handles it as an exception.
 *            
 * Xref:      STL11/85. Substantially copied from David Wheeler, 
 *            "Secure Programming for Linux and Unix HOWTO", 
 *            http://www.dwheeler.com/secure-programs/Secure-Programs-HOWTO/introduction.html.
 *            Copyright (C) 1999-2001 David A. Wheeler.
 *            Licensed under the MIT license; see Appendix C of the HOWTO.
 *            Thanks, David, for the clearest explanation of the issues 
 *            that I've seen.
 *            
 *            I also referred to H. Chen, D. Dean, and D. Wagner,
 *            "Model checking one million lines of C code", 
 *            In: Network and Distributed System Security Symposium, pp 171-185,
 *            San Diego, CA, February 2004;
 *            http://www.cs.ucdavis.edu/~hchen/paper/ndss04.pdf.
 *            Wheeler's implementation obeys Chen et al's "Property 5", 
 *            governing secure use of tempfiles.
 */
int
esl_tmpfile(char *basename6X, FILE **ret_fp)
{
  char *tmpdir = NULL;
  char *path   = NULL;
  FILE *fp     = NULL;
  int   fd;
  int   status;
  mode_t old_mode;

  /* Determine what tmp directory to use, and construct the
   * file name.
   */
  if (getuid() == geteuid() && getgid() == getegid()) 
    {
      tmpdir = getenv("TMPDIR");
      if (tmpdir == NULL) tmpdir = getenv("TMP");
    }
  if (tmpdir == NULL) tmpdir = "/tmp";
  if ((status = esl_FileConcat(tmpdir, basename6X, &path)) != eslOK) goto ERROR; 

  old_mode = umask(077);
  if ((fd = mkstemp(path)) <  0)        ESL_XEXCEPTION(eslESYS, "mkstemp() failed.");
  umask(old_mode);
  if ((fp = fdopen(fd, "w+b")) == NULL) ESL_XEXCEPTION(eslESYS, "fdopen() failed.");
  if (unlink(path) < 0)                 ESL_XEXCEPTION(eslESYS, "unlink() failed.");

  *ret_fp = fp;
  free(path);
  return eslOK;

 ERROR:
  if (path != NULL) free(path);
  if (fp   != NULL) fclose(fp);
  *ret_fp = NULL;
  return status;
}

/* Function:  esl_tmpfile_named()
 * Incept:    SRE, Sat Nov 11 09:13:25 2006 [Janelia]
 *
 * Purpose:   Open a persistent temporary file relative to the current
 *            working directory. The file name is constructed from the
 *            <basename6X> argument, which must be a modifiable string
 *            ending in the six characters "XXXXXX".  These are
 *            replaced by a unique character string by a call to POSIX
 *            <mkstemp()>. For example, <basename6X> might be
 *            <esltmpXXXXXX> on input, and <esltmp12ab34> on return; or, to
 *            put the tmp file in a subdirectory under the current
 *            working directory, something like <my_subdir/esltmpXXXXXX>
 *            on input resulting in something like
 *            <my_subdir/esltmp12ab34> on return.  The tmpfile is opened
 *            for reading and writing (in mode <w+b> with permissions
 *            0600) and the opened <FILE *> handle is returned through
 *            <ret_fp>.
 *            
 *            The created tmpfile is persistent: it will be visible in
 *            a directory listing, and will remain after program
 *            termination unless the caller explicitly removes it by a
 *            <remove()> or <unlink()> call.
 *
 *            To use this function securely, if you reopen the
 *            tmpfile, you must only reopen it for reading, not
 *            writing, and you must not trust the contents.
 *            
 *            Because the <basename6X> will be modified, it cannot be
 *            a string constant (especially on a picky compiler like
 *            gcc). You have to declare it with something like
 *               <char tmpfile[32] = "esltmpXXXXXX";> 
 *            not 
 *               <char *tmpfile    = "esltmpXXXXXX";> 
 *            because a compiler is allowed to make the <*tmpfile> version
 *            a constant.
 *
 * Returns:   <eslOK> on success, <basename6X> contains the name of the
 *            tmpfile, and <ret_fp> contains a new <FILE *> stream for the
 *            opened file. 
 *             
 *            <eslFAIL> on failure, and <ret_fp> is returned NULL and
 *            the contents of <basename6X> are undefined. The most
 *            common reason for a failure will be that the caller does
 *            not have write permission for the directory that
 *            <basename6X> is in. Easel handles this as a normal (user)
 *            failure, not an exception, because these permissions are
 *            most likely in the user's control (in contrast to
 *            <esl_tmpfile()>, which always uses a system <TMPDIR>
 *            that should always be user-writable on a properly
 *            configured POSIX system).
 *
 * Xref:      STL11/85.
 */
int
esl_tmpfile_named(char *basename6X, FILE **ret_fp)
{
  FILE  *fp;
  mode_t old_mode;
  int    fd;

  *ret_fp = NULL;
  old_mode = umask(077);
  if ((fd = mkstemp(basename6X)) <  0)    return eslFAIL;
  umask(old_mode);
  if ((fp = fdopen(fd, "w+b")) == NULL) return eslFAIL;

  *ret_fp = fp;
  return eslOK;
}


/*----------------- end of file path/name functions ------------------------*/




/*****************************************************************
 * 6. Typed comparison routines.
 *****************************************************************/

/* Function:  esl_DCompare()
 * Incept:    SRE, Mon Nov  6 10:11:47 2006 [Janelia]
 *
 * Purpose:   Compare two floating point scalars <a> and <b> for approximate equality.
 *            Return <eslOK> if equal, <eslFAIL> if not.
 *            
 *            Equality is defined by being within a relative
 *            epsilon <tol>, as <2*fabs(a-b)/(a+b)> $\leq$ <tol>.
 *            Additionally, we catch the special cases where <a>
 *            and/or <b> are 0 or -0. If both are, return <eslOK>; if
 *            one is, check that the absolute value of the other is
 *            $\leq$ <tol>.
 *            
 *            <esl_DCompare()> and <esl_FCompare()> work on <double> and <float>
 *            scalars, respectively.
 */
int
esl_DCompare(double a, double b, double tol)
{
  if (isinf(a) && isinf(b))                 return eslOK;
  if (isnan(a) && isnan(b))                 return eslOK;
  if (!isfinite(a) || !isfinite(b))         return eslFAIL;
  if (a == b)                               return eslOK;
  if (fabs(a) == 0. && fabs(b) <= tol)      return eslOK;
  if (fabs(b) == 0. && fabs(a) <= tol)      return eslOK;
  if (2.*fabs(a-b) / fabs(a+b) <= tol)      return eslOK;
  return eslFAIL;
}
int
esl_FCompare(float a, float b, float tol)
{ 
  if (isinf(a) && isinf(b))                 return eslOK;
  if (isnan(a) && isnan(b))                 return eslOK;
  if (!isfinite(a) || !isfinite(b))         return eslFAIL;
  if (a == b)                               return eslOK;
  if (fabs(a) == 0. && fabs(b) <= tol)      return eslOK;
  if (fabs(b) == 0. && fabs(a) <= tol)      return eslOK;
  if (2.*fabs(a-b) / fabs(a+b) <= tol)      return eslOK;
  return eslFAIL;
}

/* Function:  esl_CCompare()
 * Synopsis:  Compare two optional strings for equality.
 * Incept:    SRE, Wed Jun 13 10:25:06 2007 [Janelia]
 *
 * Purpose:   Compare two optional strings <s1> and <s2>
 *            for equality. 
 *            
 *            If they're non-<NULL> and identical up to their
 *            <NUL>-terminator, return <eslOK>.
 *            
 *            If they're both <NULL> (unset), return <eslOK>.
 *            
 *            Otherwise, they're not identical; return <eslFAIL>.
 */
int
esl_CCompare(char *s1, char *s2)
{
  if (s1 == NULL && s2 == NULL) return eslOK;
  if (s1 == NULL || s2 == NULL) return eslFAIL;
  if (strcmp(s1, s2) != 0)      return eslFAIL;
  return eslOK;
}


/*-------------- end, typed comparison routines --------------------*/




/*****************************************************************
 * 7. Commonly used background composition (iid) frequencies. 
 *****************************************************************/

/* Function:  esl_composition_BL62()
 * Incept:    SRE, Fri Apr 13 16:00:34 2007 [Janelia]
 *
 * Purpose:   Sets <f> to the background frequencies used in
 *            \citep{Henikoff92} to calculate the BLOSUM62
 *            substitution matrix. Caller provides space in <f>
 *            allocated for at least 20 doubles.  The entries are in
 *            alphabetic order A..Y, same as the standard Easel amino
 *            acid alphabet order.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_composition_BL62(double *f)
{
  f[0]  = 0.074;
  f[1]  = 0.025;
  f[2]  = 0.054;
  f[3]  = 0.054;
  f[4]  = 0.047;
  f[5]  = 0.074;
  f[6]  = 0.026;
  f[7]  = 0.068;
  f[8]  = 0.058;
  f[9]  = 0.099;
  f[10] = 0.025;
  f[11] = 0.045;
  f[12] = 0.039;
  f[13] = 0.034;
  f[14] = 0.052;
  f[15] = 0.057;
  f[16] = 0.051;
  f[17] = 0.073;
  f[18] = 0.013;
  f[19] = 0.032;
  return eslOK;
}

/* Function:  esl_composition_WAG()
 * Incept:    SRE, Fri Apr 13 16:02:48 2007 [Janelia]
 *
 * Purpose:   Sets <f> to the background frequencies used in
 *            \citep{WhelanGoldman01} to calculate the WAG rate
 *            matrix. Caller provides space in <f> allocated for at
 *            least 20 doubles.  The entries are in alphabetic order
 *            A..Y, same as the standard Easel amino acid alphabet
 *            order.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_composition_WAG(double *f)
{
  f[0]  = 0.086628;                     /* A */
  f[1]  = 0.019308;	                /* C */
  f[2]  = 0.057045;	                /* D */
  f[3]  = 0.058059;	                /* E */
  f[4]  = 0.038432;	                /* F */
  f[5]  = 0.083252;	                /* G */
  f[6]  = 0.024431;	                /* H */
  f[7]  = 0.048466;	                /* I */
  f[8]  = 0.062029;	                /* K */
  f[9]  = 0.086209;	                /* L */
  f[10] = 0.019503;	                /* M */
  f[11] = 0.039089;	                /* N */
  f[12] = 0.045763;	                /* P */
  f[13] = 0.036728;	                /* Q */
  f[14] = 0.043972;	                /* R */
  f[15] = 0.069518;	                /* S */
  f[16] = 0.061013;	                /* T */
  f[17] = 0.070896;	                /* V */
  f[18] = 0.014386;	                /* W */
  f[19] = 0.035274;	                /* Y */
  return eslOK;
}

/* Function:  esl_composition_SW34()
 * Incept:    SRE, Fri Apr 13 16:03:46 2007 [Janelia]
 *
 * Purpose:   Sets <f> to the background frequencies observed in
 *            Swissprot release 34 (21.2M residues).  Caller provides
 *            space in <f> allocated for at least 20 doubles.  The
 *            entries are in alphabetic order A..Y, same as the
 *            standard Easel amino acid alphabet order.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_composition_SW34(double *f)
{
  f[0]  = 0.075520;                     /* A */
  f[1]  = 0.016973;                     /* C */
  f[2]  = 0.053029;                     /* D */
  f[3]  = 0.063204;                     /* E */
  f[4]  = 0.040762;                     /* F */
  f[5]  = 0.068448;                     /* G */
  f[6]  = 0.022406;                     /* H */
  f[7]  = 0.057284;                     /* I */
  f[8]  = 0.059398;                     /* K */
  f[9]  = 0.093399;                     /* L */
  f[10] = 0.023569;                     /* M */
  f[11] = 0.045293;                     /* N */
  f[12] = 0.049262;                     /* P */
  f[13] = 0.040231;                     /* Q */
  f[14] = 0.051573;                     /* R */
  f[15] = 0.072214;                     /* S */
  f[16] = 0.057454;                     /* T */
  f[17] = 0.065252;                     /* V */
  f[18] = 0.012513;                     /* W */
  f[19] = 0.031985;                     /* Y */
  return eslOK;
}


/* Function:  esl_composition_SW50()
 * Incept:    SRE, Tue Aug 26 08:42:04 2008 [Janelia]
 *
 * Purpose:   Sets <f> to the background frequencies observed in
 *            Swissprot release 50.8 (86.0M residues; Oct 2006).
 *
 * Returns:   <eslOK> on success.
 */
int
esl_composition_SW50(double *f)
{
  f[0] = 0.0787945;		/* A */
  f[1] = 0.0151600;		/* C */
  f[2] = 0.0535222;		/* D */
  f[3] = 0.0668298;		/* E */
  f[4] = 0.0397062;		/* F */
  f[5] = 0.0695071;		/* G */
  f[6] = 0.0229198;		/* H */
  f[7] = 0.0590092;		/* I */
  f[8] = 0.0594422;		/* K */
  f[9] = 0.0963728;		/* L */
  f[10]= 0.0237718;		/* M */
  f[11]= 0.0414386;		/* N */
  f[12]= 0.0482904;		/* P */
  f[13]= 0.0395639;		/* Q */
  f[14]= 0.0540978;		/* R */
  f[15]= 0.0683364;		/* S */
  f[16]= 0.0540687;		/* T */
  f[17]= 0.0673417;		/* V */
  f[18]= 0.0114135;		/* W */
  f[19]= 0.0304133;		/* Y */
  return eslOK;
}
/*-------------- end, background compositions -------------------*/



/*****************************************************************
 * 8. Unit tests.
 *****************************************************************/
#ifdef eslEASEL_TESTDRIVE

static void
utest_strtok(void)
{
  char *msg         = "esl_strtok() unit test failed";
  char *teststring;
  char *s;
  char *tok;
  int   toklen;
  char  endc;

  if (esl_strdup("This is\t a sentence.", -1, &teststring) != eslOK) esl_fatal(msg);

  s = teststring;
  if (esl_strtok(&s, " ", &tok) != eslOK)                            esl_fatal(msg);
  if (strcmp(tok, "This")       != 0)                                esl_fatal(msg);
  if (*s != 'i')                                                     esl_fatal(msg);
  
  if (esl_strtok_adv(&s, " \t", &tok, &toklen, &endc) != eslOK)      esl_fatal(msg);
  if (strcmp(tok, "is") != 0)                                        esl_fatal(msg);
  if (*s != ' ')                                                     esl_fatal(msg);
  if (toklen != 2)                                                   esl_fatal(msg);
  if (endc != '\t')                                                  esl_fatal(msg);
  
  if (esl_strtok_adv(&s, "\n", &tok, NULL, NULL) != eslOK)           esl_fatal(msg);
  if (strcmp(tok, " a sentence.") != 0)                              esl_fatal(msg);
  if (*s != '\0')                                                    esl_fatal(msg);

  free(teststring);
  return;
}
  
static void
utest_sprintf(void)
{
  char *msg  = "unit tests for esl_[v]sprintf() failed";
  int   num  = 99;
  char *what = "beer";
  char *s    = NULL;

  if (esl_sprintf(&s, "%d bottles of %s", num, what) != eslOK) esl_fatal(msg);
  if (strcmp(s, "99 bottles of beer")                != 0)     esl_fatal(msg);
  free(s); 

  if (esl_sprintf(&s, NULL)                          != eslOK) esl_fatal(msg);
  if (s                                              != NULL)  esl_fatal(msg);
}


static void
utest_tmpfile_named(void)
{
  char *msg          = "tmpfile_named unit test failed";
  char  tmpfile[32]  = "esltmpXXXXXX";
  FILE *fp           = NULL;
  char  buf[256];

  if (esl_tmpfile_named(tmpfile, &fp) != eslOK) esl_fatal(msg);
  fprintf(fp, "Unit test.\n");
  fclose(fp);
  if ((fp = fopen(tmpfile, "r"))   == NULL)  esl_fatal(msg);
  if (fgets(buf, 256, fp)          == NULL)  esl_fatal(msg);
  if (strcmp(buf, "Unit test.\n")  != 0)     esl_fatal(msg);
  fclose(fp);
  remove(tmpfile);
  return;
}

#endif /*eslEASEL_TESTDRIVE*/


/*****************************************************************
 * 9. Test driver.
 *****************************************************************/

#ifdef eslEASEL_TESTDRIVE
/* gcc -g -Wall -o test -I. -L. -DeslEASEL_TESTDRIVE easel.c -leasel -lm
 * ./test
 */
#include "easel.h"

int main(void)
{

#ifdef eslTEST_THROWING
  esl_exception_SetHandler(&esl_nonfatal_handler);
#endif

  utest_strtok();
  utest_sprintf();
  utest_tmpfile_named();
  return eslOK;
}
#endif /*eslEASEL_TESTDRIVE*/

/*****************************************************************
 * 10. Examples.
 *****************************************************************/

#ifdef eslEASEL_EXAMPLE2
/*::cexcerpt::easel_example_tmpfiles::begin::*/
/* gcc -g -Wall -o example -I. -L. -DeslEASEL_EXAMPLE_TMPFILES easel.c -leasel -lm
 * ./example
 */
#include "easel.h"

int main(void)
{
  char  tmpfile1[32]  = "esltmpXXXXXX"; /* a transient, secure tmpfile: 6 X's are important */
  char  tmpfile2[32]  = "esltmpXXXXXX"; /* a named tmpfile                                  */
  FILE *fp            = NULL;
  char  buf[256];

  /* Example of using a secure, unnamed tmpfile. 
   * Note, the new tmpfile is automatically deleted, so to cleanup, just fclose() the FILE */
  esl_tmpfile(tmpfile1, &fp);
  fprintf(fp, "Hello world!\n");
  rewind(fp);
  fgets(buf, 256, fp);
  printf("first temp file says: %s\n", buf);
  fclose(fp);

  /* Example of reasonably securely using a named tmpfile. 
   * To cleanup, must both fclose() the FILE and remove() the file by name */
  esl_tmpfile_named(tmpfile2, &fp);
  fprintf(fp, "Hello insecure world!\n");
  fclose(fp);		/* tmpfile2 now exists on disk and can be closed/reopened */

  fp = fopen(tmpfile2, "r");
  fgets(buf, 256, fp);
  printf("second temp file says: %s\n", buf);
  fclose(fp);
  remove(tmpfile2);	/* disk file cleanup necessary with this version. */

  return eslOK;
}
/*::cexcerpt::easel_example_tmpfiles::end::*/
#endif /*eslEASEL_EXAMPLE2*/


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/  
