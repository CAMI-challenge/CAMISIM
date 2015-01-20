/* esl_stopwatch.c
 * Tracking cpu/system/elapsed time used by a process.
 *
 * Thanks to Warren Gish for assistance.
 * 
 * SRE, Wed Feb 22 20:12:55 2006 [St. Louis] [moved to Easel]
 * SRE, Thu Aug  3 08:11:52 2000 [St. Louis] [moved to SQUID]
 * SRE, Fri Nov 26 14:54:21 1999 [St. Louis] [HMMER]
 * SVN $Id: esl_stopwatch.c 385 2009-08-23 14:17:56Z eddys $
 */
#include "esl_config.h"

#include "easel.h"
#include "esl_stopwatch.h"

/*****************************************************************
 * ESL_STOPWATCH object maintenance
 *****************************************************************/

/* Function:  esl_stopwatch_Create()
 * Incept:    SRE, Wed Feb 22 20:15:05 2006 [St. Louis]
 *
 * Purpose:   Creates a new stopwatch.
 *
 * Returns:   ptr to a new <ESL_STOPWATCH> object; caller is
 *            responsible for free'ing it with 
 *            <esl_stopwatch_Destroy()>.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_STOPWATCH *
esl_stopwatch_Create(void)
{
  int status;
  ESL_STOPWATCH *w = NULL;

  ESL_ALLOC(w, sizeof(ESL_STOPWATCH));
  w->elapsed = 0.;
  w->user    = 0.;
  w->sys     = 0.;
  return w;

 ERROR:
  return NULL;
}

/* Function:  esl_stopwatch_Destroy()
 * Incept:    SRE, Thu Feb 23 07:09:23 2006 [St. Louis]
 *
 * Purpose:   Frees an <ESL_STOPWATCH>.
 */
void
esl_stopwatch_Destroy(ESL_STOPWATCH *w)
{
  free(w);
}




/* Function:  esl_stopwatch_Start()
 * Incept:    SRE, Sat Feb 25 10:41:00 2006 [St. Louis]
 *
 * Purpose:   Start a stopwatch. This sets the base 
 *            for elapsed, cpu, and system time difference
 *            calculations by subsequent calls to
 *            <esl_stopwatch_Stop()>.
 *
 * Returns:   <eslOK> on success.
 */
int 
esl_stopwatch_Start(ESL_STOPWATCH *w)
{
#ifdef HAVE_TIMES /* POSIX */
  w->t0 = times(&(w->cpu0));
#else             /* fallback to ANSI C */
  w->t0   = time(NULL);
  w->cpu0 = clock();
#endif
  w->elapsed = 0.;
  w->user    = 0.;
  w->sys     = 0.;
  return eslOK;
}

/* Function:  esl_stopwatch_Stop()
 * Incept:    SRE, Sat Feb 25 10:42:26 2006 [St. Louis]
 *
 * Purpose:   Stop a stopwatch. Record and store elapsed,
 *            cpu, and system time difference relative to the
 *            last call to <esl_stopwatch_Start()>.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stopwatch_Stop(ESL_STOPWATCH *w)
{
#ifdef HAVE_TIMES
  struct tms cpu1;
  clock_t    t1;
  double     clk_tck;
#else
  time_t  t1;
  clock_t cpu1;
#endif


#ifdef HAVE_TIMES /* POSIX */
  t1         = times(&cpu1);
  
  clk_tck    = (double) sysconf(_SC_CLK_TCK);
  w->elapsed = (double) (t1 - w->t0) / clk_tck;
  w->user    = (double) (cpu1.tms_utime + cpu1.tms_cutime -
			 w->cpu0.tms_utime - w->cpu0.tms_cutime) / clk_tck;

  w->sys     = (double) (cpu1.tms_stime + cpu1.tms_cstime -
			 w->cpu0.tms_stime - w->cpu0.tms_cstime) / clk_tck;
#else /* fallback to ANSI C */
  t1         = time(NULL);
  cpu1       = clock();
  w->elapsed = difftime(t1, w->t0);
  w->user    = (double) (cpu1- w->cpu0) / (double) CLOCKS_PER_SEC;
  w->sys     = 0.;		/* no way to portably get system time in ANSI C */

#endif
  return eslOK;
}

/* format_time_string()
 * Date:     SRE, Fri Nov 26 15:06:28 1999 [St. Louis]
 *
 * Purpose:  Given a number of seconds, format into
 *           hh:mm:ss.xx in a provided buffer.
 *
 * Args:     buf     - allocated space (128 is plenty!)
 *           sec     - number of seconds
 *           do_frac - TRUE (1) to include hundredths of a sec
 */
static void
format_time_string(char *buf, double sec, int do_frac)
{
  int h, m, s, hs;
  
  h  = (int) (sec / 3600.);
  m  = (int) (sec / 60.) - h * 60;
  s  = (int) (sec) - h * 3600 - m * 60;
  if (do_frac) {
    hs = (int) (sec * 100.) - h * 360000 - m * 6000 - s * 100;
    sprintf(buf, "%02d:%02d:%02d.%02d", h,m,s,hs);
  } else {
    sprintf(buf, "%02d:%02d:%02d", h,m,s);
  }
}

/* Function:  esl_stopwatch_Display()
 * Incept:    SRE, Sat Feb 25 10:51:09 2006 [St. Louis]
 *
 * Purpose:   Output a usage summary line from a stopped
 *            stopwatch, showing elapsed, cpu, and system time
 *            between the last calls to 
 *            <esl_stopwatch_Start()> and <esl_stopwatch_Stop()>.
 *            
 *            The string <prefix> will be prepended to the output
 *            line. Use <""> to prepend nothing. If <prefix> is NULL,
 *            a default <"CPU Time: "> prefix is used.
 *           
 *            For <prefix> = <"CPU Time: "> an example output line is:\\
 *            <CPU Time: 142.55u 7.17s 00:02:29.72 Elapsed: 00:02:35>
 *
 * Args:      fp      - output stream
 *            w       - stopped stopwatch
 *            prefix  - output line prefix ("" for nothing)
 *
 * Returns:   <eslOK> on success.
 */
int 
esl_stopwatch_Display(FILE *fp, ESL_STOPWATCH *w, char *prefix)
{
  char buf[128];	/* (safely holds up to 10^14 years) */
  
  if (prefix == NULL)
    fputs("CPU Time: ", fp);
  else 
    fputs(prefix, fp);

  format_time_string(buf, w->user+w->sys, TRUE);
#ifdef HAVE_TIMES
  fprintf(fp, "%.2fu %.2fs %s ", w->user, w->sys, buf);
#else
  fprintf(fp, "%.2fu %s ", w->user, buf);
#endif

  format_time_string(buf, w->elapsed, TRUE);
  fprintf(fp, "Elapsed: %s\n", buf);
  return eslOK;
}
  

/* Function:  esl_stopwatch_Include()
 * Incept:    SRE, Sat Feb 25 10:47:17 2006 [St. Louis]
 *
 * Purpose:   Merge the cpu and system times from a slave into
 *            a master stopwatch. Both watches must be
 *            stopped, and should not be stopped again unless
 *            You Know What You're Doing.
 *           
 *            Elapsed time is not merged. Master is assumed
 *            to be keeping track of the wall clock (real) time,
 *            and the slave/worker watch is ignored.
 *           
 *            Useful in at least two cases. One is in 
 *            PVM, where we merge in the stopwatch(es) from separate
 *            process(es) in a cluster. A second is in 
 *            threads, for broken pthreads/times() implementations
 *            that lose track of cpu times used by spawned
 *            threads.
 *
 * Args:      master  - stopwatch that's aggregating times
 *            w       - watch to add to the master.
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stopwatch_Include(ESL_STOPWATCH *master, ESL_STOPWATCH *w)
{
  master->user    += w->user;
  master->sys     += w->sys;
  return eslOK;
}



/*****************************************************************
 * Example of using the stopwatch module
 *****************************************************************/
#ifdef eslSTOPWATCH_EXAMPLE
/*::cexcerpt::stopwatch_example::begin::*/
/* compile: gcc -g -Wall -I. -o example -DeslSTOPWATCH_EXAMPLE esl_stopwatch.c easel.c -lm
 * run:     ./example
 */
#include "easel.h"
#include "esl_stopwatch.h"

int 
main(void)
{
  ESL_STOPWATCH *w;
  
  w = esl_stopwatch_Create(); 

  esl_stopwatch_Start(w);
  sleep(5);
  esl_stopwatch_Stop(w);

  esl_stopwatch_Display(stdout, w, "CPU Time: ");
  esl_stopwatch_Destroy(w);
  return 0;
}
/*::cexcerpt::stopwatch_example::end::*/
#endif /*ESL_STOPWATCH_EXAMPLE*/
