/* esl_exponential.h
 * Exponential distributions.
 * 
 * SRE, Wed Aug 10 08:32:45 2005 [St. Louis]
 * SVN $Id: esl_exponential.h 172 2007-04-26 00:54:53Z eddys $
 */
#ifndef ESL_EXP_INCLUDED
#define ESL_EXP_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif
#ifdef eslAUGMENT_HISTOGRAM
#include <esl_histogram.h>
#endif

extern double esl_exp_pdf    (double x, double mu, double lambda);
extern double esl_exp_logpdf (double x, double mu, double lambda);
extern double esl_exp_cdf    (double x, double mu, double lambda);
extern double esl_exp_logcdf (double x, double mu, double lambda);
extern double esl_exp_surv   (double x, double mu, double lambda);
extern double esl_exp_logsurv(double x, double mu, double lambda);
extern double esl_exp_invcdf (double p, double mu, double lambda);

extern double esl_exp_generic_pdf   (double x, void *params);
extern double esl_exp_generic_cdf   (double x, void *params);
extern double esl_exp_generic_surv  (double x, void *params);
extern double esl_exp_generic_invcdf(double p, void *params);

extern int    esl_exp_Plot(FILE *fp, double mu, double lambda, 
			   double (*func)(double x, double mu, double lambda), 
			   double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
extern double esl_exp_Sample(ESL_RANDOMNESS *r, double mu, double lambda);
#endif

extern int esl_exp_FitComplete     (double *x, int n, double *ret_mu, double *ret_lambda);
extern int esl_exp_FitCompleteScale(double *x, int n, double      mu, double *ret_lambda);

#ifdef eslAUGMENT_HISTOGRAM
extern int esl_exp_FitCompleteBinned(ESL_HISTOGRAM *h, 
				     double *ret_mu, double *ret_lambda);
#endif


#endif /*ESL_EXP_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
