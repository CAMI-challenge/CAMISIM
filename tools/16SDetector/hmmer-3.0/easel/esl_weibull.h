/* esl_weibull.h
 * Weibull distributions.
 * 
 * SRE, Tue Aug  9 10:48:35 2005 [St. Louis]
 * SVN $Id: esl_weibull.h 83 2005-12-13 20:54:07Z eddy $
 */
#ifndef ESL_WEI_INCLUDED
#define ESL_WEI_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

#ifdef eslAUGMENT_HISTOGRAM
#include <esl_histogram.h>
#endif

extern double esl_wei_pdf    (double x, double mu, double lambda, double tau);
extern double esl_wei_logpdf (double x, double mu, double lambda, double tau);
extern double esl_wei_cdf    (double x, double mu, double lambda, double tau);
extern double esl_wei_logcdf (double x, double mu, double lambda, double tau);
extern double esl_wei_surv   (double x, double mu, double lambda, double tau);
extern double esl_wei_logsurv(double x, double mu, double lambda, double tau);
extern double esl_wei_invcdf (double p, double mu, double lambda, double tau);

extern double esl_wei_generic_pdf   (double x, void *params);
extern double esl_wei_generic_cdf   (double x, void *params);
extern double esl_wei_generic_surv  (double x, void *params);
extern double esl_wei_generic_invcdf(double p, void *params);

extern int esl_wei_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau), 
			double xmin, double xmax, double xstep);


#ifdef eslAUGMENT_RANDOM
extern double esl_wei_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

#ifdef eslAUGMENT_MINIMIZER
extern int esl_wei_FitComplete(double *x, int n, double *ret_mu,
			       double *ret_lambda, double *ret_tau);
#ifdef eslAUGMENT_HISTOGRAM
extern int esl_wei_FitCompleteBinned(ESL_HISTOGRAM *h, double *ret_mu,
				     double *ret_lambda, double *ret_tau);
#endif /*eslAUGMENT_HISTOGRAM*/
#endif /*eslAUGMENT_MINIMIZER*/


#endif /*ESL_WEI_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
