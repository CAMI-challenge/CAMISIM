/* esl_gamma.h
 * Gamma distributions.
 * 
 * SRE, Wed Nov 16 19:15:33 2005 [St. Louis]
 * SVN $Id: esl_gamma.h 83 2005-12-13 20:54:07Z eddy $
 */
#ifndef ESL_GAMMA_INCLUDED
#define ESL_GAMMA_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

double esl_gam_pdf    (double x, double mu, double lambda, double tau);
double esl_gam_logpdf (double x, double mu, double lambda, double tau);
double esl_gam_cdf    (double x, double mu, double lambda, double tau);
double esl_gam_logcdf (double x, double mu, double lambda, double tau);
double esl_gam_surv   (double x, double mu, double lambda, double tau);
double esl_gam_logsurv(double x, double mu, double lambda, double tau);
double esl_gam_invcdf (double p, double mu, double lambda, double tau);

double esl_gam_generic_pdf   (double x, void *params);
double esl_gam_generic_cdf   (double x, void *params);
double esl_gam_generic_surv  (double x, void *params);
double esl_gam_generic_invcdf(double x, void *params);

extern int esl_gam_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau), 
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
extern double esl_gam_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

extern int esl_gam_FitComplete(double *x, int n, double mu, double *ret_lambda, double *ret_tau);

#endif /*ESL_GAMMA_INCLUDED*/
