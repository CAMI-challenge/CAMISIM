/* Statistical routines for normal distributions
 * 
 * SRE, Tue Nov 21 14:29:02 2006 [Janelia]
 * SVN $Id: esl_normal.h 269 2008-06-19 13:47:41Z eddys $
 */

#ifndef ESL_NORMAL_INCLUDED
#define ESL_NORMAL_INCLUDED

extern double esl_normal_pdf   (double x, double mu, double sigma);
extern double esl_normal_logpdf(double x, double mu, double sigma);
extern double esl_normal_cdf   (double x, double mu, double sigma);
extern double esl_normal_surv  (double x, double mu, double sigma);

#endif /*ESL_NORMAL_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
