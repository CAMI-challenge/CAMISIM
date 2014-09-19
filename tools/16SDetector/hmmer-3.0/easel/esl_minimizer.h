/* Multidimensional optimization by conjugate gradient descent.
 * 
 * SRE, Wed Jun 22 09:53:05 2005
 * SVN $Id: esl_minimizer.h 162 2007-04-10 23:50:12Z eddys $
 */
#ifndef ESL_MINIMIZER_INCLUDED
#define ESL_MINIMIZER_INCLUDED

#define MAXITERATIONS 100

extern int esl_min_Bracket(double *a, double *d, double *u, int n, 
			   double (*func)(double *, int, void *), void *prm, 
			   double *ret_fa,
			   double *b, double *ret_bx, double *ret_fb,
			   double *c, double *ret_cx, double *ret_fc);
extern int esl_min_LineSearch(double *ori, double *d, double *u, int n,
			      double (*func)(double *, int, void *), void *prm,
			      double tol, double *b, 
			      double *x, double *ret_xx, double *ret_fx);
extern int esl_min_ConjugateGradientDescent(double *x, double *u, int n, 
					    double (*func)(double *, int, void *),
					    void (*dfunc)(double *, int, void *, double *),
					    void *prm, double tol, double *wrk, double *ret_fx);

#endif /*ESL_MINIMIZER_INCLUDED*/

/*****************************************************************  
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
