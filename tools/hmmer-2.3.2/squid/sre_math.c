/*****************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 *****************************************************************/

/* sre_math.c
 * 
 * Portability for and extensions to C math library.
 * RCS $Id: sre_math.c,v 1.13 2003/04/14 16:00:16 eddy Exp $
 */

#include "squidconf.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "squid.h"

  
/* Function: Linefit()
 * 
 * Purpose:  Given points x[0..N-1] and y[0..N-1], fit to
 *           a straight line y = a + bx.
 *           a, b, and the linear correlation coefficient r
 *           are filled in for return.
 *           
 * Args:     x     - x values of data
 *           y     - y values of data               
 *           N     - number of data points
 *           ret_a - RETURN: intercept
 *           ret_b - RETURN: slope
 *           ret_r - RETURN: correlation coefficient  
 *           
 * Return:   1 on success, 0 on failure.
 */
int          
Linefit(float *x, float *y, int N, float *ret_a, float *ret_b, float *ret_r) 
{				
  float xavg, yavg;
  float sxx, syy, sxy;
  int   i;
  
  /* Calculate averages, xavg and yavg
   */
  xavg = yavg = 0.0;
  for (i = 0; i < N; i++)
    {
      xavg += x[i];
      yavg += y[i];
    }
  xavg /= (float) N;
  yavg /= (float) N;

  sxx = syy = sxy = 0.0;
  for (i = 0; i < N; i++)
    {
      sxx    += (x[i] - xavg) * (x[i] - xavg);
      syy    += (y[i] - yavg) * (y[i] - xavg);
      sxy    += (x[i] - xavg) * (y[i] - yavg);
    }
  *ret_b = sxy / sxx;
  *ret_a = yavg - xavg*(*ret_b);
  *ret_r = sxy / (sqrt(sxx) * sqrt(syy));
  return 1;
}


/* Function: WeightedLinefit()
 * 
 * Purpose:  Given points x[0..N-1] and y[0..N-1] with
 *           variances (measurement errors) var[0..N-1],  
 *           fit to a straight line y = mx + b.
 *           
 * Method:   Algorithm from Numerical Recipes in C, [Press88].
 *           
 * Return:   (void)
 *           ret_m contains slope; ret_b contains intercept 
 */                
void
WeightedLinefit(float *x, float *y, float *var, int N, float *ret_m, float *ret_b) 
{
  int    i;
  double s;
  double sx, sy;
  double sxx, sxy;
  double delta;
  double m, b;
  
  s = sx = sy = sxx = sxy = 0.;
  for (i = 0; i < N; i++)
    {
      s   += 1./var[i];
      sx  += x[i] / var[i];
      sy  += y[i] / var[i];
      sxx += x[i] * x[i] / var[i];
      sxy += x[i] * y[i] / var[i];
    }

  delta = s * sxx - (sx * sx);
  b = (sxx * sy - sx * sxy) / delta;
  m = (s * sxy - sx * sy) / delta;

  *ret_m = m;
  *ret_b = b;
}
  

/* Function: Gammln()
 *
 * Returns the natural log of the gamma function of x.
 * x is > 0.0.  
 *
 * Adapted from a public domain implementation in the
 * NCBI core math library. Thanks to John Spouge and
 * the NCBI. (According to the NCBI, that's Dr. John
 * "Gammas Galore" Spouge to you, pal.)
 */
double
Gammln(double x)
{
  int i;
  double xx, tx;
  double tmp, value;
  static double cof[11] = {
    4.694580336184385e+04,
    -1.560605207784446e+05,
    2.065049568014106e+05,
    -1.388934775095388e+05,
    5.031796415085709e+04,
    -9.601592329182778e+03,
    8.785855930895250e+02,
    -3.155153906098611e+01,
    2.908143421162229e-01,
    -2.319827630494973e-04,
    1.251639670050933e-10
  };
  
  /* Protect against x=0. We see this in Dirichlet code,
   * for terms alpha = 0. This is a severe hack but it is effective
   * and (we think?) safe. (due to GJM)
   */ 
  if (x <= 0.0) return 999999.; 

  xx       = x - 1.0;
  tx = tmp = xx + 11.0;
  value    = 1.0;
  for (i = 10; i >= 0; i--)	/* sum least significant terms first */
    {
      value += cof[i] / tmp;
      tmp   -= 1.0;
    }
  value  = log(value);
  tx    += 0.5;
  value += 0.918938533 + (xx+0.5)*log(tx) - tx;
  return value;
}


/* 2D matrix operations
 */
float **
FMX2Alloc(int rows, int cols)
{
  float **mx;
  int     r;
  
  mx    = (float **) MallocOrDie(sizeof(float *) * rows);
  mx[0] = (float *)  MallocOrDie(sizeof(float) * rows * cols);
  for (r = 1; r < rows; r++)
    mx[r] = mx[0] + r*cols;
  return mx;
}
void
FMX2Free(float **mx)
{
  free(mx[0]);
  free(mx);
}
double **
DMX2Alloc(int rows, int cols)
{
  double **mx;
  int      r;
  
  mx    = (double **) MallocOrDie(sizeof(double *) * rows);
  mx[0] = (double *)  MallocOrDie(sizeof(double) * rows * cols);
  for (r = 1; r < rows; r++)
    mx[r] = mx[0] + r*cols;
  return mx;
}
void
DMX2Free(double **mx)
{
  free(mx[0]);
  free(mx);
}
/* Function: FMX2Multiply()
 * 
 * Purpose:  Matrix multiplication.
 *           Multiply an m x p matrix A by a p x n matrix B,
 *           giving an m x n matrix C.
 *           Matrix C must be a preallocated matrix of the right
 *           size.
 */
void
FMX2Multiply(float **A, float **B, float **C, int m, int p, int n)
{
  int i, j, k;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      {
	C[i][j] = 0.;
	for (k = 0; k < p; k++)
	  C[i][j] += A[i][p] * B[p][j];
      }
}


/* Function: IncompleteGamma()
 * 
 * Purpose:  Returns 1 - P(a,x) where:
 *           P(a,x) = \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt
 *                  = \frac{\gamma(a,x)}{\Gamma(a)}
 *                  = 1 - \frac{\Gamma(a,x)}{\Gamma(a)}
 *                  
 *           Used in a chi-squared test: for a X^2 statistic x
 *           with v degrees of freedom, call:
 *                  p = IncompleteGamma(v/2., x/2.) 
 *           to get the probability p that a chi-squared value
 *           greater than x could be obtained by chance even for
 *           a correct model. (i.e. p should be large, say 
 *           0.95 or more).
 *           
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988. 
 *           
 * Args:     a  - for instance, degrees of freedom / 2     [a > 0]
 *           x  - for instance, chi-squared statistic / 2  [x >= 0] 
 *           
 * Return:   1 - P(a,x).
 */          
double
IncompleteGamma(double a, double x)
{
  int iter;			/* iteration counter */

  if (a <= 0.) Die("IncompleteGamma(): a must be > 0");
  if (x <  0.) Die("IncompleteGamma(): x must be >= 0");

  /* For x > a + 1 the following gives rapid convergence;
   * calculate 1 - P(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)}:
   *     use a continued fraction development for \Gamma(a,x).
   */
  if (x > a+1) 
    {
      double oldp;		/* previous value of p    */
      double nu0, nu1;		/* numerators for continued fraction calc   */
      double de0, de1;		/* denominators for continued fraction calc */

      nu0 = 0.;			/* A_0 = 0       */
      de0 = 1.;			/* B_0 = 1       */
      nu1 = 1.;			/* A_1 = 1       */
      de1 = x;			/* B_1 = x       */

      oldp = nu1;
      for (iter = 1; iter < 100; iter++)
	{
	  /* Continued fraction development:
	   * set A_j = b_j A_j-1 + a_j A_j-2
	   *     B_j = b_j B_j-1 + a_j B_j-2
           * We start with A_2, B_2.
	   */
				/* j = even: a_j = iter-a, b_j = 1 */
				/* A,B_j-2 are in nu0, de0; A,B_j-1 are in nu1,de1 */
	  nu0 = nu1 + ((double)iter - a) * nu0;
	  de0 = de1 + ((double)iter - a) * de0;

				/* j = odd: a_j = iter, b_j = x */
				/* A,B_j-2 are in nu1, de1; A,B_j-1 in nu0,de0 */
	  nu1 = x * nu0 + (double) iter * nu1;
	  de1 = x * de0 + (double) iter * de1;

				/* rescale */
	  if (de1 != 0.) 
	    { 
	      nu0 /= de1; 
	      de0 /= de1;
	      nu1 /= de1;
	      de1 =  1.;
	    }
				/* check for convergence */
	  if (fabs((nu1-oldp)/nu1) < 1.e-7)
	    return nu1 * exp(a * log(x) - x - Gammln(a));

	  oldp = nu1;
	}
      Die("IncompleteGamma(): failed to converge using continued fraction approx");
    }
  else /* x <= a+1 */
    {
      double p;			/* current sum               */
      double val;		/* current value used in sum */

      /* For x <= a+1 we use a convergent series instead:
       *   P(a,x) = \frac{\gamma(a,x)}{\Gamma(a)},
       * where
       *   \gamma(a,x) = e^{-x}x^a \sum_{n=0}{\infty} \frac{\Gamma{a}}{\Gamma{a+1+n}} x^n
       * which looks appalling but the sum is in fact rearrangeable to
       * a simple series without the \Gamma functions:
       *   = \frac{1}{a} + \frac{x}{a(a+1)} + \frac{x^2}{a(a+1)(a+2)} ...
       * and it's obvious that this should converge nicely for x <= a+1.
       */
      
      p = val = 1. / a;
      for (iter = 1; iter < 10000; iter++)
	{
	  val *= x / (a+(double)iter);
	  p   += val;
	  
	  if (fabs(val/p) < 1.e-7)
	    return 1. - p * exp(a * log(x) - x - Gammln(a));
	}
      Die("IncompleteGamma(): failed to converge using series approx");
    }
  /*NOTREACHED*/
  return 0.;
}
  
