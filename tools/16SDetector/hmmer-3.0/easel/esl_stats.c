/* Foundation for the statistics modules.
 * 
 * Contents:
 *   1. The stats API.
 *   2. Unit tests.
 *   3. Test driver.
 *   4. Example.
 *   5. License and copyright information.
 * 
 * SRE, Tue Jul 19 10:57:44 2005
 * SVN $Id: esl_stats.c 341 2009-06-01 12:21:15Z eddys $
 */
#include "esl_config.h"

#include <math.h>

#include "easel.h"
#include "esl_stats.h"


/* Function:  esl_stats_DMean()
 * Synopsis:  Calculates mean and $\sigma^2$ for samples $x_i$.
 * Incept:    SRE, Tue Jul 19 11:04:00 2005 [St. Louis]
 *
 * Purpose:   Calculates the sample mean and $s^2$, the unbiased
 *            estimator of the population variance, for a
 *            sample of <n> numbers <x[0]..x[n-1]>, and optionally
 *            returns either or both through <ret_mean> and
 *            <ret_var>.
 *            
 *            <esl_stats_FMean()> and <esl_stats_IMean()> do the same,
 *            for float and integer vectors.
 *
 * Args:      x        - samples x[0]..x[n-1]
 *            n        - number of samples
 *            opt_mean - optRETURN: mean
 *            opt_var  - optRETURN: estimate of population variance       
 *
 * Returns:   <eslOK> on success.
 */
int
esl_stats_DMean(const double *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}
int
esl_stats_FMean(const float *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}
int
esl_stats_IMean(const int *x, int n, double *opt_mean, double *opt_var)
{
  double sum   = 0.;
  double sqsum = 0.;
  int i;

  for (i = 0; i < n; i++) 
    { 
      sum   += x[i];
      sqsum += x[i]*x[i];
    }
  if (opt_mean != NULL)  *opt_mean = sum / (double) n;
  if (opt_var  != NULL)  *opt_var  = (sqsum - sum*sum/(double)n) / ((double)n-1);
  return eslOK;
}


/* Function:  esl_stats_LogGamma()
 * Synopsis:  Calculates $\log \Gamma(x)$.
 * Incept:    SRE, Tue Nov  2 13:47:01 2004 [St. Louis]
 *
 * Purpose:   Returns natural log of $\Gamma(x)$, for $x > 0$.
 * 
 * Credit:    Adapted from a public domain implementation in the
 *            NCBI core math library. Thanks to John Spouge and
 *            the NCBI. (According to NCBI, that's Dr. John
 *            "Gammas Galore" Spouge to you, pal.)
 *
 * Args:      x          : argument, x > 0.0
 *            ret_answer : RETURN: the answer
 *
 * Returns:   Put the answer in <ret_answer>; returns <eslOK>.
 *            
 * Throws:    <eslERANGE> if $x <= 0$.
 */
int
esl_stats_LogGamma(double x, double *ret_answer)
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
  
  /* Protect against invalid x<=0  */
  if (x <= 0.0)  ESL_EXCEPTION(eslERANGE, "invalid x <= 0 in esl_stats_LogGamma()");

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
  *ret_answer = value;
  return eslOK;
}


/* Function:  esl_stats_Psi()
 * Synopsis:  Calculates $\Psi(x)$ (the digamma function).
 * Incept:    SRE, Tue Nov 15 13:57:59 2005 [St. Louis]
 *
 * Purpose:   Computes $\Psi(x)$ (the "digamma" function), which is
 *            the derivative of log of the Gamma function:
 *            $d/dx \log \Gamma(x) = \frac{\Gamma'(x)}{\Gamma(x)} = \Psi(x)$.
 *            Argument $x$ is $> 0$. 
 * 
 *            This is J.M. Bernardo's "Algorithm AS103",
 *            Appl. Stat. 25:315-317 (1976).  
 */
int
esl_stats_Psi(double x, double *ret_answer)
{
  double answer = 0.;
  double x2;

  if (x <= 0.0) ESL_EXCEPTION(eslERANGE, "invalid x <= 0 in esl_stats_Psi()");
  
  /* For small x, Psi(x) ~= -0.5772 - 1/x + O(x), we're done.
   */
  if (x <= 1e-5) {
    *ret_answer = -eslCONST_EULER - 1./x;
    return eslOK;
  }

  /* For medium x, use Psi(1+x) = \Psi(x) + 1/x to c.o.v. x,
   * big enough for Stirling approximation to work...
   */
  while (x < 8.5) {
    answer = answer - 1./x;
    x += 1.;
  }
  
  /* For large X, use Stirling approximation
   */
  x2 = 1./x;
  answer += log(x) - 0.5 * x2;
  x2 = x2*x2;
  answer -= (1./12.)*x2;
  answer += (1./120.)*x2*x2;
  answer -= (1./252.)*x2*x2*x2;

  *ret_answer = answer;
  return eslOK;
}



/* Function: esl_stats_IncompleteGamma()
 * Synopsis: Calculates the incomplete Gamma function.
 * 
 * Purpose:  Returns $P(a,x)$ and $Q(a,x)$ where:
 *
 *           \begin{eqnarray*}
 *             P(a,x) & = & \frac{1}{\Gamma(a)} \int_{0}^{x} t^{a-1} e^{-t} dt \\
 *                    & = & \frac{\gamma(a,x)}{\Gamma(a)} \\
 *             Q(a,x) & = & \frac{1}{\Gamma(a)} \int_{x}^{\infty} t^{a-1} e^{-t} dt\\
 *                    & = & 1 - P(a,x) \\
 *           \end{eqnarray*}
 *
 *           $P(a,x)$ is the CDF of a gamma density with $\lambda = 1$,
 *           and $Q(a,x)$ is the survival function.
 *           
 *           For $x \simeq 0$, $P(a,x) \simeq 0$ and $Q(a,x) \simeq 1$; and
 *           $P(a,x)$ is less prone to roundoff error. 
 *           
 *           The opposite is the case for large $x >> a$, where
 *           $P(a,x) \simeq 1$ and $Q(a,x) \simeq 0$; there, $Q(a,x)$ is
 *           less prone to roundoff error.
 *
 * Method:   Based on ideas from Numerical Recipes in C, Press et al.,
 *           Cambridge University Press, 1988. 
 *           
 * Args:     a          - for instance, degrees of freedom / 2     [a > 0]
 *           x          - for instance, chi-squared statistic / 2  [x >= 0] 
 *           ret_pax    - RETURN: P(a,x)
 *           ret_qax    - RETURN: Q(a,x)
 *
 * Return:   <eslOK> on success.
 *
 * Throws:   <eslERANGE> if <a> or <x> is out of accepted range.
 *           <eslENOHALT> if approximation fails to converge.
 */          
int
esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax)
{
  int    iter;			/* iteration counter */
  double pax;			/* P(a,x) */
  double qax;			/* Q(a,x) */

  if (a <= 0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): a must be > 0");
  if (x <  0.) ESL_EXCEPTION(eslERANGE, "esl_stats_IncompleteGamma(): x must be >= 0");

  /* For x > a + 1 the following gives rapid convergence;
   * calculate Q(a,x) = \frac{\Gamma(a,x)}{\Gamma(a)},
   * using a continued fraction development for \Gamma(a,x).
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
	    {
	      esl_stats_LogGamma(a, &qax);	      
	      qax = nu1 * exp(a * log(x) - x - qax);

	      if (ret_pax != NULL) *ret_pax = 1 - qax;
	      if (ret_qax != NULL) *ret_qax = qax;
	      return eslOK;
	    }

	  oldp = nu1;
	}
      ESL_EXCEPTION(eslENOHALT,
		"esl_stats_IncompleteGamma(): fraction failed to converge");
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
	    {
	      esl_stats_LogGamma(a, &pax);
	      pax = p * exp(a * log(x) - x - pax);

	      if (ret_pax != NULL) *ret_pax = pax;
	      if (ret_qax != NULL) *ret_qax = 1. - pax;
	      return eslOK;
	    }
	}
      ESL_EXCEPTION(eslENOHALT,
		"esl_stats_IncompleteGamma(): series failed to converge");
    }
  /*NOTREACHED*/
  return eslOK;
}


/* Function:  esl_stats_ChiSquaredTest()
 * Synopsis:  Calculates a $\chi^2$ P-value.
 * Incept:    SRE, Tue Jul 19 11:39:32 2005 [St. Louis]
 *
 * Purpose:   Calculate the probability that a chi-squared statistic
 *            with <v> degrees of freedom would exceed the observed
 *            chi-squared value <x>; return it in <ret_answer>. If
 *            this probability is less than some small threshold (say,
 *            0.05 or 0.01), then we may reject the hypothesis we're
 *            testing.
 *
 * Args:      v          - degrees of freedom
 *            x          - observed chi-squared value
 *            ret_answer - RETURN: P(\chi^2 > x)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslERANGE> if <v> or <x> are out of valid range.
 *            <eslENOHALT> if iterative calculation fails.
 */
int
esl_stats_ChiSquaredTest(int v, double x, double *ret_answer)
{
  return esl_stats_IncompleteGamma((double)v/2, x/2, NULL, ret_answer);
}


/* Function:  esl_stats_LinearRegression()
 * Synopsis:  Fit data to a straight line.
 * Incept:    SRE, Sat May 26 11:33:46 2007 [Janelia]
 *
 * Purpose:   Fit <n> points <x[i]>, <y[i]> to a straight line
 *            $y = a + bx$ by linear regression. 
 *            
 *            The $x_i$ are taken to be known, and the $y_i$ are taken
 *            to be observed quantities associated with a sampling
 *            error $\sigma_i$. If known, the standard deviations
 *            $\sigma_i$ for $y_i$ are provided in the <sigma> array.
 *            If they are unknown, pass <sigma = NULL>, and the
 *            routine will proceed with the assumption that $\sigma_i
 *            = 1$ for all $i$.
 *            
 *            The maximum likelihood estimates for $a$ and $b$ are
 *            optionally returned in <opt_a> and <opt_b>.
 *            
 *            The estimated standard deviations of $a$ and $b$ and
 *            their estimated covariance are optionally returned in
 *            <opt_sigma_a>, <opt_sigma_b>, and <opt_cov_ab>.
 *            
 *            The Pearson correlation coefficient is optionally
 *            returned in <opt_cc>. 
 *            
 *            The $\chi^2$ P-value for the regression fit is
 *            optionally returned in <opt_Q>. This P-value may only be
 *            obtained when the $\sigma_i$ are known. If <sigma> is
 *            passed as <NULL> and <opt_Q> is requested, <*opt_Q> is
 *            set to 1.0.
 *            
 *            This routine follows the description and algorithm in
 *            \citep[pp.661-666]{Press93}.
 *
 *            <n> must be greater than 2; at least two x[i] must
 *            differ; and if <sigma> is provided, all <sigma[i]> must
 *            be $>0$. If any of these conditions isn't met, the
 *            routine throws <eslEINVAL>.
 *
 * Args:      x            - x[0..n-1]
 *            y            - y[0..n-1]
 *            sigma        - sample error in observed y_i
 *            n            - number of data points
 *            opt_a        - optRETURN: intercept estimate		
 *            opt_b        - optRETURN: slope estimate
 *            opt_sigma_a  - optRETURN: error in estimate of a
 *            opt_sigma_b  - optRETURN: error in estimate of b
 *            opt_cov_ab   - optRETURN: covariance of a,b estimates
 *            opt_cc       - optRETURN: Pearson correlation coefficient for x,y
 *            opt_Q        - optRETURN: X^2 P-value for linear fit
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error;
 *            <eslEINVAL> if a contract condition isn't met;
 *            <eslENORESULT> if the chi-squared test fails.
 *            In these cases, all optional return values are set to 0.
 */
int
esl_stats_LinearRegression(const double *x, const double *y, const double *sigma, int n,
			   double *opt_a,       double *opt_b,
			   double *opt_sigma_a, double *opt_sigma_b, double *opt_cov_ab,
			   double *opt_cc,      double *opt_Q)
{
  int     status;
  double *t      = NULL;
  double  S, Sx, Sy, Stt;
  double  Sxy, Sxx, Syy;
  double  a, b, sigma_a, sigma_b, cov_ab, cc, X2, Q;
  double  xdev, ydev;
  double  tmp;
  int     i;

  /* Contract checks. */
  if (n <= 2) ESL_XEXCEPTION(eslEINVAL, "n must be > 2 for linear regression fitting");
  if (sigma != NULL) 
    for (i = 0; i < n; i++) if (sigma[i] <= 0.) ESL_XEXCEPTION(eslEINVAL, "sigma[%d] <= 0", i);
  status = eslEINVAL;
  for (i = 0; i < n; i++) if (x[i] != 0.) { status = eslOK; break; }
  if (status != eslOK) ESL_XEXCEPTION(eslEINVAL, "all x[i] are 0.");

  /* Allocations */
  ESL_ALLOC(t, sizeof(double) * n);

  /* S = \sum_{i=1}{n} \frac{1}{\sigma_i^2}.  (S > 0.) */
  if (sigma != NULL) { for (S = 0., i = 0; i < n; i++) S += 1./ (sigma[i] * sigma[i]);  }
  else S = (double) n;

  /* S_x = \sum_{i=1}{n} \frac{x[i]}{ \sigma_i^2}  (Sx real.) */
  for (Sx = 0., i = 0; i < n; i++) { 
    if (sigma == NULL) Sx += x[i];
    else               Sx += x[i] / (sigma[i] * sigma[i]);
  }

  /* S_y = \sum_{i=1}{n} \frac{y[i]}{\sigma_i^2}  (Sy real.) */
  for (Sy = 0., i = 0; i < n; i++) { 
    if (sigma == NULL) Sy += y[i];
    else               Sy += y[i] / (sigma[i] * sigma[i]);
  }

  /* t_i = \frac{1}{\sigma_i} \left( x_i - \frac{S_x}{S} \right)   (t_i real) */
  for (i = 0; i < n; i++) {
    t[i] = x[i] - Sx/S;
    if (sigma != NULL) t[i] /= sigma[i];
  }

  /* S_{tt} = \sum_{i=1}^n t_i^2  (if at least one x is != 0, Stt > 0) */
  for (Stt = 0., i = 0; i < n; i++) { Stt += t[i] * t[i]; }

  /* b = \frac{1}{S_{tt}} \sum_{i=1}^{N} \frac{t_i y_i}{\sigma_i}  */
  for (b = 0., i = 0; i < n; i++) {
    if (sigma != NULL) { b += t[i]*y[i] / sigma[i]; }
    else               { b += t[i]*y[i]; }
  }
  b /= Stt;

  /* a = \frac{ S_y - S_x b } {S}   */
  a = (Sy - Sx * b) / S;
  
  /* \sigma_a^2 = \frac{1}{S} \left( 1 + \frac{ S_x^2 }{S S_{tt}} \right) */
  sigma_a = sqrt ((1. + (Sx*Sx) / (S*Stt)) / S);

  /* \sigma_b = \frac{1}{S_{tt}} */
  sigma_b = sqrt (1. / Stt);

  /* Cov(a,b) = - \frac{S_x}{S S_{tt}}    */
  cov_ab = -Sx / (S * Stt);
  
  /* Pearson correlation coefficient */
  Sxy = Sxx = Syy = 0.;
  for (i = 0; i < n; i++) {
    if (sigma != NULL) { 
      xdev = (x[i] / (sigma[i] * sigma[i])) - (Sx / n);
      ydev = (y[i] / (sigma[i] * sigma[i])) - (Sy / n);
    } else {
      xdev = x[i] - (Sx / n);
      ydev = y[i] - (Sy / n);
    }
    Sxy += xdev * ydev;
    Sxx += xdev * xdev;
    Syy += ydev * ydev;
  }
  cc = Sxy / (sqrt(Sxx) * sqrt(Syy));

  /* \chi^2 */
  for (X2 = 0., i = 0; i < n; i++) {
    tmp =  y[i] - a - b*x[i];
    if (sigma != NULL) tmp /= sigma[i];
    X2 += tmp*tmp;
  }
  
  /* We can calculate a goodness of fit if we know the \sigma_i */
  if (sigma != NULL) {
    if (esl_stats_ChiSquaredTest(n-2, X2, &Q) != eslOK) { status = eslENORESULT; goto ERROR; }
  } else Q = 1.0;

  /* If we didn't use \sigma_i, adjust the sigmas for a,b */
  if (sigma == NULL) {
    tmp = sqrt(X2 / (double)(n-2));
    sigma_a *= tmp;
    sigma_b *= tmp;
  }
    
  /* Done. Set up for normal return.
   */
  free(t);
  if (opt_a       != NULL) *opt_a       = a;
  if (opt_b       != NULL) *opt_b       = b;
  if (opt_sigma_a != NULL) *opt_sigma_a = sigma_a;
  if (opt_sigma_b != NULL) *opt_sigma_b = sigma_b;
  if (opt_cov_ab  != NULL) *opt_cov_ab  = cov_ab;
  if (opt_cc      != NULL) *opt_cc      = cc;
  if (opt_Q       != NULL) *opt_Q       = Q;
  return eslOK;
  
 ERROR:
  if (t != NULL) free(t);
  if (opt_a       != NULL) *opt_a       = 0.;
  if (opt_b       != NULL) *opt_b       = 0.;
  if (opt_sigma_a != NULL) *opt_sigma_a = 0.;
  if (opt_sigma_b != NULL) *opt_sigma_b = 0.;
  if (opt_cov_ab  != NULL) *opt_cov_ab  = 0.;
  if (opt_cc      != NULL) *opt_cc      = 0.;
  if (opt_Q       != NULL) *opt_Q       = 0.;
  return status;
}
/*---------------- end of API implementation --------------------*/




/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef eslSTATS_TESTDRIVE
#include "esl_random.h"
#include "esl_stopwatch.h"
#ifdef HAVE_LIBGSL 
#include <gsl/gsl_sf_gamma.h>
#endif

/* The LogGamma() function is rate-limiting in hmmbuild, because it is
 * used so heavily in mixture Dirichlet calculations.
 *    ./configure --with-gsl; [compile test driver]
 *    ./stats_utest -v
 * runs a comparison of time/precision against GSL.
 * SRE, Sat May 23 10:04:41 2009, on home Mac:
 *     LogGamma       = 1.29u  / N=1e8  =  13 nsec/call
 *     gsl_sf_lngamma = 1.43u  / N=1e8  =  14 nsec/call
 */
static void
utest_LogGamma(ESL_RANDOMNESS *r, int N, int be_verbose)
{
  char          *msg = "esl_stats_LogGamma() unit test failed";
  ESL_STOPWATCH *w   = esl_stopwatch_Create();
  double        *x   = malloc(sizeof(double) * N);
  double        *lg  = malloc(sizeof(double) * N);
  double        *lg2 = malloc(sizeof(double) * N);
  int            i;

  for (i = 0; i < N; i++) 
    x[i] = esl_random(r) * 100.;
  
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) 
    if (esl_stats_LogGamma(x[i], &(lg[i])) != eslOK) esl_fatal(msg);
  esl_stopwatch_Stop(w);

  if (be_verbose) esl_stopwatch_Display(stdout, w, "esl_stats_LogGamma() timing: ");

#ifdef HAVE_LIBGSL
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) lg2[i] = gsl_sf_lngamma(x[i]);
  esl_stopwatch_Stop(w);

  if (be_verbose) esl_stopwatch_Display(stdout, w, "gsl_sf_lngamma() timing:     ");
  
  for (i = 0; i < N; i++)
    if (esl_DCompare(lg[i], lg2[i], 1e-2) != eslOK) esl_fatal(msg);
#endif
  
  free(lg2);
  free(lg);
  free(x);
  esl_stopwatch_Destroy(w);
}


/* The test of esl_stats_LinearRegression() is a statistical test,
 * so we can't be too aggressive about testing results. 
 * 
 * Args:
 *    r          - a source of randomness
 *    use_sigma  - TRUE to pass sigma to the regression fit.
 *    be_verbose - TRUE to print results (manual, not automated test mode)
 */
static void
utest_LinearRegression(ESL_RANDOMNESS *r, int use_sigma, int be_verbose)
{
  char msg[] = "linear regression unit test failed";
  double a     = -3.;
  double b     = 1.;
  int    n     = 100;
  double xori  = -20.;
  double xstep = 1.0;
  double setsigma = 1.0;		/* sigma on all points */
  int    i;
  double *x     = NULL;
  double *y     = NULL;
  double *sigma = NULL;
  double  ae, be, siga, sigb, cov_ab, cc, Q;
  
  if ((x     = malloc(sizeof(double) * n)) == NULL) esl_fatal(msg);
  if ((y     = malloc(sizeof(double) * n)) == NULL) esl_fatal(msg);
  if ((sigma = malloc(sizeof(double) * n)) == NULL) esl_fatal(msg);
  
  /* Simulate some linear data */
  for (i = 0; i < n; i++)
    {
      sigma[i] = setsigma;
      x[i]     = xori + i*xstep;
      y[i]     = esl_rnd_Gaussian(r, a + b*x[i], sigma[i]);
    }
  
  if (use_sigma) {
    if (esl_stats_LinearRegression(x, y, sigma, n, &ae, &be, &siga, &sigb, &cov_ab, &cc, &Q) != eslOK) esl_fatal(msg);
  } else {
    if (esl_stats_LinearRegression(x, y,  NULL, n, &ae, &be, &siga, &sigb, &cov_ab, &cc, &Q) != eslOK) esl_fatal(msg);
  }

  if (be_verbose) {
    printf("Linear regression test:\n");
    printf("estimated intercept a = %8.4f   [true = %8.4f]\n", ae, a);
    printf("estimated slope b     = %8.4f   [true = %8.4f]\n", be, b);
    printf("estimated sigma on a  = %8.4f\n",                  siga);
    printf("estimated sigma on b  = %8.4f\n",                  sigb);
    printf("estimated cov(a,b)    = %8.4f\n",                  cov_ab);
    printf("correlation coeff     = %8.4f\n",                  cc);
    printf("P-value               = %8.4f\n",                  Q);
  }

  /* The following tests are statistical.
   */
  if ( fabs(ae-a) > 2*siga ) esl_fatal(msg);
  if ( fabs(be-b) > 2*sigb ) esl_fatal(msg);
  if ( cc < 0.95)            esl_fatal(msg);
  if (use_sigma) {
    if (Q < 0.001)           esl_fatal(msg);
  } else {
    if (Q != 1.0)            esl_fatal(msg);
  }

  free(x);
  free(y);
  free(sigma);
}
#endif /*eslSTATS_TESTDRIVE*/
/*-------------------- end of unit tests ------------------------*/




/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef eslSTATS_TESTDRIVE
/* gcc -g -Wall -o stats_utest  -L. -I. -DeslSTATS_TESTDRIVE esl_stats.c -leasel -lm
 * gcc -DHAVE_LIBGSL -O2 -o stats_utest -L. -I. -DeslSTATS_TESTDRIVE esl_stats.c -leasel -lgsl -lm
 */
#include <stdio.h>
#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stats.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                   0},
  {"-s",  eslARG_INT,      "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",         0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "verbose: show verbose output",          0},
  {"-N",  eslARG_INT,"10000000", NULL, NULL, NULL, NULL, NULL, "number of trials in LogGamma test",     0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for stats special functions";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r          = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  int             be_verbose = esl_opt_GetBoolean(go, "-v");
  int             N          = esl_opt_GetInteger(go, "-N");

  if (be_verbose) printf("seed = %" PRIu32 "\n", esl_randomness_GetSeed(r));

  utest_LogGamma(r, N, be_verbose);
  utest_LinearRegression(r, TRUE,  be_verbose);
  utest_LinearRegression(r, FALSE, be_verbose);
  
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  exit(0);
}
#endif /*eslSTATS_TESTDRIVE*/
/*------------------- end of test driver ------------------------*/




/*****************************************************************
 * 4. Example.
 *****************************************************************/

/* Compile:  gcc -g -Wall -o example -I. -DeslSTATS_EXAMPLE esl_stats.c esl_random.c easel.c -lm  
 * or        gcc -g -Wall -o example -I. -L. -DeslSTATS_EXAMPLE esl_stats.c -leasel -lm  
 */
#ifdef eslSTATS_EXAMPLE
/*::cexcerpt::stats_example::begin::*/
/* gcc -g -Wall -o example -I. -DeslSTATS_EXAMPLE esl_stats.c esl_random.c easel.c -lm  */
#include <stdio.h>
#include "easel.h"
#include "esl_random.h"
#include "esl_stats.h"

int main(void)
{
  ESL_RANDOMNESS *r   = esl_randomness_Create(0);
  double a            = -3.;
  double b            = 1.;
  double xori         = -20.;
  double xstep        = 1.0;
  double setsigma     = 1.0;		/* sigma on all points */
  int    n            = 100;
  double *x           = malloc(sizeof(double) * n);
  double *y           = malloc(sizeof(double) * n);
  double *sigma       = malloc(sizeof(double) * n);
  int    i;
  double  ae, be, siga, sigb, cov_ab, cc, Q;
  
  /* Simulate some linear data, with Gaussian noise added to y_i */
  for (i = 0; i < n; i++) {
    sigma[i] = setsigma;
    x[i]     = xori + i*xstep;
    y[i]     = esl_rnd_Gaussian(r, a + b*x[i], sigma[i]);
  }
  
  if (esl_stats_LinearRegression(x, y, sigma, n, &ae, &be, &siga, &sigb, &cov_ab, &cc, &Q) != eslOK)
    esl_fatal("linear regression failed");

  printf("estimated intercept a = %8.4f   [true = %8.4f]\n", ae, a);
  printf("estimated slope b     = %8.4f   [true = %8.4f]\n", be, b);
  printf("estimated sigma on a  = %8.4f\n",                  siga);
  printf("estimated sigma on b  = %8.4f\n",                  sigb);
  printf("estimated cov(a,b)    = %8.4f\n",                  cov_ab);
  printf("correlation coeff     = %8.4f\n",                  cc);
  printf("P-value               = %8.4f\n",                  Q);

  free(x);  free(y);  free(sigma); 
  esl_randomness_Destroy(r);
  exit(0);
}
/*::cexcerpt::stats_example::end::*/
#endif

/*--------------------- end of example --------------------------*/




/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
