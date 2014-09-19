/* Portable, threadsafe random number generators.
 * Provides both a fast generator and a strong generator.
 *
 *  1. The ESL_RANDOMNESS object.
 *  2. The generators and esl_random().
 *  3. Other fundamental sampling (including Gaussian, gamma).
 *  4. Multinomial sampling from discrete probability n-vectors.
 *  5. Benchmark driver
 *  6. Unit tests.
 *  7. Test driver.
 *  8. Example.
 *  9. Copyright and license information.
 *  
 * See http://csrc.nist.gov/rng/ for the NIST random number
 * generation test suite.
 * 
 * SRE, Wed Jul 14 10:54:46 2004 [St. Louis]
 * SVN $Id: esl_random.c 519 2010-02-19 14:04:51Z eddys $
 * 
 * SRE, 30 May 2009: replaced with the Mersenne Twister and Knuth LCG.
 */
#include "esl_config.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "easel.h"
#include "esl_random.h"

static uint32_t choose_arbitrary_seed(void);
static uint32_t jenkins_mix3(uint32_t a, uint32_t b, uint32_t c);
static uint32_t knuth              (ESL_RANDOMNESS *r);
static uint32_t mersenne_twister   (ESL_RANDOMNESS *r);
static void     mersenne_seed_table(ESL_RANDOMNESS *r, uint32_t seed);
static void     mersenne_fill_table(ESL_RANDOMNESS *r);

/*****************************************************************
 *# 1. The <ESL_RANDOMNESS> object.
 *****************************************************************/

/* Function:  esl_randomness_Create()
 * Synopsis:  Create the default strong random number generator.
 * Incept:    SRE, Wed Jul 14 13:02:18 2004 [St. Louis]
 *
 * Purpose:   Create a random number generator using
 *            a given random seed. The <seed> must be $\geq 0$.
 *            
 *            The default random number generator uses the Mersenne
 *            Twister MT19937 algorithm \citep{Matsumoto98}.  It has a
 *            period of $2^{19937}-1$, and equidistribution over
 *            $2^{32}$ values.
 *
 *            If <seed> is $>0$, the random number generator is
 *            reproducibly initialized with that seed.  Two RNGs
 *            created with the same nonzero seed will give exactly the
 *            same stream of pseudorandom numbers. This allows you to
 *            make reproducible stochastic simulations, for example.
 *            
 *            If <seed> is 0, an arbitrary seed is chosen.
 *            Internally, the arbitrary seed is produced by a
 *            combination of the current <time()> and the process id
 *            (if available; POSIX only). Two RNGs created with
 *            <seed>=0 will very probably (but not assuredly) give
 *            different streams of pseudorandom numbers. The true seed
 *            can be retrieved from the <ESL_RANDOMNESS> object using
 *            <esl_randomness_GetSeed()>.  The strategy used for
 *            choosing the arbitrary seed is predictable, so it is
 *            not secure in any sense, especially in the cryptographic
 *            sense.
 *            
 * Args:      seed $>= 0$.
 *
 * Returns:   an initialized <ESL_RANDOMNESS *> on success.
 *            Caller free's with <esl_randomness_Destroy()>.
 *              
 * Throws:    <NULL> on failure.
 * 
 * Xref:      STL8/p57.
 *            J5/21:    Mersenne Twister.
 */
ESL_RANDOMNESS *
esl_randomness_Create(uint32_t seed)
{
  ESL_RANDOMNESS *r      = NULL;
  int             status;

  ESL_ALLOC(r, sizeof(ESL_RANDOMNESS));
  r->type = eslRND_MERSENNE;
  r->mti  = 0;
  r->x    = 0;
  r->seed = 0;
  esl_randomness_Init(r, seed);
  return r;

 ERROR:
  return NULL;
}

/* Function:  esl_randomness_CreateFast()
 * Synopsis:  Create the alternative fast generator.
 * Incept:    SRE, Sat May 30 06:35:23 2009 [Stockholm]
 *
 * Purpose:   Same as <esl_randomness_Create()>, except that a simple
 *            linear congruential generator will be used.
 *            
 *            This is a $(a=69069, c=1)$ LCG, with a period of
 *            $2^{32}$. Because of the relatively short period, this
 *            generator should not be used for serious simulations
 *            involving large samples.
 *
 *            The properties of this generator are not as good as the
 *            default Mersenne Twister, but it is faster, especially
 *            if you only need a small number of samples from the
 *            generator; it is about 20x faster to initialize the
 *            generator, and about 25\% faster to sample a number, 
 *            compared to the default.
 *
 * Args:      seed $>= 0$.
 *
 * Returns:   an initialized <ESL_RANDOMNESS *> on success.
 *            Caller free's with <esl_randomness_Destroy()>.
 *              
 * Throws:    <NULL> on failure.
 *
 * Xref:      J5/44: for accidental proof that the period is
 *                   indeed 2^32.
 */
ESL_RANDOMNESS *
esl_randomness_CreateFast(uint32_t seed)
{
  ESL_RANDOMNESS *r      = NULL;
  int             status;

  ESL_ALLOC(r, sizeof(ESL_RANDOMNESS));
  r->type = eslRND_FAST;
  r->mti  = 0;
  r->x    = 0;
  r->seed = 0;
  esl_randomness_Init(r, seed);
  return r;

 ERROR:
  return NULL;
}


/* Function:  esl_randomness_CreateTimeseeded()
 * Synopsis:  Create an RNG with a quasirandom seed.
 * Incept:    SRE, Wed Jul 14 11:22:54 2004 [St. Louis]
 *
 * Purpose:   Like <esl_randomness_Create()>, but it initializes the
 *            the random number generator using a POSIX <time()> call 
 *            (number of seconds since the POSIX epoch).
 *            
 *            This function is deprecated. Use 
 *            <esl_randomness_Create(0)> instead.
 *
 * Returns:   an initialized <ESL_RANDOMNESS *> on success.
 *            Caller free's with <esl_randomness_Destroy()>.
 *              
 * Throws:    <NULL> on failure.
 * 
 * Xref:      STL8/p57.
 */
ESL_RANDOMNESS *
esl_randomness_CreateTimeseeded(void)
{
  return esl_randomness_Create(0);
}


/* Function:  esl_randomness_Init()
 * Synopsis:  Reinitialize a RNG.           
 * Incept:    SRE, Wed Jul 14 13:13:05 2004 [St. Louis]
 *
 * Purpose:   Reset and reinitialize an existing <ESL_RANDOMNESS>
 *            object with a new seed. 
 *            
 *            Not generally recommended. This does not make a
 *            sequence of numbers more random, and may make it less
 *            so. Sometimes, though, it's useful to reseed an RNG
 *            to guarantee a particular reproducible series of
 *            pseudorandom numbers at an arbitrary point in a program;
 *            HMMER does this, for example, to guarantee the same
 *            results from the same HMM/sequence comparison regardless
 *            of where in a search the HMM or sequence occurs.
 *
 * Args:      r     - randomness object
 *            seed  - new seed to use; >0.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if seed is $<= 0$.
 *
 * Xref:      STL8/p57.
 */
int
esl_randomness_Init(ESL_RANDOMNESS *r, uint32_t seed)
{
  if (seed == 0) seed = choose_arbitrary_seed();
  if (r->type == eslRND_MERSENNE)
    {
      mersenne_seed_table(r, seed);
      mersenne_fill_table(r);
    }
  else 
    {
      r->seed = seed;
      r->x    = jenkins_mix3(seed, 87654321, 12345678);	/* arbitrary dispersion of the seed */
      if (r->x == 0) r->x = 42;                         /* make sure we don't have a zero */
    }
  return eslOK;
}

/* Function:  esl_randomness_GetSeed()
 * Synopsis:  Returns the value of RNG's seed.
 * Incept:    SRE, Wed May 23 17:02:59 2007 [Janelia]
 *
 * Purpose:   Return the value of the seed. 
 */
uint32_t
esl_randomness_GetSeed(const ESL_RANDOMNESS *r)
{
  return r->seed;
}


/* Function:  esl_randomness_Destroy()
 * Synopsis:  Free an RNG.            
 * Incept:    SRE, Wed Jul 14 13:19:08 2004 [St. Louis]
 *
 * Purpose:   Frees an <ESL_RANDOMNESS> object.
 */
void
esl_randomness_Destroy(ESL_RANDOMNESS *r)
{
  free(r);
  return;
}

/*----------- end of ESL_RANDOMNESS object functions --------------*/



/*****************************************************************
 *# 2. The generators and <esl_random()>
 *****************************************************************/  

/* Function: esl_random()  
 * Synopsis: Generate a uniform random deviate on [0,1)
 * Incept:   SRE, Sat May 30 05:01:45 2009 [Stockholm]
 *
 * Purpose:  Returns a uniform deviate x, $0.0 <= x < 1.0$, given
 *           RNG <r>.
 *           
 *           Uses the original Mersenne Twister algorithm, MT19937
 *           [Matsumoto98]. This generator has a period of $2^19937 -
 *           1$. It generates uniformly distributed variates on the
 *           interval $0..2^32-1$. 
 *           
 * Notes:    Easel previously used a reimplementation of ran2() from
 *           Numerical Recipes in C, which uses L'Ecuyer's algorithm
 *           for combining output of two linear congruential
 *           generators, plus a Bays-Durham shuffle \citep{Press93}.
 *           MT is about 10x faster.
 *
 * Returns:  uniformly distribute random deviate on interval
 *           $0.0 \leq x < 1.0$
 *
 * Throws:   (no abnormal error conditions)
 */
double
esl_random(ESL_RANDOMNESS *r)
{
  uint32_t x = (r->type == eslRND_MERSENNE) ? mersenne_twister(r) : knuth(r);
  return ((double) x / 4294967296.0); /* 2^32: normalizes to [0,1) */
}


static uint32_t 
knuth(ESL_RANDOMNESS *r)
{
  r->x *= 69069;
  r->x += 1;
  return r->x;
}

/* mersenne_twister() and other mersenne_*() functions below:
 * A simple serial implementation of the original Mersenne Twister
 * algorithm [Matsumoto98]. 
 * 
 * There are more sophisticated and faster implementations of MT, using
 * vector instructions and/or directly generating IEEE754 doubles
 * bitwise rather than doing an expensive normalization. We can
 * improve the implementation later if necessary, but even the basic
 * MT offers ~10x speed improvement over Easel's previous RNG.
 * [SRE, 30 May 09, Stockholm]
 */
static uint32_t
mersenne_twister(ESL_RANDOMNESS *r)
{
  uint32_t x;
  if (r->mti >= 624) mersenne_fill_table(r);

  x = r->mt[r->mti++];
  x ^= (x>>11);
  x ^= (x<< 7) & 0x9d2c5680;
  x ^= (x<<15) & 0xefc60000;
  x ^= (x>>18);
  return x;
}

/* mersenne_seed_table()
 * Initialize the state of the RNG from a seed.
 * Uses the knuth linear congruential generator.
 */
static void
mersenne_seed_table(ESL_RANDOMNESS *r, uint32_t seed)
{
  int z;

  r->seed  = seed;
  r->mt[0] = seed;
  for (z = 1; z < 624; z++)
    r->mt[z] = 69069 * r->mt[z-1];
  return;
}

/* mersenne_fill_table()
 * Refill the table with 624 new random numbers.
 * We do this whenever we've reseeded, or when we 
 * run out of numbers.
 */
static void
mersenne_fill_table(ESL_RANDOMNESS *r)
{
  uint32_t y;
  int      z;
  static uint32_t mag01[2] = { 0x0, 0x9908b0df };

  for (z = 0; z < 227; z++)	/* 227 = N-M = 624-397 */
    {
      y = (r->mt[z] & 0x80000000) | (r->mt[z+1] & 0x7fffffff);
      r->mt[z] = r->mt[z+397] ^ (y>>1) ^ mag01[y & 0x1];
    }
  for (; z < 623; z++)
    {
      y = (r->mt[z] & 0x80000000) | (r->mt[z+1] & 0x7fffffff);
      r->mt[z] = r->mt[z-227] ^ (y>>1) ^ mag01[y & 0x1];
    }
  y = (r->mt[623] & 0x80000000) | (r->mt[0] & 0x7fffffff);
  r->mt[623] = r->mt[396] ^ (y>>1) ^ mag01[y & 0x1];
  r->mti = 0;
  return;
}


/* choose_arbitrary_seed()
 * Return a 'quasirandom' seed > 0.
 * This could be a *lot* better than it is now; see RFC1750
 * for a discussion of securely seeding RNGs.
 */
static uint32_t
choose_arbitrary_seed(void)
{
  uint32_t a = (uint32_t) time ((time_t *) NULL);
  uint32_t b = 87654321;	/* anything nonzero */
  uint32_t c = 12345678;	/* anything nonzero. jenkins' mix3 needs 3 numbers; add an arbitrary one */
  uint32_t seed;
#ifdef HAVE_GETPID
  b  = (uint32_t) getpid();	  /* preferable b choice, if we have POSIX getpid(); else both b,c are arbitrary */
#endif
  seed = jenkins_mix3(a,b,c);	  /* try to decorrelate closely spaced choices of pid/time */

  return (seed == 0) ? 42 : seed; /* 42 is entirely arbitrary, just to avoid seed==0. */
}

/* jenkins_mix3()
 * 
 * from Bob Jenkins: given a,b,c, generate a number that's distributed
 * reasonably uniformly on the interval 0..2^32-1 even for closely
 * spaced choices of a,b,c.
 */
static uint32_t 
jenkins_mix3(uint32_t a, uint32_t b, uint32_t c)
{
  a -= b; a -= c; a ^= (c>>13);		
  b -= c; b -= a; b ^= (a<<8); 
  c -= a; c -= b; c ^= (b>>13);
  a -= b; a -= c; a ^= (c>>12);
  b -= c; b -= a; b ^= (a<<16);
  c -= a; c -= b; c ^= (b>>5); 
  a -= b; a -= c; a ^= (c>>3); 
  b -= c; b -= a; b ^= (a<<10);
  c -= a; c -= b; c ^= (b>>15);
  return c;
}
/*----------- end of esl_random() --------------*/



/*****************************************************************
 *# 3. Other fundamental sampling (including Gaussian, gamma)
 *****************************************************************/ 

/* Function: esl_rnd_UniformPositive()
 * Synopsis: Generate a uniform positive random deviate on interval (0,1).
 * Incept:   SRE, Wed Jul 14 13:31:23 2004 [St. Louis]
 *
 * Purpose:  Same as <esl_random()>, but assure $0 < x < 1$;
 *           (positive uniform deviate).
 */
double
esl_rnd_UniformPositive(ESL_RANDOMNESS *r)
{
  double x;
  do { x = esl_random(r); } while (x == 0.0);
  return x;
}


/* Function:  esl_rnd_Gaussian()
 * Synopsis:  Generate a Gaussian-distributed sample.
 * Incept:    SRE, Wed Jul 14 13:50:36 2004 [St. Louis]
 *
 * Purpose:   Pick a Gaussian-distributed random variable
 *            with a given <mean> and standard deviation <stddev>, and
 *            return it. 
 *            
 *            Implementation is derived from the public domain
 *            RANLIB.c <gennor()> function, written by Barry W. Brown
 *            and James Lovato (M.D. Anderson Cancer Center, Texas
 *            USA) using the method described in
 *            \citep{AhrensDieter73}.
 * 
 * Method:    Impenetrability of the code is to be blamed on 
 *            FORTRAN/f2c lineage.
 *
 * Args:      r      - ESL_RANDOMNESS object
 *            mean   - mean of the Gaussian we're sampling from
 *            stddev - standard deviation of the Gaussian     
 */
double
esl_rnd_Gaussian(ESL_RANDOMNESS *r, double mean, double stddev)
{
  long   i;
  double snorm,u,s,ustar,aa,w,y,tt;

  /* These static's are threadsafe: they are magic constants
   * we will not touch.
   */
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,    
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };

  u = esl_random(r);
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
   * START CENTER
   */
  ustar = u-(double)i;
  aa = a[i-1];
S40:
  if (ustar <= t[i-1]) goto S60;
  w = (ustar - t[i-1]) * h[i-1];
S50:
  /*
   * EXIT   (BOTH CASES)
   */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return (stddev*snorm + mean);
S60:
  /*
   * CENTER CONTINUED
   */
  u = esl_random(r);
  w = u*(a[i]-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
S70:
  tt = u;
  ustar = esl_random(r);
S80:
  if(ustar > tt) goto S50;
  u = esl_random(r);
  if(ustar >= u) goto S70;
  ustar = esl_random(r);
  goto S40;
S100:
  /*
   * START TAIL
   */
  i = 6;
  aa = a[31];
  goto S120;
S110:
  aa += d[i-1];
  i += 1;
S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
S140:
  w = u*d[i-1];
  tt = (0.5*w+aa)*w;
  goto S160;
S150:
  tt = u;
S160:
  ustar = esl_random(r);
  if(ustar > tt) goto S50;
  u = esl_random(r);
  if(ustar >= u) goto S150;
  u = esl_random(r);
  goto S140;
}



/* subfunctions that esl_rnd_Gamma() is going to call:
 */
static double
gamma_ahrens(ESL_RANDOMNESS *r, double a)	/* for a >= 3 */
{
  double V;			/* uniform deviates */
  double X,Y;
  double test;
  
  do {
    do {				/* generate candidate X */
      Y = tan(eslCONST_PI * esl_random(r)); 
      X = Y * sqrt(2.*a -1.) + a - 1.;
    } while (X <= 0.);
				/* accept/reject X */
    V    = esl_random(r);
    test = (1+Y*Y) * exp( (a-1.)* log(X/(a-1.)) - Y*sqrt(2.*a-1.));
  } while (V > test);
  return X;
}
static double
gamma_integer(ESL_RANDOMNESS *r, unsigned int a)	/* for small integer a, a < 12 */
{
  int    i;
  double U,X;

  U = 1.;
  for (i = 0; i < a; i++) 
    U *= esl_rnd_UniformPositive(r);
  X = -log(U);

  return X;
}
static double
gamma_fraction(ESL_RANDOMNESS *r, double a)	/* for fractional a, 0 < a < 1 */
{				/* Knuth 3.4.1, exercise 16, pp. 586-587 */
  double p, U, V, X, q;
  
  p = eslCONST_E / (a + eslCONST_E);
  do {
    U = esl_random(r);
    V = esl_rnd_UniformPositive(r);
    if (U < p) {
      X = pow(V, 1./a);
      q = exp(-X);
    } else {
      X = 1. - log(V);
      q = pow(X, a-1.);
    }
    U = esl_random(r);
  } while (U >= q);
  return X;
}


/* Function: esl_rnd_Gamma()
 * Synopsis: Returns a random deviate from a Gamma(a, 1) distribution.
 * Incept:   SRE, Wed Apr 17 13:10:03 2002 [St. Louis]
 *
 * Purpose:  Return a random deviate distributed as Gamma(a, 1.)
 *           \citep[pp. 133--134]{Knu-81a}.
 *           
 *           The implementation follows not only Knuth \citep{Knu-81a},
 *           but also relied on examination of the implementation in
 *           the GNU Scientific Library (libgsl) \citep{Galassi06}.
 *
 * Args:     r      - random number generation seed
 *           a      - order of the gamma function; a > 0
 *
 * Throws:   <eslEINVAL> for $a <= 0$.
 */
double
esl_rnd_Gamma(ESL_RANDOMNESS *r, double a)
{
  double aint;

  aint = floor(a);
  if (a == aint && a < 12.) 
    return gamma_integer(r, (unsigned int) a);
  else if (a > 3.) 
    return gamma_ahrens(r, a);
  else if (a < 1.) 
    return gamma_fraction(r, a);
  else 
    return gamma_integer(r, aint) + gamma_fraction(r, a-aint);
  /*NOTREACHED*/
  return eslOK;
}


/*****************************************************************
 *# 4. Multinomial sampling from discrete probability n-vectors
 *****************************************************************/ 

/* Function:  esl_rnd_DChoose()
 * Synopsis:  Return random choice from discrete multinomial distribution.          
 *
 * Purpose:   Make a random choice from a normalized discrete
 *            distribution <p> of <N> elements, where <p>
 *            is double-precision. Returns the index of the
 *            selected element, $0..N-1$.
 *            
 *            <p> must be a normalized probability distribution
 *            (i.e. must sum to one). Sampling distribution is
 *            undefined otherwise: that is, a choice will always
 *            be returned, but it might be an arbitrary one.
 *
 *            All $p_i$ must be $>>$ <DBL_EPSILON> in order to 
 *            have a non-zero probability of being sampled.
 *
 *            <esl_rnd_FChoose()> is the same, but for floats in <p>.
 *
 * Note:      Why the while (1) loop? Very rarely, because of machine
 *            floating point representation, our roll is "impossibly" 
 *            >= total sum, even though any roll of esl_random() is 
 *            < 1.0 and the total sum is supposed to be 1.0 by
 *            definition. This can happen when the total_sum is not
 *            really 1.0, but something just less than that in the 
 *            machine representation, and the roll happens to also be 
 *            very very close to 1. I have not examined this analytically, 
 *            but empirically, it occurs at a frequency of about 1/10^8
 *            as measured for bug #sq5... which suggests it is on the
 *            order of machine epsilon (not surprisingly). The while 
 *            loop makes you go around and try again; it must eventually
 *            succeed.
 *            
 *            The while() loop then makes the function vulnerable to
 *            an infinite loop if <p> sums to <=0 -- which shouldn't
 *            happen, but we shouldn't infinite loop if it does,
 *            either.  That's why there's a check on the sum of
 *            <p>. We return -1 in this case, a non-standard error code
 *            for Easel.
 * 
 * Throws:    -1 on failure. (This is a non-standard error code for Easel,
 *            but the only way an error can happen is if <p> isn't a 
 *            normalized probability distribution.)
 */
int
esl_rnd_DChoose(ESL_RANDOMNESS *r, const double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;

  while (1) {	/* see note in header about this while() */
    for (i = 0; i < N; i++)
      {
	sum += p[i];
	if (roll < sum) return i;  /* success! */
      }
    if (sum < 0.99) ESL_EXCEPTION(-1, "unnormalized distribution");    /* avoid inf loop */
  }
  /*NOTREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}
int
esl_rnd_FChoose(ESL_RANDOMNESS *r, const float *p, int N)
{
  float  roll;                  /* random fraction */
  float  sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = esl_random(r);
  sum     = 0.0;

  while (1) {	/* see note in header about this while() */
    for (i = 0; i < N; i++)
      {
	sum += p[i];
	if (roll < sum) return i; /* success */
      }
    if (sum < 0.99) ESL_EXCEPTION(-1, "unnormalized distribution");    /* avoid inf loop */
  }
  /*NOTREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}


/*****************************************************************
 * 5. Benchmark driver
 *****************************************************************/
#ifdef eslRANDOM_BENCHMARK
/*
   gcc -O3 -malign-double -o esl_random_benchmark -I. -L. -DeslRANDOM_BENCHMARK esl_random.c -leasel -lm
   ./esl_random_benchmark -N 1000000000
   ./esl_random_benchmark -f -N 1000000000
   ./esl_random_benchmark -r -N1000000
   ./esl_random_benchmark -fr -N 1000000000
                               esl_random()            esl_randomness_Init()
                           iter  cpu time  per call   iter  cpu time  per call  
                           ----  --------  --------   ---- ---------- ---------
   27 Dec 08 on wanderoo:  1e7    0.78s    78 nsec     1e6   2.08s     2.1 usec   ran2() from NR
   30 May 09 on wanderoo:  1e9    8.39s     8 nsec     1e6   5.98s     6.0 usec   Mersenne Twister
                           1e9    5.73s     6 nsec     1e8   2.51s     0.03 usec  Knuth

 */
#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS options[] = {
  /* name     type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",  eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-f",  eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "run fast version instead of MT19937",              0 },
  { "-r",  eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "benchmark _Init(), not just random()",             0 },
  { "-N",  eslARG_INT, "10000000",NULL, NULL,  NULL,  NULL, NULL, "number of trials",                                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "benchmarking speed of random number generator";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r       = (esl_opt_GetBoolean(go, "-f") == TRUE ? esl_randomness_CreateFast(42) : esl_randomness_Create(42));
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  int             N       = esl_opt_GetInteger(go, "-N");

  esl_stopwatch_Start(w);
  if (esl_opt_GetBoolean(go, "-r")) {
    while (N--) esl_randomness_Init(r, 42);
  } else {
    while (N--) esl_random(r);
  }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU Time: ");

  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*eslRANDOM_BENCHMARK*/
/*----------------- end, benchmark driver -----------------------*/




/*****************************************************************
 * 6. Unit tests.
 *****************************************************************/

#ifdef eslRANDOM_TESTDRIVE
#include "esl_vectorops.h"
#include "esl_stats.h"
#include "esl_dirichlet.h"
    
  
/* The esl_random() unit test:
 * a binned frequency test.
 */
static void
utest_random(ESL_RANDOMNESS *r, int n, int nbins, int be_verbose)
{
  char            msg[]  = "esl_random() unit test failed";
  int            *counts = NULL;
  double          X2p    = 0.;
  int             i;
  int             sample;
  double          X2, exp, diff;

  if ((counts = malloc(sizeof(int) * nbins)) == NULL) esl_fatal(msg);
  esl_vec_ISet(counts, nbins, 0);

  for (i = 0; i < n; i++)
    { 
      sample = esl_rnd_Roll(r, nbins);
      if (sample < 0 || sample >= nbins) esl_fatal(msg);
      counts[sample]++;
    }

  /* X^2 value: \sum (o_i - e_i)^2 / e_i */
  for (X2 = 0., i = 0; i < nbins; i++) {
    exp  = (double) n / (double) nbins;
    diff = (double) counts[i] - exp;
    X2 +=  diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal(msg);
  if (be_verbose) printf("random():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal(msg);

  free(counts);
  return;
}

/* The DChoose() and FChoose() unit tests.
 */
static void
utest_choose(ESL_RANDOMNESS *r, int n, int nbins, int be_verbose)
{
  double *pd = NULL;
  float  *pf = NULL;
  int    *ct = NULL;
  int     i;
  double  X2, diff, exp, X2p;

  if ((pd = malloc(sizeof(double) * nbins)) == NULL) esl_fatal("malloc failed"); 
  if ((pf = malloc(sizeof(float)  * nbins)) == NULL) esl_fatal("malloc failed");
  if ((ct = malloc(sizeof(int)    * nbins)) == NULL) esl_fatal("malloc failed");

  /* Sample a random multinomial probability vector.  */
  if (esl_dirichlet_DSampleUniform(r, nbins, pd) != eslOK) esl_fatal("dirichlet sample failed");
  esl_vec_D2F(pd, nbins, pf);

  /* Sample observed counts using DChoose(). */
  esl_vec_ISet(ct, nbins, 0);
  for (i = 0; i < n; i++)
    ct[esl_rnd_DChoose(r, pd, nbins)]++;

  /* X^2 test on those observed counts. */
  for (X2 = 0., i=0; i < nbins; i++) {
    exp = (double) n * pd[i];
    diff = (double) ct[i] - exp;
    X2 += diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal("chi square eval failed");
  if (be_verbose) printf("DChoose():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal("chi squared test failed");

  /* Repeat above for FChoose(). */
  esl_vec_ISet(ct, nbins, 0);
  for (i = 0; i < n; i++)
    ct[esl_rnd_FChoose(r, pf, nbins)]++;
  for (X2 = 0., i=0; i < nbins; i++) {
    exp = (double) n * pd[i];
    diff = (double) ct[i] - exp;
    X2 += diff*diff/exp;
  }
  if (esl_stats_ChiSquaredTest(nbins, X2, &X2p) != eslOK) esl_fatal("chi square eval failed");
  if (be_verbose) printf("FChoose():  \t%g\n", X2p);
  if (X2p < 0.01) esl_fatal("chi squared test failed");
  
  free(pd);
  free(pf);
  free(ct);
  return;
}
 

#endif /*eslRANDOM_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/


/*****************************************************************
 * 7. Test driver.
 *****************************************************************/
#ifdef eslRANDOM_TESTDRIVE
/* gcc -g -Wall -o esl_random_utest -L. -I. -DeslRANDOM_TESTDRIVE esl_random.c -leasel -lm
 */
#include "esl_config.h"

#include <stdio.h>

#include "easel.h"
#include "esl_dirichlet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  {"-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",               0},
  {"-b",  eslARG_INT,      "20", NULL, "n>0",NULL, NULL, NULL, "number of test bins",               0},
  {"-n",  eslARG_INT, "1000000", NULL, "n>0",NULL, NULL, NULL, "number of samples",                 0},
  {"-s",  eslARG_INT,      "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",     0},
  {"-v",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show verbose output",               0},
  {"--mtbits",eslARG_STRING,NULL,NULL, NULL, NULL, NULL, NULL, "save MT bit file for NIST benchmark",0},
  {"--kbits", eslARG_STRING,NULL,NULL, NULL, NULL, NULL, NULL, "save Knuth bit file for NIST benchmark",0},
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for random module";

static int save_bitfile(char *bitfile, ESL_RANDOMNESS *r, int n);

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go         = esl_getopts_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r1         = esl_randomness_Create(esl_opt_GetInteger(go, "-s"));
  ESL_RANDOMNESS *r2         = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  char           *mtbitfile  = esl_opt_GetString (go, "--mtbits");
  char           *kbitfile   = esl_opt_GetString (go, "--kbits");
  int             nbins      = esl_opt_GetInteger(go, "-b");
  int             n          = esl_opt_GetInteger(go, "-n");
  int             be_verbose = esl_opt_GetBoolean(go, "-v");

  utest_random(r1, n, nbins, be_verbose);
  utest_choose(r1, n, nbins, be_verbose);
  utest_random(r2, n, nbins, be_verbose);
  utest_choose(r2, n, nbins, be_verbose);

  if (mtbitfile) save_bitfile(mtbitfile, r1, n);
  if (kbitfile)  save_bitfile(kbitfile,  r2, n);

  esl_randomness_Destroy(r1);
  esl_randomness_Destroy(r2);
  esl_getopts_Destroy(go);
  return 0;
}

static int
save_bitfile(char *bitfile, ESL_RANDOMNESS *r, int n)
{
  FILE *fp = NULL;
  int b,i;
  long x;

  /* Open the file. 
   */
  if ((fp = fopen(bitfile, "w")) == NULL) 
    esl_fatal("failed to open %s for writing", bitfile);

  /* Sample <n> random numbers, output 32n random bits to the file.
   */
  for (i = 0; i < n; i++)
    {
      x = (r->type == eslRND_FAST ? knuth(r) : mersenne_twister(r)); /* generate a 32 bit random variate by MT19937 */
      for (b = 0; b < 32; b++) 
	{
	  if (x & 01) fprintf(fp, "1");
	  else        fprintf(fp, "0");
	  x >>= 1;
	}
      fprintf(fp, "\n");
    }
  fclose(fp);
  return eslOK;
}
#endif /*eslRANDOM_TESTDRIVE*/



/*****************************************************************
 * 8. Example.
 *****************************************************************/
#ifdef eslRANDOM_EXAMPLE
/*::cexcerpt::random_example::begin::*/
/* compile: gcc -g -Wall -I. -o random_example -DeslRANDOM_EXAMPLE esl_random.c easel.c -lm
 * run:     ./random_example 42
 */
#include <stdio.h>
#include "easel.h"
#include "esl_random.h"

int 
main(int argc, char **argv)
{
  long            seed = atoi(argv[1]);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(seed); 
  int             n    = 10;

  printf("RNG seed: %" PRIu32 "\n", esl_randomness_GetSeed(r));
  printf("A sequence of %d pseudorandom numbers:\n", n);
  while (n--)  printf("%f\n", esl_random(r));

  esl_randomness_Destroy(r);
  return 0;
}
/*::cexcerpt::random_example::end::*/
#endif /*eslRANDOM_EXAMPLE*/


/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/


