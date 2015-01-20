/* Shuffling or bootstrapping multiple sequence alignments.
 * 
 * SRE, Tue Jan 22 09:18:09 2008 [Market Street Cafe, Leesburg]
 * SVN $Id: esl_msashuffle.h 249 2008-04-24 19:19:50Z eddys $
 */
#ifndef ESL_MSASHUFFLE_INCLUDED
#define ESL_MSASHUFFLE_INCLUDED

#include "esl_random.h"
#ifdef eslAUGMENT_ALPHABET
#include "esl_alphabet.h"
#endif

extern int esl_msashuffle_Shuffle  (ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *shuf);
extern int esl_msashuffle_Bootstrap(ESL_RANDOMNESS *r, ESL_MSA *msa, ESL_MSA *bootsample);

#ifdef eslAUGMENT_ALPHABET
extern int esl_msashuffle_CQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, char    *x, char    *y, char    *xs, char    *ys);
extern int esl_msashuffle_XQRNA(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, ESL_DSQ *x, ESL_DSQ *y, ESL_DSQ *xs, ESL_DSQ *ys);
#endif /*eslAUGMENT_ALPHABET*/

#endif /*ESL_MSASHUFFLE_INCLUDED*/
