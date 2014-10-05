/* esl_msaweight.h
 * Sequence weighting algorithms.
 * 
 * SVN $Id: esl_msaweight.h 302 2008-10-30 14:26:46Z eddys $
 * SRE, Sun Nov  5 09:11:13 2006 [Janelia]
 */
#ifndef ESL_MSAWEIGHT_INCLUDED
#define ESL_MSAWEIGHT_INCLUDED

#include <esl_msa.h>

extern int esl_msaweight_GSC(ESL_MSA *msa);
extern int esl_msaweight_PB(ESL_MSA *msa);
extern int esl_msaweight_BLOSUM(ESL_MSA *msa, double maxid);
extern int esl_msaweight_IDFilter(const ESL_MSA *msa, double maxid, ESL_MSA **ret_newmsa);


#endif /*ESL_MSAWEIGHT_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
