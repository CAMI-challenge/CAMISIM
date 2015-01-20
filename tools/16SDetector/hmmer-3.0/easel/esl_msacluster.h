/* Clustering sequences in an MSA by % identity.
 * 
 * SRE, Sun Nov  5 10:08:14 2006 [Janelia]
 * SVN $Id: esl_msacluster.h 238 2008-03-28 11:53:19Z eddys $
 */
#ifndef ESL_MSACLUSTER_INCLUDED
#define ESL_MSACLUSTER_INCLUDED

extern int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, 
					int **opt_c, int **opt_nin, int *opt_nc);

#endif /*ESL_MSACLUSTER_INCLUDED*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
