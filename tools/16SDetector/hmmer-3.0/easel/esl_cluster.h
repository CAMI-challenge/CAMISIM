/* Generalized single linkage clustering.
 * 
 * SRE, Mon Jan  7 09:40:06 2008 [Janelia]
 * SVN $Id: esl_cluster.h 231 2008-03-25 14:43:57Z eddys $
 */
#ifndef ESL_CLUSTER_INCLUDED
#define ESL_CLUSTER_INCLUDED

extern int esl_cluster_SingleLinkage(void *base, size_t n, size_t size, 
				     int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
				     int *workspace, int *assignments, int *ret_C);
#endif /*ESL_CLUSTER_INCLUDED*/
