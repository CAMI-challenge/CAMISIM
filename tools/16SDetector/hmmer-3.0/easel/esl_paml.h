/* PAML interface
 *
 *   "Phylogenetic Analysis by Maximum Likelihood"
 *   Ziheng Yang
 *   http://abacus.gene.ucl.ac.uk/software/paml.html
 *   [Yang97]
 * 
 *           incept: SRE, Tue Jul 13 13:20:08 2004 [St. Louis]
 * upgrade to Easel: SRE, Thu Mar  8 13:26:20 2007 [Janelia]
 * SVN $Id: esl_paml.h 158 2007-03-15 20:03:05Z eddys $
 */

#ifndef ESL_PAML_INCLUDED
#define ESL_PAML_INCLUDED

#include <stdio.h>
#include <esl_dmatrix.h>

extern int esl_paml_ReadE(FILE *fp, ESL_DMATRIX *E, double *pi);


#endif /*ESL_PAML_INCLUDED*/
