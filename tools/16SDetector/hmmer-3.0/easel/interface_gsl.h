#ifdef ESL_WITH_GSL
/* interface_gsl.h
 * Easel's interfaces to the GNU Scientific Library
 * 
 * SRE, Tue Jul 13 15:36:48 2004
 * SVN $Id: interface_gsl.h 11 2005-01-06 11:44:17Z eddy $
 */
#ifndef ESL_INTERFACE_GSL_INCLUDED
#define ESL_INTERFACE_GSL_INCLUDED

#include <stdlib.h>
#include <easel/easel.h>
#include <easel/dmatrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_eigen.h>

extern int esl_GSL_MatrixInversion(ESL_DMATRIX *A, ESL_DMATRIX **ret_Ai);


#endif /*ESL_INTERFACE_GSL_INCLUDED*/
#endif /*ESL_WITH_GSL*/
