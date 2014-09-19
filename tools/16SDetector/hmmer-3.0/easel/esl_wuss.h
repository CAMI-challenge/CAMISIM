/* wuss.h
 * RNA secondary structure markup in WUSS notation.
 * 
 * SVN $Id: esl_wuss.h 83 2005-12-13 20:54:07Z eddy $
 * SRE, Tue Feb 15 10:15:28 2005
 */
#ifndef eslWUSS_INCLUDED
#define eslWUSS_INCLUDED


extern int esl_wuss2ct(char *ss, int len, int *ct);
extern int esl_ct2wuss(int *ct, int n, char *ss);
extern int esl_wuss2kh(char *ss, char *kh);
extern int esl_kh2wuss(char *kh, char *ss);
extern int esl_wuss_full(char *oldss, char *newss);
extern int esl_wuss_nopseudo(char *ss1, char *ss2);


#endif /*eslWUSS_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
