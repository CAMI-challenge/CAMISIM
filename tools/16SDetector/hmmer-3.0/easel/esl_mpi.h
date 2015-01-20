/* Support for MPI parallelization.
 * 
 * SRE, Sat Jun  2 09:07:25 2007 [Janelia]
 * SVN $Id: esl_mpi.h 293 2008-09-19 19:08:30Z eddys $
 */

#if defined(HAVE_MPI) && defined(eslLIBRARY)
#ifndef eslMPI_INCLUDED
#define eslMPI_INCLUDED
#include "mpi.h"

#include "esl_alphabet.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_stopwatch.h"

/* 1. Communicating optional arrays */
extern int esl_mpi_PackOpt(void *inbuf, int incount, MPI_Datatype type, void *pack_buf, 
			   int pack_buf_size, int *position, MPI_Comm comm);
extern int esl_mpi_PackOptSize(void *inbuf, int incount, MPI_Datatype type, MPI_Comm comm, int *ret_n);
extern int esl_mpi_UnpackOpt(void *pack_buf, int pack_buf_size, int *pos, void **outbuf, 
			     int *opt_n, MPI_Datatype type, MPI_Comm comm);

/* 2. Communicating ESL_SQ (single sequences) */
extern int esl_sq_MPISend(ESL_SQ *sq, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int esl_sq_MPIPackSize(ESL_SQ *sq, MPI_Comm comm, int *ret_n);
extern int esl_sq_MPIPack(ESL_SQ *sq, char *buf, int n, int *pos, MPI_Comm comm);
extern int esl_sq_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_SQ **ret_sq);
extern int esl_sq_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, 
			  char **buf, int *nalloc, ESL_SQ **ret_sq);

/* 3. Communicating ESL_MSA (multiple sequence alignments) */
extern int esl_msa_MPISend(const ESL_MSA *msa, int dest, int tag, MPI_Comm comm, char **buf, int *nalloc);
extern int esl_msa_MPIPackSize(const ESL_MSA *msa, MPI_Comm comm, int *ret_n);
extern int esl_msa_MPIPack(const ESL_MSA *msa, char *buf, int n, int *position, MPI_Comm comm);
extern int esl_msa_MPIUnpack(const ESL_ALPHABET *abc, char *buf, int n, int *pos, MPI_Comm comm, ESL_MSA **ret_msa);
extern int esl_msa_MPIRecv(int source, int tag, MPI_Comm comm, const ESL_ALPHABET *abc, char **buf, int *nalloc, ESL_MSA **ret_msa);

/* 4. Communicating ESL_STOPWATCH (process timing) */
extern int esl_stopwatch_MPIReduce(ESL_STOPWATCH *w, int root, MPI_Comm comm);


#endif /*eslMPI_INCLUDED*/
#endif /*HAVE_MPI && eslLIBRARY*/

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
