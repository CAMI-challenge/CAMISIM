/* easel/esl_config.h.  Generated from esl_config.h.in by configure.  */
/* esl_config.h.in  [input to configure]
 * 
 * System-dependent configuration of Easel, by autoconf.
 * 
 * This file should be included in all Easel .c files before
 * anything else, because it may set #define's that control
 * behaviour of system includes and system libraries. An example
 * is large file support.
 * 
 * SVN $Id$
 * SRE, Fri Mar  3 08:03:32 2006 [St. Louis]
 */
#ifndef ESL_CONFIG_INCLUDED
#define ESL_CONFIG_INCLUDED

/* Version info.
 */
#define EASEL_VERSION "h3.0"
#define EASEL_DATE "March 2010"
#define EASEL_COPYRIGHT "Copyright (C) 2010 Howard Hughes Medical Institute."
#define EASEL_LICENSE "Freely distributed under the Janelia Farm Software License."

/* Large file support
 * Must precede any header file inclusion.
 */
/* #undef _FILE_OFFSET_BITS */
/* #undef _LARGE_FILES */
/* #undef _LARGEFILE_SOURCE */

/* Debugging verbosity (0=none;3=most verbose)
 */
#define eslDEBUGLEVEL 0

/* System headers
 */
#define HAVE_UNISTD_H 1
#define HAVE_STDINT_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_SYS_SYSCTL_H 1
#define HAVE_ENDIAN_H 1

/* Types
 */
/* #undef WORDS_BIGENDIAN */
/* #undef int8_t */
/* #undef int16_t */
/* #undef int32_t */
/* #undef int64_t */
/* #undef uint8_t */
/* #undef uint16_t */
/* #undef uint32_t */
/* #undef uint64_t */
/* #undef off_t */

/* Optional packages
 */
/* #undef HAVE_LIBGSL */

/* Optional parallel implementation support
 */
#define HAVE_SSE2 1
/* #undef HAVE_VMX */
/* #undef HAVE_MPI */

#define HAVE_SSE2_CAST 1

/* Functions
 */
#define HAVE_MKSTEMP 1
#define HAVE_POPEN 1
#define HAVE_STRCASECMP 1
#define HAVE_TIMES 1
#define HAVE_GETPID 1
#define HAVE_SYSCTL 1
#define HAVE_SYSCONF 1
#define HAVE_FSEEKO 1

#define HAVE_PTHREADS 1

/*****************************************************************
 * Available augmentations.
 * 
 * If you grab a single module from Easel to use it by itself,
 * leave all these #undef'd; you have no augmentations.
 * 
 * If you grab additional Easel .c files, you can enable any
 * augmentations they provide to other modules by #defining the
 * modules you have below. Alternatively, you can -D them on
 * the compile line, as in cc -DeslAUGMENT_SSI -DeslAUGMENT_MSA.
 * 
 * If you compile and install the complete Easel library, all of these
 * get #defined automatically by ./configure, plus the eslLIBRARY flag
 * which means the full library with all augmentations is
 * available. So, if you steal files from an installed library, just
 * set these all back to #undef (depending on which files you have).
 *****************************************************************/
#define eslLIBRARY 1

#ifndef eslLIBRARY
/* #undef eslAUGMENT_ALPHABET */
/* #undef eslAUGMENT_NCBI */
/* #undef eslAUGMENT_DMATRIX */
/* #undef eslAUGMENT_FILEPARSER */
/* #undef eslAUGMENT_GEV */
/* #undef eslAUGMENT_GUMBEL */
/* #undef eslAUGMENT_HISTOGRAM */
/* #undef eslAUGMENT_KEYHASH */
/* #undef eslAUGMENT_MINIMIZER */
/* #undef eslAUGMENT_MSA */
/* #undef eslAUGMENT_RANDOM */
/* #undef eslAUGMENT_SSI */
/* #undef eslAUGMENT_STATS */
#endif

#ifdef eslLIBRARY
#define eslAUGMENT_ALPHABET
#define eslAUGMENT_NCBI
#define eslAUGMENT_DMATRIX
#define eslAUGMENT_FILEPARSER
#define eslAUGMENT_GEV
#define eslAUGMENT_GUMBEL
#define eslAUGMENT_HISTOGRAM
#define eslAUGMENT_KEYHASH 
#define eslAUGMENT_MINIMIZER
#define eslAUGMENT_MSA		
#define eslAUGMENT_RANDOM
#define eslAUGMENT_SSI
#define eslAUGMENT_STATS
#endif


#endif /*ESL_CONFIG_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version h3.0; March 2010
 * Copyright (C) 2010 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
