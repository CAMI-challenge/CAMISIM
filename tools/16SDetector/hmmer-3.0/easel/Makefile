# Top level Makefile for Easel
#
# SVN $Id: Makefile.in 561 2010-03-16 18:27:06Z eddys $

top_srcdir = ..
srcdir     = .
VPATH      =  ..

SHELL       = /bin/sh

# Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
bindir      = ${exec_prefix}/bin
mandir      = ${datarootdir}/man
pdfdir      = ${docdir}

# Compiler configuration
#
CC       = gcc -std=gnu99
CFLAGS   = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing -pthread
LDFLAGS  = 
SIMDFLAGS= -msse2
CPPFLAGS = 
LIBGSL   = 
LIBS     =  

# Other tools
#
COMPRESS = gzip
AR       = /usr/bin/ar rcv
LN       = ln
RANLIB   = ranlib

#######################################################################
## You should not need to modify below this line.
#######

PACKAGE     = HMMER
BASENAME    = hmmer
RELEASE     = h3.0
RELEASEDATE = "March 2010"
COPYRIGHT   = "Copyright (C) 2010 Howard Hughes Medical Institute."
LICENSELINE = "Freely distributed under the Janelia Farm Software License."


# Code manifests
#
HDRS = \
	easel.h\
	esl_alphabet.h\
	esl_cluster.h\
	esl_dirichlet.h\
	esl_distance.h\
	esl_dmatrix.h\
	esl_exponential.h\
	esl_fileparser.h\
	esl_gamma.h\
	esl_getopts.h\
	esl_gev.h\
	esl_gumbel.h\
	esl_histogram.h\
	esl_hmm.h\
	esl_hyperexp.h\
	esl_keyhash.h\
	esl_minimizer.h\
	esl_mixgev.h\
	esl_mpi.h\
	esl_msa.h\
	esl_msacluster.h\
	esl_msashuffle.h\
	esl_msaweight.h\
	esl_normal.h\
	esl_paml.h\
	esl_random.h\
	esl_randomseq.h\
	esl_ratematrix.h\
	esl_regexp.h\
	esl_rootfinder.h\
	esl_scorematrix.h\
	esl_sq.h\
	esl_sqio.h\
	esl_sqio_ascii.h\
	esl_sqio_ncbi.h\
	esl_sse.h\
	esl_ssi.h\
	esl_stack.h\
	esl_stats.h\
	esl_stopwatch.h\
	esl_stretchexp.h\
	esl_threads.h\
	esl_tree.h\
	esl_vectorops.h\
	esl_vmx.h\
	esl_weibull.h\
	esl_workqueue.h\
	esl_wuss.h

OBJS    = \
	easel.o\
	esl_alphabet.o\
	esl_cluster.o\
	esl_dirichlet.o\
	esl_distance.o\
	esl_dmatrix.o\
	esl_exponential.o\
	esl_fileparser.o\
	esl_gamma.o\
	esl_getopts.o\
	esl_gev.o\
	esl_gumbel.o\
	esl_histogram.o\
	esl_hmm.o\
	esl_hyperexp.o\
	esl_keyhash.o\
	esl_minimizer.o\
	esl_mixgev.o\
	esl_mpi.o\
	esl_msa.o\
	esl_msacluster.o\
	esl_msashuffle.o\
	esl_msaweight.o\
	esl_normal.o\
	esl_paml.o\
	esl_random.o\
	esl_randomseq.o\
	esl_ratematrix.o\
	esl_regexp.o\
	esl_rootfinder.o\
	esl_scorematrix.o\
	esl_sq.o\
	esl_sqio.o\
	esl_sqio_ascii.o\
	esl_sqio_ncbi.o\
	esl_sse.o\
	esl_ssi.o\
	esl_stack.o\
	esl_stats.o\
	esl_stopwatch.o\
	esl_stretchexp.o\
	esl_threads.o\
	esl_tree.o\
	esl_vectorops.o\
	esl_vmx.o\
	esl_weibull.o\
	esl_workqueue.o\
	esl_wuss.o
#	esl_swat.o


UTESTS =\
	easel_utest\
	esl_alphabet_utest\
	esl_cluster_utest\
	esl_dirichlet_utest\
	esl_distance_utest\
	esl_dmatrix_utest\
	esl_exponential_utest\
	esl_fileparser_utest\
	esl_gamma_utest\
	esl_getopts_utest\
	esl_gumbel_utest\
	esl_histogram_utest\
	esl_hyperexp_utest\
	esl_keyhash_utest\
	esl_msa_utest\
	esl_msacluster_utest\
	esl_msaweight_utest\
	esl_random_utest\
	esl_randomseq_utest\
	esl_ratematrix_utest\
	esl_regexp_utest\
	esl_rootfinder_utest\
	esl_scorematrix_utest\
	esl_sq_utest\
	esl_sqio_utest\
	esl_sse_utest\
	esl_ssi_utest\
	esl_stack_utest\
	esl_stats_utest\
	esl_stretchexp_utest\
	esl_tree_utest\
	esl_vectorops_utest\
	esl_weibull_utest\
	esl_wuss_utest
#	gev_utest\
#	minimizer_utest\
#	mixgev_utest\
#	mpi_utest\
#	normal_utest\
#	paml_utest\
#	stopwatch_utest\
#	swat_utest\


all:
	${MAKE} easel_lib
	${MAKE} easel_progs

check:
	${MAKE} easel_lib
	${MAKE} easel_progs
	${MAKE} easel_utests
	cd testsuite; ${MAKE} run_sqc

dev:
	${MAKE} easel_lib
	${MAKE} easel_progs
	${MAKE} easel_utests

easel_lib:    libeasel.a
easel_utests: ${UTESTS}

easel_progs:
	cd miniapps; ${MAKE} easel_progs

run_sqc:
	cd testsuite; ${MAKE} run_sqc

libeasel.a: $(OBJS)
	$(AR) libeasel.a $(OBJS)
	$(RANLIB) libeasel.a
	chmod 644 libeasel.a

${OBJS}: ${HDRS}

.c.o: 
	${CC} -I. -I${srcdir} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -o $@ -c $<		


# Unit test compilation:
# 
# To compile unit tests, we need to construct a command line
# from module basename, module's testdriver flag, and module's
# C file. This is accomplished via a complex-looking shell scripture;
# xref J1/81.
#   $@           =  driver name            esl_msa_utest     easel_utest
#   $${BASENAME} =  base module name       msa               easel
#   $${DFLAG}    =  driver flag            eslMSA_TESTDRIVE  eslEASEL_TESTDRIVE
#   $${DFILE}    =  driver's C file        esl_msa.c         easel.c
#
${UTESTS}:  libeasel.a
	@BASENAME=`echo $@ | sed -e 's/_utest//'| sed -e 's/^esl_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=esl$${DFLAG}_TESTDRIVE ;\
	if test $@ = "easel_utest"; then \
           DFILE=${srcdir}/easel.c ;\
        else \
           DFILE=${srcdir}/esl_$${BASENAME}.c ;\
        fi;\
	echo ${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} -o $@ -I. -I${srcdir} -L. -D$${DFLAG} $${DFILE} -leasel -lm;\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} -o $@ -I. -I${srcdir} -L. -D$${DFLAG} $${DFILE} -leasel -lm


# install:   installs libeasel.a in ${libdir}/
#            installs Easel headers in ${includedir}/
#            installs miniapp binaries in ${bindir}/
#
# All install paths are prefaced with DESTDIR/ to facilitate
# building RPMS in Linux distros. Usually DESTDIR is empty.
#
install:
	mkdir -p ${DESTDIR}${libdir}
	mkdir -p ${DESTDIR}${includedir}
	mkdir -p ${DESTDIR}${bindir}
	cp libeasel.a ${DESTDIR}${libdir}/
	for file in ${HDRS}; do\
	   cp $$file ${DESTDIR}${includedir}/;\
	done
	(cd miniapps; make install)

uninstall:
	-rm -f ${DESTDIR}${libdir}/libeasel.a
	for file in ${HDRS}; do\
	   rm -f ${DESTDIR}${includedir}/$$file;\
	done
	(cd miniapps; make uninstall)

clean:
	if test -d documentation; then  (cd documentation; ${MAKE} clean); fi
	if test -d testsuite;     then  (cd testsuite;     ${MAKE} clean); fi
	if test -d miniapps;      then  (cd miniapps;      ${MAKE} clean); fi
	-rm -f ${OBJS} libeasel.a
	-rm -f ${UTESTS}
	-rm -f *~ TAGS
	-rm -f *.gcno *.gcda *.gcov
	-rm -f core.[0-9]*
	-rm -f config.log config.status

distclean:
	${MAKE} clean
	-rm esl_config.h
	-rm -rf autom4te.cache
	if test -d documentation; then (cd documentation; ${MAKE} distclean); fi
	if test -d testsuite;     then (cd testsuite;     ${MAKE} distclean); fi
	if test -d miniapps;      then (cd miniapps;      ${MAKE} distclean); fi
	-rm LICENSE.csh
	-rm Makefile

#######
## Targets below this line are meant for developers, not users.
## Some may only work in context of a developer's Subversion working directory.
#######

# 'make REPOSITORY=... dist' : create a new distribution tarball
#
dist:
	# Make sure we pointed the make dist to a subversion repository
	@if test x${REPOSITORY} = 'x'; then\
	  echo "REPOSITORY isn't set; use make REPOSITORY=... dist";\
	  exit 1;\
	fi
	# Delete any old versions of the same release
	@if test -d ../${BASENAME}-${RELEASE};        then rm -rf ../${BASENAME}-${RELEASE};        fi
	@if test -e ../${BASENAME}-${RELEASE}.tar;    then rm -f  ../${BASENAME}-${RELEASE}.tar;    fi
	@if test -e ../${BASENAME}-${RELEASE}.tar.Z;  then rm -f  ../${BASENAME}-${RELEASE}.tar.Z;  fi
	@if test -e ../${BASENAME}-${RELEASE}.tar.gz; then rm -f  ../${BASENAME}-${RELEASE}.tar.gz; fi
	# Export the release from Subversion control
	svn export ${REPOSITORY} ../${BASENAME}-${RELEASE}
	# Make the configure script from configure.ac
	(cd ../${BASENAME}-${RELEASE};  autoconf)
	# We'll need some basic Makefiles
	(cd ../${BASENAME}-${RELEASE}/; ./configure)
	# Set file permissions.
	find ../${BASENAME}-${RELEASE} -type f -exec chmod +r {} \;
	# Compile the PDF documentation and move it to top level dir
	(cd ../${BASENAME}-${RELEASE}/documentation/; make pdf)
	mv ../${BASENAME}-${RELEASE}/documentation/current.pdf  ../${BASENAME}-${RELEASE}/Userguide.pdf
	# Attach license stamps by replacing LICENSE keyword everywhere
	find ../${BASENAME}-${RELEASE} -type f -exec devkit/sedition-pp LICENSE LICENSE.tag {} \;
	# Other keyword replacements
	devkit/sedition RELEASE ${RELEASE} RELEASEDATE ${RELEASEDATE} COPYRIGHT ${COPYRIGHT} ../${BASENAME}-${RELEASE}/00README
	devkit/sedition RELEASE ${RELEASE}                                                   ../${BASENAME}-${RELEASE}/INSTALL
	# Cleanup
	(cd ../${BASENAME}-${RELEASE}/; make implode)
	(cd ../${BASENAME}-${RELEASE}/; make distclean)
	# Packaging
	(cd ..;	tar cvf ${BASENAME}-${RELEASE}.tar ${BASENAME}-${RELEASE})
	(cd ..; ${COMPRESS} ${BASENAME}-${RELEASE}.tar)

# `make dcheck`: run developer testsuite, save output in testsuite/dcheck.report
dcheck:
	(cd testsuite; make dcheck)

# The "implode" target is used when we've got Easel as a subdirectory
# of a larger project, like Infernal or HMMER, and we're building a
# distribution for that larger project. We remove any "personal
# information" from Easel. Eventually that'd include our bug tracking,
# for instance, but right now there's nothing to do.
implode:
	-rm -f 00CHECKLIST
	-rm -f LICENSE.tag

TAGS:
	etags *.c *.h *.tex documentation/*.tex miniapps/*.c miniapps/*.man Makefile.in

# "supertags": a tags table for all Easel, HMMER, Infernal source.
#              This lets me make changes in Easel, then propagate 
#              necessary modifications throughout code that depends
#              on Easel. 
#              The target depends on how I have my source directories
#              laid out (~src/hmmer, ~src/infernal)
supertags:
	etags *.c *.h *.tex documentation/*.tex 
	etags -a miniapps/*.c miniapps/*.man Makefile.in
	etags -a ~/src/hmmer/src/*.c ~/src/hmmer/src/*.h 
	etags -a ~/src/hmmer/src/impl_sse/*.c ~/src/hmmer/src/impl_sse/*.h 
	etags -a ~/src/hmmer/src/impl_vmx/*.c ~/src/hmmer/src/impl_vmx/*.h 
	etags -a ~/src/hmmer/Makefile.in ~/src/hmmer/src/Makefile.in 
	etags -a ~/src/infernal/src/*.c ~/src/infernal/src/*.h 
	etags -a ~/src/infernal/hmmer/src/*.c ~/src/infernal/hmmer/src/*.h 
	etags -a ~/src/infernal/testsuite/*.c 
	etags -a ~/src/infernal/Makefile.in ~/src/infernal/src/Makefile.in 

# magic SVN for setting keyword ID replacement on a new module foo:
# svn propset svn:keywords "Id" foo.[ch]

#################################################################
## For building as a library with an SSU-ALIGN distribution.
##
ssu:
	cp miniapps/esl-alimanip miniapps/ssu-esl-alimanip
	cp miniapps/esl-alimask miniapps/ssu-esl-alimask
	cp miniapps/esl-alimerge miniapps/ssu-esl-alimerge
	cp miniapps/esl-alistat miniapps/ssu-esl-alistat
	cp miniapps/esl-reformat miniapps/ssu-esl-reformat
	cp miniapps/esl-seqstat miniapps/ssu-esl-seqstat
	cp miniapps/esl-sfetch miniapps/ssu-esl-sfetch
	cp miniapps/esl-ssdraw miniapps/ssu-esl-ssdraw
	cp miniapps/esl-weight miniapps/ssu-esl-weight

################################################################
# Easel - a library of C functions for biological sequence analysis
# Version h3.0; March 2010
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# Easel is distributed under the Janelia Farm Software License, a BSD
# license. See the LICENSE file for more details.
################################################################
