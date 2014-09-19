###############################################################
# Makefile for SQUID library
# CVS $Id: Makefile.in,v 1.53 2003/06/13 20:05:31 eddy Exp $
#
# Note: The autoconf variables in this file must be coordinated
#       with HMMER. HMMER creates a Makefile from this 
#       Makefile.in using its own configure script, not SQUID's.
#################################################################
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2003 Washington University School of Medicine
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
#################################################################


### Installation points
###
prefix      = /usr/local
exec_prefix = ${prefix}
BINDIR      = ${exec_prefix}/bin
MANDIR      = ${prefix}/man
INCLUDEDIR  = ${prefix}/include
LIBDIR      = ${exec_prefix}/lib
SCRIPTDIR   = ${exec_prefix}/bin

## your compiler and compiler flags
#
CC       = gcc
CFLAGS   = -g -O2
CPPFLAGS = 
LDFLAGS  = 
DEFS     = -DHAVE_CONFIG_H
LIBS     =  -lm

## Archiver command
#
AR     = /usr/bin/ar rcv
RANLIB = ranlib

## instructions for installing man pages
#
INSTMAN   = cp
MANSUFFIX = 1

# Configuration for compiling in optional PVM support
#
PVMLIBDIR = 
PVMINCDIR = 
PVMLIBS   = 

#######
## You should not need to modify below this line
#######
SHELL       = /bin/sh
BASENAME    = squid
PACKAGE     = SQUID
RELEASE     = 1.9g
RELCODE     = squid1_9g
RELEASEDATE = "January 2003"
SEEALSO     = "http:\/\/www.genetics.wustl.edu\/eddy\/software\/\#squid"
FTPDIR      = /nfs/ftp/eddy/software/
LICENSETAG  = gnu
COMPRESS    = gzip

PROGS =	afetch\
	alistat\
	compalign\
	compstruct\
	revcomp\
	seqsplit\
	seqstat\
	sfetch\
	shuffle\
	sindex\
	sreformat\
	translate\
	weight

MANS =  afetch\
	alistat\
	seqstat\
	sfetch\
	shuffle\
	sreformat\

READMES = 00README INSTALL Makefile.in 

SCRIPTS = 

PRECONFHDRS = \
	squid.h.in\
	squidconf.h.in

POSTCONFHDRS = \
	squid.h\
	squidconf.h

HDRS =	rk.h\
	sqfuncs.h\
	gki.h\
	gsi.h\
	msa.h\
	sre_random.h\
	ssi.h\
	stopwatch.h\
	vectorops.h

OBJS =	a2m.o\
	aligneval.o\
	alignio.o\
	clustal.o\
	cluster.o\
	dayhoff.o\
	eps.o\
	file.o\
	getopt.o\
	gki.o\
	gsi.o\
	hsregex.o\
	iupac.o\
	msa.o\
	msf.o\
	phylip.o\
	revcomp.o\
	rk.o\
	selex.o\
	seqencode.o\
	shuffle.o\
	sqerror.o\
	sqio.o\
	squidcore.o\
	sre_ctype.o\
	sre_math.o\
	sre_random.o\
	sre_string.o\
	ssi.o\
	stack.o\
	stockholm.o\
	stopwatch.o\
	translate.o\
	types.o\
	vectorops.o\
	weight.o

################################################################
# Targets that actually build the squid executables
all: $(PROGS) libsquid.a

$(PROGS): %: %_main.o ${OBJS} 
	${CC} ${CFLAGS} ${LDFLAGS} ${PVMLIBDIR} ${DEFS} -o $@ $@_main.o ${OBJS} ${PVMLIBS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${CPPFLAGS} ${PVMINCDIR} ${DEFS} -c $<		

################################################################


################################################################
# Targets expected by packages (e.g. HMMER) that
# include SQUID as a module. 
#
module: libsquid.a

libsquid.a: $(OBJS)
	$(AR) libsquid.a $(OBJS)
	$(RANLIB) libsquid.a
	chmod 644 libsquid.a
#################################################################


install:  $(PROGS) libsquid.a
	test -d $(LIBDIR)    || mkdir -p $(LIBDIR)
	test -d $(BINDIR)    || mkdir -p $(BINDIR)
	test -d $(SCRIPTDIR) || mkdir -p $(SCRIPTDIR)
	test -d $(INCLUDEDIR)|| mkdir -p $(INCLUDEDIR)
	test -d $(MANDIR)/man$(MANSUFFIX) || mkdir -p $(MANDIR)/man$(MANSUFFIX)
	cp libsquid.a $(LIBDIR)/
	cp $(HDRS)  $(INCLUDEDIR)/
	cp $(POSTCONFHDRS)  $(INCLUDEDIR)/
	cp $(PROGS) $(BINDIR)/
#	for scriptfile in $(SCRIPTS); do\
#	  cp Scripts/$$scriptfile $(SCRIPTDIR)/;\
#	done
	@for manpage in $(MANS); do\
	  $(INSTMAN) Man/$$manpage.man $(MANDIR)/man$(MANSUFFIX)/$$manpage.$(MANSUFFIX);\
	done

uninstall:
	rm $(LIBDIR)/libsquid.a
	for file in $(HDRS); do\
	   rm $(INCLUDEDIR)/$$file;\
	done
	for file in $(PROGS); do\
	   rm $(BINDIR)/$$file;\
	done
#	for file in $(SCRIPTS); do\
#	   rm $(SCRIPTDIR)/$$file;\
#	done
	for file in $(MANS); do\
	   rm $(MANDIR)/man$(MANSUFFIX)/$$file.$(MANSUFFIX);\
	done

check:	libsquid.a 
	(cd Testsuite; make CC="$(CC)" CFLAGS="$(CFLAGS)")
	(cd Testsuite; make check)

distclean:
	-rm -f *.o *~ core TAGS llib-lsquid.ln ccmalloc.log $(PROGS) 
	-rm -f Makefile libsquid.a ${POSTCONFHDRS}
	-rm -f config.log config.status
	-rm -rf autom4te.cache
	(cd Testsuite; make distclean)

clean:
	-rm -f *.o *~ core TAGS llib-lsquid.ln ccmalloc.log $(PROGS) 
	(cd Testsuite; make clean)

binclean:
	(cd Testsuite; make binclean)
	-rm -f *.o *~ core TAGS llib-lsquid.ln ccmalloc.log 
	-rm -f libsquid.a ${POSTCONFHDRS}
	-rm -f config.log config.status
	-rm -rf autom4te.cache

# dist: build a new distribution directory in squid-$RELEASE, and make a tarball.
#       Extracts straight from the CVS repository, so you must first do
#       a "cvs commit" (it checks to be sure you do, at least for the current
#       working directory). 
dist:
# Delete old versions of the same release
#
	@if test -d ${BASENAME}-$(RELEASE);        then rm -rf ${BASENAME}-$(RELEASE);        fi
	@if test -e ${BASENAME}-$(RELEASE).tar;    then rm -f  ${BASENAME}-$(RELEASE).tar;    fi
	@if test -e ${BASENAME}-$(RELEASE).tar.Z;  then rm -f  ${BASENAME}-$(RELEASE).tar.Z;  fi
	@if test -e ${BASENAME}-$(RELEASE).tar.gz; then rm -f  ${BASENAME}-$(RELEASE).tar.gz; fi
#
# CVS tag and extract. -c: make sure we committed;
#                      -F: allow more than one "make dist" per rel
# prep: must have done "cvs commit", and CVSROOT must be set
#
	cvs tag -c -F ${RELCODE}
	cvs export -r ${RELCODE} -d ${BASENAME}-${RELEASE} ${BASENAME}
# 
# Make the configure script from configure.ac
#
	(cd ${BASENAME}-${RELEASE}; autoconf)
#
# Include the appropriate license files
#
	cp Licenses/LICENSE.${LICENSETAG}     ${BASENAME}-${RELEASE}/LICENSE
	cp Licenses/COPYRIGHT.${LICENSETAG}   ${BASENAME}-${RELEASE}/COPYRIGHT
#
# Put license tags (short licenses) on files that need 'em (replace LICENSE keyword)
#
	for file in $(READMES) *.c ${HDRS} ${PRECONFHDRS}; do\
	   licenseadd.pl Licenses/$(LICENSETAG) ${BASENAME}-${RELEASE}/$$file;\
	done;
#
# Remove files/directories that aren't supposed to go out in the distro.
# Do this last, so other steps (license adding, etc.) have simple loops.
#
	-rm -rf ${BASENAME}-${RELEASE}/Licenses
	-rm -rf ${BASENAME}-${RELEASE}/Docs
	-rm ${BASENAME}-${RELEASE}/LOG
	-rm ${BASENAME}-${RELEASE}/configure.ac
	-rm ${BASENAME}-${RELEASE}/test_main.c
# 
# Do replacements
# (move this to configure)
#	for file in ${MANS}; do\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT} HMMER - Biological sequence analysis with profile HMMs ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT} Copyright (C) 1992-2003 Washington University School of Medicine ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT} All Rights Reserved ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT}  ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT}     This source code is distributed under the terms of the ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT}     GNU General Public License. See the files COPYING and LICENSE ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#	   sedition2 @RELEASEDATE@ ${RELEASEDATE} @PACKAGE@ ${PACKAGE} @RELEASE@ ${RELEASE} @COPYRIGHT@ ${COPYRIGHT}     for details. ${LICENSE} @SEEALSO@ ${SEEALSO} ${BASENAME}-${RELEASE}/Man/$$file.man;\
#
	done	
# 
# Set permissions.
#
	chmod -R ugo+rX ${BASENAME}-${RELEASE}
	chmod +x ${BASENAME}-${RELEASE}/install-sh
	chmod +x ${BASENAME}-${RELEASE}/Testsuite/bug-*
	chmod +x ${BASENAME}-${RELEASE}/Testsuite/x-base-*
#
# pack it up!
#
	tar cvf ${BASENAME}-${RELEASE}.tar ${BASENAME}-${RELEASE}
	${COMPRESS} ${BASENAME}-$(RELEASE).tar


ftpdist:
	cp -f ${BASENAME}-${RELEASE}.tar.gz ${FTPDIR}/
	rm -f ${FTPDIR}/${BASENAME}.tar.gz 
	(cd ${FTPDIR}; 	ln -s ${BASENAME}-${RELEASE}.tar.gz ${BASENAME}.tar.gz )


# implode:
#   Collapse a distro version of the squid library, removing all StL-specific
#   development material.
#
implode:
	-rm -rf 00CHECKLIST Bugs INSTALL LOG configure.ac
	-rm -rf Docs Formats Licenses Man 

TAGS:
	etags -t *.h *.c Makefile.in


