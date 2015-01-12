################################################################
# Makefile for the Easel testsuite
# Makefile is generated from Makefile.in by autoconf
# DO NOT EDIT Makefile: only edit Makefile.in
# 
# SRE, Mon Jun 11 08:19:32 2007 [Janelia]
# SVN $Id: Makefile.in 509 2010-02-07 22:56:55Z eddys $

top_srcdir = ../..
srcdir     = .
VPATH      =  ../..

CC       = gcc -std=gnu99
CFLAGS   = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing
CPPFLAGS = 
LDFLAGS  = 
SIMDFLAGS= -msse2
DEFS     = -DHAVE_CONFIG_H
LIBS     =   -lm
MPILIBS  = 

AR      = /usr/bin/ar rcv
RANLIB  = ranlib

ESLDIR   = ..
SQCLEVEL = 2

PROGS = parse-t1\
	dmatrix-t1\
	dmatrix-t2\

all: 
	cd ${ESLDIR}; ${MAKE} easel_lib
	cd ${ESLDIR}; ${MAKE} easel_progs
	cd ${ESLDIR}; ${MAKE} easel_utests

check:  
	cd ${ESLDIR}; ${MAKE} easel_lib
	cd ${ESLDIR}; ${MAKE} easel_progs
	cd ${ESLDIR}; ${MAKE} easel_utests
	${MAKE} run_sqc

dev:
	cd ${ESLDIR}; ${MAKE} easel_lib
	cd ${ESLDIR}; ${MAKE} easel_progs
	cd ${ESLDIR}; ${MAKE} easel_utests

testsuite_progs: 

run_sqc:
	@echo 
	@echo Running Easel test suite...
	@echo
	@${srcdir}/../devkit/sqc 2 ${srcdir}/testsuite.sqc .. ${srcdir}/..

${PROGS}: %: %.o ${ESLDIR}/libeasel.a
	${CC} ${CFLAGS} ${SIMDFLAGS} -o $@ -L ${ESLDIR} $@.o -leasel -lm

.c.o:
	${CC} -I ${ESLDIR} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -o $@ -c $<		

dcheck: 
	@echo    Running developer testsuite...
	@echo -n Running unit tests...
	@${srcdir}/../devkit/sqc -D ${srcdir} 2 ${srcdir}/testsuite.sqc . ${ESLDIR} ${ESLDIR}/miniapps/ ${srcdir}  >  utest.report
	@echo    done. [See utest.report]
	@echo -n Testing that all drivers compile...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./driver_report.pl       >  driver.report
	@echo    done. [See driver.report]
	@echo -n Measuring test code coverage...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./coverage_report.pl     >  coverage.report
	@echo    done. [See coverage.report]
	@echo -n Looking for memory leaks...
	@env CC="${CC}" CFLAGS="${CFLAGS}" ./valgrind_report.pl     >  valgrind.report
	@echo    done. [See valgrind.report]
	@cat driver.report coverage.report valgrind.report utest.report > dcheck.report
	@echo
	@echo    Finished developer testsuite.
	@echo    For a complete report, see dcheck.report.

clean:
	-rm -f esl[^_]?????
	-rm -f *.gcno *.gcda *.gcov
	-rm -f *.o *~ 
	-rm -f driver.report coverage.report valgrind.report utest.report
	-rm -f ${PROGS} 

distclean:
	${MAKE} clean
	rm Makefile
