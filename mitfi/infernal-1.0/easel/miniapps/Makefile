# Makefile.in for Easel's miniapps
# SVN $Id: Makefile.in 267 2008-06-18 12:19:44Z eddys $
SHELL = /bin/sh

## Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin


## your compiler and compiler flags.
#
CC       = gcc -std=gnu99
CFLAGS   = -g -O2 -msse2
DEFS     = -DHAVE_CONFIG_H
LDFLAGS  = 
LIBS     = -leasel   -lm

## where is Easel?
#
ESLLIBDIR = -L..
ESLINCDIR = -I..


## list of the miniapps to compile.
#
PROGS = esl-afetch\
	esl-alimanip\
	esl-alistat\
	esl-compstruct\
	esl-histplot\
	esl-reformat\
	esl-selectn\
	esl-seqstat\
	esl-sfetch\
	esl-shuffle\
	esl-weight

## Targets defining how to make miniapp executables.
#
all: 	$(PROGS) 

${PROGS}: %: %.o  ../libeasel.a
	${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} ${ESLLIBDIR} -o $@ $@.o ${OBJS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${DEFS} ${ESLINCDIR} -c $<

install:
	cp ${PROGS} ${DESTDIR}${bindir}/

clean:
	-rm ${PROGS}
	-rm *.o *~
	-rm TAGS

distclean:
	make clean
	-rm -f ${PROGS}
	-rm Makefile

TAGS:
	etags *.[ch] Makefile.in


################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
################################################################
