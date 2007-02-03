#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#

LOCAL_WORK = $(shell cd ../../; pwd)

LIN_ALG_LIB = libAS_LIN.a

#ifeq ($(OSTYPE), osf3.2)
#LIN_ALG_LIB = libdxml.so
#endif

LIBSOURCES = ConsistencyChecksREZ.c \
             UtilsREZ.c \
             GapFillREZ.c \
             UpdateREZ.c \
             GapWalkerREZ.c \
             StatisticsREZ.c \
             BccREZ.c \
             GWDriversREZ.c \
             SubgraphREZ.c \
             ConfirmStones.c \
             GreedyOverlapREZ.c \
             MicroHetREZ_test3.c \
             MicroHetScoreREZ_test3.c \
             MicroHetInterfaceREZ_test3.c \
             MicroHetPartitionsREZ_test3.c \
             FbacREZ.c

SOURCES = MicroHetIUMREZ_test3.c MicroHetIUMdiffsPrint.c batchify.c chunk_analyze.c OverlapChunks.c $(LIBSOURCES)

SCRIPTS = get-repchunk print-ali extract-uoms.pl extract-lsf-uoms.pl concat-uoms.pl

OBJECTS    = $(SOURCES:.c=.o) 
LIBOBJECTS = $(LIBSOURCES:.c=.o) 

LIBRARIES  = libAS_REZ.a libCA.a

PROGS = batchify chunk_analyze microhet-test3 MicroHetIUMdiffsPrint

include $(LOCAL_WORK)/src/c_make.as

install: $(SCRIPTS)
	@chmod ugo+x $(SCRIPTS)
	@cp -f $(SCRIPTS) $(LOCAL_BIN)

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

libAS_REZ.a: $(LIBOBJECTS) 

libCA.a: $(LIBOBJECTS)


batchify:  batchify.o libAS_MSG.a libAS_UTL.a

chunk_analyze:  chunk_analyze.o libAS_MSG.a libAS_UTL.a

microhet-test3 : MicroHetIUMREZ_test3.o libCA.a

MicroHetIUMdiffsPrint : MicroHetIUMdiffsPrint.o libCA.a
