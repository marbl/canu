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
# gmake Makefile for AS_CNS
VERBOSE=1
REQD_LIBS = libAS_PER.a libAS_UTL.a libAS_MSG.a libAS_ALN.a libAS_CNS.a libAS_REZ.a libAS_CGW.a

LIB_SOURCES = MultiAlignment_CNS.c MultiAlignStore_CNS.c \
	Array_CNS.c Utilities_CNS.c colCorr_CNS.c
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)
LIBRARIES = libAS_CNS.a
LOCAL_WORK = $(shell cd ../..; pwd)


AS_CNS_SRCS = $(LIB_SOURCES) Consensus_CNS.c null_scaffolder.c quality_assess.c \
              placed_contigs.c TestBaseCall_CNS.c TestAlign_CNS.c \
		RebuildMultiAlign_CNS.c RebuildScaffolds_CNS.c ScaffoldUnitigProfile_CNS.c
# post_analysis.c AlignStats_CNS.c

AS_CNS_OBJS = $(AS_CNS_SRCS:.c=.o)

SOURCES = $(AS_CNS_SRCS)
OBJECTS = $(AS_CNS_OBJS)

PROGS = consensus null_scaffolder quality_assess placed_contigs testbasecall testalignment rebuildalignment rebuildscaffolds scaffoldunitigs
# post_analysis 

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(LIB_OBJECTS) $(LIBRARIES) $(OBJECTS) $(PROGS)

install: all
	echo "Nothing to be done for AS_CNS"

libAS_CNS.a: $(LIB_OBJECTS)

consensus: Consensus_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

testbasecall: TestBaseCall_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

testalignment: TestAlign_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

rebuildalignment: RebuildMultiAlign_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

rebuildscaffolds: RebuildScaffolds_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o colCorr_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

scaffoldunitigs: ScaffoldUnitigProfile_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o \
	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CGW.a libAS_LIN.a -lm

#alignstats: AlignStats_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
#	Array_CNS.o Utilities_CNS.o \
#	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libm.a

#post_analysis: post_analysis.o \
#        libm.a libAS_CNS.a libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libm.a

placed_contigs: placed_contigs.o \
        libAS_CNS.a libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libm.a

quality_assess: quality_assess.o \
        libAS_CNS.a libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libm.a

null_scaffolder: null_scaffolder.o libAS_MSG.a libAS_UTL.a
