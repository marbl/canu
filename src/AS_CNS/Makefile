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
#VERBOSE=1
REQD_LIBS = libAS_PER.a libAS_UTL.a libAS_MSG.a libAS_ALN.a libAS_CNS.a libAS_REZ.a libAS_CGW.a

LIB_SOURCES = MultiAlignment_CNS.c MultiAlignStore_CNS.c \
	Array_CNS.c Utilities_CNS.c
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)
LIBRARIES = libAS_CNS.a
LOCAL_WORK = $(shell cd ../..; pwd)


AS_CNS_SRCS = $(LIB_SOURCES) Consensus_CNS.c null_scaffolder.c quality_assess.c \
              placed_contigs.c TestBaseCall_CNS.c TestAlign_CNS.c \
		RebuildMultiAlign_CNS.c RebuildScaffolds_CNS.c ScaffoldUnitigProfile_CNS.c \
	      colCorr_CNS.c
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

LIBS = -lAS_CGW -lAS_REZ -lAS_CGW -lAS_PER -lAS_MSG \
       -lAS_SDB -lAS_ALN -lAS_OVL -lAS_UTL -lAS_LIN

consensus: Consensus_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o $(LIBS)

testbasecall: TestBaseCall_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o $(LIBS)

testalignment: TestAlign_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o $(LIBS)

rebuildalignment: RebuildMultiAlign_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o $(LIBS)

rebuildscaffolds: RebuildScaffolds_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o colCorr_CNS.o $(LIBS)

scaffoldunitigs: ScaffoldUnitigProfile_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
	Array_CNS.o Utilities_CNS.o $(LIBS)

#alignstats: AlignStats_CNS.o  MultiAlignment_CNS.o MultiAlignStore_CNS.o \
#	Array_CNS.o Utilities_CNS.o \
#	libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a 

#post_analysis: post_analysis.o \
#         libAS_CNS.a libAS_MSG.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libm.a

placed_contigs: placed_contigs.o $(LIBS)

quality_assess: quality_assess.o $(LIBS)

null_scaffolder: null_scaffolder.o libAS_MSG.a libAS_UTL.a
