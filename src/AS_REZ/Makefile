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
#VERBOSE=1

#     SOURCES = 
#       Lists all sources compiled by this makefile.
#
#     OBJECTS = 
#       Lists all objects compiled by this makefile (i.e., one-to-one
#       with SOURCES). OBJECTS can be listed explicitly or can be
#       specified with the following definition:
#
#          OBJECTS = $(SOURCES:.c=.o)
#
#     LIBRARIES = 
#       Lists simple name (in "libx.a" form) of each library that is
#       built by this makefile.
#
#     PROGS = 
#       Lists each executable created by this makefile.
#
#     LOCAL_WORK =
#       Identifies the <subsystem_level> directory (parent of
#       the inc, src, lib, obj, and bin subdirectories) with
#       an absolute path.
#       For a makefile located in
#         <subsystem_level>/src/<component_level_1>
#       the following definition would work:
#          LOCAL_WORK = $(shell cd ../..; pwd)

#
# Also include dependencies of programs and libraries on objects and
# libraries, e.g.:
#
#    my_program: $(OBJECTS) libxxx.a
#

# Also include as first target, an "all" target, e.g.:
#
#    all: $(PROGS)
#

LOCAL_WORK = $(shell cd ../../; pwd)

LIN_ALG_LIB = libAS_LIN.a
#ifeq ($(OSTYPE), osf3.2)
#LIN_ALG_LIB = libdxml.so
#endif

ARTSOURCES = batchify.c chunk_analyze.c OverlapChunks.c

#MHETSOURCES =  MicroHetREZ.c MicroHetSimulatorREZ.c CheckDistributionREZ. PartitionsREZ.c ScoreREZ.c MicroHetPartitionsREZ.c MicroHetScoreREZ.c MicroHetIUMREZ.c MicroHetInterfaceREZ.c
MHETSOURCES =   MicroHetSimulatorREZ.c CheckDistributionREZ. PartitionsREZ.c ScoreREZ.c \
                MicroHetPartitionsREZ.c MicroHetScoreREZ.c MicroHetIUMREZ.c  MicroHetPairSimREZ.c \
                MicroHetIUMREZ_test3.c MicroHetIUMdiffsPrint.c
# MicroHetInterfaceREZ.c MicroHetREZ.c
#OVERLAPSOURCES = GreedyOverlapMainREZ.c LengthAdaptionREZ.c
OVERLAPSOURCES = 
LIBSOURCES =  ConsistencyChecksREZ.c UtilsREZ.c GapFillREZ.c UpdateREZ.c GapWalkerREZ.c StatisticsREZ.c \
        BccREZ.c GWDriversREZ.c SubgraphREZ.c ConfirmStones.c GreedyOverlapREZ.c \
	MicroHetREZ_test3.c MicroHetScoreREZ_test3.c MicroHetInterfaceREZ_test3.c MicroHetPartitionsREZ_test3.c FbacREZ.c

#ContainmentREZ.c

TMPSOURCES = 

SOURCES =  $(ARTSOURCES) $(LIBSOURCES) $(MHETSOURCES) $(OVERLAPSOURCES) $(TMPSOURCES)

:
SCRIPTS=get-repchunk print-ali extract-uoms.pl extract-lsf-uoms.pl concat-uoms.pl

OBJECTS = $(SOURCES:.c=.o) 
LIBOBJECTS = $(LIBSOURCES:.c=.o) 

LIBRARIES = libAS_REZ.a libCA.a

PROGS = batchify chunk_analyze microhet-test3 MicroHetIUMdiffsPrint

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

#-include make.conf

install: $(SCRIPTS)
	@chmod ugo+x $(SCRIPTS)
	@cp -f $(SCRIPTS) $(LOCAL_BIN)


all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

libAS_REZ.a: $(LIBOBJECTS) 

libCA.a: $(LIBOBJECTS)

batchify:  batchify.o libAS_MSG.a libAS_UTL.a

chunk_analyze:  chunk_analyze.o libAS_MSG.a libAS_UTL.a

greedy-overlap : GreedyOverlapMainREZ.o libCA.a

length-adaption : LengthAdaptionREZ.o libCA.a

getrep : GetRepeats.o libAS_MSG.a

overlap-chunks: OverlapChunks.o libCA.a

CGW_OBJECTS =  AS_CGW_dataTypes.o GraphCGW_T.o AS_CGW_histo.o DiagnosticsCGW.o Celamy_CGW.o  ScaffoldGraph_CGW.o CIEdgeT_CGW.o CIScaffoldT_CGW.o SEdgeT_CGW.o\
MergeEdges_CGW.o ContigT_CGW.o LeastSquaresGaps_CGW.o\
Output_CGW.o Input_CGW.o TransitiveReduction_CGW.o\
ChunkOverlap_CGW.o SplitScaffolds_CGW.o

extract-quality-uoms : ExtractQualityREZ.o libCA.a

microhet : MicroHetIUMREZ.o MicroHetREZ.o MicroHetPartitionsREZ.o MicroHetScoreREZ.o MicroHetInterfaceREZ.o  libCA.a

microhet-simulator : MicroHetInterfaceREZ.o MicroHetSimulatorREZ.o MicroHetREZ.o MicroHetPartitionsREZ.o MicroHetScoreREZ.o  libCA.a

microhet-pairsim : MicroHetInterfaceREZ.o MicroHetPairSimREZ.o MicroHetREZ.o MicroHetPartitionsREZ.o MicroHetScoreREZ.o  libCA.a

check-distributions: CheckDistributionREZ.o PartitionsREZ.o ScoreREZ.o UtilsREZ.o libAS_MSG.a libAS_SDB.a libAS_UTL.a

microhet-clean :
	rm  MicroHet*.o ${LOCAL_OBJ}/MicroHet*.o


microhet-test3 :  MicroHetIUMREZ_test3.o libCA.a

MicroHetIUMdiffsPrint : MicroHetIUMdiffsPrint.o libCA.a
