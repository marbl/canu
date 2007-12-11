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

CGW_SOURCES =     AS_CGW_main.c \
                  eCR.c \
                  eCR-examineGap.c \
                  eCR-diagnostic.c \
                  resolveSurrogates.c \
                  dumpCloneMiddles.c \
                  exploreGapStructure.c \
                  frgs2clones.c \
                  createFrgDeletes.c \
                  dumpSingletons.c

#  Not all of these are external, most are probably private to cgw itself.
CGW_LIB_SOURCES = AS_CGW_dataTypes.c \
                  AS_CGW_EdgeDiagnostics.c \
                  CIEdgeT_CGW.c \
                  CIScaffoldT_Biconnected_CGW.c \
                  CIScaffoldT_CGW.c \
                  CIScaffoldT_Cleanup_CGW.c \
                  CIScaffoldT_Merge_CGW.c \
                  Celamy_CGW.c \
                  ChunkOverlap_CGW.c \
                  ContigT_CGW.c \
                  DemoteUnitigsWithRBP_CGW.c \
                  DiagnosticsCGW.c \
                  fragmentPlacement.c \
                  GraphCGW_T.c \
                  Input_CGW.c \
                  Instrument_CGW.c \
                  InterleavedMerging.c \
                  LeastSquaresGaps_CGW.c \
                  MergeEdges_CGW.c \
                  Output_CGW.c \
                  SEdgeT_CGW.c \
                  ScaffoldGraph_CGW.c \
                  SplitChunks_CGW.c \
                  SplitScaffolds_CGW.c \
                  Stats_CGW.c \
                  TransitiveReduction_CGW.c

SOURCES         = $(CGW_SOURCES) $(CGW_LIB_SOURCES)
SCRIPTS		= 

OBJECTS         = $(SOURCES:.c=.o) 
CGW_LIB_OBJECTS = $(CGW_LIB_SOURCES:.c=.o)

LIBRARIES     = libAS_CGW.a libCA.a
LIBS          = libCA.a

PROGS = cgw extendClearRanges resolveSurrogates dumpCloneMiddles dumpSingletons exploreGapStructure frgs2clones createFrgDeletes


include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop


libAS_CGW.a: $(CGW_LIB_OBJECTS)

libCA.a: $(CGW_LIB_OBJECTS)


cgw: AS_CGW_main.o $(LIBS)

extendClearRanges: eCR.o eCR-examineGap.o eCR-diagnostic.o $(LIBS)

resolveSurrogates: resolveSurrogates.o $(LIBS)

dumpCloneMiddles: dumpCloneMiddles.o $(LIBS)

dumpSingletons: dumpSingletons.o $(LIBS)

exploreGapStructure: exploreGapStructure.o $(LIBS)

frgs2clones: frgs2clones.o $(LIBS)

createFrgDeletes: createFrgDeletes.o $(LIBS)
