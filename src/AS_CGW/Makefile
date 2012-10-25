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

CGW_SOURCES =     AS_CGW_main.C \
                  cgwDump.C \
                  analyzeScaffolds.C \
                  eCR.C \
                  eCR-examineGap.C \
                  eCR-diagnostic.C \
                  eCR-partition.C \
                  resolveSurrogates.C \
                  dumpCloneMiddles.C \
                  frgs2clones.C \
                  dumpSingletons.C

#  Not all of these are external, most are probably private to cgw itself.
CGW_LIB_SOURCES = Globals_CGW.C \
                  CIEdgeT_CGW.C \
                  CIScaffoldT_Analysis.C \
                  CIScaffoldT_Biconnected_CGW.C \
                  CIScaffoldT_CGW.C \
                  CIScaffoldT_Cleanup_CGW.C \
                  CIScaffoldT_Merge_CGW.C \
                  CIScaffoldT_Merge_Interleaved.C \
                  CIScaffoldT_Merge_AlignScaffold.C \
                  CIScaffoldT_MergeScaffolds.C \
                  Celamy_CGW.C \
                  ChunkOverlap_CGW.C \
                  ContigT_CGW.C \
                  DemoteUnitigsWithRBP_CGW.C \
                  fragmentPlacement.C \
                  GraphCGW_T.C \
                  GraphEdgeIterator.C \
                  Input_CGW.C \
                  Instrument_CGW.C \
                  InterleavedMerging.C \
                  LeastSquaresGaps_CGW.C \
                  MarkInternalEdgeStatus.C \
                  MergeEdges_CGW.C \
                  Output_CGW.C \
                  SEdgeT_CGW.C \
                  ScaffoldGraph_CGW.C \
                  SplitScaffolds_CGW.C \
                  Stats_CGW.C \
                  TransitiveReduction_CGW.C \
		  ShatterScaffolds_CGW.C

SOURCES         = $(CGW_SOURCES) $(CGW_LIB_SOURCES)
SCRIPTS		=

OBJECTS         = $(SOURCES:.C=.o)
CGW_LIB_OBJECTS = $(CGW_LIB_SOURCES:.C=.o)

LIBRARIES     = libAS_CGW.a libCA.a
LIBS          = libCA.a

CXX_PROGS = cgw cgwDump analyzeScaffolds extendClearRanges extendClearRangesPartition resolveSurrogates dumpCloneMiddles dumpSingletons frgs2clones

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop


libAS_CGW.a: $(CGW_LIB_OBJECTS)

libCA.a: $(CGW_LIB_OBJECTS)

cgw: AS_CGW_main.o $(LIBS)

analyzeScaffolds: analyzeScaffolds.o $(LIBS)

extendClearRanges: eCR.o eCR-examineGap.o eCR-diagnostic.o $(LIBS)

extendClearRangesPartition: eCR-partition.o $(LIBS)

resolveSurrogates: resolveSurrogates.o $(LIBS)

dumpCloneMiddles: dumpCloneMiddles.o $(LIBS)

dumpSingletons: dumpSingletons.o $(LIBS)

frgs2clones: frgs2clones.o $(LIBS)

cgwDump: cgwDump.o $(LIBS)
