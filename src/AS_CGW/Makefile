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
                  fragmentPlacement.c \
                  dumpCloneMiddles.c \
                  dumpDistanceEstimates.c \
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
                  DiagnosticsCGW.c \
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

CRUD_SOURCES = \
                  AS_CGW_LoadCheckpoint.c \
                  CheckILK.c \
                  Checkpoint.c \
                  InstrumentCheckpoint.c \
                  ProcessScaffolds_CGW.c \
                  ScaffoldMap.c \
                  TraceNodes.c \
                  asmInfoForFrg.c \
                  assemblyStructure.c \
                  ckp_scflen.c \
                  combineMates.c \
                  createFrgDeletes.c \
                  dfsFragmentTiling.c \
                  dumpBinaryUnitigMAs.c \
                  dumpCheckpointAsFasta.c \
                  examineFragStore.c \
                  exploreGapStructure.c \
                  exploreMates.c \
                  exploreScaffoldEnds.c \
                  findAnEdge.c \
                  findMissedOverlaps.c \
                  fixECRCheckpoint.c \
                  fixZLFContigs.c \
                  frgs2clones.c \
                  greedyFragmentTiling.c \
                  localUnitigging_CGW.c \
                  smallLargeScaffolds.c \
                  testSplitChunks.c \
                  whereAreMyOverlaps.c

SOURCES         = $(CGW_SOURCES) $(CGW_LIB_SOURCES)
OBJECTS         = $(SOURCES:.c=.o) 
CGW_LIB_OBJECTS = $(CGW_LIB_SOURCES:.c=.o)

LIBRARIES     = libAS_CGW.a libCA.a
LIBS          = libCA.a

PROGS = cgw extendClearRanges resolveSurrogates dumpCloneMiddles dumpDistanceEstimates dumpSingletons

CRUD =  CheckILK Checkpoint process_scaffolds loadcgw \
        ScaffoldMap testSplitChunks instrumentCheckpoint \
	traceNodes findMissedOverlaps localUnitigging_CGW \
	dumpCheckpointAsFasta examineFragStore exploreMates fixECRCheckpoint \
	dumpBinaryUnitigMAs assemblyStructure \
	exploreGapStructure exploreScaffoldEnds asmInfoForFrg \
	whereAreMyOverlaps findAnEdge \
	combineMates ckp_scflen frgs2clones createFrgDeletes \
	greedyFragmentTiling dfsFragmentTiling


include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop


libAS_CGW.a: $(CGW_LIB_OBJECTS)

libCA.a: $(CGW_LIB_OBJECTS)


cgw: AS_CGW_main.o $(LIBS)

extendClearRanges: eCR.o eCR-examineGap.o eCR-diagnostic.o $(LIBS)

resolveSurrogates: resolveSurrogates.o fragmentPlacement.o $(LIBS)

dumpCloneMiddles: dumpCloneMiddles.o $(LIBS)

dumpDistanceEstimates: dumpDistanceEstimates.o $(LIBS) 

dumpSingletons: dumpSingletons.o $(LIBS)



splitutgs: UnitigSplitter.o $(LIBS)

loadcgw: AS_CGW_LoadCheckpoint.o $(LIBS)

dumpBinaryUnitigMAs: dumpBinaryUnitigMAs.o $(LIBS)

fixECRCheckpoint: fixECRCheckpoint.o $(LIBS)

exploreMates: exploreMates.o $(LIBS)

dumpCheckpointAsFasta: dumpCheckpointAsFasta.o $(LIBS)

examineFragStore: examineFragStore.o $(LIBS)

findMissedOverlaps: findMissedOverlaps.o $(LIBS)

localUnitigging_CGW: localUnitigging_CGW.o $(LIBS)

GetDregsStats: GetDregsStats.o $(LIBS)

findFragments: FindFragments.o $(LIBS)

CheckILK: CheckILK.o $(LIBS)

Checkpoint: Checkpoint.o $(LIBS)

process_scaffolds: ProcessScaffolds_CGW.o $(LIBS)

ScaffoldMap: ScaffoldMap.o $(LIBS)

testSplitChunks: testSplitChunks.o $(LIBS)

instrumentCheckpoint: InstrumentCheckpoint.o $(LIBS)

traceNodes: TraceNodes.o $(LIBS)

asmInfoForFrg: asmInfoForFrg.o $(LIBS)

assemblyStructure: assemblyStructure.o $(LIBS)

exploreGapStructure: exploreGapStructure.o $(LIBS)

greedyFragmentTiling: greedyFragmentTiling.o $(LIBS)

dfsFragmentTiling: dfsFragmentTiling.o $(LIBS)

exploreScaffoldEnds: exploreScaffoldEnds.o $(LIBS)

ckp_scflen: ckp_scflen.o $(LIBS)

findAnEdge: findAnEdge.o $(LIBS)

whereAreMyOverlaps: whereAreMyOverlaps.o $(LIBS)

combineMates: combineMates.o $(LIBS)

createFrgDeletes: createFrgDeletes.o $(LIBS)

frgs2clones: frgs2clones.o $(LIBS)
