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
#######################################################################

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

LOCAL_WORK = $(shell cd ../../; pwd)

CGW_SOURCES = AS_CGW_LoadCheckpoint.c AS_CGW_main.c AS_CGW_dataTypes.c \
	GraphCGW_T.c Stats_CGW.c\
        DiagnosticsCGW.c Celamy_CGW.c \
	ScaffoldGraph_CGW.c CIEdgeT_CGW.c CIScaffoldT_CGW.c \
	CIScaffoldT_Merge_CGW.c CIScaffoldT_Cleanup_CGW.c SEdgeT_CGW.c \
	MergeEdges_CGW.c ContigT_CGW.c LeastSquaresGaps_CGW.c \
	Output_CGW.c Input_CGW.c TransitiveReduction_CGW.c \
	ChunkOverlap_CGW.c SplitScaffolds_CGW.c \
	CIScaffoldT_Biconnected_CGW.c SplitChunks_CGW.c Instrument_CGW.c \
	eCR.c eCR-examineGap.c eCR-diagnostic.c localUnitigging_CGW.c \
	AS_CGW_EdgeDiagnostics.c findMissedOverlaps.c dumpCheckpointAsFasta.c \
	examineFragStore.c InterleavedMerging.c exploreMates.c \
	fixZLFContigs.c fixECRCheckpoint.c dumpBinaryUnitigMAs.c \
	assemblyStructure.c exploreGapStructure.c \
	asmInfoForFrg.c exploreScaffoldEnds.c dumpCloneMiddles.c \
	whereAreMyOverlaps.c findAnEdge.c dumpSingletons.c combineMates.c \
	resolveSurrogates.c fragmentPlacement.c ckp_scflen.c frgs2clones.c createFrgDeletes.c \
        greedyFragmentTiling.c dfsFragmentTiling.c dumpDistanceEstimates.c

CGW_LIB_SOURCES = AS_CGW_dataTypes.c GraphCGW_T.c Stats_CGW.c\
        DiagnosticsCGW.c Celamy_CGW.c \
	ScaffoldGraph_CGW.c CIEdgeT_CGW.c CIScaffoldT_CGW.c \
	CIScaffoldT_Merge_CGW.c CIScaffoldT_Cleanup_CGW.c SEdgeT_CGW.c \
	MergeEdges_CGW.c ContigT_CGW.c LeastSquaresGaps_CGW.c \
	Output_CGW.c Input_CGW.c TransitiveReduction_CGW.c \
	ChunkOverlap_CGW.c SplitScaffolds_CGW.c CIScaffoldT_Biconnected_CGW.c \
	SplitChunks_CGW.c Instrument_CGW.c AS_CGW_EdgeDiagnostics.c \
	InterleavedMerging.c fixZLFContigs.c fragmentPlacement.c

TSC_SOURCES = testSplitChunks.c
TSC_OBJECTS = $(TSC_SOURCES:.c=.o)

INST_SOURCES = InstrumentCheckpoint.c
INST_OBJECTS = $(INST_SOURCES:.c=.o)

TRAC_SOURCES = TraceNodes.c
TRAC_OBJECTS = $(TRAC_SOURCES:.c=.o)

EXPL_SOURCES = exploreMates.c
EXPL_OBJECTS = $(EXPL_SOURCES:.c=.o)

#SLS_SOURCES = smallLargeScaffolds.c
#SLS_OBJECTS = $(SLS_SOURCES:.c=.o)

SOURCES = $(CGW_SOURCES) CheckILK.c Checkpoint.c ProcessScaffolds_CGW.c \
          ScaffoldMap.c $(RPT_SOURCES) $(TSC_SOURCES) $(INST_SOURCES) \
	  $(TRAC_SOURCES)
CGW_OBJECTS = $(CGW_SOURCES:.c=.o) 
OBJECTS = $(CGW_OBJECTS) 
CGW_LIB_OBJECTS = $(CGW_LIB_SOURCES:.c=.o)

LIBRARIES = libAS_CGW.a

PROGS = cgw extendClearRanges CheckILK Checkpoint process_scaffolds loadcgw \
        ScaffoldMap testSplitChunks instrumentCheckpoint \
	traceNodes findMissedOverlaps localUnitigging_CGW \
	dumpCheckpointAsFasta examineFragStore exploreMates fixECRCheckpoint \
	dumpBinaryUnitigMAs assemblyStructure \
	exploreGapStructure exploreScaffoldEnds asmInfoForFrg \
	dumpCloneMiddles whereAreMyOverlaps findAnEdge \
	dumpSingletons combineMates resolveSurrogates ckp_scflen frgs2clones createFrgDeletes \
	greedyFragmentTiling dfsFragmentTiling dumpDistanceEstimates

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)

ALL_LIBRARIES = libCA.a

cgw: AS_CGW_main.o $(ALL_LIBRARIES)

splitutgs: UnitigSplitter.o $(ALL_LIBRARIES)

loadcgw: AS_CGW_LoadCheckpoint.o $(ALL_LIBRARIES)

dumpBinaryUnitigMAs: dumpBinaryUnitigMAs.o $(ALL_LIBRARIES)

fixECRCheckpoint: fixECRCheckpoint.o $(ALL_LIBRARIES)

exploreMates: exploreMates.o $(ALL_LIBRARIES)

extendClearRanges: eCR.o eCR-examineGap.o eCR-diagnostic.o $(ALL_LIBRARIES)

dumpCheckpointAsFasta: dumpCheckpointAsFasta.o $(ALL_LIBRARIES)

examineFragStore: examineFragStore.o $(ALL_LIBRARIES)

findMissedOverlaps: findMissedOverlaps.o $(ALL_LIBRARIES)

localUnitigging_CGW: localUnitigging_CGW.o $(ALL_LIBRARIES)

GetDregsStats: GetDregsStats.o $(ALL_LIBRARIES)

findFragments: FindFragments.o $(ALL_LIBRARIES)

CheckILK: CheckILK.o $(ALL_LIBRARIES)

Checkpoint: Checkpoint.o $(ALL_LIBRARIES)

process_scaffolds: ProcessScaffolds_CGW.o $(ALL_LIBRARIES)

ScaffoldMap: ScaffoldMap.o $(ALL_LIBRARIES)

testSplitChunks: testSplitChunks.o $(ALL_LIBRARIES)

instrumentCheckpoint: InstrumentCheckpoint.o $(ALL_LIBRARIES)

traceNodes: TraceNodes.o $(ALL_LIBRARIES)

asmInfoForFrg: asmInfoForFrg.o $(ALL_LIBRARIES)

dumpSingletons: dumpSingletons.o $(ALL_LIBRARIES)

assemblyStructure: assemblyStructure.o $(ALL_LIBRARIES)

exploreGapStructure: exploreGapStructure.o $(ALL_LIBRARIES)

greedyFragmentTiling: greedyFragmentTiling.o $(ALL_LIBRARIES)

dfsFragmentTiling: dfsFragmentTiling.o $(ALL_LIBRARIES)

exploreScaffoldEnds: exploreScaffoldEnds.o $(ALL_LIBRARIES)

ckp_scflen: ckp_scflen.o $(ALL_LIBRARIES)

dumpCloneMiddles: dumpCloneMiddles.o $(ALL_LIBRARIES)

findAnEdge: findAnEdge.o $(ALL_LIBRARIES)

whereAreMyOverlaps: whereAreMyOverlaps.o $(ALL_LIBRARIES)

combineMates: combineMates.o $(ALL_LIBRARIES)

createFrgDeletes: createFrgDeletes.o $(ALL_LIBRARIES)

frgs2clones: frgs2clones.o $(ALL_LIBRARIES)

dumpDistanceEstimates: dumpDistanceEstimates.o $(ALL_LIBRARIES) 

resolveSurrogates: resolveSurrogates.o $(ALL_LIBRARIES)

libAS_CGW.a: $(CGW_LIB_OBJECTS)
