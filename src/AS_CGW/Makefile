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


LIN_ALG_LIB = libAS_LIN.a
#ifeq ($(OSTYPE), osf3.2)
#LIN_ALG_LIB = libdxml.so
#endif

LOCAL_WORK = $(shell cd ../../; pwd)

OLD_CGW_SOURCES = AS_CGW_MateEdges.c AS_CGW_Scaffold.c \
	          AS_CGW_Scaffolds.c AS_CGW_ScaffoldsPrimitives.c\
	         ChunkOfScaffoldsCGW.c COSMateEdges_CGW.c \
	         ChunkGraph_CGW.c AS_CGW_PartialOrder.c ScaffoldOverlapsCGW.c

CGW_SOURCES = AS_CGW_LoadCheckpoint.c AS_CGW_main.c AS_CGW_dataTypes.c \
	GraphCGW_T.c Stats_CGW.c\
        DiagnosticsCGW.c Celamy_CGW.c \
	ScaffoldGraph_CGW.c CIEdgeT_CGW.c CIScaffoldT_CGW.c \
	CIScaffoldT_Merge_CGW.c CIScaffoldT_Cleanup_CGW.c SEdgeT_CGW.c \
	MergeEdges_CGW.c ContigT_CGW.c LeastSquaresGaps_CGW.c \
	Output_CGW.c Input_CGW.c TransitiveReduction_CGW.c \
	ChunkOverlap_CGW.c SplitScaffolds_CGW.c \
	CIScaffoldT_Biconnected_CGW.c SplitChunks_CGW.c Instrument_CGW.c \
	extendClearRanges.c localUnitigging_CGW.c \
	AS_CGW_EdgeDiagnostics.c findMissedOverlaps.c scaffold2fasta.c \
	examineFragStore.c InterleavedMerging.c exploreMates.c \
	fixZLFContigs.c fixECRCheckpoint.c dumpBinaryUnitigMAs.c
#	smallLargeScaffolds.c FixContigMultiAligns.c fixZLFContigs.c
#	ChunkOverlap_CGW.c SplitScaffolds_CGW.c CIScaffoldT_Biconnected_CGW.c UnitigSplitter.c FindFragments.c

CGW_LIB_SOURCES = AS_CGW_dataTypes.c GraphCGW_T.c Stats_CGW.c\
        DiagnosticsCGW.c Celamy_CGW.c \
	ScaffoldGraph_CGW.c CIEdgeT_CGW.c CIScaffoldT_CGW.c \
	CIScaffoldT_Merge_CGW.c CIScaffoldT_Cleanup_CGW.c SEdgeT_CGW.c \
	MergeEdges_CGW.c ContigT_CGW.c LeastSquaresGaps_CGW.c \
	Output_CGW.c Input_CGW.c TransitiveReduction_CGW.c \
	ChunkOverlap_CGW.c SplitScaffolds_CGW.c CIScaffoldT_Biconnected_CGW.c \
	SplitChunks_CGW.c Instrument_CGW.c AS_CGW_EdgeDiagnostics.c \
	InterleavedMerging.c fixZLFContigs.c

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

#CHUNKSTORE_SOURCES = CreateChunkStore.c LoadChunkStore.c testChunkOverlap.c
#CHUNKSTORE_OBJECTS = $(CHUNKSTORE_SOURCES:.c=.o) 
SOURCES = $(CGW_SOURCES) CheckILK.c Checkpoint.c ProcessScaffolds_CGW.c \
          ScaffoldMap.c $(RPT_SOURCES) $(TSC_SOURCES) $(INST_SOURCES) \
	  $(TRAC_SOURCES)
CGW_OBJECTS = $(CGW_SOURCES:.c=.o) 
OBJECTS = $(CGW_OBJECTS) 
CGW_LIB_OBJECTS = $(CGW_LIB_SOURCES:.c=.o)


LIBRARIES = libAS_CGW.a

#PROGS = cgw CreateChunkStore LoadChunkStore testChunkOverlap
PROGS = cgw CheckILK Checkpoint process_scaffolds loadcgw \
        ScaffoldMap testSplitChunks instrumentCheckpoint \
	traceNodes extendClearRanges findMissedOverlaps localUnitigging_CGW \
	scaffold2fasta examineFragStore exploreMates fixECRCheckpoint \
	dumpBinaryUnitigMAs
#	smallLargeScaffolds fixContigMultiAligns
#        ScaffoldMap splitutgs findFragments 

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

#-include make.conf


all: $(OBJECTS) $(LIBRARIES) $(PROGS)

cgw: AS_CGW_main.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CGW.a libAS_CNS.a libAS_UTL.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

splitutgs: UnitigSplitter.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CNS.a libAS_UTL.a libAS_SDB.a libAS_CGW.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

loadcgw: AS_CGW_LoadCheckpoint.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

dumpBinaryUnitigMAs: dumpBinaryUnitigMAs.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

fixECRCheckpoint: fixECRCheckpoint.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

exploreMates: exploreMates.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

#smallLargeScaffolds: smallLargeScaffolds.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

extendClearRanges: extendClearRanges.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

#fixContigMultiAligns: FixContigMultiAligns.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

scaffold2fasta: scaffold2fasta.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

examineFragStore: examineFragStore.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_SDB.a  libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

findMissedOverlaps: findMissedOverlaps.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_SDB.a  libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

localUnitigging_CGW: localUnitigging_CGW.o libAS_CGB.a libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CNS.a libAS_OVL.a libAS_CGW.a libAS_ALN.a libAS_SDB.a  libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_UTL.a libAS_CGW.a libAS_UTL.a libm.a

GetDregsStats: GetDregsStats.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_SDB.a  libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

findFragments: FindFragments.o libAS_CGW.a libAS_MSG.a  libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_SDB.a  libAS_CGW.a $(LIN_ALG_LIB) libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a

CheckILK: CheckILK.o libAS_CNS.a libAS_UTL.a libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a 

Checkpoint: Checkpoint.o                  libAS_CGW.a libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a libAS_ALN.a libAS_MSG.a libAS_SDB.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

process_scaffolds: ProcessScaffolds_CGW.o libAS_MSG.a libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_SDB.a libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_CGW.a $(LIN_ALG_LIB) libAS_UTL.a libm.a

ScaffoldMap: ScaffoldMap.o libAS_MSG.a  libAS_UTL.a libAS_CGW.a libAS_UTL.a libm.a

testSplitChunks: testSplitChunks.o libAS_CGW.a libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CGW.a libAS_ALN.a libAS_CNS.a libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

instrumentCheckpoint: InstrumentCheckpoint.o  libAS_CGW.a libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CGW.a libAS_ALN.a libAS_CNS.a libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

traceNodes: TraceNodes.o  libAS_CGW.a libAS_CNS.a libAS_PER.a libAS_ALN.a libAS_REZ.a libAS_CGW.a libAS_ALN.a libAS_CNS.a libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_CGW.a $(LIN_ALG_LIB) libAS_CGW.a libAS_UTL.a libm.a

libAS_CGW.a: $(CGW_LIB_OBJECTS)

#LoadChunkStore: LoadChunkStore.o  libAS_CGW.a libAS_PER.a libAS_UTL.a libAS_MSG.a libAS_CGW.a libAS_UTL.a libm.a 

#CreateChunkStore: CreateChunkStore.o  libAS_CGW.a libAS_UTL.a libAS_MSG.a libAS_PER.a  libAS_CGW.a libAS_UTL.a libm.a 

#testChunkOverlap: testChunkOverlap.o ChunkStore_AS.o libAS_UTL.a libAS_MSG.a  \
#	libAS_ALN.a libAS_PER.a  libAS_CGW.a libAS_UTL.a libm.a
