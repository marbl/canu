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
# gmake Makefile for AS_MPA

LOCAL_WORK = $(shell cd ../..; pwd)

C_PROG_SOURCES = getUniqueMatchesFromCGM.c 
C_PROG_OBJECTS = $(C_PROG_SOURCES:.c=.o)

C_SOURCES = $(C_PROG_SOURCES)
C_OBJECTS = $(C_SOURCES:.c=.o)

CXX_SOURCES = mapContigsScaffolds.cc \
	getIntervalIntersections.cc \
	pullFragmentLocations.cc \
	matchCounts.cc \
	processIntra.cc \
	processInter.cc \
	processInterCG.cc \
	testConvexHull.cc \
	findCliques.cc \
	findMPSets.cc \
	findTTIs.cc \
	quads2mbrs.cc \
	rotatePoints.cc \
	clusterMPs.cc
#	compareAssemblyMPs.cc
CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)

SOURCES = $(C_SOURCES) $(CXX_SOURCES)
OBJECTS = $(C_OBJECTS) $(CXX_OBJECTS)

PROGS = getUniqueMatchesFromCGM \
	mapContigsScaffolds \
	getIntervalIntersections \
	pullFragmentLocations \
	matchCounts \
	processIntra \
	processInter \
	processInterCG \
	testConvexHull \
	findCliques \
	findMPSets \
	findTTIs \
	quads2mbrs \
	rotatePoints
#	compareAssemblyMPs

SCRIPTS = atai2Lengths.sh \
	breakpoints2StartLenNs.sh \
	getDoubleCountedInTranpositions.sh \
	analyzeGaps.sh \
	getInterFromResults.sh \
	getInterFromResultsATA.sh \
	getInterToResults.sh \
	getIntraResults.sh \
	morphInter.sh \
	processMPs.sh \
	processAssemblyMPs.sh \
	buildMapStoreFromMatches.sh \
	getAllAssemblyMPs.sh \
	findBacDeletions.sh \
	getLibSpecifics.sh \
	processIntraMPs.sh \
	processInterMPs.sh \
	processInterMPsATA.sh \
	processClumps.sh \
	generateMatchRunFiles.sh \
	processInterMatches.sh \
	processInterRuns.sh \
	explainIntervals.sh \
	explainBreakpoints.sh \
	explainIntervalsInterChrom.sh \
	explainIntervalsInterChrom2.sh \
	fixOffsets.sh \
	ivals2Intervals.sh \
	dumpChromFiles.sh \
	filterExplain.sh \
	isolateIntraScaffoldUnmapped.sh \
	convertICRuns.sh \
	processInterClumpRunsFile.sh \
	consolidateInterClumpRunsResults.sh \
	combineCSVs.sh \
	filterIntersWithSatisfieds.sh \
	processHiLoInterChromIntervalFiles.sh

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS) $(SCRIPTS)
	cp -f atai2Lengths.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/atai2Lengths.sh
	cp -f breakpoints2StartLenNs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/breakpoints2StartLenNs.sh
	cp -f getDoubleCountedInTranpositions.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getDoubleCountedInTranpositions.sh
	cp -f analyzeGaps.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/analyzeGaps.sh
	cp -f getInterFromResults.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getInterFromResults.sh
	cp -f getInterFromResultsATA.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getInterFromResultsATA.sh
	cp -f getInterToResults.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getInterToResults.sh
	cp -f getIntraResults.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getIntraResults.sh
	cp -f morphInter.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/morphInter.sh
	cp -f processMPs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processMPs.sh
	cp -f processAssemblyMPs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processAssemblyMPs.sh
	cp -f buildMapStoreFromMatches.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/buildMapStoreFromMatches.sh
	cp -f getAllAssemblyMPs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getAllAssemblyMPs.sh
	cp -f findBacDeletions.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/findBacDeletions.sh
	cp -f getLibSpecifics.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/getLibSpecifics.sh
	cp -f processIntraMPs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processIntraMPs.sh
	cp -f processInterMPs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processInterMPs.sh
	cp -f processInterMPsATA.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processInterMPsATA.sh
	cp -f processClumps.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processClumps.sh
	cp -f generateMatchRunFiles.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/generateMatchRunFiles.sh
	cp -f processInterMatches.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processInterMatches.sh
	cp -f processInterRuns.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processInterRuns.sh
	cp -f explainIntervals.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/explainIntervals.sh
	cp -f explainBreakpoints.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/explainBreakpoints.sh
	cp -f explainIntervalsInterChrom.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/explainIntervalsInterChrom.sh
	cp -f explainIntervalsInterChrom2.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/explainIntervalsInterChrom2.sh
	cp -f fixOffsets.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/fixOffsets.sh
	cp -f ivals2Intervals.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/ivals2Intervals.sh
	cp -f dumpChromFiles.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/dumpChromFiles.sh
	cp -f filterExplain.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/filterExplain.sh
	cp -f isolateIntraScaffoldUnmapped.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/isolateIntraScaffoldUnmapped.sh
	cp -f convertICRuns.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/convertICRuns.sh
	cp -f processInterClumpRunsFile.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processInterClumpRunsFile.sh
	cp -f consolidateInterClumpRunsResults.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/consolidateInterClumpRunsResults.sh
	cp -f combineCSVs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/combineCSVs.sh
	cp -f filterIntersWithSatisfieds.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/filterIntersWithSatisfieds.sh
	cp -f processHiLoInterChromIntervalFiles.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/processHiLoInterChromIntervalFiles.sh

install: all

mapContigsScaffolds: mapContigsScaffolds.o 

getIntervalIntersections: getIntervalIntersections.o

pullFragmentLocations: pullFragmentLocations.o

getUniqueMatchesFromCGM: getUniqueMatchesFromCGM.o 

matchCounts: matchCounts.o

processIntra: processIntra.o clusterMPs.o

compareAssemblyMPs: compareAssemblyMPs.o clusterMPs.o

processInter: processInter.o clusterMPs.o

processInterCG: processInterCG.o clusterMPs.o

testConvexHull: testConvexHull.o

findCliques: findCliques.o

findMPSets: findMPSets.o

findTTIs: findTTIs.o

quads2mbrs: quads2mbrs.o

rotatePoints: rotatePoints.o

