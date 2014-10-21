#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# Copyright (C) 2005, J Craig Venter Institute. All rights reserved.
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

LOCAL_WORK = $(shell cd ../..; pwd)

BAT_SOURCES = AS_BAT_BestOverlapGraph.C \
              AS_BAT_Breaking.C \
              AS_BAT_ChunkGraph.C \
              AS_BAT_EvaluateMates.C \
              AS_BAT_ExtendByMates.C \
              AS_BAT_FragmentInfo.C \
              AS_BAT_Instrumentation.C \
              AS_BAT_IntersectBubble.C \
              AS_BAT_IntersectSplit.C \
              AS_BAT_InsertSizes.C \
              AS_BAT_Joining.C \
              AS_BAT_Logging.C \
              AS_BAT_MateBubble.C \
              AS_BAT_MateChecker.C \
              AS_BAT_MateLocation.C \
              AS_BAT_MergeSplitJoin.C \
              AS_BAT_MoveContains.C \
              AS_BAT_Outputs.C \
              AS_BAT_OverlapCache.C \
              AS_BAT_PlaceContains.C \
              AS_BAT_PlaceFragUsingOverlaps.C \
              AS_BAT_PlaceZombies.C \
              AS_BAT_PopulateUnitig.C \
              AS_BAT_PromoteToSingleton.C \
              AS_BAT_ReconstructRepeats.C \
              AS_BAT_SetParentAndHang.C \
              AS_BAT_SplitDiscontinuous.C \
              AS_BAT_Unitig_AddAndPlaceFrag.C \
              AS_BAT_Unitig_AddFrag.C \
              AS_BAT_Unitig_PlaceFragUsingEdges.C \
              AS_BAT_Unitig.C
BAT_OBJECTS = $(BAT_SOURCES:.C=.o)

DNC_SOURCES = classifyMates-BFS.C \
              classifyMates-DFS.C \
              classifyMates-RFS.C \
              classifyMates-Suspicious.C \
              classifyMates-globalData.C
DNC_OBJECTS = $(DNC_SOURCES:.C=.o)


SOURCES     = bogart.C bogus.C bogusness.C AS_BAT_bogusUtil.C classifyMates.C classifyMatesApply.C classifyMatesPairwise.C splitUnitigs.C computeCoverageStat.C markRepeatUnique.C rewriteCache.C erate-estimate.C $(BAT_SOURCES) $(DNC_SOURCES)
OBJECTS     = bogart.o bogus.o bogusness.o AS_BAT_bogusUtil.o classifyMates.o classifyMatesApply.o classifyMatesPairwise.o splitUnitigs.o computeCoverageStat.o markRepeatUnique.o rewriteCache.o erate-estimate.o $(BAT_OBJECTS) $(DNC_OBJECTS)

CXX_PROGS   = bogart bogus bogusness classifyMates classifyMatesApply classifyMatesPairwise splitUnitigs computeCoverageStat markRepeatUnique rewriteCache erate-estimate



# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

bogart:                 bogart.o                $(BAT_OBJECTS) libCA.a
#petey:                  petey.o                 $(BAT_OBJECTS) libCA.a
bogus:                  bogus.o                                libCA.a AS_BAT_bogusUtil.o
bogusness:              bogusness.o                            libCA.a AS_BAT_bogusUtil.o
classifyMates:          classifyMates.o         $(DNC_OBJECTS) libCA.a libutil.a
classifyMatesApply:     classifyMatesApply.o    $(DNC_OBJECTS) libCA.a libutil.a
classifyMatesPairwise:  classifyMatesPairwise.o $(DNC_OBJECTS) libCA.a libutil.a
splitUnitigs:           splitUnitigs.o                         libCA.a
computeCoverageStat:    computeCoverageStat.o                  libCA.a
markRepeatUnique:       markRepeatUnique.o                     libCA.a
rewriteCache:           rewriteCache.o AS_BAT_FragmentInfo.o AS_BAT_Logging.o libCA.a
erate-estimate:         erate-estimate.o                       libCA.a
