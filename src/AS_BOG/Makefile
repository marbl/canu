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

BOG_SOURCES = AS_BOG_BestOverlapGraph.C \
              AS_BOG_Breaking.C \
              AS_BOG_ChunkGraph.C \
              AS_BOG_EvaluateMates.C \
              AS_BOG_FragmentInfo.C \
              AS_BOG_Instrumentation.C \
              AS_BOG_IntersectBubble.C \
              AS_BOG_IntersectSplit.C \
              AS_BOG_InsertSizes.C \
              AS_BOG_Joining.C \
              AS_BOG_MateBubble.C \
              AS_BOG_MateChecker.C \
              AS_BOG_MateLocation.C \
              AS_BOG_MoveContains.C \
              AS_BOG_Outputs.C \
              AS_BOG_OverlapBubble.C \
              AS_BOG_PlaceContains.C \
              AS_BOG_PlaceFragUsingOverlaps.C \
              AS_BOG_PlaceZombies.C \
              AS_BOG_PopulateUnitig.C \
              AS_BOG_SetParentAndHang.C \
              AS_BOG_SplitDiscontinuous.C \
              AS_BOG_Unitig_AddAndPlaceFrag.C \
              AS_BOG_Unitig_AddFrag.C \
              AS_BOG_Unitig_PlaceFragUsingEdges.C \
              AS_BOG_Unitig.C \
              AS_BOG_UnitigGraph.C
BOG_OBJECTS = $(BOG_SOURCES:.C=.o)

SOURCES     = BuildUnitigs.C  fixUnitigs.C analyzeBest.C $(BOG_SOURCES)
OBJECTS     = BuildUnitigs.o  fixUnitigs.o analyzeBest.o $(BOG_OBJECTS)

CXX_PROGS = buildUnitigs fixUnitigs analyzeBest

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

buildUnitigs:    BuildUnitigs.o $(BOG_OBJECTS) libCA.a
fixUnitigs:      fixUnitigs.o                  libCA.a
analyzeBest:     analyzeBest.o                 libCA.a
