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

BOG_SOURCES = AS_BOG_ArrivalRate.cc \
              AS_BOG_BestOverlapGraph.cc \
              AS_BOG_Breaking.cc \
              AS_BOG_ChunkGraph.cc \
              AS_BOG_Instrumentation.cc \
              AS_BOG_IntersectBubble.cc \
              AS_BOG_IntersectSplit.cc \
              AS_BOG_Joining.cc \
              AS_BOG_MateBubble.cc \
              AS_BOG_MateChecker.cc \
              AS_BOG_Outputs.cc \
              AS_BOG_PlaceContains.cc \
              AS_BOG_PlaceZombies.cc \
              AS_BOG_PopulateUnitig.cc \
              AS_BOG_SetParentAndHang.cc \
              AS_BOG_Unitig.cc \
              AS_BOG_UnitigGraph.cc
BOG_OBJECTS = $(BOG_SOURCES:.cc=.o)

SOURCES     = BuildUnitigs.cc fixUnitigs.C scoreOverlaps.C $(BOG_SOURCES)
OBJECTS     = BuildUnitigs.o  fixUnitigs.o scoreOverlaps.o $(BOG_OBJECTS)

CXX_PROGS = buildUnitigs fixUnitigs scoreOverlaps

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

buildUnitigs:    BuildUnitigs.o $(BOG_OBJECTS) libCA.a
fixUnitigs:      fixUnitigs.o                  libCA.a
scoreOverlaps:   scoreOverlaps.o               libCA.a
