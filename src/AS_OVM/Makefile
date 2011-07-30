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

OVM_SOURCES = overlapInCore-Build_Hash_Index.C \
              overlapInCore-Extend_Alignment.C \
              overlapInCore-Find_Overlaps.C \
              overlapInCore-Output.C \
              overlapInCore-Process_Overlaps.C \
              overlapInCore-Process_String_Overlaps.C \
              overlapInCore-Read_Next_Frag.C \
              overlapInCore.C
OVM_OBJECTS = $(OVM_SOURCES:.C=.o)

SOURCES     = $(OVM_SOURCES)
OBJECTS     = $(OVM_OBJECTS)

CXX_PROGS   = overlapInCore

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

overlapInCore: $(OVM_OBJECTS) libCA.a
