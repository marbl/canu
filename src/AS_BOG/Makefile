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
# gmake Makefile for AS_OVL

LOCAL_WORK = $(shell cd ../..; pwd)

TOOL_SRCS = BuildUnitigs.cc FindContainedFrags.cc

HEADERS = AS_BOG_BestOverlapGraph.hh AS_BOG_ChunkGraph.hh AS_BOG_UnitigGraph.hh
LIB_SOURCES = AS_BOG_BestOverlapGraph.cc AS_BOG_ChunkGraph.cc AS_BOG_UnitigGraph.cc 

SOURCES = $(LIB_SOURCES) $(TOOL_SRCS)

TOOL_OBJS = $(TOOL_SRCS:.cc=.o)
LIB_OBJS = $(LIB_SOURCES:.cc=.o) 
OBJECTS =  $(LIB_OBJS)  $(TOOL_OBJS)

#LIBRARIES = libAS_OVL.a

AS_LIBRARIES = libAS_PER.a libAS_UTL.a libAS_MSG.a libAS_OVL.a
#libAS_ALN.a libAS_CNS.a \
#	libAS_CGW.a libAS_REZ.a libAS_SDB.a libAS_CNS.a libAS_CGW.a \
#	libAS_REZ.a libAS_SDB.a libAS_ALN.a libAS_UTL.a libAS_PER.a \
#	libAS_OVL.a libAS_MSG.a libAS_LIN.a
SYS_LIBRARIES = 

TOOL_PROGS = buildUnitigs findContainedFrags

PROGS = $(TOOL_PROGS)


# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(PROGS) 

buildUnitigs: BuildUnitigs.o $(LIB_OBJS) $(AS_LIBRARIES) $(SYS_LIBRARIES) 

findContainedFrags:  FindContainedFrags.o $(LIB_OBJS) $(AS_LIBRARIES) $(SYS_LIBRARIES)

