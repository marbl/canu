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

LIB_SOURCES = Array_CNS.c \
              MultiAlign.c \
              MultiAlignFromCCO.c \
              MultiAlignFromICM.c \
              MultiAlignFromIUM.c \
              MultiAlignPrint.c \
              MultiAlignStore_CNS.c \
              MultiAlignment_CNS.c

LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

LOCAL_WORK  = $(shell cd ../..; pwd)

AS_CNS_SRCS = $(LIB_SOURCES) \
              Consensus_CNS.c \
              AS_CNS_asmReBaseCall.c \
              dumpUnitigConsensusFromIUM.c

AS_CNS_OBJS = $(AS_CNS_SRCS:.c=.o)

SOURCES   = $(AS_CNS_SRCS)
OBJECTS   = $(AS_CNS_OBJS)
LIBRARIES = libAS_CNS.a libCA.a

PROGS = consensus \
        asmReBaseCall \
        dumpUnitigConsensusFromIUM

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(LIB_OBJECTS) $(LIBRARIES) $(OBJECTS) $(PROGS)
	@test -n nop

libAS_CNS.a: $(LIB_OBJECTS)
libCA.a:     $(LIB_OBJECTS)

LIBS =  libCA.a

consensus:                    Consensus_CNS.o                $(LIBS)
asmReBaseCall:                AS_CNS_asmReBaseCall.o         $(LIBS)
dumpUnitigConsensusFromIUM:   dumpUnitigConsensusFromIUM.o   $(LIBS)
