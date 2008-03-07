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
LOCAL_WORK = $(shell cd ../..; pwd)

#
# gmake Makefile for AS_CNS, Consensus Modules

#  SEQ_AN should point to a configured, compiled and installed copy of
#  http://www.seqan.de/
#
#  We attempt to find it in the same place that our src/ directory is.


#  We'd like to include the global include (c_make.as) to get these
#  definitions, but we cannot without defining rules.

OSTYPE      = $(shell echo `uname`)
MACHINETYPE = $(shell echo `uname -m`)

ifeq ($(MACHINETYPE), x86_64)
  MACHINETYPE = amd64
endif
ifeq ($(MACHINETYPE), Power Macintosh)
  MACHINETYPE = ppc
endif

#  This isn't perfect; if we're building debug here, we _usually_ want
#  to use the debug SeqAn.  But checking if SeqAn exists, and failing if
#  we don't find the debug SeqAn, is a little tricky.
ifneq "$(origin SEQ_AN)" "environment"

ifeq ($(shell ls -d $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE)-debug 2> /dev/null), $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE)-debug)
  SEQ_AN = $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE)-debug
endif

ifeq ($(shell ls -d $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE) 2> /dev/null), $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE))
  SEQ_AN = $(LOCAL_WORK)/seqan/$(OSTYPE)-$(MACHINETYPE)
endif

endif

LIB_SOURCES = Array_CNS.c \
              MultiAlign.c \
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
              SeqAn_CNS.c

AS_CNS_OBJS = $(AS_CNS_SRCS:.c=.o)

SOURCES   = $(AS_CNS_SRCS)
OBJECTS   = $(AS_CNS_OBJS)
LIBRARIES = libAS_CNS.a libCA.a

PROGS = consensus \
        asmReBaseCall \
	SeqAn_CNS

SCRIPTS =

ifdef SEQ_AN
  SCRIPTS += graph_consensus 
endif

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(LIB_OBJECTS) $(LIBRARIES) $(OBJECTS) $(PROGS) $(SCRIPTS)
	@test -n nop

libAS_CNS.a: $(LIB_OBJECTS)
libCA.a:     $(LIB_OBJECTS)

LIBS =  libCA.a

consensus:                    Consensus_CNS.o                $(LIBS)
asmReBaseCall:                AS_CNS_asmReBaseCall.o         $(LIBS)
SeqAn_CNS:                    SeqAn_CNS.o                    $(LIBS)
graph_consensus: 	      $(SEQ_AN)/demos/graph_consensus
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

