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

LOCAL_WORK  = $(shell cd ../..; pwd)

LIB_SOURCES = gkStore.C
LIB_OBJECTS = $(LIB_SOURCES:.C=.o)

SOURCES   = $(LIB_SOURCES) gatekeeperCreate.C gatekeeperDumpFASTQ.C
OBJECTS   = $(LIB_OBJECTS) gatekeeperCreate.o gatekeeperDumpFASTQ.o
LIBRARIES = libCA.a

CXX_PROGS = gatekeeperCreate gatekeeperDumpFASTQ

SCRIPTS =

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(LIB_OBJECTS) $(LIBRARIES) $(OBJECTS) $(CXX_PROGS) $(SCRIPTS)
	@test -n nop

libCA.a:     $(LIB_OBJECTS)


gatekeeperCreate:       gatekeeperCreate.o             libCA.a
gatekeeperDumpFASTQ:    gatekeeperDumpFASTQ.o          libCA.a
