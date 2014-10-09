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

LOCAL_WORK = $(shell cd ../..; pwd)

LIB_SOURCES = AS_UTL_Hash.C \
              AS_UTL_heap.C \
              AS_UTL_Var.C \
              AS_UTL_rand.C \
              AS_UTL_interval.C \
              UnionFind_AS.C \
              AS_UTL_skiplist.C \
              AS_UTL_alloc.C \
              AS_UTL_fileIO.C \
              AS_UTL_fasta.C \
              AS_UTL_UID.C \
              AS_UTL_reverseComplement.C \
              AS_UTL_decodeRange.C \
              AS_UTL_stackTrace.C

LIB_OBJECTS = $(LIB_SOURCES:.C=.o)

SOURCES     = $(LIB_SOURCES) 

OBJECTS     = $(SOURCES:.C=.o)
LIBRARIES   = libAS_UTL.a libCA.a
SCRIPTS     =
CXX_PROGS   = 

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

libAS_UTL.a: $(LIB_OBJECTS)
libCA.a: $(LIB_OBJECTS)

.PHONY: test
test:
	cc -O3 -o testHashTable -I.. -I. testHashTable.C AS_UTL_Hash.C AS_UTL_heap.C AS_UTL_alloc.C AS_UTL_fileIO.C -lm
	cc -O3 -o testRand      -I.. -I. testRand.C      AS_UTL_rand.C                                              -lm
	cc -O3 -o testVar       -I.. -I. testVar.C       AS_UTL_Var.C                AS_UTL_alloc.C AS_UTL_fileIO.C -lm
