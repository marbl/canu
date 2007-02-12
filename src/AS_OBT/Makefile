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
# gmake Makefile for AS_OBT, Overlap Based Trimming

LOCAL_WORK = $(shell cd ../..; pwd)

SOURCES = chimera.C consolidate.C dumpFragStoreAsFasta.C dumpFragStoreAsOFG.c \
          dumpRandomMatedFrags.C \
          findDuplicateFragsInStore.C getNumScaffolds.c initialTrim.C \
          merge-trimming.C sort-overlaps.C \
          overlap-convert.C trim.C util++.C \
          overlap-filter.C
OBJECTS = $(SOURCES:.C=.o) $(SOURCES:.c=.o)
PROGS   = chimera consolidate initialTrim merge-trimming sort-overlaps \
          overlap-convert \
          overlap-filter \
          dumpRandomMatedFrags \
          dumpFragStoreAsFasta dumpFragStoreAsOFG findDuplicateFragsInStore getNumScaffolds
SCRIPTS = 

all: $(PROGS) $(SCRIPTS)

initialTrim: initialTrim.o trim.o -lAS_PER -lAS_UTL

#sort-overlap-trim: sort-overlap-trim.o util++.o

merge-trimming: merge-trimming.o trim.o -lAS_PER -lAS_UTL

sort-overlaps: sort-overlaps.o util++.o

overlap-convert: overlap-convert.o util++.o

overlap-filter: overlap-filter.o util++.o

chimera: chimera.o util++.o -lAS_PER -lAS_UTL

consolidate: consolidate.o util++.o

dumpFragStoreAsFasta: dumpFragStoreAsFasta.o trim.o -lAS_PER -lAS_UTL

dumpFragStoreAsOFG: dumpFragStoreAsOFG.o -lAS_PER -lAS_UTL

dumpRandomMatedFrags: dumpRandomMatedFrags.o -lAS_PER -lAS_UTL

findDuplicateFragsInStore: findDuplicateFragsInStore.o -lAS_PER -lAS_UTL

getNumScaffolds: getNumScaffolds.o libCA.a

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

