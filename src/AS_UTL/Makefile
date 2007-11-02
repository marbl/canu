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

LIB_SOURCES = AS_UTL_Hash.c AS_UTL_heap.c AS_UTL_Var.c AS_UTL_rand.c AS_UTL_version.c AS_UTL_interval.c UnionFind_AS.c AS_UTL_skiplist.c AS_UTL_SequenceBucket.c getopt.c getopt1.c AS_UTL_alloc.c AS_UTL_fileIO.c AS_UTL_qsort_mt.c AS_UTL_fasta.c
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

SOURCES     = AS_UTL_skiplist_test.c testHashTable.c testVar.c testRand.c testSequenceBucket.c $(LIB_SOURCES)
OBJECTS     = $(SOURCES:.c=.o)
LIBRARIES   = libAS_UTL.a libCA.a
SCRIPTS     = AS_UTL_extractMSG.awk extractMSG 
PROGS       = testVar testRand testSkiplist testBuckets testHashTable

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install
	@test -n nop

install: $(SCRIPTS)
	@cp -f $(SCRIPTS) $(LOCAL_BIN)

libAS_UTL.a: $(LIB_OBJECTS)
libCA.a: $(LIB_OBJECTS)

testBuckets:    testSequenceBucket.o libAS_UTL.a
testSkiplist:   AS_UTL_skiplist_test.o libAS_UTL.a
testHashTable:  testHashTable.o AS_UTL_Hash.o libAS_UTL.a
testRand:       testRand.o libAS_UTL.a
testVar:        testVar.o libAS_UTL.a

extractMSG: AS_UTL_extractMSG
	cp -f AS_UTL_extractMSG extractMSG
