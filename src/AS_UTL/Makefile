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

# gmake Makefile for AS_UTL

LOCAL_WORK = $(shell cd ../..; pwd)

LIBRARIES = libAS_UTL.a

REQ_LIB = libAS_MSG.a

SCRIPTS = AS_UTL_extractMSG.awk extractMSG 

PROGS = testVar \
	queryPHashTable \
	testRand \
	testSkiplist \
	testBuckets \
	testHashTable
# testhisto
# testPHashTable
# testfname

LIB_SOURCES = AS_UTL_HashCommon.c AS_UTL_Hash.c AS_UTL_PHash.c \
	AS_UTL_heap.c AS_UTL_Var.c AS_UTL_rand.c AS_UTL_version.c \
	AS_UTL_interval.c UnionFind_AS.c AS_UTL_skiplist.c \
        AS_UTL_ID_store.c AS_UTL_SequenceBucket.c AS_UTL_param_proc.c \
        getopt.c getopt1.c AS_UTL_alloc.c
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

SL_SOURCES = AS_UTL_skiplist_test.c
SL_OBJECTS = $(SL_SOURCES:.c=.o)

HASH_SOURCES = testHashTable.c 
HASH_OBJECTS = $(HASH_SOURCES:.c=.o)

VAR_SOURCES = testVar.c
VAR_OBJECTS = $(VAR_SOURCES:.c=.o)

PHASH_SOURCES = queryPHashTable.c
PHASH_OBJECTS = $(PHASH_SOURCES:.c=.o)

RAND_SOURCES = testRand.c
RAND_OBJECTS = $(RAND_SOURCES:.c=.o)

#HISTO_SOURCES = AS_UTL_testhisto.c
#HISTO_OBJECTS = $(HISTO_SOURCES:.c=.o)

BUCKET_SOURCES = testSequenceBucket.c
BUCKET_OBJECTS = $(BUCKET_SOURCES:.c=.o)

FNAME_SOURCES = AS_UTL_fileNames.c
FNAME_OBJECTS = $(FNAME_SOURCES:.c=.o)

SOURCES = $(LIB_SOURCES) $(HASH_SOURCES) $(PHASH_SOURCES) \
	$(RAND_SOURCES) $(SL_SOURCES) $(VAR_SOURCES) $(BUCKET_SOURCES) 
#$(FNAME_SOURCES)

OBJECTS = $(SOURCES:.c=.o)

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

install: $(SCRIPTS)
	cp -f $(SCRIPTS) $(LOCAL_BIN)

testBuckets: $(BUCKET_OBJECTS) $(LIB_OBJECTS) libm.a $(REQ_LIB)

#testfname: $(FNAME_OBJECTS) $(LIB_OBJECTS) libm.a $(REQ_LIB)

testSkiplist: $(SL_OBJECTS) $(LIB_OBJECTS) libm.a $(REQ_LIB)

testHashTable: $(HASH_OBJECTS) $(LIB_OBJECTS) libm.a  $(REQ_LIB)

# testPHashTable: $(PHASH_OBJECTS) $(LIB_OBJECTS) libm.a  $(REQ_LIB)

queryPHashTable: $(PHASH_OBJECTS) $(LIB_OBJECTS) libm.a  $(REQ_LIB)

testRand: $(RAND_OBJECTS) $(LIB_OBJECTS) libm.a  $(REQ_LIB)

testVar: $(VAR_OBJECTS) $(LIB_OBJECTS) libm.a  $(REQ_LIB)

##testhisto: $(HISTO_OBJECTS) $(LIB_OBJECTS) libm.a $(REQ_LIB)

libAS_UTL.a: $(LIB_OBJECTS)

extractMSG: AS_UTL_extractMSG
	cp -f AS_UTL_extractMSG extractMSG

