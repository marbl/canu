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
# gmake Makefile for AS_PER

LOCAL_WORK = $(shell cd ../..; pwd)

AS_PER_SOURCES = AS_PER_ReadStruct.c AS_PER_fragStore.c \
	AS_PER_genericStore.c AS_PER_SafeIO.c AS_PER_gkpStore.c \
	AS_PER_encodeSequenceQuality.c AS_PER_fragStorePartition.c \
	AS_PER_asmStore.c

AS_PER_TEST_SOURCES = testFragStore.c dumpFragStore.c testLoad.c \
	FragStoreStats.c binaryDumpFragStore.c loadBinaryDump.c \
	AS_PER_fragDB.c dumpFragDB.c PopulateFragStore.c \
	partitionFragStore.c dumpFrag.c dumpFragFromPartition.c\
	testFrgStoreSizes.c test_include.c resetFrgStoreToCNSClr.c \
        dumpFragClearLens.c

#stressFragStore.c 

SOURCES = $(AS_PER_SOURCES) $(AS_PER_TEST_SOURCES) 

AS_PER_OBJS = $(AS_PER_SOURCES:.c=.o)

AS_PER_TEST_OBJS = $(AS_PER_TEST_SOURCES:.c=.o)

OBJECTS = $(AS_PER_OBJS) $(AS_PER_TEST_OBJS)

LIBRARIES = libAS_PER.a libCA.a

PROGS = testFragStore dumpFragStore testLoad FragStoreStats \
	dumpFragDB PopulateFragStore partitionFragStore \
	dumpFrag dumpFragFromPartition testFrgStoreSizes \
	resetFrgStoreToCNSClr         dumpFragClearLens

#binaryDumpFragStore loadBinaryDump 
#stressFragStore 

INC_IMPORT_DIRS = \
	$(ORACLE_HOME)/rdbms/demo $(ORACLE_HOME)/network/public

#LLIBCLNTSH = /home/flanigmj/sandbox/oracle/lib/libclntsh.so.8.0
OCISHAREDLIBS=$(LLIBCLNTSH)

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

# FragStore_version.
# Use this DEF for backward compatibility with old fragStore databases.
# If no version is defined here, code will assume 5.
# If you change the version, you should `make clean all` in AS_PER.
# This line must occur after the c_make.as include.
#CDEFS += -DFRAGSTORE_VERSION=2    # added branch-point info (1999)
#CDEFS += -DFRAGSTORE_VERSION=3    # added screen match info (1999)
#CDEFS += -DFRAGSTORE_VERSION=4    # added indexed store features (2000)
#CDEFS += -DFRAGSTORE_VERSION=5    # added multiple clear ranges (Oct 2001)
CDEFS += -DFRAGSTORE_VERSION=6     # remove time_t and old cruft 2005/7/19

all: $(OBJECTS) $(LIBRARIES) $(PROGS) 

install: all
	echo "Nothing to do for AS_PER"

objects: $(AS_PER_OBJS)

libAS_PER.a: $(AS_PER_OBJS)

libCA.a: $(AS_PER_OBJS)

# Test for Frag Store
# Should have additional tests for generic stores
testFragStore: $(AS_PER_OBJS) testFragStore.o libAS_UTL.a 

#stressFragStore: $(AS_PER_OBJS) stressFragStore.o AS_PER_fragDB.o libAS_UTL.a  $(OCISHAREDLIBS) 

testLoad: $(AS_PER_OBJS) testLoad.o libAS_UTL.a 

dumpFragDB: $(AS_PER_OBJS) AS_PER_fragDB.o dumpFragDB.o libAS_UTL.a $(OCISHAREDLIBS) 

dumpFragStore: $(AS_PER_OBJS)  dumpFragStore.o libAS_UTL.a 

dumpFragClearLens: $(AS_PER_OBJS) dumpFragClearLens.o libAS_UTL.a

dumpFrag: $(AS_PER_OBJS)  dumpFrag.o libAS_UTL.a 

dumpFragFromPartition: $(AS_PER_OBJS) dumpFragFromPartition.o libAS_UTL.a 

PopulateFragStore: $(AS_PER_OBJS)  PopulateFragStore.o libAS_MSG.a libAS_PER.a libAS_UTL.a

binaryDumpFragStore: $(AS_PER_OBJS) binaryDumpFragStore.o libAS_UTL.a 

loadBinaryDump: $(AS_PER_OBJS) loadBinaryDump.o libAS_UTL.a 

FragStoreStats: $(AS_PER_OBJS) FragStoreStats.o libAS_UTL.a 

partitionFragStore: $(AS_PER_OBJS) partitionFragStore.o libAS_UTL.a 

testFrgStoreSizes: $(AS_PER_OBJS) testFrgStoreSizes.o libAS_UTL.a

resetFrgStoreToCNSClr: $(AS_PER_OBJS) resetFrgStoreToCNSClr.o libAS_UTL.a
