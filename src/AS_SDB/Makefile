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
VERBOSE=1
#######################################################################

#     SOURCES = 
#       Lists all sources compiled by this makefile.
#
#     OBJECTS = 
#       Lists all objects compiled by this makefile (i.e., one-to-one
#       with SOURCES). OBJECTS can be listed explicitly or can be
#       specified with the following definition:
#
#          OBJECTS = $(SOURCES:.c=.o)
#
#     LIBRARIES = 
#       Lists simple name (in "libx.a" form) of each library that is
#       built by this makefile.
#
#     PROGS = 
#       Lists each executable created by this makefile.
#
#     LOCAL_WORK =
#       Identifies the <subsystem_level> directory (parent of
#       the inc, src, lib, obj, and bin subdirectories) with
#       an absolute path.
#       For a makefile located in
#         <subsystem_level>/src/<component_level_1>
#       the following definition would work:
#          LOCAL_WORK = $(shell cd ../..; pwd)

#
# Also include dependencies of programs and libraries on objects and
# libraries, e.g.:
#
#    my_program: $(OBJECTS) libxxx.a
#

# Also include as first target, an "all" target, e.g.:
#
#    all: $(PROGS)
#


LIN_ALG_LIB = libAS_LIN.a
#ifeq ($(OSTYPE), osf3.2)
#LIN_ALG_LIB = libdxml.so
#endif

LOCAL_WORK = $(shell cd ../../; pwd)

SDB_SOURCES = BuildSequenceDB.c
SDB_OBJECTS = $(SDB_SOURCES:.c=.o) 

SDB_LIB_SOURCES = AS_SDB_SequenceDBPartition.c AS_SDB_SequenceDB.c 
LIB_SOURCES = $(SDB_LIB_SOURCES) 
SDB_LIB_OBJECTS = $(SDB_LIB_SOURCES:.c=.o) 

TST_SOURCES = LookAtSDB.c ViewMultiAlignment.c AlignContigs.c BeatOnSDB.c FetchContigs.c PartitionSequenceDB1.c  PartitionSequenceDB2.c \
	TestSequenceDBPartition.c 

OBJECTS = $(SDB_OBJECTS) $(SDB_LIB_OBJECTS) $(TST_SOURCES)
SOURCES = $(SDB_SOURCES) $(SDB_LIB_SOURCES) $(TST_SOURCES)


LIBRARIES = libAS_SDB.a

#PROGS = AlignContigs 
PROGS = TestSDBP BuildSequenceDB LookAtSDB ViewMultiAlignT AlignContigs BeatOnSDB FetchContigs PartitionSDB1 PartitionSDB2

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

#-include make.conf


all: $(OBJECTS) $(LIBRARIES) $(PROGS)

libAS_SDB.a: $(SDB_LIB_OBJECTS)

TestSDBP: TestSequenceDBPartition.o libAS_SDB.a libAS_CNS.a libAS_SDB.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
BeatOnSDB: BeatOnSDB.o libAS_SDB.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
AlignContigs: AlignContigs.o libAS_SDB.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_REZ.a libAS_CNS.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
FetchContigs: FetchContigs.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
BuildSequenceDB: BuildSequenceDB.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
LookAtSDB: LookAtSDB.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
ViewMultiAlignT: ViewMultiAlignment.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_REZ.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
PartitionSDB1: PartitionSequenceDB1.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm
PartitionSDB2: PartitionSequenceDB2.o libAS_SDB.a libAS_CNS.a libAS_REZ.a libAS_SDB.a libAS_CNS.a libAS_PER.a  libAS_UTL.a libAS_SDB.a libAS_MSG.a libAS_ALN.a libAS_CGW.a $(LIN_ALG_LIB) -lm

