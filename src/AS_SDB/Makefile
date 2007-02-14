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

LOCAL_WORK = $(shell cd ../../; pwd)

LIB_SOURCES     = AS_SDB_SequenceDBPartition.c AS_SDB_SequenceDB.c
LIB_OBJECTS     = $(LIB_SOURCES:.c=.o)

SOURCES = $(LIB_SOURCES) \
          SDButil.c \
          PartitionSDB.c \
          PartitionSequenceDB1.c \
          PartitionSequenceDB2.c
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_SDB.a libCA.a

PROGS = SDButil \
        PartitionSDB \
        PartitionSDB1 \
        PartitionSDB2

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop

libAS_SDB.a: $(LIB_OBJECTS)

libCA.a: $(LIB_OBJECTS)

SDButil:                      SDButil.o                       libCA.a
PartitionSDB:                 PartitionSDB.o                  libCA.a
PartitionSDB1:                PartitionSequenceDB1.o          libCA.a
PartitionSDB2:                PartitionSequenceDB2.o          libCA.a
