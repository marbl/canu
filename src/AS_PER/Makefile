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

SOURCES = AS_PER_genericStore.C \
          AS_PER_gkLibrary.C \
          AS_PER_gkStore.C \
          AS_PER_gkStore_IID.C \
          AS_PER_gkStore_PLC.C \
          AS_PER_gkStore_UID.C \
          AS_PER_gkStore_clearRange.C \
          AS_PER_gkStore_fragments.C \
          AS_PER_gkStore_load.C \
          AS_PER_gkStore_partition.C \
          AS_PER_gkStore_stats.C \
          AS_PER_gkStream.C \
          AS_PER_encodeSequenceQuality.C

OBJECTS = $(SOURCES:.C=.o)

LIBRARIES = libAS_PER.a libCA.a

CXX_PROGS =

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

objects: $(OBJECTS)

libAS_PER.a: $(OBJECTS)

libCA.a: $(OBJECTS)

