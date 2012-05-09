
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

OVS_LIB_SRC = AS_OVS_overlap.C AS_OVS_overlapFile.C AS_OVS_overlapStore.C
OVS_STR_SRC = overlapStore.C overlapStore_dump.C overlapStore_erates.C overlapStore_genomeLength.C
OVS_STA_SRC = overlapStats.C
OVS_CVT_SRC = convertOverlap.C filterOverlap.C

OVS_LIB_OBJ = $(OVS_LIB_SRC:.C=.o)
OVS_STR_OBJ = $(OVS_STR_SRC:.C=.o)
OVS_STA_OBJ = $(OVS_STA_SRC:.C=.o)
OVS_CVT_OBJ = $(OVS_CVT_SRC:.C=.o)

SOURCES = $(OVS_LIB_SRC) $(OVS_STR_SRC) $(OVS_STA_SRC) $(OVS_CVT_SRC) overlapStoreBucketizer.C overlapStoreSorter.C overlapStoreIndexer.C overlapStoreBuild.C
OBJECTS = $(SOURCES:.C=.o)

LIBRARIES = libAS_OVS.a libCA.a

CXX_PROGS = overlapStore overlapStats convertOverlap filterOverlap overlapStoreBucketizer overlapStoreSorter overlapStoreIndexer overlapStoreBuild

include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

libAS_OVS.a: $(OVS_LIB_OBJ)
libCA.a:     $(OVS_LIB_OBJ)

overlapStore:   $(OVS_STR_OBJ)   libCA.a
overlapStoreBucketizer:   overlapStoreBucketizer.o   libCA.a
overlapStoreSorter:       overlapStoreSorter.o       libCA.a
overlapStoreIndexer:      overlapStoreIndexer.o      libCA.a
overlapStoreBuild:        overlapStoreBuild.o        libCA.a

overlapStats:   overlapStats.o   libCA.a
convertOverlap: convertOverlap.o libCA.a
filterOverlap:  filterOverlap.o  libCA.a
