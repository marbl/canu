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

OVS_LIB_SRC = AS_OVS_overlapFile.c AS_OVS_overlapStore.c
OVS_STR_SRC = overlapStore.c overlapStore_build.c overlapStore_merge.c overlapStore_dump.c overlapStore_stats.c overlapStore_erates.c
OVS_CVT_SRC = convertOverlap.c 

OVS_LIB_OBJ = $(OVS_LIB_SRC:.c=.o)
OVS_STR_OBJ = $(OVS_STR_SRC:.c=.o)
OVS_CVT_OBJ = $(OVS_CVT_SRC:.c=.o)

SOURCES = $(OVS_LIB_SRC) $(OVS_STR_SRC) $(OVS_CVT_SRC)
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_OVS.a libAS.a

PROGS = convertOverlap overlapStore

include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop

libAS_OVS.a: $(OVS_LIB_OBJ)
libCA.a:     $(OVS_LIB_OBJ)

overlapStore:   $(OVS_STR_OBJ) libCA.a
convertOverlap: $(OVS_CVT_OBJ) libCA.a
