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

LOCAL_WORK = $(shell cd ../../; pwd)

LIB_SOURCES = AS_ALN_dpaligner.c \
              AS_ALN_qvaligner.c \
	      AS_ALN_loverlapper.c \
              AS_ALN_pieceOlap.c AS_ALN_forcns.c \
	      CA_ALN_local.c \
              CA_ALN_overlap.c \
              CA_ALN_scafcomp.c \
              AS_ALN_bruteforcedp.c

LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

SOURCES = $(LIB_SOURCES)
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_ALN.a libCA.a

PROGS = 

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop

libAS_ALN.a: $(LIB_OBJECTS)
libCA.a: $(LIB_OBJECTS)
