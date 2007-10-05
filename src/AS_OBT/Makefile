#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2005-2007, J. Craig Venter Institute.
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

SOURCES = acceptableOBToverlap.c \
          chimera.C \
          consolidate.C \
          prefixDelete.C \
          getNumScaffolds.c \
          initialTrim.C \
          merge-trimming.C \
          trim.C \
          readOverlap.C \
          util++.C

OBJECTS = $(SOURCES:.C=.o) $(SOURCES:.c=.o)
PROGS   = acceptableOBToverlap \
          initialTrim \
          merge-trimming \
          consolidate \
          chimera \
          prefixDelete \
          getNumScaffolds
SCRIPTS = 

include $(LOCAL_WORK)/src/c_make.as

all: $(PROGS) $(SCRIPTS)
	@test -n nop

acceptableOBToverlap: acceptableOBToverlap.o                    -lAS_OVS -lAS_PER -lAS_UTL
initialTrim:          initialTrim.o trim.o                               -lAS_PER -lAS_UTL
consolidate:          consolidate.o util++.o readOverlap.o      -lAS_OVS -lAS_PER -lAS_UTL
merge-trimming:       merge-trimming.o trim.o                            -lAS_PER -lAS_UTL
chimera:              chimera.o util++.o readOverlap.o          -lAS_OVS -lAS_PER -lAS_UTL

prefixDelete:         prefixDelete.o                                     -lAS_PER -lAS_UTL

getNumScaffolds: getNumScaffolds.o libCA.a
