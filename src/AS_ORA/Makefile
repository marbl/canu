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
#VERBOSE=1
# Makefile for AS_ORA

LOCAL_WORK = $(shell cd ../..; pwd)

AS_ORA_SRCS = 	AS_ORA_statistics.c \
		AS_ORA_repeats.c \
		AS_ORA_overlaps.c \
		AS_ORA_fragments.c

AS_ORA_OBJS = $(AS_ORA_SRCS:.c=.o)

ORA_MAIN_SRCS = AS_ORA_main.c  
ORA_MAIN_OBJS = $(ORA_MAIN_SRCS:.c=.o)

SOURCES = $(AS_ORA_SRCS) $(ORA_MAIN_SRCS)  get-olaps.c get-perfect-olaps.c nodups.c fraguid2iid.c fraguid2iid_byhash.c make-canonical.c compare-canonical.c fragiid2uid.c dump_fragment_as_fasta.c computeGCcontent.c dumpMatePairing.c

OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = 

PROGS = overlap_regressor get-olaps get-perfect-olaps nodups fraguid2iid fraguid2iid_byhash make-canonical compare-canonical fragiid2uid dump_fragment_as_fasta computeGCcontent dumpMatePairing

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

all: $(PROGS) $(LIBRARIES) $(OBJECTS)

get-olaps: $(AS_ORA_OBJS) get-olaps.o libAS_MSG.a libAS_PER.a libAS_UTL.a

get-perfect-olaps: $(AS_ORA_OBJS) get-perfect-olaps.o libAS_MSG.a libAS_PER.a libAS_UTL.a -lAS_SDB -lAS_CNS -lAS_REZ -lAS_ALN -lAS_MSG

nodups: $(AS_ORA_OBJS) nodups.o libAS_MSG.a libAS_PER.a libAS_UTL.a

overlap_regressor: $(AS_ORA_OBJS) $(ORA_MAIN_OBJS) libAS_MSG.a libAS_PER.a libAS_UTL.a 

fraguid2iid: $(AS_ORA_OBJS) fraguid2iid.o libAS_MSG.a libAS_PER.a libAS_UTL.a libAS_CGW.a libAS_LIN.a

fraguid2iid_byhash: $(AS_ORA_OBJS) fraguid2iid_byhash.o libAS_MSG.a libAS_PER.a libAS_UTL.a libAS_CGW.a libAS_LIN.a

fragiid2uid : fragiid2uid.o libAS_PER.a  libAS_UTL.a libAS_CGW.a libAS_REZ.a libAS_CGW.a libAS_REZ.a libAS_LIN.a -lAS_CNS -lAS_SDB -lAS_REZ -lAS_ALN -lAS_MSG -lAS_PER -lAS_OVL -lAS_UTL

dump_fragment_as_fasta : dump_fragment_as_fasta.o libAS_CGW.a libAS_REZ.a libAS_CGW.a libAS_REZ.a libAS_LIN.a -lAS_CNS -lAS_SDB -lAS_REZ -lAS_ALN -lAS_MSG -lAS_PER -lAS_OVL -lAS_UTL

dumpMatePairing : dumpMatePairing.o libAS_CGW.a libAS_REZ.a libAS_CGW.a libAS_REZ.a libAS_LIN.a -lAS_CNS -lAS_SDB -lAS_REZ -lAS_ALN -lAS_MSG -lAS_PER -lAS_OVL -lAS_UTL

computeGCcontent : computeGCcontent.o libAS_CGW.a libAS_REZ.a libAS_CGW.a libAS_REZ.a libAS_LIN.a -lAS_CNS -lAS_SDB -lAS_REZ -lAS_ALN -lAS_MSG -lAS_PER -lAS_OVL -lAS_UTL 

make-canonical: $(AS_ORA_OBJS) make-canonical.o libAS_MSG.a libAS_PER.a libAS_UTL.a

compare-canonical: $(AS_ORA_OBJS) compare-canonical.o libAS_MSG.a libAS_PER.a libAS_UTL.a

install: $(PROGS) $(LIBRARIES) $(OBJECTS)
	@echo "NO INSTALL"
