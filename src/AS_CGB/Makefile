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
# Makefile for AS_CGB

LOCAL_WORK = $(shell cd ../..; pwd)

AS_FGB_SRCS = AS_FGB_main.c AS_FGB_io.c 

AS_OGB_SRCS = AS_CGB_histo.c \
              AS_CGB_fga.c AS_CGB_fgb.c AS_FGB_hanging_fragment.c AS_FGB_contained.c \
              AS_CGB_store.c AS_FGB_FragmentHash.c AS_FGB_buildFragmentHash.c  \
              AS_CGB_traversal.c AS_CGB_walk.c \
              AS_CGB_count_fragment_and_edge_labels.c \
              AS_CGB_edgemate.c

AS_CGB_SRCS = AS_CGB_miniunitigger.c AS_CGB_main.c AS_CGB_io.c \
              AS_CGB_cgb.c AS_CGB_classify.c AS_CGB_chimeras.c

AS_CGA_SRCS = AS_CGB_cga.c

AS_CGB_BUB_SRCS = AS_CGB_Bubble.c AS_CGB_Bubble_Graph.c \
                  AS_CGB_Bubble_VertexSet.c AS_CGB_Bubble_top.c \
	          AS_CGB_Bubble_dfs.c AS_CGB_Bubble_Popper.c \
	          AS_CGB_Bubble_PopperMethods.c

UNITIG_SRCS   = $(AS_OGB_SRCS) $(AS_CGB_SRCS) $(AS_CGA_SRCS) $(AS_CGB_BUB_SRCS)

AS_FGB_OBJS     = $(AS_FGB_SRCS:.c=.o)
AS_OGB_OBJS     = $(AS_OGB_SRCS:.c=.o)
AS_CGB_OBJS     = $(AS_CGB_SRCS:.c=.o)
AS_CGA_OBJS     = $(AS_CGA_SRCS:.c=.o)
AS_CGB_BUB_OBJS = $(AS_CGB_BUB_SRCS:.c=.o)

UNITIG_OBJS     = $(AS_OGB_OBJS) $(AS_CGB_OBJS) $(AS_CGA_OBJS) $(AS_CGB_BUB_OBJS)

LIB_SOURCES = $(UNITIG_SRCS) $(AS_FGB_SRCS) AS_CGB_breakers.c AS_CGB_util.c
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

SOURCES = $(LIB_SOURCES) \
          AS_CGB_unitigger.c \
          AS_CGB_iid2cid.c \
          AS_CGB_fom2uom.c \
          check_breakers.c \
          make_OFG_from_FragStore.c \
          miniunitigger_test.c
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_CGB.a

PROGS = unitigger \
        AS_CGB_fom2uom \
        AS_CGB_iid2cid \
        check_breakers \
        make_OFG_from_FragStore \
        miniunitigger_test

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

libAS_CGB.a: $(LIB_OBJECTS)

make_OFG_from_FragStore: make_OFG_from_FragStore.o libAS_PER.a libAS_UTL.a 

miniunitigger_test: miniunitigger_test.o $(LIB_OBJECTS) libAS_OVL.a libAS_PER.a libAS_ALN.a libAS_MSG.a libAS_UTL.a 

unitigger: AS_CGB_unitigger.o $(LIB_OBJECTS) libAS_OVL.a libAS_PER.a libAS_ALN.a libAS_MSG.a libAS_UTL.a 

AS_CGB_fom2uom:        AS_CGB_fom2uom.o AS_CGB_histo.o libAS_PER.a libAS_OVL.a libAS_MSG.a libAS_UTL.a 

AS_CGB_iid2cid:        AS_CGB_iid2cid.o AS_CGB_histo.o libAS_PER.a libAS_OVL.a libAS_MSG.a libAS_UTL.a 

check_breakers:        check_breakers.o  libAS_ALN.a libAS_MSG.a libAS_PER.a libAS_UTL.a libAS_UID.a

install: 

regression:
	@echo "NO REGRESSION"
tester:
	@echo "NO TESTER"




