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

AS_FGB_INCS = AS_FGB_hanging_fragment.h
AS_OGB_INCS = AS_CGB_all.h AS_CGB_histo.h \
AS_CGB_fga.h AS_CGB_fgb.h AS_CGB_io.h AS_CGB_store.h \
AS_CGB_count_fragment_and_edge_labels.h \
AS_CGB_traversal.h AS_CGB_walk.h AS_CGB_edgemate.h
# AS_FGB_buddy.h
AS_CGB_INCS = \
AS_CGB_cgb.h AS_CGB_chimeras.h AS_CGB_classify.h \
AS_FGB_FragmentHash.h AS_FGB_buildFragmentHash.h \
AS_CGB_miniunitigger.h AS_CGB_unitigger.h
# AS_CGB_dfs.h
# AS_CGB_edgetrimmer.h
# AS_CGB_branchpts.h

AS_CGA_INCS = AS_CGB_cga.h

AS_FGB_SRCS = AS_FGB_main.c AS_FGB_io.c 

AS_OGB_SRCS = \
AS_CGB_histo.c \
AS_CGB_fga.c AS_CGB_fgb.c AS_FGB_hanging_fragment.c AS_FGB_contained.c \
AS_CGB_store.c AS_FGB_FragmentHash.c AS_FGB_buildFragmentHash.c  \
AS_CGB_traversal.c AS_CGB_walk.c \
AS_CGB_count_fragment_and_edge_labels.c \
AS_CGB_edgemate.c
# AS_FGB_buddy.c

# AS_CGB_unitigger.c 
AS_CGB_SRCS = \
AS_CGB_miniunitigger.c AS_CGB_main.c AS_CGB_io.c \
AS_CGB_cgb.c AS_CGB_classify.c AS_CGB_chimeras.c
# AS_CGB_branchpts.c
# AS_CGB_dfs.c
# AS_CGB_edgetrimmer.c

AS_CGA_SRCS = AS_CGB_cga.c

AS_CGB_BUB_INCS = AS_CGB_Bubble.h  \
                  AS_CGB_Bubble_Graph.h AS_CGB_Bubble_VertexSet.h \
                  AS_CGB_Bubble_GraphMethods.h AS_CGB_Bubble_Popper.h \
	          AS_CGB_Bubble_PopperMethods.h AS_CGB_Bubble_log.h
AS_CGB_BUB_SRCS = AS_CGB_Bubble.c AS_CGB_Bubble_Graph.c \
                  AS_CGB_Bubble_VertexSet.c AS_CGB_Bubble_top.c \
	          AS_CGB_Bubble_dfs.c AS_CGB_Bubble_Popper.c \
	          AS_CGB_Bubble_PopperMethods.c

UNITIG_SRCS = $(AS_OGB_SRCS) $(AS_CGB_SRCS) $(AS_CGA_SRCS) $(AS_CGB_BUB_SRCS)
VIEWER_SRCS = aug-cam.c delcher.c
UTILITY1_SRCS = strip-out.c
UTILITY2_SRCS = get-iid-from-frag-src.c
UTILITY3_SRCS = cgi2fasta.c
UTILITY4_SRCS = AS_CGB_util.c
FOM2UOM_SRCS = AS_CGB_fom2uom.c
IID2CID_SRCS = AS_CGB_iid2cid.c

AS_BRK_INCS = AS_CGB_breakers.h
AS_BRK_SRCS = AS_CGB_breakers.c
CHKBRK_SRCS = check_breakers.c
# RPRBRKS_SRCS = AS_CGB_repair_breakers.c 
#NO_OVL_SRCS = AS_CGB_no_overlaps.c

AS_FGB_OBJS = $(AS_FGB_SRCS:.c=.o)
AS_OGB_OBJS = $(AS_OGB_SRCS:.c=.o)
AS_CGB_OBJS = $(AS_CGB_SRCS:.c=.o)
AS_CGA_OBJS = $(AS_CGA_SRCS:.c=.o)
AS_CGB_BUB_OBJS = $(AS_CGB_BUB_SRCS:.c=.o)
UNITIG_OBJS = $(AS_OGB_OBJS) $(AS_CGB_OBJS) $(AS_CGA_OBJS) $(AS_CGB_BUB_OBJS)
VIEWER_OBJS = $(VIEWER_SRCS:.c=.o)
UTILITY1_OBJS = $(UTILITY1_SRCS:.c=.o)
UTILITY2_OBJS = $(UTILITY2_SRCS:.c=.o)
UTILITY3_OBJS = $(UTILITY3_SRCS:.c=.o)
UTILITY4_OBJS = $(UTILITY4_SRCS:.c=.o)
IID2CID_OBJS = $(IID2CID_SRCS:.c=.o)
FOM2UOM_OBJS = $(FOM2UOM_SRCS:.c=.o)

AS_BRK_OBJS = $(AS_BRK_SRCS:.c=.o)
CHKBRK_OBJS = $(CHKBRK_SRCS:.c=.o)
# RPRBRKS_OBJS = $(RPRBRKS_SRCS:.c=.o)
#NO_OVL_OBJS = $(NO_OVL_SRCS:.c=.o)

LIB_SOURCES = $(UNITIG_SRCS) $(AS_FGB_SRCS) $(AS_BRK_SRCS) $(UTILITY4_SRCS)
# $(RPRBRKS_SRCS) 
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)


SOURCES = $(VIEWER_SRCS) $(UTILITY1_SRCS) $(UTILITY2_SRCS) $(UTILITY3_SRCS) $(IID2CID_SRCS) $(FOM2UOM_SRCS) $(CHKBRK_SRCS) make_OFG_from_FragStore.c miniunitigger_test.c AS_CGB_unitigger.c WriteUnitigsFromCheckPoint.c $(LIB_SOURCES)
# $(NO_OVL_SRCS)
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_CGB.a


PROGS = unitigger aug-cam strip-out AS_CGB_fom2uom get-iid-from-frag-src cgi2fasta AS_CGB_iid2cid check_breakers make_OFG_from_FragStore miniunitigger_test WriteUnitigsFromCheckPoint
# repair_breakers 
# no_overlaps

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as
#CFLAGS += -g -fno-inline
#CFLAGS += -g

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)
# lint -std -I../../inc -I../AS_MSG AS_CGB_cgb.c

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

libAS_CGB.a: $(LIB_OBJECTS)

WriteUnitigsFromCheckPoint : WriteUnitigsFromCheckPoint.o ${LIB_OBJECTS} libAS_OVL.a libAS_ALN.a \
 libAS_MSG.a libAS_PER.a libAS_UTL.a libm.a

make_OFG_from_FragStore: make_OFG_from_FragStore.o libAS_PER.a libAS_UTL.a libm.a

miniunitigger_test: miniunitigger_test.o $(LIB_OBJECTS) libAS_OVL.a libAS_PER.a libAS_ALN.a libAS_MSG.a libAS_UTL.a libm.a

unitigger: AS_CGB_unitigger.o $(LIB_OBJECTS) libAS_OVL.a libAS_PER.a libAS_ALN.a libAS_MSG.a libAS_UTL.a libm.a


# fgb: $(AS_FGB_OBJS) $(AS_OGB_OBJS) $(UTILITY4_OBJS) libAS_OVL.a libAS_PER.a libAS_MSG.a libAS_UTL.a libm.a

# cgb: $(UNITIG_OBJS) libAS_PER.a libAS_ALN.a libAS_MSG.a libAS_UTL.a libm.a
# libAS_CNS.a

aug-cam: $(VIEWER_OBJS) libm.a

strip-out: $(UTILITY1_OBJS) libAS_MSG.a libAS_UTL.a

get-iid-from-frag-src: $(UTILITY2_OBJS) AS_CGB_histo.o libAS_MSG.a libAS_UTL.a libm.a

cgi2fasta:  $(UTILITY3_OBJS) AS_CGB_histo.o libAS_MSG.a libAS_UTL.a libm.a

AS_CGB_fom2uom: $(FOM2UOM_OBJS) AS_CGB_histo.o libAS_MSG.a libAS_UTL.a libm.a

AS_CGB_iid2cid: $(IID2CID_OBJS) AS_CGB_histo.o libAS_MSG.a libAS_UTL.a libm.a

check_breakers: $(CHKBRK_OBJS) libAS_ALN.a libAS_MSG.a libAS_PER.a libAS_UTL.a libm.a

#repair_breakers: $(RPRBRKS_OBJS) $(AS_BRK_OBJS) $(UTILITY4_OBJS) AS_FGB_io.o $(AS_OGB_OBJS) libAS_OVL.a libAS_ALN.a libAS_PER.a libAS_MSG.a libAS_UTL.a libm.a

#no_overlaps: $(NO_OVL_OBJS) $(AS_OGB_OBJS) libAS_ALN.a libAS_PER.a libAS_MSG.a libAS_UTL.a libm.a

install: 
	cp -f display-repeats show-feature $(LOCAL_WORK)/bin
regression:
	@echo "NO REGRESSION"
tester:
	@echo "NO TESTER"




