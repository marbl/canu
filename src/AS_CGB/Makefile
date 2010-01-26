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

AS_FGB_SRCS = AS_FGB_main.C AS_FGB_io.C
AS_OGB_SRCS = AS_CGB_histo.C \
              AS_CGB_fga.C AS_CGB_fgb.C AS_FGB_hanging_fragment.C AS_FGB_contained.C \
              AS_CGB_traversal.C AS_CGB_walk.C \
              AS_CGB_edgemate.C
AS_CGB_SRCS = AS_CGB_main.C AS_CGB_cgb.C AS_CGB_classify.C AS_CGB_chimeras.C
AS_CGA_SRCS = AS_CGB_cga.C
AS_BUB_SRCS = AS_CGB_Bubble.C AS_CGB_Bubble_Graph.C \
              AS_CGB_Bubble_VertexSet.C \
              AS_CGB_Bubble_dfs.C AS_CGB_Bubble_Popper.C \
              AS_CGB_Bubble_PopperMethods.C
AS_EDT_OBJECTS = $(AS_EDT_SRC:.C=.o)

AS_UTG_SRC = AS_CGB_unitigger.C $(AS_FGB_SRCS) $(AS_OGB_SRCS) $(AS_CGB_SRCS) $(AS_CGA_SRCS) $(AS_BUB_SRCS)
AS_UTG_OBJECTS = $(AS_UTG_SRC:.C=.o)

SOURCES = $(AS_UTG_SRC) $(AS_EDT_SRC)
OBJECTS = $(AS_UTG_SRC:.C=.o)

#  This only for BOG
LIB_OBJECTS   = AS_CGB_histo.o

LIBRARIES = libAS_CGB.a libCA.a

CXX_PROGS = unitigger

include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.C *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(CXX_PROGS)

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

libAS_CGB.a:     $(LIB_OBJECTS)
libCA.a:         $(LIB_OBJECTS)

unitigger:       $(AS_UTG_OBJECTS) libCA.a

