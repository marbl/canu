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

AS_FGB_SRCS = AS_FGB_main.c AS_FGB_io.c 
AS_OGB_SRCS = AS_CGB_histo.c \
              AS_CGB_fga.c AS_CGB_fgb.c AS_FGB_hanging_fragment.c AS_FGB_contained.c \
              AS_CGB_traversal.c AS_CGB_walk.c \
              AS_CGB_edgemate.c
AS_CGB_SRCS = AS_CGB_main.c AS_CGB_cgb.c AS_CGB_classify.c AS_CGB_chimeras.c
AS_CGA_SRCS = AS_CGB_cga.c
AS_BUB_SRCS = AS_CGB_Bubble.c AS_CGB_Bubble_Graph.c \
              AS_CGB_Bubble_VertexSet.c \
              AS_CGB_Bubble_dfs.c AS_CGB_Bubble_Popper.c \
              AS_CGB_Bubble_PopperMethods.c
AS_EDT_SRC = adjustAStat.c
AS_EDT_OBJECTS = $(AS_EDT_SRC:.c=.o)

AS_UTG_SRC = AS_CGB_unitigger.c $(AS_FGB_SRCS) $(AS_OGB_SRCS) $(AS_CGB_SRCS) $(AS_CGA_SRCS) $(AS_BUB_SRCS)
AS_UTG_OBJECTS = $(AS_UTG_SRC:.c=.o)

SOURCES = $(AS_UTG_SRC) $(AS_EDT_SRC)
OBJECTS = $(AS_UTG_SRC:.c=.o)

#  This only for BOG
LIB_OBJECTS   = AS_CGB_histo.o

LIBRARIES = libAS_CGB.a libCA.a

PROGS = unitigger adjustAStat

include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@cp -f stripCGBReads.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/stripCGBReads.pl

libAS_CGB.a:     $(LIB_OBJECTS)
libCA.a:         $(LIB_OBJECTS)

unitigger:       $(AS_UTG_OBJECTS) libCA.a

adjustAStat:     $(AS_EDT_OBJECTS) libCA.a

