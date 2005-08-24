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
# Makefile for alignment library

LOCAL_WORK = $(shell cd ../../; pwd)

LIB_SOURCES = AS_ALN_dpaligner.c AS_ALN_qvaligner.c \
              AS_ALN_sataligner.c AS_ALN_satfinder.c \
              AS_ALN_trim_vector.c \
	      AS_ALN_loverlapper.c AS_ALN_pieceOlap.c AS_ALN_forcns.c \
	      CA_ALN_local.c CA_ALN_overlap.c CA_ALN_scafcomp.c\
	      AS_ALN_boxfilling.c
	      
LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

TEST1_SOURCES = AS_ALN_testdriver1.c
TEST1_OBJECTS = $(TEST1_SOURCES:.c=.o)

TEST2_SOURCES = AS_ALN_testdriver2.c
TEST2_OBJECTS = $(TEST2_SOURCES:.c=.o)

TEST4_SOURCES = AS_ALN_testdriver4.c
TEST4_OBJECTS = $(TEST4_SOURCES:.c=.o)

TEST5_SOURCES = AS_ALN_testdriver5.c
TEST5_OBJECTS = $(TEST5_SOURCES:.c=.o)

SATTEST_SOURCES = AS_ALN_testsat.c
SATTEST_OBJECTS = $(SATTEST_SOURCES:.c=.o)

VT_TEST_SOURCES = AS_ALN_test_vt.c
VT_TEST_OBJECTS = $(VT_TEST_SOURCES:.c=.o)

RC_SOURCES = AS_ALN_revcompl.c
RC_OBJECTS = $(RC_SOURCES:.c=.o)

BOXFILL_SOURCES = AS_ALN_boxfiller.c
BOXFILL_OBJECTS = $(BOXFILL_SOURCES:.c=.o)

SOURCES = $(TEST1_SOURCES) $(TEST2_SOURCES)\
	  $(TEST4_SOURCES) $(TEST5_SOURCES)\
          $(SATTEST_SOURCES) $(VT_TEST_SOURCES) \
	  $(BOXFILL_SOURCES)\
          $(LIB_SOURCES) $(RC_SOURCES) 
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = libAS_ALN.a

PROGS = AS_ALN_testdriver1 AS_ALN_testdriver2 AS_ALN_testdriver4\
	AS_ALN_testdriver5 AS_ALN_testsat AS_ALN_test_vt AS_ALN_revcompl\
	AS_ALN_boxfiller

include $(LOCAL_WORK)/src/c_make.as
#CFLAGS = -g

default: all

all: $(OBJECTS) $(LIBRARIES) $(PROGS)

libAS_ALN.a: $(LIB_OBJECTS)

AS_ALN_testdriver1: $(TEST1_OBJECTS) $(LIBRARIES) libAS_MSG.a libAS_UTL.a libAS_PER.a libAS_OVL.a

AS_ALN_testdriver2: $(TEST2_OBJECTS) $(LIB_OBJECTS) libAS_UTL.a

AS_ALN_testdriver4: $(TEST4_OBJECTS) $(LIBRARIES) libAS_MSG.a libAS_UTL.a libAS_PER.a libAS_OVL.a

AS_ALN_testdriver5: $(TEST5_OBJECTS) $(LIBRARIES) libAS_MSG.a libAS_UTL.a libAS_PER.a libAS_OVL.a

AS_ALN_boxfiller: $(BOXFILL_OBJECTS) $(LIBRARIES) libAS_MSG.a libAS_UTL.a libAS_PER.a libAS_OVL.a

AS_ALN_revcompl:  $(RC_OBJECTS) $(LIBRARIES) libAS_UTL.a

AS_ALN_testsat: $(SATTEST_OBJECTS) $(LIB_OBJECTS) libAS_UTL.a

AS_ALN_test_vt: $(VT_TEST_OBJECTS) $(LIB_OBJECTS) libAS_MSG.a libAS_UTL.a

AS_ALN_tester:
	time AS_ALN_testdriver1 <dataset1 >out1
	time AS_ALN_testdriver1 <dataset2 >out2
	time AS_ALN_testdriver4 <dataset4 >out4
	time AS_ALN_testdriver5 <dataset5 >out5

package:
	shar *.c *.h >package
