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
# *********************************************************************
# $Id: Makefile,v 1.3 2005-03-22 19:04:58 jason_miller Exp $
# Module: AS_MER
# Description: The merCounter.
# Assumptions:
# *********************************************************************
#
#######################################################################
#VERBOSE=1
LOCAL_WORK = $(shell cd ../..; pwd)

AS_MER_GENERAL_SOURCES = \
	AS_MER_stream.cc

AS_MER_MAIN_SOURCES = \
	AS_MER_meryl.cc

CXX_SOURCES = $(AS_MER_GENERAL_SOURCES) $(AS_MER_MAIN_SOURCES)
CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)

TEST_SOURCES = 

SOURCES = $(AS_MER_GENERAL_SOURCES) $(AS_MER_MAIN_SOURCES) $(TEST_SOURCES)

AS_MER_GENERAL_OBJECTS = $(AS_MER_GENERAL_SOURCES:.cc=.o)
AS_MER_MAIN_OBJECTS = $(AS_MER_MAIN_SOURCES:.cc=.o)
TEST_OBJECTS = $(TEST_SOURCES:.cc=.o)
OBJECTS = $(SOURCES:.cc=.o)

LIBRARIES = 

PROGS = meryl 

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

# These variables influence cds/SEP/makefiles/c_make.gen
#VERBOSE = 1
#CFLAGS += -O
LDFLAGS += -lm

CXXDEFS = 
#CXXFLAGS += -g -pg -qcpluscmt -qstrict 
# DEC cxx requires tlocal for template linking. Not recognized by g++.
ifeq "$(OSTYPE)" "osf1"
CXXFLAGS += -tlocal
endif

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

# Careful! Library order matters.
STORE_LIBS = libAS_PER.a libAS_UTL.a 

all: $(OBJECTS) $(LIBRARIES) $(PROGS) 

meryl: AS_MER_meryl.o \
	$(AS_MER_GENERAL_OBJECTS) \
	$(STORE_LIBS)
