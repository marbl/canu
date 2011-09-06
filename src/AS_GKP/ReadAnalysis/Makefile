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

LOCAL_WORK = $(shell cd ../../..; pwd)

SOURCES_LIB = \
	AlignmentError.C \
	AlignmentDataFilter.C \
	AlignmentDataReader.C \
	AlignmentDataStats.C \
	BaseAlignment.C \
	ReadAlignment.C \
	ReadAlignmentProfiler.C \
	ReadLengthAlignmentDataFilter.C \
	SnapperAlignmentDataReader.C

SOURCES_PROFILE_READ_ALIGNMENT = \
	profileReadAlignment.C

OBJECTS_LIB = $(SOURCES_LIB:.C=.o)
OBJECTS_PROFILE_READ_ALIGNMENT = $(SOURCES_PROFILE_READ_ALIGNMENT:.C=.o)

SOURCES = $(SOURCES_LIB) $(SOURCES_PROFILE_READ_ALIGNMENT)
OBJECTS = $(SOURCES:.C=.o)
LIBRARIES = libReadAnalysis.a libCA.a
CXX_PROGS = profileReadAlignment
SCRIPT_PROGS = plotReadAlignmentProfile.sh

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS) $(SCRIPT_PROGS)
	@test -n nop

libReadAnalysis.a: $(OBJECTS_LIB)
libCA.a: $(OBJECTS_LIB)

profileReadAlignment: $(OBJECTS_PROFILE_READ_ALIGNMENT) libCA.a
