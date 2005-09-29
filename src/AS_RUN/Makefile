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
# gmake Makefile for AS_RUN

LOCAL_WORK = $(shell cd ../..; pwd)

C_PROG_SOURCES = 
C_PROG_OBJECTS = $(C_PROG_SOURCES:.c=.o)

C_SOURCES = $(C_PROG_SOURCES)
C_OBJECTS = $(C_SOURCES:.c=.o)

CXX_SOURCES = 
CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)

SOURCES = $(C_SOURCES) $(CXX_SOURCES)
OBJECTS = $(C_OBJECTS) $(CXX_OBJECTS)

PROGS = 

SCRIPTS = runCA-OBT.pl

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

#  all: builds the usual binaries, then explicitly copies scripts to
#  the bin directory.
#	cp -f getCloneLengths.sh $(LOCAL_BIN)
#	chmod 775 $(LOCAL_BIN)/getCloneLengths.sh


all: $(OBJECTS) $(LIBRARIES) $(PROGS) $(SCRIPTS)

install: all

RUNCAOBTUTIL = runCA-OBT/runCAutil/applyOverlapCorrection.pl \
               runCA-OBT/runCAutil/checkOverlap.pl \
               runCA-OBT/runCAutil/checkPostUnitiggerConsensus.pl \
               runCA-OBT/runCAutil/createConsensusJobs.pl \
               runCA-OBT/runCAutil/createFragmentCorrectionJobs.pl \
               runCA-OBT/runCAutil/createOverlapCorrectionJobs.pl \
               runCA-OBT/runCAutil/createOverlapJobs.pl \
               runCA-OBT/runCAutil/createOverlapStore.pl \
               runCA-OBT/runCAutil/createPostUnitiggerConsensus.pl \
               runCA-OBT/runCAutil/mergeFragmentCorrection.pl \
               runCA-OBT/runCAutil/meryl.pl \
               runCA-OBT/runCAutil/overlapTrim.pl \
               runCA-OBT/runCAutil/preoverlap.pl \
               runCA-OBT/runCAutil/scaffolder.pl \
               runCA-OBT/runCAutil/terminate.pl \
               runCA-OBT/runCAutil/unitigger.pl \
               runCA-OBT/runCAutil/util.pl \
               runCA-OBT/runCAutil/scheduler.pm

runCA-OBT.pl: runCA-OBT/runCA.pl $(RUNCAOBTUTIL)
	cat runCA-OBT/runCA.pl $(RUNCAOBTUTIL) > $(LOCAL_BIN)/runCA-OBT.pl
	chmod 775 $(LOCAL_BIN)/runCA-OBT.pl
