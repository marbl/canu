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


LOCAL_WORK = $(shell cd ../../; pwd)

include $(LOCAL_WORK)/src/site_name.as

ANALYZE_SOURCES = analyzePosMap.C \
                  analyzePosMap-libraryFate.C \
                  analyzePosMap-load.C \
                  analyzePosMap-gapFillProbability.C
ANALYZE_OBJECTS = $(ANALYZE_SOURCES:.C=.o)

SOURCES   = terminator.C asmOutputFasta.C asmOutputStatistics.C fragsInVars.C buildPosMap.C $(ANALYZE_SOURCES) demotePosMap.C
OBJECTS   = (SOURCES:.C=.o)
CXX_PROGS = terminator asmOutputFasta asmOutputStatistics fragsInVars buildPosMap analyzePosMap demotePosMap

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

terminator:              terminator.o               libCA.a
asmOutputFasta:          asmOutputFasta.o           libCA.a
asmOutputStatistics:     asmOutputStatistics.o      libCA.a
fragsInVars:             fragsInVars.o              libCA.a
buildPosMap:             buildPosMap.o              libCA.a
analyzePosMap:           $(ANALYZE_OBJECTS)         libCA.a
demotePosMap:            demotePosMap.o             libCA.a