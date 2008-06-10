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

SOURCES = AS_TER_terminator.c asmOutputFasta.c fragsInVars.c buildPosMap.c
OBJECTS = (SOURCES:.c=.o)
PROGS   = terminator asmOutputFasta fragsInVars buildPosMap

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@cp -f buildFragContigPosMap.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/buildFragContigPosMap.pl

terminator:              AS_TER_terminator.o        libCA.a
asmOutputFasta:          asmOutputFasta.o           libCA.a
fragsInVars:             fragsInVars.o              libCA.a
buildPosMap:             buildPosMap.o              libCA.a