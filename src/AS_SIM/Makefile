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
#VERBOSE=1

# gmake Makefile for AS_SIM

LIBRARIES =
SOURCES = AS_SIM_frag.c AS_SIM_poly.c AS_SIM_massage.c AS_SIM_outputADT.c AS_SIM_outputBAC.c AS_SIM_outputDST.c AS_SIM_tobatches.c
OBJECTS = $(SOURCES:.c=.o)
PROGS = massage frag poly outputADT outputDST tobatches outputBAC
SCRIPTS = AS_SIM_celsim.awk AS_SIM_labelfrag.pl AS_SIM_extractLength.pl celsim

# Include for AS project rules
LOCAL_WORK = $(shell cd ../..; pwd)
include $(LOCAL_WORK)/src/c_make.as

install: $(SCRIPTS) 
	cp -f $(SCRIPTS) $(LOCAL_BIN)

objects: $(OBJECTS)

all: $(OBJECTS) $(LIBRARIES) $(PROGS) install

massage: AS_SIM_massage.o libAS_MSG.a libAS_UTL.a libm.a

outputADT: AS_SIM_outputADT.o libAS_MSG.a libAS_UTL.a

outputDST: AS_SIM_outputDST.o libAS_MSG.a libAS_UTL.a

outputBAC: AS_SIM_outputBAC.o libAS_MSG.a libAS_UTL.a

frag: AS_SIM_frag.o libAS_UTL.a libm.a

poly: AS_SIM_poly.o libm.a

tobatches: AS_SIM_tobatches.o libAS_MSG.a libAS_UTL.a

celsim: AS_SIM_celsim
	cp -f AS_SIM_celsim celsim
	chmod ug+x celsim
