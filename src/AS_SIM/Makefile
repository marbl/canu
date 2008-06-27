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

LIBRARIES =
SOURCES = AS_SIM_frag.c \
          AS_SIM_poly.c \
          AS_SIM_massage.c \
          AS_SIM_outputADT.c \
          AS_SIM_outputBAC.c \
          AS_SIM_outputDST.c \
          AS_SIM_tobatches.c
OBJECTS = $(SOURCES:.c=.o)
PROGS = AS_SIM_massage \
        AS_SIM_frag \
        AS_SIM_poly \
        AS_SIM_outputADT \
        AS_SIM_outputDST \
        AS_SIM_outputBAC \
        AS_SIM_tobatches

# Include for AS project rules
LOCAL_WORK = $(shell cd ../..; pwd)
include $(LOCAL_WORK)/src/c_make.as

objects: $(OBJECTS)

all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@cp -f AS_SIM_celsim.pl $(LOCAL_BIN)/AS_SIM_celsim.pl
	@chmod ugo+x $(LOCAL_BIN)/AS_SIM_celsim.pl
	@cp -f AS_SIM_celsim.awk $(LOCAL_BIN)/AS_SIM_celsim.awk
	@chmod ugo+x $(LOCAL_BIN)/AS_SIM_celsim.awk
	@cp -f AS_SIM_labelfrag.pl $(LOCAL_BIN)/AS_SIM_labelfrag.pl
	@chmod ugo+x $(LOCAL_BIN)/AS_SIM_labelfrag.pl
	@cp -f AS_SIM_extractLength.pl $(LOCAL_BIN)/AS_SIM_extractLength.pl
	@chmod ugo+x $(LOCAL_BIN)/AS_SIM_extractLength.pl

AS_SIM_massage: AS_SIM_massage.o libAS_MSG.a libAS_UTL.a

AS_SIM_frag: AS_SIM_frag.o libAS_UTL.a

AS_SIM_poly: AS_SIM_poly.o

AS_SIM_outputADT: AS_SIM_outputADT.o libAS_MSG.a libAS_UTL.a

AS_SIM_outputDST: AS_SIM_outputDST.o libAS_MSG.a libAS_UTL.a

AS_SIM_outputBAC: AS_SIM_outputBAC.o libAS_MSG.a libAS_UTL.a

AS_SIM_tobatches: AS_SIM_tobatches.o libAS_MSG.a libAS_UTL.a
