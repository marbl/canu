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
# gmake Makefile for AS_GKP

LOCAL_WORK = $(shell cd ../..; pwd)


GKPSOURCES = AS_GKP_main.c AS_GKP_check.c AS_GKP_checkFrag.c  AS_GKP_checkBac.c
GKPOBJECTS = AS_GKP_main.o AS_GKP_check.o AS_GKP_checkFrag.o AS_GKP_checkBac.o
OTHERSOURCES = dumpGatekeeperRecords.c dumpGatekeeper.c \
	GateKeeperNextName.c GateKeeperCurrName.c debugGatekeeperStore.c
OTHEROBJECTS = dumpGatekeeperRecords.o dumpGatekeeper.o \
	GateKeeperNextName.o GateKeeperCurrName.o debugGatekeeperStore.o

SOURCES = $(GKPSOURCES) $(OTHERSOURCES)
OBJECTS = $(GKPOBJECTS) $(OTHEROBJECTS)

LIBRARIES =

PROGS = dumpGatekeeperRecords gatekeeper dumpGatekeeper \
	GateKeeperNextName  GateKeeperCurrName  debugGatekeeperStore

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(PROGS)

install: all
	echo "Nothing to be done for AS_GKP"

gatekeeper: $(GKPOBJECTS) \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 

debugGatekeeperStore: debugGatekeeperStore.o \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 

dumpGatekeeper: dumpGatekeeper.o \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 

dumpGatekeeperRecords: dumpGatekeeperRecords.o \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 

GateKeeperNextName: GateKeeperNextName.o \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 

GateKeeperCurrName: GateKeeperCurrName.o \
	libAS_PER.a libAS_UTL.a libAS_MSG.a 
