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

#
#  These executables are for interacting with AS_PER_asmStore.
#  Please do not add unrelated code here!  Thanks.
#

SOURCES = ASMData.c \
          asm2asmStore.c \
          buildMapStore.c \
          dumpClones.c \
          dumpFragScaffoldCoordinates.c \
	  dumpContigLengths.c \
          dumpSurrogateFrags.c \
          dumpMultiInstanceFrags.c \
	  dumpSurrogateCoords.c \
          dumpSurrogateSeqCoords.c \
	  dumpUnreferencedFrags.c \
          dumpMissingFrags.c \
          dumpDeletedFrags.c \
	  dumpScaffoldLengths.c \
          dumpBadMPs.c \
	  dumpMapStore.c \
          dumpLibs.c \
          dumpMappedMatePairs.c \
	  dumpMappedMatePairsBinary.c \
	  dumpMappedElsewheres.c \
          dumpElsewheres.c \
          dumpForTampa.c \
          generateReadsFiles.c
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = 

PROGS = asm2asmStore \
        buildMapStore \
        dumpClones \
	dumpFragScaffCoords \
	dumpContigLengths \
        dumpSurrogateFrags \
        dumpMultiInstanceFrags \
	dumpSurrogateCoords \
        dumpSurrogateSeqCoords \
	dumpUnreferencedFrags \
        dumpMissingFrags \
        dumpDeletedFrags \
	dumpScaffoldLengths \
        dumpBadMPs \
	dumpMapStore \
        dumpLibs \
        dumpMappedMatePairs \
	dumpMappedMatePairsBinary \
	dumpMappedElsewheres \
        dumpElsewheres \
        dumpForTampa \
        generateReadsFiles

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

all: $(PROGS) $(LIBRARIES) $(OBJECTS) $(SCRIPTS)
	@test -n nop


asm2asmStore: asm2asmStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

buildMapStore: buildMapStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a


dumpMappedMatePairsBinary: dumpMappedMatePairsBinary.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMappedMatePairs: dumpMappedMatePairs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpElsewheres: dumpElsewheres.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpForTampa: dumpForTampa.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMappedElsewheres: dumpMappedElsewheres.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpBadMPs: dumpBadMPs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpLibs: dumpLibs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMapStore: dumpMapStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpDeletedFrags: dumpDeletedFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMissingFrags: dumpMissingFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpUnreferencedFrags: dumpUnreferencedFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpSurrogateCoords: dumpSurrogateCoords.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpSurrogateSeqCoords: dumpSurrogateSeqCoords.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpClones: dumpClones.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpFragScaffCoords: dumpFragScaffoldCoordinates.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpContigLengths: dumpContigLengths.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpScaffoldLengths: dumpScaffoldLengths.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpSurrogateFrags: dumpSurrogateFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMultiInstanceFrags: dumpMultiInstanceFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a


generateReadsFiles: generateReadsFiles.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a
