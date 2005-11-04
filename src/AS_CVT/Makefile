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
# Makefile for AS_CVT

LOCAL_WORK = $(shell cd ../..; pwd)

SOURCES = AS_CVT_hashtable.c AS_CVT_uid_interface.c \
          extract_component.c make_cns_list.c make_range_file.c \
          extract_from_uids.c extract_multiple.c extract_fraction.c \
	  ASMData.c asm2asmStore.c dumpClones.c \
	  getMDIDataFromASM.c dumpFragScaffoldCoordinates.c loadAsmStore.c \
	  dumpContigLengths.c dumpSurrogateFrags.c dumpMultiInstanceFrags.c \
	  dumpSurrogateCoords.c dumpSurrogateSeqCoords.c \
	  dumpUnreferencedFrags.c dumpMissingFrags.c dumpDeletedFrags.c \
	  dumpScaffoldLengths.c dumpBadMPs.c buildMapStore.c \
	  dumpMapStore.c dumpLibs.c dumpMappedMatePairs.c \
	  dumpMappedMatePairsBinary.c scf2chrSimple.c \
	  dumpMappedElsewheres.c dumpMappedBACFragIntervals.c \
	  concatRegions.c dumpElsewheres.c dumpForTampa.c \
          generateReadsFiles.c
OBJECTS = $(SOURCES:.c=.o)

LIBRARIES = 

PROGS = extract_component make_cns_list make_range_file \
	extract_from_uids extract_multiple extract_fraction \
	asm2asmStore dumpClones \
	getMDIDataFromASM dumpFragScaffCoords loadAsmStore \
	dumpContigLengths dumpSurrogateFrags dumpMultiInstanceFrags \
	dumpSurrogateCoords dumpSurrogateSeqCoords \
	dumpUnreferencedFrags dumpMissingFrags dumpDeletedFrags \
	dumpScaffoldLengths dumpBadMPs buildMapStore \
	dumpMapStore dumpLibs dumpMappedMatePairs \
	dumpMappedMatePairsBinary scf2chrSimple \
	dumpMappedElsewheres dumpMappedBACFragIntervals \
	concatRegions dumpElsewheres dumpForTampa \
        generateReadsFiles


SCRIPTS = run_components_thru_grande \
          run_grande_on_component \
          run_grande_on_frgs \
          setup_status_file \
	  mdi2ps.sh

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

DEPEND_FILES = *.c *.h
CLEANABLE_FILES = *.o *~

INST_SET= $(PROGS)

all: $(PROGS) $(LIBRARIES) $(OBJECTS) $(SCRIPTS)
install: $(PROGS) $(LIBRARIES) $(OBJECTS) $(SCRIPTS)

scf2chrSimple: scf2chrSimple.o libAS_UTL.a

concatRegions: concatRegions.o

dumpMappedBACFragIntervals: dumpMappedBACFragIntervals.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMappedMatePairsBinary: dumpMappedMatePairsBinary.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMappedMatePairs: dumpMappedMatePairs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpElsewheres: dumpElsewheres.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpForTampa: dumpForTampa.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

generateReadsFiles: generateReadsFiles.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMappedElsewheres: dumpMappedElsewheres.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpBadMPs: dumpBadMPs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpLibs: dumpLibs.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMapStore: dumpMapStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

buildMapStore: buildMapStore.o ASMData.o \
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

loadAsmStore: loadAsmStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

asm2asmStore: asm2asmStore.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpClones: dumpClones.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpFragScaffCoords: dumpFragScaffoldCoordinates.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

getMDIDataFromASM: getMDIDataFromASM.o libAS_MSG.a libAS_UTL.a

dumpContigLengths: dumpContigLengths.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpScaffoldLengths: dumpScaffoldLengths.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpSurrogateFrags: dumpSurrogateFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

dumpMultiInstanceFrags: dumpMultiInstanceFrags.o ASMData.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

extract_fraction: extract_fraction.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

extract_multiple: extract_multiple.o \
	AS_CVT_hashtable.o AS_CVT_uid_interface.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

extract_component: extract_component.o \
	AS_CVT_hashtable.o AS_CVT_uid_interface.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

extract_from_uids: extract_from_uids.o \
	AS_CVT_hashtable.o AS_CVT_uid_interface.o \
	libAS_MSG.a libAS_PER.a libAS_UTL.a

make_cns_list: make_cns_list.o libAS_MSG.a

make_range_file: make_range_file.o libAS_UTL.a

run_components_thru_grande:
	cp -f run_components_thru_grande.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/run_components_thru_grande.sh

run_grande_on_component:
	cp -f run_grande_on_component.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/run_grande_on_component.sh

run_grande_on_frgs:
	cp -f run_grande_on_frgs.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/run_grande_on_frgs.sh

setup_status_file:
	cp -f setup_status_file.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/setup_status_file.sh

mdi2ps.sh:
	cp -f mdi2ps.sh $(LOCAL_BIN)
	chmod 775 $(LOCAL_BIN)/mdi2ps.sh

