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
VERBOSE=1
# gmake Makefile for AS_MSG

LOCAL_WORK = $(shell cd ../..; pwd)

LIBRARIES = libAS_MSG.a

SOURCES = AS_MSG_pmesg.c AS_MSG_bmesg.c  \
          proto2binary.c binary2proto.c  proto2proto.c \
          proto_append.c binary_append.c remove_fragment.c \
          AS_MSG_pmesg_test.c CountMessages.c ExtractMessages.c \
	  track_fragment.c AS_MSG_Utility.c asm2bac.c
OBJECTS = $(SOURCES:.c=.o)

PROGS = AS_MSG_pmesg_test proto2binary binary2proto proto2proto \
	proto_append binary_append remove_fragment countmessages \
	extractmessages track_fragment asm2bac 

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)

install: all
	echo "Nothing to be done for AS_MSG"

Clean:

default: clean all regression tester install 

libAS_MSG.a:   AS_MSG_pmesg.o AS_MSG_bmesg.o AS_MSG_Utility.o

AS_MSG_pmesg_test: AS_MSG_pmesg_test.o AS_MSG_pmesg.o AS_MSG_bmesg.o \
	libAS_UTL.a

proto2binary: proto2binary.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

binary2proto: binary2proto.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

proto2proto: proto2proto.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

proto_append: proto_append.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

binary_append: binary_append.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

remove_fragment: remove_fragment.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

track_fragment: track_fragment.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

countmessages: CountMessages.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

extractmessages: ExtractMessages.o AS_MSG_bmesg.o AS_MSG_pmesg.o \
	libAS_UTL.a

extractxchunkmessages: ExtractXChunkMessages.o AS_MSG_bmesg.o AS_MSG_pmesg.o

asm2bac: asm2bac.o libAS_CNS.a libAS_REZ.a libAS_SDB.a libAS_MSG.a libAS_TER.a libAS_PER.a libAS_UTL.a libAS_ALN.a libAS_CGW.a libAS_LIN.a -lm

# asm_tables: asm_tables.o  AS_MSG_bmesg.o AS_MSG_pmesg.o


tester: AS_MSG_pmesg_test AS_MSG_testfile.in
	AS_MSG_pmesg_test -d < AS_MSG_testfile.in > testfile.out
	diff -w -n AS_MSG_testfile.in testfile.out
	proto2binary < AS_MSG_testfile.in | AS_MSG_pmesg_test -d >testfile.out
	diff -w -n AS_MSG_testfile.in testfile.out
	proto2binary < AS_MSG_testfile.in | binary2proto > testfile.out
	diff -w -n AS_MSG_testfile.in testfile.out
