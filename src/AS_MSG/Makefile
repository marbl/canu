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

LIBRARIES = libAS_MSG.a libCA.a

SOURCES = AS_MSG_pmesg.c  \
          proto_append.c \
          remove_fragment.c \
          AS_MSG_pmesg_test.c \
          ExtractMessages.c

OBJECTS = $(SOURCES:.c=.o)

PROGS = AS_MSG_pmesg_test \
	proto_append \
        remove_fragment \
	extractmessages

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS)

install: all
	echo "Nothing to be done for AS_MSG"

default: clean all regression tester install 

libAS_MSG.a:   AS_MSG_pmesg.o

libCA.a:       AS_MSG_pmesg.o

AS_MSG_pmesg_test: AS_MSG_pmesg_test.o AS_MSG_pmesg.o libAS_UTL.a

proto_append: proto_append.o AS_MSG_pmesg.o libAS_UTL.a

remove_fragment: remove_fragment.o AS_MSG_pmesg.o libAS_UTL.a

extractmessages: ExtractMessages.o AS_MSG_pmesg.o libAS_UTL.a

tester: AS_MSG_pmesg_test AS_MSG_testfile.in
	AS_MSG_pmesg_test -d < AS_MSG_testfile.in > testfile.out
	diff -w -n AS_MSG_testfile.in testfile.out
