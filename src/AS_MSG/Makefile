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

SOURCES    = AS_MSG_pmesg.c AS_MSG_pmesg1.c AS_MSG_pmesg2.c remove_fragment.c ExtractMessages.c
OBJECTS    = $(SOURCES:.c=.o)
PROGS      = AS_MSG_pmesg_test remove_fragment extractmessages
LIBRARIES  = libAS_MSG.a libCA.a
SCRIPTS    = convert-fasta-to-v2.pl convert-v1-to-v2.pl tracedb-to-frg.pl


# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS) $(SCRIPTS)
	@cp -f convert-fasta-to-v2.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/convert-fasta-to-v2.pl
	@cp convert-v1-to-v2.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/convert-v1-to-v2.pl
	@cp tracedb-to-frg.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/tracedb-to-frg.pl

libAS_MSG.a:       AS_MSG_pmesg.o AS_MSG_pmesg1.o AS_MSG_pmesg2.o
libCA.a:           AS_MSG_pmesg.o AS_MSG_pmesg1.o AS_MSG_pmesg2.o

remove_fragment:   remove_fragment.o   libCA.a
extractmessages:   ExtractMessages.o   libCA.a

perlmodule:
	cd p5-AS-MSG-Parser && perl ./Makefile.PL AS_BASE=$(LOCAL_WORK) INSTALL_BASE=$(LOCAL_OS)
	cd p5-AS-MSG-Parser && make install
