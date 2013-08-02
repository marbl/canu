#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2005-2007, J. Craig Venter Institute.
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
LOCAL_WORK = $(shell cd ../..; pwd)

#
# gmake Makefile for AS_OBT, Overlap Based Trimming

#  FIGARO should point to a configured, compiled and installed copy of
#  http://amos.sourceforge.net/Figaro/Figaro.html
#
#  We attempt to find it in the same place that our src/ directory is.


#  We'd like to include the global include (c_make.as) to get these
#  definitions, but we cannot without defining rules.

OSTYPE      = $(shell echo `uname`)
MACHINETYPE = $(shell echo `uname -m`)

ifeq ($(MACHINETYPE), x86_64)
  MACHINETYPE = amd64
endif
ifeq ($(MACHINETYPE), Power Macintosh)
  MACHINETYPE = ppc
endif

#  This isn't perfect; if we're building debug here, we _usually_ want
#  to use the debug figaro.  But checking if figaro exists, and failing if
#  we don't find the debug figaro, is a little tricky.
ifneq "$(origin FIGARO)" "environment"

ifeq ($(shell ls -d $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE)-debug 2> /dev/null), $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE)-debug)
  FIGARO = $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE)-debug
endif

ifeq ($(shell ls -d $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE) 2> /dev/null), $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE))
  FIGARO = $(LOCAL_WORK)/figaro/$(OSTYPE)-$(MACHINETYPE)
endif

endif

EXE_SOURCES = initialTrim.C \
              deduplicate.C \
              finalTrim.C \
              finalTrim-largestCovered.C \
              finalTrim-evidenceBased.C \
              finalTrim-bestEdge.C \
              chimera.C \
              trim.C \
              removeMateOverlap.C \

LIB_SOURCES = AS_OBT_acceptableOverlap.C \
              AS_OBT_overlaps.C

SOURCES     = $(EXE_SOURCES) $(LIB_SOURCES)
OBJECTS     = $(SOURCES:.C=.o) $(LIB_SOURCES:.C=.o)
LIBRARIES   = libAS_OBT.a libCA.a

CXX_PROGS   = initialTrim \
              deduplicate \
              finalTrim \
              chimera \
              removeMateOverlap
SCRIPTS =
ifdef FIGARO
  SCRIPTS += figaro
endif

include $(LOCAL_WORK)/src/c_make.as

all: $(CXX_PROGS) $(SCRIPTS)
	@test -n nop

libAS_OBT.a: $(LIB_SOURCES:.C=.o)
libCA.a:     $(LIB_SOURCES:.C=.o)

initialTrim:          initialTrim.o trim.o libCA.a
deduplicate:          deduplicate.o libCA.a
finalTrim:            finalTrim.o trim.o finalTrim-largestCovered.o finalTrim-evidenceBased.o finalTrim-bestEdge.o libCA.a
chimera:              chimera.o libCA.a
removeMateOverlap:    removeMateOverlap.o libCA.a

figaro: $(FIGARO)/figaro
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

