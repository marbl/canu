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

#  KMER should point to a configured, compiled and installed copy of
#  http://sourceforge.net/projects/kmer/
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
#  to use the debug kmer.  But checking if kmer exists, and failing if
#  we don't find the debug kmer, is a little tricky.

ifneq "$(origin KMER)" "environment"

ifeq ($(shell ls -d $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE)-debug 2> /dev/null), $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE)-debug)
  KMER = $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE)-debug
endif

ifeq ($(shell ls -d $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE) 2> /dev/null), $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE))
  KMER = $(LOCAL_WORK)/kmer/$(OSTYPE)-$(MACHINETYPE)
endif

endif

ifdef KMER
  MERYLSOURCE      = meryl.C
  MERYLOBJECTS     = meryl.o AS_MER_gkpStore_to_FastABase.o
  SOURCES          = $(MERYLSOURCE) mercy.C mercy-regions.C overmerry.C AS_MER_gkpStore_to_FastABase.C
  OBJECTS          = $(SOURCES:.C=.o)
  CXX_PROGS        = meryl mercy overmerry
  INC_IMPORT_DIRS += $(KMER)/include
  LIB_IMPORT_DIRS += $(KMER)/lib
  KMERLIBS         = libmerylguts.a libkmer.a libmeryl.a libseq.a libbio.a libutil.a
else
  MERYLSOURCE      = AS_MER_meryl.cc
  MERYLOBJECTS     = AS_MER_meryl.o
  SOURCES          = $(MERYLSOURCE)
  OBJECTS          = $(SOURCES:.cc=.o)
  CXX_PROGS        = meryl
endif


#  include again to get everything else updated
include $(LOCAL_WORK)/src/c_make.as


#  This must come after c_make.as, grrr.
ifdef KMER
  #  See also AS_OVL/Makefile
  ifeq ($(OSTYPE), Linux)
    LDFLAGS  += -lpthread
    CFLAGS   += -pthread
    CPPFLAGS += -pthread
  endif

  ifeq ($(OSTYPE), OSF1)
    LDFLAGS  += -lpthread
    CFLAGS   += -pthread
    CPPFLAGS += -pthread
  endif

  ifeq ($(OSTYPE), FreeBSD)
  endif

  ifeq ($(OSTYPE), Darwin)
  endif
endif

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

meryl:     $(MERYLOBJECTS)                            libCA.a $(KMERLIBS)
mercy:     mercy.o                                    libCA.a $(KMERLIBS)
overmerry: overmerry.o AS_MER_gkpStore_to_FastABase.o libCA.a $(KMERLIBS)
