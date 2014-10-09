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

MERYLSOURCE      = meryl.C
MERYLOBJECTS     = meryl.o AS_MER_gkpStore_to_FastABase.o AS_MER_gkpStoreChain.o
SOURCES          = $(MERYLSOURCE) mercy.C mercy-regions.C overmerry.C merTrim.C merTrimAdapter.C merTrimApply.C estimate-mer-threshold.C AS_MER_gkpStore_to_FastABase.C AS_MER_gkpStoreChain.C
OBJECTS          = $(SOURCES:.C=.o)
CXX_PROGS        = meryl mercy overmerry merTrim merTrimApply estimate-mer-threshold
INC_IMPORT_DIRS += $(KMER)/include
LIB_IMPORT_DIRS += $(KMER)/lib
KMERLIBS         = libmerylguts.a libkmer.a libmeryl.a libseq.a libbio.a libutil.a

#  include again to get everything else updated
include $(LOCAL_WORK)/src/c_make.as


#  This must come after c_make.as, grrr.
#
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


all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

meryl:                  $(MERYLOBJECTS)                                           libCA.a $(KMERLIBS)
mercy:                  mercy.o                                                   libCA.a $(KMERLIBS)
overmerry:              overmerry.o                AS_MER_gkpStore_to_FastABase.o libCA.a $(KMERLIBS)
merTrim:                merTrim.o merTrimAdapter.o AS_MER_gkpStore_to_FastABase.o libCA.a $(KMERLIBS)
merTrimApply:           merTrimApply.o                                            libCA.a $(KMERLIBS)
estimate-mer-threshold: estimate-mer-threshold.o   AS_MER_gkpStore_to_FastABase.o libCA.a $(KMERLIBS)
