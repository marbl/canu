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

#  KMER should point to a configured and compiled copy of
#  http://sourceforge.net/projects/kmer/

ifdef KMER
  MERYLSOURCE      = meryl.C
  MERYLOBJECTS     = meryl.o
  SOURCES          = $(MERYLSOURCE) mercy.C mercy-regions.C
  OBJECTS          = $(SOURCES:.C=.o)
  PROGS            = meryl mercy
  INC_IMPORT_DIRS += $(KMER)/meryl $(KMER)/libutil $(KMER)/libmeryl $(KMER)/libbio $(KMER)/libkmer
  LIB_IMPORT_DIRS += $(KMER)/meryl $(KMER)/libutil $(KMER)/libmeryl $(KMER)/libbio $(KMER)/libkmer
  KMERLIBS         = libmerylguts.a libkmer.a libmeryl.a libbio.a libutil.a
else
  MERYLSOURCE      = AS_MER_meryl.cc
  MERYLOBJECTS     = AS_MER_meryl.o
  SOURCES          = $(MERYLSOURCE)
  OBJECTS          = $(SOURCES:.cc=.o)
  PROGS            = meryl
endif


include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(PROGS) 
	@test -n nop

meryl: $(MERYLOBJECTS) libAS_PER.a libAS_UTL.a $(KMERLIBS)

mercy: mercy.o libAS_PER.a libAS_UTL.a $(KMERLIBS)
