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

#  UMD_OVL should point to a configured, compiled and installed copy of
#  http://www.genome.umd.edu/overlapper.htm
#
#  We attempt to find it in the same place that our src/ directory is.
#
#  Currently, nothing is build in the UMD overlapper, it must be pre-compiled. All the necessary executables
#  are copied into the CA bin tree.
#

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

ifneq "$(origin UMD_OVL)" "environment"

ifeq ($(shell ls -d $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE)-debug 2> /dev/null), $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE)-debug)
  UMD_OVL = $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE)-debug
endif

ifeq ($(shell ls -d $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE) 2> /dev/null), $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE))
  UMD_OVL = $(LOCAL_WORK)/UMDOverlapper/$(OSTYPE)-$(MACHINETYPE)
endif

endif

ifdef UMD_OVL
  AS_SCRIPTS = runUMDOverlapper overlapPipeline dataWorkReduced
endif

AS_OVL_CA_SRCS = AS_OVL_overlap_ca.c  AS_OVL_driver_ca.c
AS_OVL_CA_OBJS = $(AS_OVL_CA_SRCS:.c=.o)

AS_OVL_SRCS = AS_OVL_overlap.c AS_OVL_driver.c
AS_OVL_OBJS = $(AS_OVL_SRCS:.c=.o)

AS_OVL_CMN_SRCS = AS_OVL_delcher.c
AS_OVL_CMN_OBJS = $(AS_OVL_CMN_SRCS:.c=.o)

SEED_OLAP_SRCS = OlapFromSeedsOVL.c
SEED_OLAP_OBJS = $(SEED_OLAP_SRCS:.c=.o)

CORRECT_SRCS = FragCorrectOVL.c ShowCorrectsOVL.c CorrectOlapsOVL.c CatCorrectsOVL.c CatEratesOVL.c
CORRECT_OBJS = $(CORRECT_SRCS:.c=.o)

LIBSOURCES =  SharedOVL.c
LIBOBJECTS = $(LIBSOURCES:.c=.o)

SOURCES = $(AS_OVL_SRCS) $(AS_OVL_CA_SRCS) $(AS_OVL_CMN_SRCS) $(SEED_OLAP_SRCS) $(CORRECT_SRCS) $(LIBSOURCES)
OBJECTS = $(SOURCES:.cc=.o)

AS_OVL_PROGS    = overlap overlap_ca
SEED_OLAP_PROGS = olap-from-seeds
CORRECT_PROGS   = correct-frags show-corrects correct-olaps cat-corrects cat-erates

PROGS      = $(AS_OVL_PROGS) $(SEED_OLAP_PROGS) $(CORRECT_PROGS)
SCRIPTS    = $(AS_SCRIPTS)
LIBRARIES  = libAS_OVL.a
LIBRARIES += libCA.a
LIBS       = libCA.a

include $(LOCAL_WORK)/src/c_make.as

#  See also AS_MER/Makefile
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


all:    $(OBJECTS) $(LIBRARIES) $(PROGS) $(SCRIPTS)
	@test -n nop

libAS_OVL.a: $(LIBOBJECTS) AS_OVL_delcher.o
libCA.a:     $(LIBOBJECTS) AS_OVL_delcher.o

show-corrects:  ShowCorrectsOVL.o AS_OVL_delcher.o $(LIBS)
update-erates:  UpdateEratesOVL.o AS_OVL_delcher.o $(LIBS)
force-erates:   ForceEratesOVL.o  AS_OVL_delcher.o $(LIBS)
cat-corrects:   CatCorrectsOVL.o  AS_OVL_delcher.o $(LIBS)
cat-erates:     CatEratesOVL.o    AS_OVL_delcher.o $(LIBS)
correct-frags:  FragCorrectOVL.o  AS_OVL_delcher.o $(LIBS)
correct-olaps:  CorrectOlapsOVL.o AS_OVL_delcher.o $(LIBS)

overlap:    $(AS_OVL_OBJS)    $(AS_OVL_CMN_OBJS) $(LIBS)
overlap_ca: $(AS_OVL_CA_OBJS) $(AS_OVL_CMN_OBJS) $(LIBS)

olap-from-seeds:  OlapFromSeedsOVL.o AS_OVL_delcher.o $(LIBS)


runUMDOverlapper: $(UMD_OVL)/runUMDOverlapper.perl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

overlapPipeline: $(UMD_OVL)/overlapPipeline
	@mkdir -p $(LOCAL_BIN)/overlapPipeline
	@chmod -f 775 $(LOCAL_BIN)/overlapPipeline
	@cp $</* $(LOCAL_BIN)/overlapPipeline
	@chmod 775 $(LOCAL_BIN)/overlapPipeline/

dataWorkReduced: $(UMD_OVL)/dataWorkReduced
	@mkdir -p $(LOCAL_BIN)/dataWorkReduced
	@chmod -f 775 $(LOCAL_BIN)/dataWorkReduced
	@cp $</* $(LOCAL_BIN)/dataWorkReduced
	@chmod 775 $(LOCAL_BIN)/dataWorkReduced/*
