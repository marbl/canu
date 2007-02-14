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

AS_OVL_CA_SRCS = AS_OVL_overlap_ca.c  AS_OVL_driver_ca.c 
AS_OVL_CA_OBJS = $(AS_OVL_CA_SRCS:.c=.o)

AS_OVL_SRCS = AS_OVL_overlap.c AS_OVL_driver.c
AS_OVL_OBJS = $(AS_OVL_SRCS:.c=.o)

AS_OVL_CMN_SRCS = AS_OVL_delcher.c 
AS_OVL_CMN_OBJS = $(AS_OVL_CMN_SRCS:.c=.o)

CORRECT_SRCS = FragCorrectOVL.c ShowCorrectsOVL.c CorrectOlapsOVL.c DeleteOlapsOVL.c FilterOlapsOVL.c CatCorrectsOVL.c CatEratesOVL.c
CORRECT_OBJS = $(CORRECT_SRCS:.c=.o)

LIBSOURCES =  OlapStoreOVL.c  MatchListOVL.c
LIBOBJECTS = $(LIBSOURCES:.c=.o) 

OLAPSTORE_SRCS = DumpOlapStoreOVL.c ForceEratesOVL.c GrowOlapStoreOVL.c UpdateEratesOVL.c
OLAPSTORE_OBJS = $(OLAPSTORE_SRCS:.c=.o)

SCREEN_SRCS = frag-anomaly.c AutoScreenOVL.c CombineMUMsOVL.c RemoveDupScreenOVL.c removedups.c
SCREEN_OBJS = $(SCREEN_SRCS:.c=.o)

TOOL_SRCS = get-degrees.c get-subgraph.c get-iidgraph.c get-uid-olaps.c fgbovlhits.c rev-olaps.c Hit_Limit_OVL.c
TOOL_OBJS = $(TOOL_SRCS:.c=.o)

CXX_SOURCES =  olapcoverage.cc
CXX_OBJECTSS = $(TOOL_CXX_SRCS:.cc=.o)

SOURCES = $(AS_OVL_SRCS) $(AS_OVL_CA_SRCS) $(AS_OVL_CMN_SRCS) $(CORRECT_SRCS) $(OLAPSTORE_SRCS) $(LIBSOURCES) LastFragInStoreOVL.c 
OBJECTS = $(SOURCES:.cc=.o)

AS_OVL_PROGS    = overlap overlap_ca
CORRECT_PROGS   = correct-frags show-corrects correct-olaps delete-olaps filter-olaps cat-corrects cat-erates
OLAPSTORE_PROGS = dump-olap-store force-erates grow-olap-store update-erates

OBSOLETE        = 
SCREEN_PROGS    = frag-anomaly auto-screen combine-mums remove-dup-screen removedups
TOOL_PROGS      = get-degrees get-subgraph get-iidgraph get-uid-olaps olapcoverage fgbovlhits rev-olaps hit-limit

PROGS     = $(AS_OVL_PROGS) $(OLAPSTORE_PROGS) $(CORRECT_PROGS) lastfraginstore
SCRIPTS   = 
LIBRARIES = libAS_OVL.a libCA.a

LIBS = libCA.a

include $(LOCAL_WORK)/src/c_make.as

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


all:    $(OBJECTS) $(LIBRARIES) $(PROGS)
	@test -n nop


libAS_OVL.a: $(LIBOBJECTS) AS_OVL_delcher.o

libCA.a: $(LIBOBJECTS) AS_OVL_delcher.o


auto-screen:  AutoScreenOVL.o $(LIBS)

cat-corrects:  CatCorrectsOVL.o AS_OVL_delcher.o $(LIBS)

cat-erates:  CatEratesOVL.o AS_OVL_delcher.o $(LIBS)

combine-mums:  CombineMUMsOVL.o AS_OVL_delcher.o $(LIBS)

correct-frags:  FragCorrectOVL.o AS_OVL_delcher.o $(LIBS) $(PTHREADS_LIBS)

correct-olaps:  CorrectOlapsOVL.o AS_OVL_delcher.o $(LIBS)

delete-olaps:  DeleteOlapsOVL.o AS_OVL_delcher.o $(LIBS)

dump-olap-store:  DumpOlapStoreOVL.o $(LIBS)

filter-olaps:  FilterOlapsOVL.o AS_OVL_delcher.o $(LIBS)

frag-anomaly:  frag-anomaly.o AS_OVL_delcher.o $(LIBS)

get-degrees:  get-degrees.o $(LIBS)

get-subgraph:  get-subgraph.o AS_OVL_delcher.o $(LIBS)

get-iidgraph:  get-iidgraph.o AS_OVL_delcher.o $(LIBS)

get-uid-olaps:  get-uid-olaps.o AS_OVL_delcher.o $(LIBS)

grow-olap-store:  GrowOlapStoreOVL.o $(LIBS)

hit-limit:  Hit_Limit_OVL.o

lastfraginstore:  LastFragInStoreOVL.o  $(LIBS)

olapcoverage:  olapcoverage.o

fgbovlhits:  fgbovlhits.o

remove-dup-screen:  RemoveDupScreenOVL.o AS_OVL_delcher.o $(LIBS)

removedups:  removedups.o  AS_OVL_delcher.o libAS_UTL.a

rev-olaps:  rev-olaps.o

show-corrects:  ShowCorrectsOVL.o AS_OVL_delcher.o $(LIBS)

update-erates:  UpdateEratesOVL.o AS_OVL_delcher.o $(LIBS)

force-erates:  ForceEratesOVL.o AS_OVL_delcher.o $(LIBS)

overlap: $(AS_OVL_OBJS) $(AS_OVL_CMN_OBJS) $(LIBS) $(PTHREADS_LIBS)

overlap_ca: $(AS_OVL_CA_OBJS) $(AS_OVL_CMN_OBJS) $(LIBS) $(PTHREADS_LIBS)
