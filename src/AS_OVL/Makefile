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
# gmake Makefile for AS_OVL

LOCAL_WORK = $(shell cd ../..; pwd)

AS_OVL_CA_SRCS = AS_OVL_overlap_ca.c  AS_OVL_driver_ca.c 
AS_OVL_CA_OBJS = $(AS_OVL_CA_SRCS:.c=.o)

AS_OVL_SRCS = AS_OVL_overlap.c AS_OVL_driver.c
AS_OVL_OBJS = $(AS_OVL_SRCS:.c=.o)

AS_OVL_CMN_SRCS = AS_OVL_delcher.c 
AS_OVL_CMN_OBJS = $(AS_OVL_CMN_SRCS:.c=.o)

CORRECT_SRCS = FragCorrectOVL.c ShowCorrectsOVL.c CorrectOlapsOVL.c DeleteOlapsOVL.c \
  FilterOlapsOVL.c CatCorrectsOVL.c CatEratesOVL.c
CORRECT_OBJS = $(CORRECT_SRCS:.c=.o)

LIBSOURCES =  OlapStoreOVL.c  MatchListOVL.c
LIBOBJECTS = $(LIBSOURCES:.c=.o) 

OLAPSTORE_SRCS = DumpOlapStoreOVL.c ForceEratesOVL.c GrowOlapStoreOVL.c UpdateEratesOVL.c
OLAPSTORE_OBJS = $(OLAPSTORE_SRCS:.c=.o)

SCREEN_SRCS = screen_analyze.c frag-anomaly.c screen-anomaly.c AutoScreenOVL.c \
  CombineMUMsOVL.c RemoveDupScreenOVL.c removedups.c
SCREEN_OBJS = $(SCREEN_SRCS:.c=.o)

SCRIPT_SRCS = CorrectScriptOVL.c MakeScriptOVL.c NewBactigScriptOVL.c NewFragScriptOVL.c
SCRIPT_OBJS = $(SCRIPT_SRCS:.c=.o)

TOOL_SRCS = get-degrees.c get-subgraph.c get-iidgraph.c get-uid-olaps.c \
  iid2fastaOVL.c isn2frgOVL.c LastFragInStoreOVL.c \
  fgbovlhits.c ModifyClearRangeOVL.c rev-olaps.c Hit_Limit_OVL.c
TOOL_OBJS = $(TOOL_SRCS:.c=.o)

CXX_SOURCES =  olapcoverage.cc
CXX_OBJECTSS = $(TOOL_CXX_SRCS:.cc=.o)

SOURCES = $(AS_OVL_CA_SRCS) $(AS_OVL_SRCS) $(AS_OVL_CMN_SRCS) $(CORRECT_SRCS) \
  $(SCREEN_SRCS) $(SCRIPT_SRCS) $(TOOL_SRCS) $(OLAPSTORE_SRCS) $(LIBSOURCES) $(CXX_SOURCES)
OBJECTS = $(AS_OVL_CA_OBJS) $(AS_OVL_OBJS) $(AS_OVL_CMN_OBJS) $(CORRECT_OBJS) \
  $(SCREEN_OBJS) $(SCRIPT_OBJS) $(TOOL_OBJS) $(OLAPSTORE_OBJS) $(CXX_OBJECTS)

LIBRARIES = libAS_OVL.a

AS_LIBRARIES = libCA.a
SYS_LIBRARIES = 

AS_OVL_PROGS = overlap overlap_ca
CORRECT_PROGS = correct-frags show-corrects correct-olaps delete-olaps filter-olaps \
  cat-corrects cat-erates
OLAPSTORE_PROGS = dump-olap-store force-erates grow-olap-store update-erates
SCREEN_PROGS = screen_analyze frag-anomaly screen-anomaly auto-screen combine-mums \
  remove-dup-screen removedups
SCRIPT_PROGS = correct-script make-ovl-script make-newbactigovl-script \
  make-newfragovl-script
TOOL_PROGS = get-degrees get-subgraph get-iidgraph get-uid-olaps iid2fasta \
  isn2frg lastfraginstore olapcoverage fgbovlhits modify-clear-range \
  rev-olaps hit-limit

PROGS = $(AS_OVL_PROGS) $(CORRECT_PROGS) $(SCREEN_PROGS) $(SCRIPT_PROGS) \
  $(TOOL_PROGS) $(OLAPSTORE_PROGS)

SCRIPTS = overlap_script.pl make-screen.csh


# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

# (MP) this should take care of intermittent "don't know how to make 
# libpthread.a" problems 

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
  #LDFLAGS  += -pthread
  #CFLAGS   += -pthread
  #CPPFLAGS += -pthread
endif

ifeq ($(OSTYPE), Darwin)
endif


#ifeq ($(OSTYPE), aix)
#  PTHREADS_LIBS=libc_r.so /usr/lib/libpthread.so
#  PTHREADS_LIBS=libc_r.a /usr/lib/libpthread.a
#endif

#ifeq ($(OSTYPE), Darwin)
#  PTHREADS_LIBS=libc_r.so /usr/lib/libpthread.dylib
#  PTHREADS_LIBS=libc_r.a /usr/lib/libpthread.dylib
#endif

#ifeq ($(OSTYPE), OSF1)
#  PTHREADS_LIBS=-lpthread
#  PTHREADS_LIBS=/usr/lib/libpthread.so
#endif

all:    $(OBJECTS) $(LIBRARIES) $(PROGS) install

install: $(SCRIPTS)
	chmod ugo+x $(SCRIPTS)
	cp -f $(SCRIPTS) $(LOCAL_BIN)

auto-screen:  AutoScreenOVL.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

cat-corrects:  CatCorrectsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

cat-erates:  CatEratesOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

combine-mums:  CombineMUMsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

correct-frags:  FragCorrectOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES) $(PTHREADS_LIBS)

correct-olaps:  CorrectOlapsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

correct-script:  CorrectScriptOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

delete-olaps:  DeleteOlapsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

dump-olap-store:  DumpOlapStoreOVL.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

filter-olaps:  FilterOlapsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

frag-anomaly:  frag-anomaly.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

get-degrees:  get-degrees.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

get-subgraph:  get-subgraph.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

get-iidgraph:  get-iidgraph.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

get-uid-olaps:  get-uid-olaps.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

grow-olap-store:  GrowOlapStoreOVL.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

hit-limit:  Hit_Limit_OVL.o $(SYS_LIBRARIES)

iid2fasta:  iid2fastaOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

isn2frg:  isn2frgOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

lastfraginstore:  LastFragInStoreOVL.o  $(AS_LIBRARIES) $(SYS_LIBRARIES)

libAS_OVL.a: $(LIBOBJECTS) AS_OVL_delcher.o

make-newbactigovl-script:  NewBactigScriptOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

make-newfragovl-script:  NewFragScriptOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

make-ovl-script:  MakeScriptOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

modify-clear-range:  ModifyClearRangeOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

olapcoverage:  olapcoverage.o

fgbovlhits:  fgbovlhits.o

remove-dup-screen:  RemoveDupScreenOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

removedups:  removedups.o  AS_OVL_delcher.o libAS_UTL.a

rev-olaps:  rev-olaps.o

screen_analyze:  screen_analyze.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

screen-anomaly:  screen-anomaly.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

show-corrects:  ShowCorrectsOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

update-erates:  UpdateEratesOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

force-erates:  ForceEratesOVL.o AS_OVL_delcher.o $(AS_LIBRARIES) $(SYS_LIBRARIES)

overlap: $(AS_OVL_OBJS) $(AS_OVL_CMN_OBJS) $(AS_LIBRARIES) $(SYS_LIBRARIES) $(PTHREADS_LIBS)

overlap_ca: $(AS_OVL_CA_OBJS) $(AS_OVL_CMN_OBJS) $(AS_LIBRARIES) $(SYS_LIBRARIES) $(PTHREADS_LIBS)

#build_hash: AS_OVL_FragHash.o
#	$(CC) $(LDFLAGS) -o $@ AS_OVL_FragHash.o \
#	../AS_MSG/AS_MSG_pmesg.o ../AS_UTL/AS_UTL_Hash.o $(LOADLIBES)
#
#



