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
# Tell c_make.gen to echo commands.

## Uncomment to enable
VERBOSE=1
ifdef $(VERBOSE)
export VERBOSE
endif

OSTYPE = $(shell echo `uname`)
export OSTYPE

MACHINETYPE = $(shell echo `uname -m`)
export MACHINETYPE

include $(LOCAL_WORK)/src/site_name.as

# AS Project Standards
include $(LOCAL_WORK)/src/c_make.gen

SUBDIRS = AS_MSG \
          AS_UTL \
          AS_PER \
          AS_ALN \
          AS_SIM \
          AS_OVL \
          AS_CNS \
          AS_CGB \
          AS_ORA \
          AS_GKP \
          AS_REZ \
          AS_LIN \
          AS_CGW \
          AS_VWR \
          AS_TER \
          AS_SDB \
          AS_CVT \
          AS_MER \
          AS_MPA \
          AS_BOG \
          AS_OBT \
          AS_UID \
          AS_RUN

# Compiler & flags

CC         = gcc
CXX        = g++
CFLAGS_OPT = -g
CXXDEFS    = -D__cplusplus
ARFLAGS    = rvs
USRLIB     = /usr/lib
LDFLAGS   += -O3

ifeq ($(OSTYPE), Linux)
  CC         = gcc
  CXX        = g++
  CFLAGS_OPT = -g
  CFLAGS    += -Os -DANSI_C -DX86_GCC_LINUX -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
  CXXFLAGS  += -Os -DANSI_C -DX86_GCC_LINUX -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
  USRLIB     = /usr/X11R6/lib
  CXXDEFS    = -D__cplusplus
  LDFLAGS   += --gc-sections

  ifeq ($(MACHINETYPE), i686)
    CFLAGS   += -march=i686
    CXXFLAGS += -march=i686
    LDFLAGS  += -march=i686
  endif
  ifeq ($(MACHINETYPE), x86_64)
    CC        = gcc
    CXX       = g++
    CFLAGS   += -m64 -mcmodel=medium 
    CXXFLAGS += -m64 -mcmodel=medium
    ifneq ($(SITE_NAME), TIGR)
      USRLIB    = /usr/lib64 /usr/X11R6/lib64
    endif
  endif
endif


#  Maddeningly, you absolutely must specify both LARGEFILE_SOURCE and
#  FILE_OFFSET_BITS=64 if your off_t is 8 bytes natively and your
#  off_t is 4 bytes.
#
ifeq ($(OSTYPE), FreeBSD)
  CC               = gcc
  CXX              = g++
  CFLAGS_OPT       = -g 
  CFLAGS          += -O3 -DNEEDXDRUHYPER -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 
  CXXDEFS          = -D__cplusplus
  INC_IMPORT_DIRS += /usr/local/include /usr/X11R6/include
  USRLIB          += /usr/local/lib /usr/X11R6/lib
  LDFLAGS         += -llthread -llgcc_r
endif

ifeq ($(OSTYPE), Darwin)
  CC         = gcc
  CXX        = g++
  CFLAGS_OPT = -O3 -g
  CFLAGS    += -DNEEDXDRUHYPER -D_THREAD_SAFE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 
  USRLIB    += /usr/local/lib /usr/X11R6/lib
endif

ifeq ($(OSTYPE), SunOS)
  ifeq ($(MACHINETYPE), i86pc)
    CC = gcc
    CXX = g++
    CFLAGS += -DBYTE_ORDER=LITTLE_ENDIAN
    LDFLAGS += -lm
  else
    CC = cc
    CXX = cxx
    CFLAGS += -DBYTE_ORDER=BIG_ENDIAN -DMAP_FILE=0 -g
    LDFLAGS += -lsocket -lnsl
  endif
endif

ifeq ($(OSTYPE), aix)
  CC = xlc
  CFLAGS_OPT = -O3
  # CFLAGS_OPT = -g
  # CFLAGS_OPT = -q64 -g -D_LARGE_FILES -qcpluscmt
  # THRINC  = -qthreaded -D_REENTRANT
  # THRLIB  = -lpthread
  # FLAGS  = -O3 -qarch=pwr3 -qtune=pwr3 -qthreaded -D_REENTRANT
  # FLAGS  = -O3 -qarch=pwr3 -qtune=pwr3
  # DEBUG  = -DDEBUG
  CFLAGS_OPT += -qflag=W:W \
             -qcpluscmt \
             -q64 \
             -qlonglong \
             -D_LARGE_FILES \
             -qarch=auto \
             -qtune=auto \
             -qcache=auto \
             -bhalt:8 \
             -qstaticinline \
             -qlanglvl=ext \
             -qignerrno \
             -qupconv \
             -qstrict
  CXX     = xlC_r $(CFLAGS_OPT)
  CXXDEFS =
  CXXFLAGS  = $(CFLAGS)
  # LDFLAGS = -g -q64 \
  LDFLAGS = -O3 -q64 \
          -qlonglong \
          -D_LARGE_FILES \
          -qarch=auto \
          -qtune=auto \
          -qcache=auto \
          -bhalt:8 \
          -qstaticinline \
          -qlanglvl=ext \
          -qignerrno \
          -qupconv \
          -qstrict \
          -lpthread
  # -bnoquiet -lpthread
  #AR = ar -X32
  #CPP = cc -E
  #CPPFLAGS  =  -I/scratch/hannone/incl -qcpluscmt -DAIX
  #CcOPTS  =  -O -qstrict -D_LARGE_FILES -qcpluscmt
  #ccWARNS = 
  #CC      = xlC
  #CCOPTS  = -q64 -O3 -qstrict -D_LARGE_FILES
  #CCWARNS =
endif

ifeq ($(OSTYPE), OSF1)
  CC = cc
  CXX = c++
  CFLAGS_OPT=  -pthread -fast -O4
  # To get accurate variable prints in the debuggers:
  # CFLAGS_OPT= -g 
  # CFLAGS_OPT= -g -check_bounds
  # Note that "-g" implies "-O0", which turns off all inlining.
  # Note that "-fast" implies "-inline speed"
  # CFLAGS_OPT= -g -O1 -inline manual
  # CFLAGS_OPT= -g -O1 -inline manual -check_bounds
  # CFLAGS_OPT= -O4 
  # CFLAGS_OPT= -fast -O4 -tune ev56 -arch ev56
  # CFLAGS_OPT= -fast -O4 -tune ev6 -arch ev6
  # -omp
  # -msg_enable c_to_cxx obsolescent
  # -check_bounds
  #
  CFLAGS_WARNINGS = -w0 -warnprotos -trapuv -float_const -readonly_strings \
    -msg_enable overflow \
    -msg_enable check \
    -msg_enable defunct \
    -msg_enable alignment \
    -msg_enable obsolescent \
    -msg_enable performance \
    -msg_enable preprocessor \
    -msg_enable nestedtype \
    -msg_disable valuepres \
    -msg_disable cxxcomment \
    -msg_disable c_to_cxx \
    -msg_disable unrefsdecl \
    -msg_disable strctpadding \
    -msg_disable truncintcast \
    -msg_disable unusedtop \
    -msg_disable unusedincl \
    -msg_disable unnecincl \
    -msg_disable nestincl \
    -msg_disable uncalled \
    -msg_disable ignorecallval \
    -msg_disable questcompare2 \
    -msg_disable unrefadecl \
    -msg_disable unrefsfunc \
    -msg_disable truncintasn \
    -msg_disable truncfltasn \
    -msg_disable truncfltint \
    -msg_disable extrasemi
  #  -msg_disable hexoctunsign
  #  -omp  -check_omp 
  #  -msg_enable unused 
  #  -msg_enable returnchecks
  #  -msg_enable  strctpadding
  #  -msg_enable level6
  INC_IMPORT_DIRS +=  /usr/local/include /usr/include
endif


CFLAGS          += $(CFLAGS_OPT) $(CFLAGS_WARNINGS)
INC_IMPORT_DIRS += $(patsubst %, $(LOCAL_WORK)/src/%, $(strip $(SUBDIRS))) \
                   $(LOCAL_WORK)/inc 
LIB_IMPORT_DIRS += $(LOCAL_LIB) $(USRLIB) 
OBJ_SEARCH_PATH = $(LOCAL_OBJ)

#/usr/shlib /usr/X11R6/lib /usr/X/lib /usr/shlib/X11



## Load if we are using SOAP or CURL as our UID transport
## This has to be external to this file so that AS_UID/Makefile will work
#include $(LOCAL_WORK)/src/AS_UID/uid_transport.as

ifeq ($(USE_SOAP_UID), 1)
  CFLAGS   += -DUSE_SOAP_UID
  CXXFLAGS += -DUSE_SOAP_UID
else
  LDFLAGS += -lcurl
endif

