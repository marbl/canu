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

# Compiler & flags
#
# Please, please do not put defaults here.  They tend to break
# compiles on other architectures!
#
# The following sets compiler and linker flags based on the
# architecture detected.  To avoid name conflicts, we set:
#
# ARCH_CFLAGS    any flags for the compiler
# ARCH_LDFLAGS   flags for the linker
#
# ARCH_INC       any include directory paths
# ARCH_LIB       any library directory paths
#
# The final CFLAGS, CXXFLAGS and LDFLAGS are constructed from these.


# You can enable a debugging build, disabling all optimizations, by
# setting this to 1.
#
ifneq "$(origin BUILDDEBUG)" "environment"
BUILDDEBUG  = 0
endif

# You can enable a profiling build by setting this to 1.
#
ifneq "$(origin BUILDPROFILE)" "environment"
BUILDPROFILE  = 0
endif

# You can enable a line coverage build by setting this to 1.  This
# implies a debug build with no optimization.
#
ifneq "$(origin BUILDCOVERAGE)" "environment"
BUILDCOVERAGE  = 0
endif



OSTYPE      = $(shell echo `uname`)
MACHINETYPE = $(shell echo `uname -m`)

include $(LOCAL_WORK)/src/site_name.as
include $(LOCAL_WORK)/src/c_make.gen

ifeq ($(OSTYPE), Linux)
  CC         = gcc
  CXX        = g++

  ARCH_CFLAGS = -DANSI_C -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS  += -g -Wimplicit
    ARCH_LDFLAGS  = 
  else
    ARCH_CFLAGS  += -O2 -Wimplicit
    ARCH_LDFLAGS  = -Wl,-O1
  endif

  ARCH_LIB      = /usr/X11R6/lib

  ifeq ($(MACHINETYPE), i686)
    ARCH_CFLAGS   += -march=i686 -DX86_GCC_LINUX
    ARCH_LDFLAGS  += -march=i686
  endif
  ifeq ($(MACHINETYPE), x86_64)
    ARCH_CFLAGS   += -m64 -DX86_GCC_LINUX

    #  JCVI's opteron has 32-bit libraries in /usr/lib and /usr/X11R6/lib,
    #  so we reset ARCH_LIB to the appropriate things.
    ifeq ($(SITE_NAME), JCVI)
      ARCH_LIB     = /usr/lib64 /usr/X11R6/lib64
    endif
  endif
  ifeq ($(MACHINETYPE), ia64)
    # Don't set X86_GCC_LINUX because IEEE floating point mode is not available.
  endif
endif


ifeq ($(OSTYPE), FreeBSD)
  CC               = gcc
  CXX              = g++

  ARCH_LDFLAGS    += -llthread -llgcc_r
  ARCH_CFLAGS      = -D_THREAD_SAFE -I/usr/local/include/pthread/linuxthreads 

  ifeq ($(BUILDCOVERAGE), 1)
    ARCH_CFLAGS   += -g -Wimplicit -fprofile-arcs -ftest-coverage
    ARCH_LDFLAGS  += -lgcov
  else
    ifeq ($(BUILDDEBUG), 1)
      ARCH_CFLAGS   += -g -Wimplicit -pg
    else
      ARCH_CFLAGS   += -O3 -Wimplicit
    endif
  endif

  ARCH_INC         = /usr/local/include /usr/X11R6/include
  ARCH_LIB         = /usr/local/lib     /usr/X11R6/lib
endif

ifeq ($(OSTYPE), Darwin)
  CC               = gcc
  CXX              = g++
  ARCH_CFLAGS      = -pipe -DNEEDXDRUHYPER -D_THREAD_SAFE
  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g
  else
    ARCH_CFLAGS   += -fast
  endif
  ARCH_LIB         = /usr/local/lib /usr/X11R6/lib
endif

ifeq ($(OSTYPE), SunOS)
  ifeq ($(MACHINETYPE), i86pc)
    CC             = gcc
    CXX            = g++
    ARCH_CFLAGS    = -DBYTE_ORDER=LITTLE_ENDIAN
    ARCH_LDFLAGS  += -lm
  else
    CC             = cc
    CXX            = cxx
    ARCH_CFLAGS    = -DBYTE_ORDER=BIG_ENDIAN -DMAP_FILE=0
    ARCH_LDFLAGS  += -lsocket -lnsl
  endif
  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g
  else
    ARCH_CFLAGS   += -fast
  endif
endif

ifeq ($(OSTYPE), aix)
  CC               = xlc
  CXX              = xlC_r
  ARCH_CFLAGS      = -qflag=W:W \
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
  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g
  else
    ARCH_CFLAGS   += -O3
  endif
  ARCH_LDFLAGS     = -q64 \
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
endif

ifeq ($(OSTYPE), OSF1)
  CC               = cc
  CXX              = cxx
  ARCH_CFLAGS      = -w0 \
                     -warnprotos \
                     -trapuv \
                     -float_const \
                     -readonly_strings \
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
                     -msg_disable bitnotint
  ARCH_CFLAGS     += -pthread
  ARCH_CFLAGS      = -warnprotos \
                     -trapuv \
                     -float_const \
                     -readonly_strings \
                     -msg_disable extrasemi
  ARCH_CXXFLAGS    = -tlocal

  #  We used to include -fast in optimized builds, but that results in F_pow() (the fastmath pow function)
  #  being undefined for c++ code.

  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g
  else
    ARCH_CFLAGS   += -O4

  endif
  ARCH_INC         =  /usr/local/include /usr/include
endif

ifeq ($(BUILDPROFILE), 1)
  ARCH_CFLAGS  += -pg
  ARCH_LDFLAGS += -pg
endif

# One can argue that CXXFLAGS should be separate.  For now, we only
# add to the flags.

CFLAGS          += $(ARCH_CFLAGS)
CXXFLAGS        += $(ARCH_CFLAGS) $(ARCH_CXXFLAGS)
LDFLAGS         += $(ARCH_LDFLAGS)

INC_IMPORT_DIRS += $(LOCAL_WORK)/src $(patsubst %, $(LOCAL_WORK)/src/%, $(strip $(SUBDIRS)))
INC_IMPORT_DIRS += $(ARCH_INC)

LIB_IMPORT_DIRS += $(LOCAL_LIB)
LIB_IMPORT_DIRS += $(ARCH_LIB)

OBJ_SEARCH_PATH  = $(LOCAL_OBJ)


ifeq ($(SITE_NAME), TIGR)
  CFLAGS   += -DUSE_SOAP_UID
  CXXFLAGS += -DUSE_SOAP_UID
endif

ifeq ($(SITE_NAME), JCVI)
  LDFLAGS += -lcurl
endif

#
#  AS_SIM is no longer a supported component.
#
#  The order of compilation here is very carefully chosen to be the
#  same as the order used in running an assembly.  It is extremely
#  useful if you happen to be making changes to, say, the persistent
#  stores.  Break the continuation lines after AS_GKP and you'll build
#  just gatekeeper.
#
#  It's also more-or-less telling us link relationships.  CGB doesn't
#  use REZ or CNS, etc.  It gets hairy at the end; REZ, CNS and CGW are
#  all dependent on each other.  Everything after TER isn't needed for
#  an assembly.

#CFLAGS += -I/n8/wgs/src/AS_CNS -I/n8/wgs/src/AS_CGW -I/n8/wgs/src/AS_ALN -I/n8/wgs/src/AS_REZ -I/n8/wgs/src/AS_SDB

SUBDIRS = AS_RUN \
          AS_UTL \
          AS_UID \
          AS_MSG \
          AS_PER \
          AS_GKP \
          AS_OBT \
          AS_MER \
          AS_OVL \
          AS_ALN \
          AS_CGB \
          AS_BOG \
          AS_REZ \
          AS_CNS \
          AS_SDB \
          AS_LIN \
          AS_CGW \
          AS_TER \
          AS_CVT \
          AS_VWR \
          AS_MPA

D=\
          AS_ORA \
