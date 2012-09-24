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


# Use the CXXFLAGS environment variable to enhance optimizations, e.g.:
#   -mtune=native -march=native


# You can enable a debugging build, disabling all optimizations, by
# setting BUILDDEBUG to 1.
#
# You can enable a profiling build by setting BUILDPROFILE to 1.
#
ifneq "$(origin BUILDDEBUG)" "environment"
BUILDDEBUG     = 0
endif
ifneq "$(origin BUILDPROFILE)" "environment"
BUILDPROFILE   = 0
endif


include $(LOCAL_WORK)/src/site_name.as
include $(LOCAL_WORK)/src/c_make.gen



ifeq ($(OSTYPE), Linux)
  ARCH_CFLAGS  += -D_GLIBCXX_PARALLEL -fopenmp
  ARCH_LDFLAGS += -D_GLIBCXX_PARALLEL -fopenmp

  ifeq ($(MACHINETYPE), i686)
    ARCH_LDFLAGS  += -pthread -lm
    ARCH_CFLAGS   += -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DX86_GCC_LINUX
    ARCH_LIB       = /usr/X11R6/lib
  endif
  ifeq ($(MACHINETYPE), amd64)
    ARCH_LDFLAGS  += -pthread -lm
    ARCH_CFLAGS   += -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DX86_GCC_LINUX
    ARCH_LIB       = /usr/lib64 /usr/X11R6/lib64
  endif
  ifeq ($(MACHINETYPE), ia64)
    ARCH_LDFLAGS  += -pthread -lm
    ARCH_CFLAGS   += -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
    ARCH_LIB       = /usr/X11R6/lib
  endif

  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS  += -g
    ARCH_LDFLAGS +=
  else
    #  gcc412 doesn't know these
    ARCH_CFLAGS  += -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    ARCH_LDFLAGS += -Wl,-O1
  endif
endif



#  NOTE:  using -pedantic generates many warnings "use of C99 long long integer constant" in kmer

ifeq ($(OSTYPE), FreeBSD)
  CC  = gcc46
  CXX = g++46

  ifeq ($(MACHINETYPE), i386)
    ARCH_LDFLAGS    += -pthread -lthr -lm
    ARCH_CFLAGS      = -pthread -Wall -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare
    ARCH_CFLAGS      = -pthread       -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare
  endif
  ifeq ($(MACHINETYPE), amd64)
    ARCH_LDFLAGS    += -pthread -lthr -lm
    ARCH_CFLAGS      =  -pthread               -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat
    ARCH_CFLAGS      =  -pthread -Wall -Wextra -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare -Wformat
  endif

  ifeq ($(BUILDPROFILE), 1)
    ARCH_CFLAGS  += -fopenmp
    ARCH_LDFLAGS += -fopenmp -rpath /usr/local/lib/gcc46
  else
    ARCH_CFLAGS  += -D_GLIBCXX_PARALLEL -fopenmp
    ARCH_LDFLAGS += -D_GLIBCXX_PARALLEL -fopenmp -rpath /usr/local/lib/gcc46
  endif

  ifeq ($(BUILDDEBUG), 1)
    #  Inconveniently fix a problem where gcc46 doesn't work with gdb.
    CC  = gcc44
    CXX = g++44
    ARCH_CFLAGS   += -g
  else
    ifeq ($(BUILDPROFILE), 1)
      ARCH_CFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions
    else
      ARCH_CFLAGS   += -O4 -funroll-loops -fexpensive-optimizations -finline-functions -fomit-frame-pointer
    endif
  endif

  ARCH_INC         = /usr/local/include /usr/X11R6/include
  ARCH_LIB         = /usr/local/lib     /usr/X11R6/lib
endif



ifeq ($(OSTYPE), Darwin)
  CC               = gcc
  CXX              = g++
  ARCH_CFLAGS      = -D_THREAD_SAFE

  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g -Wall
  else
    ARCH_CFLAGS   += -fast
  endif

  ARCH_CFLAGS  += -D_GLIBCXX_PARALLEL -fopenmp
  ARCH_LDFLAGS += -D_GLIBCXX_PARALLEL -fopenmp

  ARCH_CFLAGS += -fPIC -m64 -fmessage-length=0 -Wimplicit -Wno-write-strings -Wno-unused -Wno-char-subscripts -Wno-sign-compare
# ARCH_CFLAGS += -Wshorten-64-to-32  # Wow, tough
# ARCH_CFLAGS += -Wextra
# ARCH_CFLAGS += -pedantic  (see above about pedantic)

  ARCH_LDFLAGS  += -m64 -lm

  ARCH_LIB         = /opt/local/lib /usr/X11R6/lib
endif



# Use "gmake SHELL=/bin/bash".
#
ifeq ($(OSTYPE), SunOS)
  ifeq ($(MACHINETYPE), i86pc)
    ARCH_CFLAGS    = -DBYTE_ORDER=LITTLE_ENDIAN -DANSI_C -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -pthreads
    ARCH_LDFLAGS  += -lm
  endif

  ifeq ($(MACHINETYPE), sparc32)
    ARCH_CFLAGS    = -DANSI_C -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -pthreads
    ARCH_LDFLAGS  += -lm -lnsl -lsocket
  endif

  ifeq ($(MACHINETYPE), sparc64)
    ARCH_CFLAGS    = -m64 -DANSI_C -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -pthreads
    ARCH_LDFLAGS  += -m64 -lm -lnsl -lsocket
  endif

  ifeq ($(BUILDDEBUG), 1)
    ARCH_CFLAGS   += -g
  else
    ARCH_CFLAGS   += -O
  endif
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

ifeq ($(SITE_NAME), JCVI)
  LDFLAGS += -lcurl
endif

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

# Broken by BPW's string UID hack: AS_CVT, AS_MPA.  AS_CVT might work,
# but its only used by AS_MPA.

# Broken by the C++ switch: AS_VWR (won't link on Linux64)

SUBDIRS = AS_RUN \
          AS_UTL \
          AS_UID \
          AS_MSG \
          AS_PER \
          AS_GKP \
          AS_OBT \
          AS_MER \
          AS_OVL \
          AS_OVM \
          AS_OVS \
          AS_ALN \
          AS_CGB \
          AS_BOG \
          AS_BAT \
          AS_PBR \
          AS_REZ \
          AS_CNS \
          AS_LIN \
          AS_CGW \
          AS_TER \
          AS_ENV \
          AS_ARD \
          AS_REF
