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

LOCAL_WORK = $(shell cd ../..; pwd)

AS_ARD_INC =
AS_ARD_LIB =
AS_ARD_LDFLAGS =
AS_ARD_CXXFLAGS =

#currently we build everything but sources aren't defined based on -D flag
#this is because SRC must be fully defnied before including c_make.as
#but we need to include it to get OSTYPE and MACHINETIME and all the settings
COMMON_SOURCES = Sybase.cc
DB_SOURCES = 	AS_ARD_database.cc \
                DBTextOutput.cc \
		BCPOutput.cc \
	 	SQLOutput.cc
SOURCES = loadAssembly.cc getAssembly.cc $(COMMON_SOURCES) $(DB_SOURCES)
OBJECTS = $(SOURCES:.cc=.o)
LIBRARIES =
CXX_PROGS = loadAssembly getAssembly

#minimum version of sybase required to build code
SYBASE_REQUIRED=15.0

include $(LOCAL_WORK)/src/c_make.as

ifeq ($(SITE_NAME), JCVI)
   ifeq ($(OSTYPE), Linux)
		ifndef SYBASE
			SYBASE = /usr/local/packages/sybase15
		endif

                #### check for sybase directory, don't build if it doesn't exist
		SYBASE_EXISTS=$(shell if test -d $(SYBASE); then echo 1; else echo 0; fi;)

		ifeq ($(SYBASE_EXISTS), 1)
                        ### check to make sure version is compatible
			SYBASE_VERSION=$(shell $(SYBASE)/bin/bcp -v | grep -P -o "(\/(\d|\.)+\/)" | awk '{gsub("/", "", $$1); print $$1}')

			ifeq ($(SYBASE_VERSION), $(SYBASE_REQUIRED))
				AS_ARD_CXXFLAGS = -DSYBASE

				ifeq ($(MACHINETYPE), i686)
					SYB_LIB_DIR = $(SYBASE)/OCS/lib $(SYBASE)/OCS/lib3p
					LIB_PATH = $(SYBASE)/OCS/lib:$(SYBASE)/OCS/lib3p
					AS_ARD_LDFLAGS = -lct -lcs
				else
					AS_ARD_CXXFLAGS += -DSYBPLATFORM=linuxamd64
					SYBPLATFORM = linuxamd64
					SYB_LIB_DIR = $(SYBASE)/OCS/lib $(SYBASE)/OCS/lib3p64
					LIB_PATH = $(SYBASE)/OCS/lib:$(SYBASE)/OCS/lib3p64
					AS_ARD_LDFLAGS = -lct64 -lcs64
				endif

				ifdef LD_LIBRARY_PATH
					LIB_PATH := $(LIB_PATH):$(LD_LIBRARY_PATH)
				endif
				export LD_LIBRARY_PATH := $(LIB_PATH)

				AS_ARD_INC = $(SYBASE)/include
				AS_ARD_LIB = $(SYB_LIB_DIR)
			endif
		endif
	endif
endif

INC_IMPORT_DIRS += $(AS_ARD_INC)
LIB_IMPORT_DIRS += $(AS_ARD_LIB)
LDFLAGS += $(AS_ARD_LDFLAGS)
CXXFLAGS += $(AS_ARD_CXXFLAGS)

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

loadAssembly: loadAssembly.o $(COMMON_SOURCES:.cc=.o) $(DB_SOURCES:.cc=.o) libCA.a
getAssembly: getAssembly.o $(COMMON_SOURCES:.cc=.o) libCA.a
