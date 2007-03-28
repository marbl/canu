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

CLIENT_SOURCES = \
  SYS_UIDcommon.c \
  SYS_UIDerror.c \
  SYS_UIDclient_stand_alone.c \
  SYS_UIDclient.c \
  SYS_UIDclient_simple.c \
  SYS_UIDclient_test.c

SERVER_SOURCES = \
  SYS_UIDserver.c \
  SYS_UIDserver_main.c

SOURCES = $(CLIENT_SOURCES) $(SERVER_SOURCES)

#  The C++ client looks outdated.
#
#  SYS_UIDclient_cpp.C \
#  SYS_UIDclient_test_cpp.C \
#
#  uid_client_test_cpp 
#  libAS_UIDcpp.a

CLIENT_PROGS = \
  uidclient \
  uid_client_test

SERVER_PROGS = \
  uidsvrd

PROGS = $(CLIENT_PROGS) $(SERVER_PROGS)


LIBRARIES = libAS_UID.a libCA.a


LOCAL_WORK = $(shell cd ../../; pwd)

## Include to set UID_SERVER
#include uid_transport.as

include $(LOCAL_WORK)/src/site_name.as

SITEOBJS = 

## TIGR EUID server is now deprecated and is no longer used
#ifeq ($(SITE_NAME), TIGR)
#  SOURCES        += SYS_UIDclient_SOAP.c checkEuidService.c
#  PROGS          += checkEuidService
#  SOAPSOURCES     = stdsoap2.c euidC.c euidClientLib.c 
#  SOURCES        += $(SOAPSOURCES)
#  SITEOBJS        = SYS_UIDclient_SOAP.o $(SOAPSOURCES:.c=.o)
#endif

ifeq ($(SITE_NAME), JCVI)
  SOURCES        += SYS_UIDclient_JTC.c checkEuidService.c
  SITEOBJS        = SYS_UIDclient_JTC.o
  PROGS          += checkEuidService
endif

ifeq ($(SITE_NAME), LOCAL)
  SOURCES        += SYS_UIDclient_LOCAL.c
  SITEOBJS       += SYS_UIDclient_LOCAL.o
endif


## Now include c_make.as with the updated SOURCES and OBJECTS
include $(LOCAL_WORK)/src/c_make.as

OBJECTS = \
  $(SOURCES:.c=.o)

all: $(LIBRARIES) $(PROGS)
	@test -n nop

libAS_UID.a: SYS_UIDcommon.o SYS_UIDclient.o SYS_UIDclient_simple.o SYS_UIDerror.o $(SITEOBJS)

libCA.a: SYS_UIDcommon.o SYS_UIDclient.o SYS_UIDclient_simple.o SYS_UIDerror.o $(SITEOBJS)

#libAS_UIDcpp.a: SYS_UIDclient_cpp.o

uidsvrd: SYS_UIDcommon.o SYS_UIDserver.o SYS_UIDerror.o SYS_UIDserver_main.o libAS_UTL.a

uid_client_test: SYS_UIDclient_test.o libAS_UID.a libAS_UTL.a

checkEuidService: checkEuidService.o libAS_UID.a libAS_UTL.a

#uid_client_test_cpp: SYS_UIDclient_cpp.o SYS_UIDclient_test_cpp.o libAS_UID.a libAS_UTL.a

uidclient: SYS_UIDclient_stand_alone.o libAS_UID.a libAS_UTL.a

