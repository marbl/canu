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

SOURCES = uidclient.c uidsrvd.c SYS_UIDcommon.c SYS_UIDerror.c SYS_UIDclient.c
C_PROGS = uidclient uidsvrd

LIBRARIES = libAS_UID.a libCA.a

LOCAL_WORK = $(shell cd ../../; pwd)

include $(LOCAL_WORK)/src/site_name.as

SITEOBJS =

ifeq ($(SITE_NAME), JCVI)
  SOURCES        += SYS_UIDclient_JTC.c
  SITEOBJS        = SYS_UIDclient_JTC.o
endif

ifeq ($(SITE_NAME), LOCAL)
  SOURCES        += SYS_UIDclient_LOCAL.c
  SITEOBJS       += SYS_UIDclient_LOCAL.o
endif

ifeq ($(SITE_NAME), BRI)
  SOURCES        += SYS_UIDclient_LOCAL.c
  SITEOBJS       += SYS_UIDclient_LOCAL.o
endif

OBJECTS = $(SOURCES:.c=.o)

## Now include c_make.as with the updated SOURCES and OBJECTS
include $(LOCAL_WORK)/src/c_make.as

all: $(LIBRARIES) $(C_PROGS)
	@test -n nop

libAS_UID.a: SYS_UIDcommon.o SYS_UIDclient.o SYS_UIDerror.o $(SITEOBJS)
libCA.a:     SYS_UIDcommon.o SYS_UIDclient.o SYS_UIDerror.o $(SITEOBJS)

uidsvrd:   uidsrvd.o   libAS_UID.a libAS_UTL.a
uidclient: uidclient.o libAS_UID.a libAS_UTL.a
