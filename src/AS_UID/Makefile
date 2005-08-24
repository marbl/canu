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

######################################################################
# $Source: /work/NIGHTLY/wgs-assembler-cvs/src/AS_UID/Makefile,v $
# $Revision: 1.5 $
# $Date: 2005-08-24 10:57:43 $
# $Name: not supported by cvs2svn $
# $Author: brianwalenz $
# $Log: not supported by cvs2svn $
# Revision 1.4  2005/08/17 16:57:11  brianwalenz
# Rework the way the UID server is selected; now define either
# UID_SERVER=TIGR or UID_SERVER=JCVI in AS_UID/uid_transport.as.
# c_make.as then sets CFLAGS,CXXFLAGS,LDFLAGS appropriately.  Makefiles
# that need to access the UID server via curl needed to be modified
# to remove the explicit link to libcurl.a; everyone now links
# against it.  Hopefully, cleaner, and easier to extent to other
# sites.
#
# Revision 1.3  2005/03/22 19:49:28  jason_miller
# The TIGR tip as of March 22 2005. Commit by Jason Miller at TIGR.
#
# Revision 1.3  2004/09/10 12:31:43  mschatz
# Add standard copyright notice
#
# Revision 1.2  2004/09/09 22:34:12  mschatz
# USE_SOAP_UID support, cleanup
#
# Revision 1.1  2004/06/24 12:51:06  mpop
# Added AS_UID
#
# Revision 1.3  2003/05/13 19:20:08  mpop
# Fixed the UID client code, now the terminator can terminate :)
#
# Revision 1.2  2003/05/09 21:04:02  mpop
# Dos2unixed all files.
# Modified c_make.as to set SEP_PATH relative to LOCAL_WORK
#
# Revision 1.1.1.1  2003/05/08 18:40:11  aaronhalpern
# versions from TIGR
#
# Revision 1.2  2001/09/25 23:03:20  mpop
# Dos2Unixed
#
# Revision 1.1.1.1  2001/09/25 20:21:05  mpop
# Celera Assembler
#
# Revision 1.8  1999/10/13 19:03:53  sdmurphy
# added standalone client
#
# Revision 1.7  1999/07/22 01:13:46  dewim
# added += to *_IMPORT_DIRS assignments
# now including $(SEP_PATH)/makefiles/c_make.gen
#
# Revision 1.6  1999/06/30 20:57:52  stine
# Changed the makefile to use the cds.h from cds/SYS as determined
# by the SEPTIC group.
#
# Revision 1.5  1999/03/04 22:41:33  sdmurphy
# added SYS_UIDerror.o to libuid.a
#
# Revision 1.4  1999/01/13 20:09:43  sdmurphy
# changed uid_server to uidsvrd
#
# Revision 1.3  1999/01/13 14:27:47  sdmurphy
# version 0 prelim
#
# Revision 1.2  1998/12/22 21:17:01  sdmurphy
# added c++ support
#
# Revision 1.1  1998/12/17 21:20:12  sdmurphy
# Makefile for UID libs and progs
#
######################################################################

######################################################################
# Module:
#
# Description:
#
# Assumptions:
#
######################################################################


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


LIBRARIES = \
  libAS_UID.a 


## Include to test if TIGR_UID_SERVER is true
include uid_transport.as

SOAPOBJS = 
ifeq ($(UID_SERVER), "TIGR")
  SOURCES     += checkEuidService.c
  PROGS       += checkEuidService
  SOAPSOURCES  = stdsoap2.c euidC.c euidClientLib.c 
  SOURCES     += $(SOAPSOURCES)
  SOAPOBJS     = $(SOAPSOURCES:.c=.o)
endif

OBJECTS = \
  $(SOURCES:.c=.o)

## Now reinclude c_make.as with the updated SOURCES and OBJECTS
LOCAL_WORK = $(shell cd ../../; pwd)
include $(LOCAL_WORK)/src/c_make.as

all: $(LIBRARIES) $(PROGS)

libAS_UID.a: SYS_UIDcommon.o SYS_UIDclient.o SYS_UIDclient_simple.o SYS_UIDerror.o $(SOAPOBJS)

libAS_UIDcpp.a: SYS_UIDclient_cpp.o

uidsvrd: SYS_UIDcommon.o SYS_UIDserver.o SYS_UIDerror.o SYS_UIDserver_main.o

## Commenting out libssl.a libcrypto.a libz.a
uid_client_test: SYS_UIDclient_test.o libAS_UID.a

checkEuidService: checkEuidService.o libAS_UID.a

uid_client_test_cpp: SYS_UIDclient_cpp.o SYS_UIDclient_test_cpp.o libAS_UID.a

uidclient: SYS_UIDclient_stand_alone.o libAS_UID.a

