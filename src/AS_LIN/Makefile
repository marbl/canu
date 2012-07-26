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
# gmake Makefile for AS_LIN

REQD_LIBS = libAS_LIN.a

LIB_SOURCES = dgbtrf.c dgbtf2.c dgbtrs.c \
              iparmq.c ieeeck.c dlaswp.c disnan.c dlaisnan.c \
              dcopy.c dswap.c dger.c \
              dpbtrf.c dpbtrs.c daxpy.c dgemm.c dgemv.c dlamch.c dpbtf2.c \
              dpotf2.c dptrfs.c dpttrs.c dsyr.c dsyrk.c dtbsv.c dtrsm.c \
              idamax.c ilaenv.c lsame.c xerbla.c s_copy.c s_cmp.c \
              ddot.c dscal.c

LIB_OBJECTS = $(LIB_SOURCES:.c=.o)

LIBRARIES = libAS_LIN.a libCA.a

AS_LIN_SRCS = dgbtrf.c dgbtf2.c dgbtrs.c \
              iparmq.c ieeeck.c dlaswp.c disnan.c dlaisnan.c \
              dcopy.c dswap.c dger.c \
              dpbtrf.c dpbtrs.c daxpy.c dgemm.c dgemv.c dlamch.c dpbtf2.c \
              dpotf2.c dptrfs.c dpttrs.c dsyr.c dsyrk.c dtbsv.c dtrsm.c \
              idamax.c ilaenv.c lsame.c xerbla.c s_copy.c s_cmp.c \
              ddot.c dscal.c

AS_LIN_OBJS = $(AS_LIN_SRCS:.c=.o)

LOCAL_WORK = $(shell cd ../..; pwd)

SOURCES = $(AS_LIN_SRCS)
OBJECTS = $(AS_LIN_OBJS)

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(LIB_OBJECTS) $(LIBRARIES)
	@test -n nop

libAS_LIN.a: $(LIB_OBJECTS)

libCA.a: $(LIB_OBJECTS)
