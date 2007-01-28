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
#VERBOSE=1
#######################################################################
#                     Standard Celera Leaf Makefile
#
# Leaf makefiles are intended to be used in a directory containing source
# code to compile that code into libraries and/or programs. They should
# have the following structure:
#   1. variable definitions
#      see the notes below under Variables
#   2. inclusion of project makefile
#      the line in this file should work as-is.
#   3. target definitions - dependencies and commands
#      see the notes below under Targets
#
# A CVS log information block should be maintained at the end of the file.
# NOTE: When adapting a template makefile or copying an existing makefile
#       for use in another source code directory, remember to clean out
#       the log information.
#
# All Celera makefiles should be compatible with gmake, the GNU make
# utility.
#######################################################################

############################## Variables ##############################
#
#   Note that variable definitions may continue onto multiple lines by using
#   the '\' character. Lines beginning with the '#' character are treated as
#   comments.
#
#     LOCAL_WORK =
#       Identifies the <subsystem_level> directory (parent of
#       the inc, src, lib, obj, and bin subdirectories) with
#       an absolute path.
#       For a makefile located in <subsystem_level>/src/<component_level_1>
#       the following definition would work:
#
#          LOCAL_WORK = $(shell cd ../..; pwd)
#
#     SOURCES = 
#       Lists all sources compiled by this makefile.
#
#     OBJECTS = 
#       Lists all objects compiled by this makefile (i.e., one-to-one
#       with SOURCES). OBJECTS should be specified with the following
#       definition to reduce typo errors & ease maintenance:
#
#          OBJECTS = $(SOURCES:.c=.o)
#
#     LIBRARIES = 
#       Lists simple name (in "libx.a" form) of each library that is
#       built by this makefile.
#
#     PROGS = 
#       Lists each executable created by this makefile.
#
# If different sources should be compiled/linked into different
# libraries/applications, define a distinct pair of variables (e.g.,
# PROG1_SOURCES & PROG1_OBJECTS, LIB4_SOURCES & LIB4_OBJECTS) to
# distinguish which source and object files are associated with which
# target. The appropriate objects variable(s) may then be used to
# identify the correct dependencies for each target. Take care in
# defining these sets of source/object files so that the target
# dependencies include everything they should & nothing they shouldn't.
#
############################### Targets ###############################
#
# For each program in PROGS, add a 'target' line identifying the program
# and its dependencies, for example,
#
#          my_program: $(MY_PROGRAM_OBJECTS) $(OBJECTS) libx.a
#
# where MY_PROGRAM_OBJECTS is the set or combination of objects that only
# my_program requires, and OBJECTS is a set of objects that other programs
# or libraries may require.
#
# Also include as the first target, an "all" 'phony' target, e.g.:
#
#          all: $(PROGS)
#
# If you want commands executed when gmake is run with clean, install,
# or regression specified, add a target line for the phony target, and
# enter the commands to be executed on subsequent lines. Note that the
# linese specifying commands MUST start with a tab. e.g.:
#
#          install:
#               cp <file1> <bin_dir>
#               cp <file2> <bin_dir>
#
#######################################################################




######################## Variable Definitions #########################
LOCAL_WORK = $(shell cd ../../; pwd)

#### RESTORED FOR TIGR UNTIL WE ELIMINATE THE SOAP INTERFACE --- JASON
include $(LOCAL_WORK)/src/site_name.as
ifeq ($(SITE_NAME), TIGR)
  USE_SOAP_UID = 1
endif


TERM_SRCS=AS_TER_terminator.c AS_TER_terminator_funcs.c
TERM_OBJS=$(TERM_SRCS:.c=.o)

PRCSCF_SRCS=asmProcessScaffolds_TER.c
PRCSCF_OBJS=$(PRCSCF_SRCS:.c=.o)

CTGOUT_SRCS=asmOutputContigsFasta.c
CTGOUT_OBJS=$(CTGOUT_SRCS:.c=.o)

ASM2BAC_SRCS=asm2bac.c
ASM2BAC_OBJS=$(ASM2BAC_SRCS:.c=.o)

SOURCES = $(TERM_SRCS) $(PRCSCF_SRCS) $(CTGOUT_SRCS) $(ASM2BAC_SRCS)
OBJECTS = $(TERM_OBJS)
PROGS   = terminator asmProcessScaffolds_TER asmOutputContigsFasta asm2bac


################### Project Rule Makefile Inclusion ###################
#include $(LOCAL_WORK)/src/AS_UID/uid_transport.as
include $(LOCAL_WORK)/src/c_make.as

######################### Target Definitions ##########################
all: $(OBJECTS) $(LIBRARIES) $(PROGS)
	@cp -f buildFragContigPosMap.pl $(LOCAL_BIN)
	@chmod 775 $(LOCAL_BIN)/buildFragContigPosMap.pl

terminator:              $(TERM_OBJS)    libCA.a
asmProcessScaffolds_TER: $(PRCSCF_OBJS)  libCA.a
asmOutputContigsFasta:   $(CTGOUT_OBJS)  libCA.a
asm2bac:                 $(ASM2BAC_OBJS) libCA.a
