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

COMMON_SOURCES = wpack.c button.c
COMMON_OBJECTS = $(COMMON_SOURCES:.c=.o)

CELAMY_SOURCES = layout.c celamy.c hash.c agrep.c
CELAMY_OBJECTS = $(CELAMY_SOURCES:.c=.o)

CELAGRAM_SOURCES = celagram.c
CELAGRAM_OBJECTS = $(CELAGRAM_SOURCES:.c=.o)

QUERY_SOURCES = query_celamy.c layout_uint.c
SOURCES = $(CELAMY_SOURCES) $(CELAGRAM_SOURCES) $(QUERY_SOURCES) $(COMMON_SOURCES)
OBJECTS = $(SOURCES:.c=.o)

CXX_PROGS = celamy celagram query_celamy

## Pick up global definitions after defining OBJECTS and C_PROGS, but before
## Using OSTYPE
LOCAL_WORK = $(shell cd ../../; pwd)
include $(LOCAL_WORK)/src/c_make.as

#  This _should_ be localized to this module, and not propagate to
#  other modules.
LDFLAGS += -lXt -lX11 -lm -lSM -lICE
ifeq ($(OSTYPE),Linux)
  LDFLAGS += -ldl
  ifeq ($(MACHINETYPE), x86_64)
    LDFLAGS += -L/usr/X11R6/lib64
  endif
else
  ifeq ($(OSTYPE),OSF1)
    LDFLAGS += -ldnet_stub
  endif
endif

all: $(OBJECTS) $(C_PROGS)
	@test -n nop

celamy: $(CELAMY_OBJECTS) $(COMMON_OBJECTS) libAS_UTL.a

query_celamy: layout_uint.o query_celamy.o hash.o agrep.o $(COMMON_OBJECTS) libAS_UTL.a

celagram: $(CELAGRAM_OBJECTS) $(COMMON_OBJECTS) libAS_UTL.a

package:
	shar *.c *.h >package
