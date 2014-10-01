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

LOCAL_WORK = $(shell cd ../..; pwd)

SOURCES = fragmentDepth.C markUniqueUnique.C
OBJECTS = fragmentDepth.o markUniqueUnique.o

CXX_PROGS = fragmentDepth markUniqueUnique

SCRIPTS = runCA \
          runCA-dedupe \
          runCA-overlapStoreBuild \
          mergeqc.pl \
          caqc.pl \
          caqc_help.ini \
          ca2ace.pl \
          asmToAGP.pl \
          replaceUIDwithName-fastq \
          replaceUIDwithName-posmap


# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS) $(SCRIPTS)
	@test -n nop

RUNCAOBTUTIL =

runCA: runCA.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

runCA-dedupe: runCA-dedupe.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

runCA-overlapStoreBuild: runCA-overlapStoreBuild.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

TIGR: asmQC/TIGR
	@mkdir $(LOCAL_BIN)/TIGR
	@cp $</AsmLib.pm $(LOCAL_BIN)/TIGR
	@cp $</Foundation.pm $(LOCAL_BIN)/TIGR

caqc_help.ini: asmQC/caqc_help.ini
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

caqc.pl: asmQC/caqc.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

ca2ace.pl: asmQC/ca2ace.pl TIGR
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

asmToAGP.pl: asmQC/asmToAGP.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

mergeqc.pl: asmQC/mergeqc.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

replaceUIDwithName-fastq: replaceUIDwithName-fastq.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

replaceUIDwithName-posmap: replaceUIDwithName-posmap.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

fragmentDepth: fragmentDepth.o libCA.a

markUniqueUnique: markUniqueUnique.o libCA.a
