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

C_PROG_SOURCES = fragmentDepth.c
C_PROG_OBJECTS = $(C_PROG_SOURCES:.c=.o)

C_SOURCES = $(C_PROG_SOURCES)
C_OBJECTS = $(C_SOURCES:.c=.o)

CXX_SOURCES =
CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)

SOURCES = $(C_SOURCES) $(CXX_SOURCES)
OBJECTS = $(C_OBJECTS) $(CXX_OBJECTS)

C_PROGS = fragmentDepth

SCRIPTS = runCA \
          assemblyCompare \
          qcCompare \
          tampaCompare \
          mummerCompare \
          testVersionOnPlatform.sh \
          apply_assembler \
          glen_est_truncadjusted.pl \
          caqc.pl \
          caqc_help.ini \
          caqc-diff.pl \
          ca2ace.pl \
          asmToAGP.pl \
          markUniqueUnique.rb \
          toggleAndRun.sh \
          spec/OBT.specFile \
          spec/noOBT.specFile \
          spec/noVec.specFile \
          spec/metaGenomic.specFile


#  JCVI likes to have both a modified (carun) and unmodified (runCA)
#  version.  Other sites will probably be happy just redefining
#  RUNCALOCAL.

RUNCALOCAL = runCA-OBT/runCAutil/site-generic.pl

ifeq ($(SITE_NAME),JCVI)
  SCRIPTS   += carun
  RUNCALOCAL = runCA-OBT/runCAutil/site-generic.pl
  CARUNLOCAL = runCA-OBT/runCAutil/site-jcvi.pl
endif

ifeq ($(SITE_NAME),BRI)
  RUNCALOCAL = runCA-OBT/runCAutil/site-bri.pl
endif

# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(C_PROGS) $(SCRIPTS)
	@test -n nop

RUNCAOBTUTIL = runCA-OBT/runCAutil/checkOverlap.pl \
               runCA-OBT/runCAutil/createConsensusJobs.pl \
               runCA-OBT/runCAutil/merOverlapper.pl \
               runCA-OBT/runCAutil/createOverlapJobs.pl \
               runCA-OBT/runCAutil/createOverlapStore.pl \
               runCA-OBT/runCAutil/createPostUnitiggerConsensus.pl \
               runCA-OBT/runCAutil/meryl.pl \
               runCA-OBT/runCAutil/overlapCorrection.pl \
               runCA-OBT/runCAutil/overlapTrim.pl \
               runCA-OBT/runCAutil/preoverlap.pl \
               runCA-OBT/runCAutil/scaffolder.pl \
               runCA-OBT/runCAutil/terminate.pl \
               runCA-OBT/runCAutil/cleaner.pl \
               runCA-OBT/runCAutil/unitigger.pl \
	       runCA-OBT/runCAutil/umdOverlapper.pl \
	       runCA-OBT/runCAutil/vectorTrim.pl

runCA: runCA-OBT/runCA.pl.head $(RUNCALOCAL) runCA-OBT/runCAutil/util.pl runCA-OBT/runCA.pl $(RUNCAOBTUTIL) runCA-OBT/runCAutil/scheduler.pm
	@cat runCA-OBT/runCA.pl.head $(RUNCALOCAL) runCA-OBT/runCAutil/util.pl $(RUNCAOBTUTIL) runCA-OBT/runCA.pl runCA-OBT/runCAutil/scheduler.pm > $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@
	@ln -fs $@ $(LOCAL_BIN)/runCA-OBT.pl

carun: runCA-OBT/runCA.pl.head $(CARUNLOCAL) runCA-OBT/runCAutil/util.pl runCA-OBT/runCA.pl $(RUNCAOBTUTIL) runCA-OBT/runCAutil/scheduler.pm
	@cat runCA-OBT/runCA.pl.head $(CARUNLOCAL) runCA-OBT/runCAutil/util.pl $(RUNCAOBTUTIL) runCA-OBT/runCA.pl runCA-OBT/runCAutil/scheduler.pm > $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

markUniqueUnique.rb: asmQC/markUniqueUnique.rb
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

toggleAndRun.sh: runCA-OBT/toggleAndRun.sh
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

TIGR: asmQC/TIGR
	@cp -rp $< $(LOCAL_BIN)

caqc_help.ini: asmQC/caqc_help.ini
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

caqc.pl: asmQC/caqc.pl
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

caqc-diff.pl: asmQC/caqc-diff.pl
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

ca2ace.pl: asmQC/ca2ace.pl TIGR
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

asmToAGP.pl: asmQC/asmToAGP.pl
	@cp -p $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

assemblyCompare: testing/assemblyCompare.pl
	@cp -p testing/assemblyCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

qcCompare: testing/qcCompare.pl
	@cp -p testing/qcCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

tampaCompare: testing/tampaCompare.pl
	@cp -p testing/tampaCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

mummerCompare: testing/mummerCompare.pl
	@cp -p testing/mummerCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

apply_assembler: testing/apply_assembler.pl
	@cp -p testing/apply_assembler.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

testVersionOnPlatform.sh: testing/$@
	@cp -p testing/$@ $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

glen_est_truncadjusted.pl: dataQC/$@
	@cp -p dataQC/$@ $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

spec:
	@mkdir $(LOCAL_BIN)/spec

spec/%.specFile: spec runCA-OBT/spec/%.specFile
	@cp -p runCA-OBT/$@ $(LOCAL_BIN)/spec

fragmentDepth: fragmentDepth.o libCA.a
