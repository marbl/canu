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
          assemblyCompare \
          qcCompare \
          tampaCompare \
          mummerCompare \
          glen_est_truncadjusted.pl \
          caqc.pl \
          caqc_help.ini \
          caqc-diff.pl \
          ca2ace.pl \
          asmToAGP.pl \
          spec/OBT.specFile \
          spec/noOBT.specFile \
          spec/noVec.specFile \
          spec/metaGenomic.specFile


# Include for AS project rules
include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS) $(SCRIPTS)
	@test -n nop

RUNCAOBTUTIL = runCA-OBT/runCAutil/checkOverlap.pl \
               runCA-OBT/runCAutil/createConsensusJobs.pl \
               runCA-OBT/runCAutil/merOverlapper.pl \
               runCA-OBT/runCAutil/merTrim.pl \
               runCA-OBT/runCAutil/createOverlapJobs.pl \
               runCA-OBT/runCAutil/createOverlapStore.pl \
               runCA-OBT/runCAutil/createPostUnitiggerConsensus.pl \
               runCA-OBT/runCAutil/meryl.pl \
               runCA-OBT/runCAutil/overlapCorrection.pl \
               runCA-OBT/runCAutil/overlapTrim.pl \
               runCA-OBT/runCAutil/preoverlap.pl \
               runCA-OBT/runCAutil/shred.pl \
               runCA-OBT/runCAutil/scaffolder.pl \
               runCA-OBT/runCAutil/terminate.pl \
               runCA-OBT/runCAutil/cleaner.pl \
               runCA-OBT/runCAutil/unitigger.pl \
	       runCA-OBT/runCAutil/umdOverlapper.pl \
	       runCA-OBT/runCAutil/vectorTrim.pl \
	       runCA-OBT/runCAutil/toggler.pl

runCA: runCA-OBT/runCA.pl.head runCA-OBT/runCAutil/util.pl runCA-OBT/runCA.pl $(RUNCAOBTUTIL) runCA-OBT/runCAutil/scheduler.pm
	@cat runCA-OBT/runCA.pl.head runCA-OBT/runCAutil/util.pl $(RUNCAOBTUTIL) runCA-OBT/runCA.pl runCA-OBT/runCAutil/scheduler.pm > $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@
	@ln -fs $@ $(LOCAL_BIN)/runCA-OBT.pl

TIGR: asmQC/TIGR
	@cp -r $< $(LOCAL_BIN)

caqc_help.ini: asmQC/caqc_help.ini
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

caqc.pl: asmQC/caqc.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

caqc-diff.pl: asmQC/caqc-diff.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

ca2ace.pl: asmQC/ca2ace.pl TIGR
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

asmToAGP.pl: asmQC/asmToAGP.pl
	@cp $< $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

assemblyCompare: testing/assemblyCompare.pl
	@cp testing/assemblyCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

qcCompare: testing/qcCompare.pl
	@cp testing/qcCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

tampaCompare: testing/tampaCompare.pl
	@cp testing/tampaCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

mummerCompare: testing/mummerCompare.pl
	@cp testing/mummerCompare.pl $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

glen_est_truncadjusted.pl: dataQC/$@
	@cp dataQC/$@ $(LOCAL_BIN)/$@
	@chmod 775 $(LOCAL_BIN)/$@

spec:
	@mkdir $(LOCAL_BIN)/spec

spec/%.specFile: spec runCA-OBT/spec/%.specFile
	@cp runCA-OBT/$@ $(LOCAL_BIN)/spec

fragmentDepth: fragmentDepth.o libCA.a

markUniqueUnique: markUniqueUnique.o libCA.a
