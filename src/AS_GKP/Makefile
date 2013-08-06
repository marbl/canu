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

GKPSRC    = AS_GKP_buildPartition.C \
            AS_GKP_checkLibrary.C \
            AS_GKP_checkFrag.C \
            AS_GKP_checkLink.C \
            AS_GKP_checkPlace.C \
            AS_GKP_dump.C \
            AS_GKP_edit.C \
            AS_GKP_errors.C \
            AS_GKP_illumina.C

SOURCES   = AS_GKP_main.C \
            AS_GKP_bench.C \
            $(GKPSRC) \
            sffToCA.C \
            fastqToCA.C \
            fastqSample.C \
            fastqSimulate.C \
            fastqSimulate-sort.C \
            fastqAnalyze.C \
            gkpStoreDumpFASTQ.C \
            gkpStoreCreate.C \
            upgrade-v8-to-v9.C \
            upgrade-v9-to-v10.C
OBJECTS   = $(SOURCES:.C=.o)
LIBRARIES =
CXX_PROGS = gatekeeper \
            gatekeeperbench \
            sffToCA \
            fastqToCA \
            fastqSample \
            fastqSimulate \
            fastqSimulate-sort \
            fastqAnalyze \
            gkpStoreDumpFASTQ \
            gkpStoreCreate \
            upgrade-v8-to-v9 \
            upgrade-v9-to-v10

include $(LOCAL_WORK)/src/c_make.as

all:    $(OBJECTS) $(LIBRARIES) $(CXX_PROGS)
	@test -n nop

gatekeeper:           AS_GKP_main.o $(GKPSRC:.C=.o) libCA.a
gatekeeperbench:      AS_GKP_bench.o                libCA.a
sffToCA:              sffToCA.o                     libCA.a
fastqToCA:            fastqToCA.o                   libCA.a
fastqSample:          fastqSample.o                 libCA.a
fastqSimulate:        fastqSimulate.o               libCA.a
fastqSimulate-sort:   fastqSimulate-sort.o          libCA.a
fastqAnalyze:         fastqAnalyze.o                libCA.a
gkpStoreDumpFASTQ:    gkpStoreDumpFASTQ.o           libCA.a
gkpStoreCreate:       gkpStoreCreate.o              libCA.a
upgrade-v8-to-v9:     upgrade-v8-to-v9.o            libCA.a
upgrade-v9-to-v10:    upgrade-v9-to-v10.o           libCA.a
