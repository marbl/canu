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

LOCAL_WORK = $(shell cd ../../; pwd)

SOURCES         = greedyFragmentTiling.C \
		  chimChe.C \
                  metagenomics_ovl_analyses.C

SCRIPTS		= greedy_layout_to_IUM.pl \
		  utg2fasta.pl 		  \
		  run_greedy.csh

OBJECTS         = $(SOURCES:.C=.o)
LIBRARIES       =
CXX_PROGS       = greedyFragmentTiling metagenomics_ovl_analyses chimChe

include $(LOCAL_WORK)/src/c_make.as

all: $(OBJECTS) $(LIBRARIES) $(CXX_PROGS) $(SCRIPTS)
	@cp -f greedy_layout_to_IUM.pl $(LOCAL_BIN)/greedy_layout_to_IUM
	@chmod 775 $(LOCAL_BIN)/greedy_layout_to_IUM
	@cp -f utg2fasta.pl $(LOCAL_BIN)/utg2fasta
	@chmod 775 $(LOCAL_BIN)/utg2fasta
	@cp -f run_greedy.csh $(LOCAL_BIN)/run_greedy.csh
	@chmod 775 $(LOCAL_BIN)/run_greedy.csh

greedyFragmentTiling: greedyFragmentTiling.o libCA.a

metagenomics_ovl_analyses: metagenomics_ovl_analyses.o libCA.a

chimChe: chimChe.o libCA.a
