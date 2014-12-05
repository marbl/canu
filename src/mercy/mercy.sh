#!/bin/sh

# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2006-2007, J. Craig Venter Institute
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

#  Test if the mers in the consensus sequence are supported by mers in
#  the fragments.

#  If we count just the clear, we get a clearer (ha, ha) picture of
#  the assembly quality, while if we count all reads we get a picture
#  of trimming.
#
onlyClear=1
onlyReal=1

mem=8192
mem=16384
mem=24576

ms=22

binroot=/bioinfo/assembly/walenz/src/genomics
asmMerQC=$binroot/meryl/asmMerQC
mapMers=$binroot/meryl/mapMers

dir=/scratch/drosnightly
asm=willi

dir=/project/huref6/redo-consensus_gennady
asm=h6tmp

dir=/project/huref6/assembly/h6
asm=h6


#  Count mers in reads
#
if [ ! -e $asm-ms$ms-clr-frags.mcidx ] ; then
  bin/dumpFragStoreAsFasta -frg $dir/$asm.frgStore | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-clr-frags -threads 4 -memory $mem -v
fi

if [ ! -e $asm-ms$ms-all-frags.mcidx ] ; then
  bin/dumpFragStoreAsFasta -allbases -allfrags -frg $dir/$asm.frgStore | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-all-frags -threads 4 -memory $mem -v
fi

echo Finding contigs.

if [ ! -e $asm.normalcontigs.fasta ] ; then
  bin/asmOutputContigsFasta    < $dir/9-terminator/$asm.asm > $asm.normalcontigs.fasta &
fi
if [ ! -e $asm.degeneratecontigs.fasta ] ; then
  bin/asmOutputContigsFasta -D < $dir/9-terminator/$asm.asm > $asm.degeneratecontigs.fasta &
fi
if [ ! -e $asm.allcontigs.fasta ] ; then
  bin/asmOutputContigsFasta -d < $dir/9-terminator/$asm.asm > $asm.allcontigs.fasta &
fi

#  Count mers in contigs
#
if [ ! -e $asm-ms$ms-normal-contigs.mcidx ] ; then
  meryl -B -C -m $ms -s $asm.normalcontigs.fasta -o $asm-ms$ms-normal-contigs -threads 4 -segments 4 -v &
fi
if [ ! -e $asm-ms$ms-degenerate-contigs.mcidx ] ; then
  meryl -B -C -m $ms -s $asm.degeneratecontigs.fasta -o $asm-ms$ms-degenerate-contigs -threads 4 -segments 4 -v &
fi
if [ ! -e $asm-ms$ms-all-contigs.mcidx ] ; then
  meryl -B -C -m $ms -s $asm.allcontigs.fasta -o $asm-ms$ms-all-contigs -threads 4 -segments 4 -v &
fi

if [ ! -e $asm-ms$ms.asmMerQC ] ; then
  $asmMerQC -af $asm-ms$ms-all-frags \
            -tf $asm-ms$ms-clr-frags \
            -co $asm-ms$ms-normal-contigs \
            -ac $asm-ms$ms-all-contigs \
            -dc $asm-ms$ms-degenerate-contigs \
  > $asm-ms$ms.asmMerQC &
fi

echo Finding badmers.

if [ ! -e $asm-ms$ms-allfrags-normalcontigs.badmers.asmMerQC ] ; then
  $asmMerQC -af $asm-ms$ms-all-frags \
            -co $asm-ms$ms-normal-contigs \
            -dump $asm-ms$ms-allfrags-normalcontigs.badmers \
  > $asm-ms$ms-allfrags-normalcontigs.badmers.asmMerQC &
fi
if [ ! -e $asm-ms$ms-allfrags-allcontigs.badmers.asmMerQC ] ; then
  $asmMerQC -af $asm-ms$ms-all-frags \
            -ac $asm-ms$ms-all-contigs \
            -dump $asm-ms$ms-allfrags-allcontigs.badmers \
  > $asm-ms$ms-allfrags-allcontigs.badmers.asmMerQC &
fi
if [ ! -e $asm-ms$ms-allfrags-degeneratecontigs.badmers.asmMerQC ] ; then
  $asmMerQC -af $asm-ms$ms-all-frags \
            -dc $asm-ms$ms-degenerate-contigs \
            -dump $asm-ms$ms-allfrags-degeneratecontigs.badmers \
  > $asm-ms$ms-allfrags-degeneratecontigs.badmers.asmMerQC &
fi

if [ ! -e $asm-ms$ms-clrfrags-normalcontigs.badmers.asmMerQC ] ; then
  $asmMerQC -tf $asm-ms$ms-clr-frags \
            -co $asm-ms$ms-normal-contigs \
            -dump $asm-ms$ms-clrfrags-normalcontigs.badmers \
  > $asm-ms$ms-clrfrags-normalcontigs.badmers.asmMerQC &
fi
if [ ! -e $asm-ms$ms-clrfrags-allcontigs.badmers.asmMerQC ] ; then
  $asmMerQC -tf $asm-ms$ms-clr-frags \
            -ac $asm-ms$ms-all-contigs \
            -dump $asm-ms$ms-clrfrags-allcontigs.badmers \
  > $asm-ms$ms-clrfrags-allcontigs.badmers.asmMerQC &
fi
if [ ! -e $asm-ms$ms-clrfrags-degeneratecontigs.badmers.asmMerQC ] ; then
  $asmMerQC -tf $asm-ms$ms-clr-frags \
            -dc $asm-ms$ms-degenerate-contigs \
            -dump $asm-ms$ms-clrfrags-degeneratecontigs.badmers \
  > $asm-ms$ms-clrfrags-degeneratecontigs.badmers.asmMerQC &
fi

echo Mapping.

if [ ! -e $asm-ms$ms-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.normalcontigs.fasta \
  > $asm-ms$ms-allfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers &
fi
if [ ! -e $asm-ms$ms-allfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-allfrags-allcontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.allcontigs.fasta \
  > $asm-ms$ms-allfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers &
fi
if [ ! -e $asm-ms$ms-allfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-allfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.degeneratecontigs.fasta \
  > $asm-ms$ms-allfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.badmers &
fi

if [ ! -e $asm-ms$ms-clrfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-clrfrags-normalcontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.normalcontigs.fasta \
  > $asm-ms$ms-clrfrags-normalcontigs.badmers.0.singlecontig.zerofrag.badmers &
fi
if [ ! -e $asm-ms$ms-clrfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-clrfrags-allcontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.allcontigs.fasta \
  > $asm-ms$ms-clrfrags-allcontigs.badmers.0.singlecontig.zerofrag.badmers &
fi
if [ ! -e $asm-ms$ms-clrfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.badmers ] ; then
  $mapMers -m 22 \
           -mers $asm-ms$ms-clrfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.fasta \
           -seq $asm.degeneratecontigs.fasta \
  > $asm-ms$ms-clrfrags-degeneratecontigs.badmers.0.singlecontig.zerofrag.badmers &
fi

if [ ! -e $asm-ms$ms-allfrags-normalcontigs.badmers.5.all.badmers ] ; then
  cat $asm-ms$ms-allfrags-normalcontigs.badmers.[01].*.fasta > $asm-ms$ms-allfrags-normalcontigs.badmers.5.allzero.fasta
  $mapMers -m 22 \
           -mers $asm-ms$ms-allfrags-normalcontigs.badmers.5.allzero.fasta \
           -seq $asm.normalcontigs.fasta \
  > $asm-ms$ms-allfrags-normalcontigs.badmers.5.allzero.badmers &
fi

date
