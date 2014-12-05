#!/bin/sh

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
