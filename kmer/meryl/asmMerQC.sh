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

asmMerQC=/bioinfo/assembly/walenz/src/genomics/meryl/asmMerQC

dir=/scratch/drosnightly
asm=willi

dir=/project/huref6/assembly/h6
asm=h6

date

#  Count mers in reads
#
if [ ! -e $asm-ms$ms-clr-frags.mcidx ] ; then
  bin/dumpFragStoreAsFasta -frg $dir/$asm.frgStore | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-clr-frags -threads 4 -memory $mem -v
fi

date

if [ ! -e $asm-ms$ms-all-frags.mcidx ] ; then
  bin/dumpFragStoreAsFasta -allbases -allfrags -frg $dir/$asm.frgStore | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-all-frags -threads 4 -memory $mem -v
fi

date

#  Count mers in contigs
#
if [ ! -e $asm-ms$ms-normal-contigs.mcidx ] ; then
  bin/asmOutputContigsFasta < $dir/9-terminator/$asm.asm | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-normal-contigs -threads 4 -segments 4 -v
fi

date

if [ ! -e $asm-ms$ms-degenerate-contigs.mcidx ] ; then
  bin/asmOutputContigsFasta -D < $dir/9-terminator/$asm.asm | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-degenerate-contigs -threads 4 -segments 4 -v
fi

date

if [ ! -e $asm-ms$ms-all-contigs.mcidx ] ; then
  bin/asmOutputContigsFasta -d < $dir/9-terminator/$asm.asm | \
  meryl -B -C -m $ms -s - -o $asm-ms$ms-all-contigs -threads 4 -segments 4 -v
fi

date

echo $asmMerQC -af $asm-ms$ms-all-frags \
          -tf $asm-ms$ms-clr-frags \
          -co $asm-ms$ms-normal-contigs \
          -ac $asm-ms$ms-all-contigs \
          -dc $asm-ms$ms-degenerate-contigs \
\> $asm-ms$ms.asmMerQC

date
