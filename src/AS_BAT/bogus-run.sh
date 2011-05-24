#!/bin/sh

#  Build ideal unitigs given a fasta or frg file.
#  Align ideal unitigs back to the reference with nucmer.
#  Generate a mummerplot.
#
#  ASSUMES all programs are in your path.
#    mummer:        nucmer, mummerplot
#    kmer:          snapper2, convertToExtent
#    wgs-assembler: bogus


#  Edit reference to make defline have only one work.  This is necessary because snapper reports the
#  whole line, but bogus truncates to the first word.

FAS="../FRAGS/porphyromonas_gingivalis_w83.flx.3200bp.0900bp.FJRUAFO0.fasta"
REF="../AE015924.fasta"

if [ $# -gt 0 ] ; then
  FAS=$1
  shift
fi

if [ $# -gt 0 ] ; then
  REF=$1
  shift
fi

OUT=`echo $FAS | tr '/' ' ' | awk '{ print $NF}' | sed s/.fasta//`

if [ ! -e $FAS ] ; then
  echo "Failed to find FASTA $FAS"
  exit
fi

if [ ! -e $REF ] ; then
  echo "Failed to find REFERENCE $REF"
  exit
fi

if [ `cat $REF | wc -l` != 2 ] ; then
  echo "REFERENCE is multi-line FASTA, sequence must be on one line."
  exit
fi

if [ ! -e $OUT.snapper ] ; then
  snapper2 \
    -queries $FAS \
    -genomic $REF \
    -minmatchidentity 94 -minmatchcoverage 1 -verbose -mersize 16 \
    -output $OUT.snapper
fi

if [ ! -e $OUT.snapper.extent ] ; then
  convertToExtent \
   < $OUT.snapper \
   | grep -v ^cDNAid \
   | sort -k7n \
   > $OUT.snapper.extent
fi

if [ ! -e $OUT.ideal.fasta ] ; then
  bogus \
    -snapper   $OUT.snapper.extent \
    -reference $REF \
    -output    $OUT.ideal
fi

if [ ! -e $OUT.ideal.fasta ] ; then
  echo "No fasta output from bogus, not mapping to reference."
  exit
fi

if [ ! -e $OUT.ideal.delta ] ; then
  nucmer --maxmatch --coords -p $OUT.ideal \
    $REF \
    $OUT.ideal.fasta
fi

if [ ! -e $OUT.ideal.png ] ; then
  mummerplot --layout --filter -p $OUT.ideal -t png $OUT.ideal.delta
fi
