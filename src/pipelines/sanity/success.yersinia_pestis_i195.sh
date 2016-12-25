#!/bin/sh

name=$1

if [ x$name = x ] ; then
  echo "usage: $0 assembly-prefix"
  exit 1
fi

if [ -e "/work/software/MUMmer3.23/mummerplot" ] ; then  PATH=$PATH:/work/software/MUMmer3.23 ; fi

dnadiff=`which dnadiff`
mummerplot=`which mummerplot`

if [ ! -e $dnadiff -o ! -e $mummerplot ] ; then
  echo "Didn't find dnadiff or mummerplot in your path."
  exit 1
fi

#  195/P, but checking against CO92 is reasonable (a few inversions exist)
#    https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=ACNR01
#
#  yersinia_pestis_CO92_AL109969.fasta - pPCP1
#  yersinia_pestis_CO92_AL117189.fasta - pCD1
#  yersinia_pestis_CO92_AL117211.fasta - pMT1
#  yersinia_pestis_CO92_AL590842.fasta - genome
#
#  yersinia_pestis_India195.ACNR01.1.contigs.fasta

$dnadiff -p $name.co92.dnadiff \
  /data/references/yersinia_pestis_CO92_AL590842.fasta \
  $name.contigs.fasta

$mummerplot --filter --fat -t png -p $name.co92.dnadiff $name.co92.dnadiff.delta

rm -f $name.co92.dnadiff.1coords
rm -f $name.co92.dnadiff.1delta
rm -f $name.co92.dnadiff.mcoords
rm -f $name.co92.dnadiff.mdelta
rm -f $name.co92.dnadiff.qdiff
rm -f $name.co92.dnadiff.rdiff

rm -f $name.co92.dnadiff.filter
rm -f $name.co92.dnadiff.rplot
rm -f $name.co92.dnadiff.fplot
rm -f $name.co92.dnadiff.gp


$dnadiff -p $name.i195.dnadiff \
  /data/references/yersinia_pestis_India195.ACNR01.1.contigs.fasta \
  $name.contigs.fasta

$mummerplot --filter --fat -t png -p $name.i195.dnadiff $name.i195.dnadiff.delta

rm -f $name.i195.dnadiff.1coords
rm -f $name.i195.dnadiff.1delta
rm -f $name.i195.dnadiff.mcoords
rm -f $name.i195.dnadiff.mdelta
rm -f $name.i195.dnadiff.qdiff
rm -f $name.i195.dnadiff.rdiff

rm -f $name.i195.dnadiff.filter
rm -f $name.i195.dnadiff.rplot
rm -f $name.i195.dnadiff.fplot
rm -f $name.i195.dnadiff.gp
