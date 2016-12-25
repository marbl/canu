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

$dnadiff -p $name.dnadiff \
  /data/references/francisella_tularensis_99a-2628-GCA_000478845.1.fasta \
  $name.contigs.fasta

$mummerplot --fat -t png -p $name.dnadiff $name.dnadiff.delta

rm -f $name.dnadiff.1coords
rm -f $name.dnadiff.1delta
rm -f $name.dnadiff.mcoords
rm -f $name.dnadiff.mdelta
rm -f $name.dnadiff.qdiff
rm -f $name.dnadiff.rdiff

rm -f $name.dnadiff.filter
rm -f $name.dnadiff.rplot
rm -f $name.dnadiff.fplot
rm -f $name.dnadiff.gp
