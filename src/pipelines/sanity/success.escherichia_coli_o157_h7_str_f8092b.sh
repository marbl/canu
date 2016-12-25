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

#$dnadiff -p $name-EC4115.dnadiff \
#  /data/references/escherichia_coli_o157_h7_EC4115_uid59091__NC_011353.fasta \
#  $name.contigs.fasta
#$mummerplot --fat -t png -p $name-EC4115.dnadiff $name-EC4115.dnadiff.delta
#
#$dnadiff -p $name-EDL933.dnadiff \
#  /data/references/escherichia_coli_o157_h7_EDL933_uid57831__NC_002655.fasta \
#  $name.contigs.fasta
#$mummerplot --fat -t png -p $name-EDL933.dnadiff $name-EDL933.dnadiff.delta
#
#$dnadiff -p $name-Sakai.dnadiff \
#  /data/references/escherichia_coli_o157_h7_Sakai_uid57781__NC_002695.fasta \
#  $name.contigs.fasta
#$mummerplot --fat -t png -p $name-Sakai.dnadiff $name-Sakai.dnadiff.delta
#
#$dnadiff -p $name-TW14359.dnadiff \
#  /data/references/escherichia_coli_o157_h7_TW14359_uid59235__NC_013008.fasta \
#  $name.contigs.fasta
#$mummerplot --fat -t png -p $name-TW14359.dnadiff $name-TW14359.dnadiff.delta

$dnadiff -p $name-F8092B.dnadiff \
  /data/references/escherichia_coli_o157_h7_F8092B__GCA_000513035.1.fasta \
  $name.contigs.fasta
$mummerplot --fat -t png -p $name-F8092B.dnadiff $name-F8092B.dnadiff.delta

rm -f $name-*.dnadiff.1coords
rm -f $name-*.dnadiff.1delta
rm -f $name-*.dnadiff.mcoords
rm -f $name-*.dnadiff.mdelta
rm -f $name-*.dnadiff.qdiff
rm -f $name-*.dnadiff.rdiff

rm -f $name-*.dnadiff.filter
rm -f $name-*.dnadiff.rplot
rm -f $name-*.dnadiff.fplot
rm -f $name-*.dnadiff.gp

