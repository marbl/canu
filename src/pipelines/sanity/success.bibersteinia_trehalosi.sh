#!/bin/sh

name=$1

dnadiff -p $name.dnadiff \
  /data/references/bibersteinia_trehalosi/bibersteinia_trehalosi_USDA-ARS-USMARC-192-CP003745.1.fasta \
  $name.contigs.fasta

mummerplot --fat -t png -p $name.dnadiff $name.dnadiff.delta

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
