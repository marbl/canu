#!/bin/sh

name=$1

/work/software/MUMmer3.23/dnadiff -p $name.dnadiff \
  /data/references/escherichia_coli_k12_substr__MG1655_uid57779__NC_000913.3.fasta \
  $name.contigs.fasta

/work/software/MUMmer3.23/mummerplot --fat -t png -p $name.dnadiff $name.dnadiff.delta

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
