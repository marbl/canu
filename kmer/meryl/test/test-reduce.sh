#!/bin/sh

#  Build the initial table
./meryl -B -f -m 20 -s test/test-seq1.fasta -o 1

#  Dump the initial table as fasta
./meryl -Dt -n 0 -s 1 > 2.reduce.fasta

#  Build a new table on the dumped fasta
./meryl -B -f -m 20 -s 2.reduce.fasta -o 2

#  Remove one copy of each mer
./meryl -M sub -s 1 -s 2 -o 3


#  Dump / test the mask
#./meryl -M max -s 1 -s 2 -mask 3 -o 4
#./meryl -Dt -s 4

#rm -f 1.mcidx 1.mcdat 1.merStream
#rm -f 2.reduce.fasta
#rm -f 2.mcidx 2.mcdat 2.merStream
#rm -f 3.mcidx 3.mcdat
