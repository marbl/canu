#!/bin/sh

#  Build the initial table
./meryl -B -H 20 -f -m 20 -t 20 -s test/test-seq1.fasta -o 1

#  Dump the initial table as fasta
./meryl -Dt -n 0 -s 1 > 2.reduce.fasta

#  Build a new table on the dumped fasta
./meryl -B -H 20 -f -m 20 -t 20 -s 2.reduce.fasta -o 2

#  Remove one copy of each mer
#./meryl -R -s 1 -s 2 -o 3
./meryl -H 20 -M sub -s 1 -s 2 -o 3

echo 

#  Dump the reduced version
#echo "1"
#./meryl -Dt -n 0 -s 1
#echo 
#
#echo "2"
#./meryl -Dt -n 0 -s 2
#echo 
#
#echo "3"
#./meryl -Dt -n 0 -s 3
#echo 

#  Dump / test the mask
./meryl -M max -s 1 -s 2 -mask 3 -o 4
#./meryl -M min -s 1 -s 2 -mask 3

./meryl -Dt -s 4
