#!/bin/sh
#
#  Green  - read is retained for assembly.
#  Purple - read is discarded before assembly starts.
#
#  Y axis is read length.
#  X axis is the score, sorted from highest to lowest.
#
#  The algorithm will keep reads from highest to lowest score until the coverage threshold is met.
#
#  Reads are scored with 
#        readLen[ii].score  = mtctx.mtRandomRealOpen() * pow(len, lengthBias);
#
#  then sorted by decreasing score.  Reads with high score are retained,
#  up until we hit the desired coverage, then all (lower scoring) reads
#  are deleted.
#
#  In the plots, this can be interpreted as a line scanning down
#  the Y axis.  Reads above the line are retained, reads below
#  are discarded.
#
#  The actual value of the score doesn't matter, just the ordering
#  of lengths.
#

rm -rf test-c.seqStore
sqStoreCreate -bias -10 -o test-c.seqStore -genomesize 135000000 -coverage 1 -pacbio-hifi LIB /data/reads/drosophila/F1-hifi-11k/1x-split/x00[01].fastq.xz
grep    REMOVED test-c.seqStore/filteredReads | sort -k2n > remove
grep -v REMOVED test-c.seqStore/filteredReads | sort -k2n > retain

gnuplot <<EOF
set title  "readSamplingBias=-10.0"
set xlabel "Read Score Ordinal"
set ylabel "Read Length"
set pointsize 0.5

set terminal png size 320,200
set output "bias=-10.png"
plot [] [0:20000] 'remove' using 2:3 with points title '', 'retain' using 2:3 with points title ''
EOF


rm -rf test-c.seqStore
sqStoreCreate -bias -1 -o test-c.seqStore -genomesize 135000000 -coverage 1 -pacbio-hifi LIB /data/reads/drosophila/F1-hifi-11k/1x-split/x00[01].fastq.xz
grep    REMOVED test-c.seqStore/filteredReads | sort -k2n > remove
grep -v REMOVED test-c.seqStore/filteredReads | sort -k2n > retain

gnuplot <<EOF
set title  "readSamplingBias=-1.0"
set xlabel "Read Score Ordinal"
set ylabel "Read Length"
set pointsize 0.5

set terminal png size 320,200
set output "bias=-01.png"
plot [] [0:20000] 'remove' using 2:3 with points title '', 'retain' using 2:3 with points title ''
EOF


rm -rf test-c.seqStore
sqStoreCreate -bias 0 -o test-c.seqStore -genomesize 135000000 -coverage 1 -pacbio-hifi LIB /data/reads/drosophila/F1-hifi-11k/1x-split/x00[01].fastq.xz
grep    REMOVED test-c.seqStore/filteredReads | sort -k2n > remove
grep -v REMOVED test-c.seqStore/filteredReads | sort -k2n > retain

gnuplot <<EOF
set title  "readSamplingBias=0.0"
set xlabel "Read Score Ordinal"
set ylabel "Read Length"
set pointsize 0.5

set terminal png size 320,200
set output "bias=+00.png"
plot [] [0:20000] 'remove' using 2:3 with points title '', 'retain' using 2:3 with points title ''
EOF


rm -rf test-c.seqStore
sqStoreCreate -bias 1 -o test-c.seqStore -genomesize 135000000 -coverage 1 -pacbio-hifi LIB /data/reads/drosophila/F1-hifi-11k/1x-split/x00[01].fastq.xz
grep    REMOVED test-c.seqStore/filteredReads | sort -k2n > remove
grep -v REMOVED test-c.seqStore/filteredReads | sort -k2n > retain

gnuplot <<EOF
set title  "readSamplingBias=1.0"
set xlabel "Read Score Ordinal"
set ylabel "Read Length"
set pointsize 0.5

set terminal png size 320,200
set output "bias=+01.png"
plot [] [0:20000] 'remove' using 2:3 with points title '', 'retain' using 2:3 with points title ''
EOF


rm -rf test-c.seqStore
sqStoreCreate -bias 10 -o test-c.seqStore -genomesize 135000000 -coverage 1 -pacbio-hifi LIB /data/reads/drosophila/F1-hifi-11k/1x-split/x00[01].fastq.xz
grep    REMOVED test-c.seqStore/filteredReads | sort -k2n > remove
grep -v REMOVED test-c.seqStore/filteredReads | sort -k2n > retain

gnuplot <<EOF
set title  "readSamplingBias=10.0"
set xlabel "Read Score Ordinal"
set ylabel "Read Length"
set pointsize 0.5

set terminal png size 320,200
set output "bias=+10.png"
plot [] [0:20000] 'remove' using 2:3 with points title '', 'retain' using 2:3 with points title ''
EOF

exit
