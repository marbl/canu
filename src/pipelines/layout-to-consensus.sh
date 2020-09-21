#!/bin/sh

#  Runs consensus, in parallel, on a ctgStore from bogart.
#
#  Magic values for partitioning:
#    partitionSize    - make partitions this big, relative to the biggest contig
#                       0.10 - each big contig gets its own partition
#    partitionScaling - expect contigs to expand by this much due to homopoly compression
#                       1.50 - 50% expansion of contig length, affects memory estimate
#    partitionTigs    - put this many reads (as fraction of total) in a partition
#                       0.01 - puts small contigs into many partitions

asmpre=asm
cnspre=asm
errorRate=0.01

if [ ! -e ./${cnspre}.ctgStore/partitioning ] ; then
  echo Partitioning reads.

  utgcns -V \
    -S ../../${asmpre}.seqStore \
    -T ./${cnspre}.ctgStore 1 \
    -partition 0.10 1.50 0.01 \
  > ./${cnspre}.ctgStore/partitioning.log 2>&1
fi

parts=`cd ./${cnspre}.ctgStore ; ls partition.???? | sed s/partition.//`

for pp in $parts ; do
  if [ ! -e "./${cnspre}.${pp}.cns" ] ; then
    echo \
    utgcns \
      -R ./${cnspre}.ctgStore/partition.${pp} \
      -T ./${cnspre}.ctgStore 1 \
      -P ${pp} \
      -A ./${cnspre}.${pp}.cns.fasta \
      -O ./${cnspre}.${pp}.cns \
      -maxcoverage 50 \
      -e ${errorRate} \
      -threads 8 \&
  fi
done

exit 0
