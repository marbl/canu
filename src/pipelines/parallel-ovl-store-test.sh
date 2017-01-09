#!/bin/sh

bin="/work/canu/FreeBSD-amd64/bin"

if [ ! -e test.gkpStore ] ; then
  echo Didn\'t find \'test.gkpStore\', can\'t make fake overlaps.
  exit
fi

rm -rf *ovb *counts test.ovlStore test.ovlStore.seq

jobs=50
size=100000000

#  Create a bunch of overlap outputs, one for each job.

echo ""
echo "Creating overlaps"
echo ""

for ii in `seq 1 $jobs` ; do
  name=`printf %03d $ii`
  $bin/overlapImport -G test.gkpStore -o $name.ovb -random $size
done

#  Configure.

$bin/ovStoreBuild -G test.gkpStore -O test.ovlStore -M 0.26 -config config *ovb > config.err 2>&1

buckets=`grep Will config.err | grep buckets | awk '{ print $9 }'`

echo ""
echo "Using $buckets buckets for sorting"
echo ""

#  Bucketize each job.

for ii in `seq 1 $jobs` ; do
  name=`printf %03d $ii`
  $bin/ovStoreBucketizer -G test.gkpStore -O test.ovlStore -C config -i $name.ovb -job $ii
done

#  Sort each bucket.

echo ""
echo "Sorting each bucket"
echo ""

for ii in `seq 1 $buckets` ; do
  echo ""
  echo "Sorting bucket $ii"
  echo ""
  $bin/ovStoreSorter -deletelate -G test.gkpStore -O test.ovlStore -F $jobs -job $ii $jobs
done

#  And build the index

$bin/ovStoreIndexer -O test.ovlStore -F $buckets

#  And a sequential store?

exit

echo ""
echo "Building sequential store"
echo ""

$bin/ovStoreBuild -G test.gkpStore -O test.ovlStore.seq *ovb > test.ovlStore.seq.err 2>&1
