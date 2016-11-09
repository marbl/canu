#!/bin/sh

#  (re)Load the sge config.
. $SGE_ROOT/$SGE_CELL/common/settings.sh

#
#  Master controller of nightly sanity checks.  Optional date on command line.
#
#  To start the next nightly, launch with something like:
#    qsub -cwd -j y -o init.err -A assembly-nightly -l fast -a 200907180001.00 -b n sanity.sh 2009-07-18-0001 grid
#
#  To start a previous nightly, the sge hold isn't needed, well, neither is sge for that matter.
#  Just run sanity.sh with the date/time you want to start at.
#
#  CAREFUL!  Submitting something in the past will automagically submit every nightly since that date!
#

date=$1
bins=/work/canu/src/pipelines/sanity
spec=/work/canu/src/pipelines/sanity

if [ x$date = x ] ; then
  date=`date +%Y-%m-%d-%H%M`
fi

echo "SANITY BEGINS for $date at `date`"

#  Remove old versions
#perl $bins/sanity-purge-old.pl purge
#rm -rf DEL

perl $bins/sanity.pl fetch                #  Update the local repo.
perl $bins/sanity.pl checkout $date       #  Checkout from that repo.
perl $bins/sanity.pl build    $date       #  Compile.
perl $bins/sanity.pl submit   $date       #  Submit the next iteration.

#  Run small stuff daily.

perl $bins/sanity.pl assemble $date $spec/small1.spec
perl $bins/sanity.pl assemble $date $spec/small2.spec
perl $bins/sanity.pl assemble $date $spec/small3.spec
perl $bins/sanity.pl assemble $date $spec/small4.spec

#  Run big stuff weekly.

#if [ `date +%u` = 6] ; then
#    sh sanity-weekly-dros.sh  $date
#    sh sanity-weekly-moore.sh $date
#fi

