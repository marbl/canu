#!/bin/sh

#  (re)Load the sge config.
. $SGE_ROOT/$SGE_CELL/common/settings.sh

#  Needed for old checkouts.  Newer checkouts default to LOCAL.
#With SITE=JCVI, every EUID is different, so every QC report is different.
#With default behavior, EUIDs always start at the same number, so QC reports can be identical.
#We want default bahavior so the summary email reports no differences.
#export SITE_NAME=JCVI

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
grid=$2

if [ x$date != x ] ; then
  echo "SANITY BEGINS for $date at `date`"
fi


#  Remove old versions
perl sanity-purge-old.pl purge
rm -rf DEL

#  Update the repository.
perl sanity.pl rsync


#  Checkout the latest version.
perl sanity.pl checkout $date


#  Figure out what was checked out.
if [ x$date = x ] ; then
  date=`ls -1d 20??-??-??-???? | tail -n 1`
fi


#  Build it.
perl sanity.pl build $date

#  Let the user pick one to run
if [ x$grid = x ] ; then
    echo "$date checked out and compiled.  Run some of:"
    echo "  sh sanity-daily-test.sh   $date"
    echo "  sh sanity-daily-pging.sh  $date"

else
    nextofft=604800  # one week
    nextofft=86400   # one day
    nextofft=14400   # four hours
    nextofft=7200    # two hours
    nextofft=21600   # six hours
    nextofft=43200   # twelve hours
    nextdate=`perl sanity-get-next-date.pl $date $nextofft next`
    nexthold=`perl sanity-get-next-date.pl $date $nextofft hold`
    nextname=`perl sanity-get-next-date.pl $date $nextofft name`

    #  Which tests should we run?  For now, just the simple p.ging

    sh sanity-daily-pging.sh  $date

    #if [ `date +%u` = 6] ; then
    #    sh sanity-weekly-dros.sh  $date
    #    sh sanity-weekly-moore.sh $date
    #fi

    echo "SUBMIT for $nextdate ($nexthold)"

    echo \
    qsub -cwd -j y -o $nextdate.err -P 334007 -A assembly-nightly -N CAsnty$nextname -a $nexthold -b n sanity.sh $nextdate grid

    qsub -cwd -j y -o $nextdate.err -P 334007 -A assembly-nightly -N CAsnty$nextname -a $nexthold -b n sanity.sh $nextdate grid
fi
