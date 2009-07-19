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
grid=$2

if [ x$date != x ] ; then
  echo "SANITY BEGINS for $date at `date`"
fi

#  Update the repository.
#perl sanity.pl rsync


#  Checkout the latest version.
#perl sanity.pl checkout $date


#  Figure out what was checked out.
if [ x$date = x ] ; then
  date=`ls -1d 200?-??-??-???? | tail -n 1`
fi


#  Build it.
#
#perl sanity.pl build $date


#  Let the user pick one to run
#
if [ x$grid = x ] ; then
    echo "$date checked out and compiled.  Run some of:"
    echo "  sh sanity-daily-small.sh  $date"
    echo "  sh sanity-daily-large.sh  $date"
    echo "  sh sanity-weekly-dros.sh  $date"
    echo "  sh sanity-weekly-moore.sh $date"

else
    nextofft=1800
    nextdate=`perl sanity-get-next-date.pl $date $nextofft next`
    nexthold=`perl sanity-get-next-date.pl $date $nextofft hold`
    nextname=`perl sanity-get-next-date.pl $date $nextofft name`

    #sh sanity-daily-test.sh   $date

    #sh sanity-daily-small.sh  $date
    #sh sanity-daily-large.sh  $date

    #if [ `date +%u` = 6] ; then
    #    sh sanity-weekly-dros.sh  $date
    #    sh sanity-weekly-moore.sh $date
    #fi

    echo "SUBMIT for $nextdate ($nexthold)"

    echo \
    qsub -cwd -j y -o $nextdate.err -A assembly-nightly -N CAsnty$nextname -l fast -a $nexthold -b n sanity.sh $nextdate grid

    qsub -cwd -j y -o $nextdate.err -A assembly-nightly -N CAsnty$nextname -l fast -a $nexthold -b n sanity.sh $nextdate grid
fi
