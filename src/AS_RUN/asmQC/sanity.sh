#!/bin/sh

#
#  Master controller of nightly sanity checks.  Optional date on command line.
#

date=$1

#  Update the repository.
perl sanity.pl rsync


#  Checkout the latest version.
perl sanity.pl checkout $date


#  Figure out what was checked out.
if [ x$date = x ] ; then
  date=`ls -1d 200?-??-??-???? | tail -1`
fi


#  Build it.
#
perl sanity.pl build $date

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org"

#  Nightly small
#
if [ x$2 != x"weekly" ] ; then
perl sanity.pl assemble $date microbes \
  SPECFILES/e.coli.higherror.spec \
  SPECFILES/e.coli.merbog.spec \
  SPECFILES/e.litoralis.utg.deep.mix.spec \
  SPECFILES/p.cnpt3.merbog.spec \
  SPECFILES/p.cnpt3.ovlutg.spec \
  SPECFILES/p.cnpt3.shred.spec \
  SPECFILES/p.gingivalis.merbog.spec \
  SPECFILES/p.gingivalis.merbogdeep.spec \
  SPECFILES/p.gingivalis.ovlutg.spec \
  SPECFILES/p.gingivalis.shreddeep.spec \
  $email
fi

#  Nightly big
#
if [ x$2 != x"weekly" ] ; then
perl sanity.pl assemble $date eukaryotes \
  SPECFILES/p.humanus.utg.higherror.spec \
  SPECFILES/p.ultimum.merbog.spec \
  SPECFILES/d.erecta.spec \
  $email
fi

#  Weekly, every Friday
#    Can't run erecta here too!
#
if [ `date +%w` -eq 5 -o x$2 = x"weekly" ] ; then
perl sanity.pl assemble $date drosophila-weekly \
  SPECFILES/d.ananassae.spec \
  SPECFILES/d.grimshawi.spec \
  SPECFILES/d.mojavensis.spec \
  SPECFILES/d.persimilis.spec \
  SPECFILES/d.pseudoobscura.spec \
  SPECFILES/d.sechellia.spec \
  SPECFILES/d.simulans.spec \
  SPECFILES/d.virilis.spec \
  SPECFILES/d.willistoni.spec \
  SPECFILES/d.yakuba.spec \
  $email
fi

#  Resubmit ourself for tomorrow.
