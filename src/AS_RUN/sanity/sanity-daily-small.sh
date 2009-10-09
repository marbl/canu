#!/bin/sh

date=$1

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org"
email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org"

#  Nightly small
#
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

#  Resubmit ourself for tomorrow.
