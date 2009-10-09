#!/bin/sh

date=$1

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org"
email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org"

#  Nightly small
#
perl sanity.pl assemble $date test \
  SPECFILES/p.cnpt3.ovlutg.spec \
  SPECFILES/p.cnpt3.shred.spec \
  SPECFILES/p.gingivalis.ovlutg.spec \
  $email

#  Resubmit ourself for tomorrow.
