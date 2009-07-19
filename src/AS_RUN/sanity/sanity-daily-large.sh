#!/bin/sh

date=$1

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org"
email="bwalenz@jcvi.org"

#  Nightly big
#
perl sanity.pl assemble $date eukaryotes \
  SPECFILES/p.humanus.utg.higherror.spec \
  SPECFILES/p.ultimum.merbog.spec \
  SPECFILES/d.erecta.spec \
  $email

#  Resubmit ourself for tomorrow.
