#!/bin/sh

date=$1

email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org gsutton@jcvi.org jjohnson@jcvi.org"
email="bwalenz@jcvi.org jmiller@jcvi.org skoren@jcvi.org"
email="bwalenz@jcvi.org"

#  Weekly, every Friday
#    Can't run erecta here too - it is run in daily-large.
#
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

#  Resubmit ourself for tomorrow.
