#!/bin/sh

date=$1

email="bwalenz@jcvi.org"

#  Nightly small
#
perl sanity.pl assemble $date microbes \
  SPECFILES/p.cnpt3.ovlutg.spec \
  $email

