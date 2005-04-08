#!/bin/sh
#
#  Mostly for SGE, since it cannot accept perl scripts.  If this
#  fails, Try submitting with "-S /bin/sh" or fixing your install of
#  SGE to use /bin/sh instead of /bin/csh

if [ -r /usr/local/sge/default/common/settings.sh ]; then
  . /usr/local/sge/default/common/settings.sh
fi

/usr/bin/perl "$@"
