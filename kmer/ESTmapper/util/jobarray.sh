#!/bin/sh
#
#  SGE can't run perl directly, so this wrapper just launches the perl
#  version.  Further, SGE copies this script somewhere else, so we
#  need to pass in the real script we want to run.

if [ -r /usr/local/sge/default/common/settings.sh ]; then
  . /usr/local/sge/default/common/settings.sh
fi

/usr/bin/perl "$@"
