#!/bin/sh

#  Attempt to install ESTmapper in a specified location ($dir).
#  It will make $dir/bin and $dir/lib.

if [ $# -eq 0 ] ; then
  echo "usage: $0 install-directory"
  exit
fi

installdir=$1

mkdir -p $installdir
mkdir -p $installdir/bin
mkdir -p $installdir/bin/util
mkdir -p $installdir/lib

#  Install binaries
#
binaries="mergeCounts terminate"
binaries="$binaries ../leaff/leaff ../sim4db/sim4db ../meryl/meryl ../libkmer/positionDB"
binaries="$binaries ../seagen/seagen ../seagen/filterEST ../seagen/filterMRNA ../seagen/filterNULL ../seagen/sortHits"
binaries="$binaries ../sim4dbutils/cleanPolishes ../sim4dbutils/sortPolishes ../sim4dbutils/parseSNP ../sim4dbutils/filterPolishes ../sim4dbutils/pickBestPolish"

for i in $binaries ; do
  cp -fp $i $installdir/bin/
done

#  Install ESTmapper itself
#
cp -fp ESTmapper.pl             $installdir/bin/
cp -fp configureESTmapper.pl    $installdir/bin/
cp -p ../scripts/scheduler.pm   $installdir/bin/util/

chmod 755 $installdir/bin/ESTmapper.pl
chmod 755 $installdir/bin/configureESTmapper.pl

#  Make sure ESTmapper.pl runs.
perl $installdir/bin/ESTmapper.pl          -justtestingifitworks || echo "WARNING:  ESTmapper.pl failed to install."
perl $installdir/bin/configureESTmapper.pl -justtestingifitworks || echo "WARNING:  configureESTmapper.pl failed to install."

echo "Unless you see warnings immediately before this line, ESTmapper"
echo "is now installed in $installdir/bin/".

exit 0
