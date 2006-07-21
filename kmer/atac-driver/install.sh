#!/bin/sh

#  Attempt to install ATAC/A2Amapper in a specified location ($dir).
#  It will make $dir/bin and $dir/lib.

if [ $# -eq 0 ] ; then
  echo "usage: $0 install-directory"
  exit
fi

installdir=$1

mkdir -p $installdir
mkdir -p $installdir/bin
mkdir -p $installdir/lib

#  Install binaries
#
binaries="../leaff/leaff ../meryl/meryl ../libkmer/existDB ../seatac/seatac"
binaries="$binaries briatac.pl matchExtender/matchExtender clumpMaker/clumpMaker"

for i in $binaries ; do
  #echo $i
  cp -fp $i $installdir/bin/
  #chmod 755 $installdir/bin/
done

#  Install seatac's filter libraries
#
cp -p ../seatac/filter-heavychains.so   $installdir/lib
cp -p ../seatac/filter-nop.so           $installdir/lib

#  Install the chainer
#
cp -fp chainer/halignmodule.so                $installdir/lib
cp -fp chainer/localAlignerInterfacemodule.so $installdir/lib
cp -fp chainer/python/*py                     $installdir/lib
cp -fp chainer/python/AtacDriver.py           $installdir/bin


#  Make sure briatac.pl runs.
perl $installdir/bin/briatac.pl -justtestingifitworks || echo "WARNING:  briatac.pl failed to install."

#  Make sure the chainer runs.
#
export PYTHONPATH=$installdir/lib
python $installdir/bin/AtacDriver.py justtestingifitworks || echo "WARNING:  Chainer failed to install."


echo "Unless you see warnings immediately before this line, ATAC/A2Amapper"
echo "is now installed in $installdir/bin/".

exit 0
