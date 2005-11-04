#!/usr/local/bin/bash


if [ "$1" == "" ] || [ ! -d $1 ] || [ "$2" == "" ]; then
  echo "Usage: makeTarball.sh  rootDir  tarballFile"
  return
fi

cwd=`pwd`

cd $1

echo "Deleting CVS subdirectories"
for cvsDir in `find . -type d -name CVS -print`; do
  rm -rf $cvsDir
done

if [ -d CVSROOT ]; then
  echo "Deleting CVSROOT directory"
  rm -rf CVSROOT
fi

echo "Deleting #-prefixed and ~-suffixed files"
for backFile in `find . -type f -name "*~" -print`; do
  echo "deleting $backFile"
  rm $backFile
done
for backFile in `find . -type f -name "#*" -print`; do
  echo "deleting $backFile"
  rm $backFile
done

# get rid of binaries and other build directories
# find should return ./PLATFORM/dep
for depDir in `find . -type d -name dep -print -maxdepth 2`; do
  platformDir=${depDir%/*}
  platformDir=${platformDir##*/}
  echo "Deleting $platformDir"
  rm -rf $platformDir
done

# tar -czvf
cd ..
tar -czvf $2 wgs-assembler

cd $cwd
