#!/bin/sh

version=$1

if [ x$version = x ] ; then
  echo usage: $0 numeric-version
  exit
fi

#  From the tarball

if [ ! -e canu-$version.tar.gz ] ; then
  echo Fetch.
  curl -L -R -o canu-$version.tar.gz https://github.com/marbl/canu/archive/v$version.tar.gz
fi
if [ ! -d canu-$versioin ] ; then
  echo Unpack.
  gzip -dc canu-$version.tar.gz | tar -xf -
fi
cd canu-$version

#  From the repo

#git clone git@github.com:marbl/canu.git
#mv canu canu-$version
#cd canu-$version
#git tag v$version
#git checkout v$version

echo Build MacOS.
cd src
gmake -j 12 > ../Darwin-amd64.out 2>&1
cd ../..

rm -f canu-$version/linux.sh

echo >> canu-$version/linux.sh  \#\!/bin/bash
#echo >> canu-$version/linux.sh  yum install -y git
echo >> canu-$version/linux.sh  cd /build/canu-$version/src
echo >> canu-$version/linux.sh  gmake -j 12 \> ../Linux-amd64.out 2\>\&1
echo >> canu-$version/linux.sh  cd ../..
echo >> canu-$version/linux.sh  rm -rf canu-$version/Darwin-amd64/obj
echo >> canu-$version/linux.sh  rm -rf canu-$version/Linux-amd64/obj
echo >> canu-$version/linux.sh  tar -cf canu-$version.Darwin-amd64.tar canu-$version/README* canu-$version/Darwin-amd64
echo >> canu-$version/linux.sh  tar -cf canu-$version.Linux-amd64.tar  canu-$version/README* canu-$version/Linux-amd64

chmod 755 canu-$version/linux.sh

echo Build Linux and make tarballs.
docker run -v `pwd`:/build -t -i --rm phusion/holy-build-box-64:latest /hbb_exe/activate-exec bash /build/canu-$version/linux.sh

echo Compress.
xz -9v canu-$version.Darwin-amd64.tar
xz -9v canu-$version.Linux-amd64.tar

exit
