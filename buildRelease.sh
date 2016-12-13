#!/bin/sh

version=$1

if [ x$version = x ] ; then
  echo usage: $0 numeric-version
  exit
fi

git clone git@github.com:marbl/canu.git

mv canu canu-$version

cd canu-$version
git tag v$version
cd src
gmake -j 12 > ../Darwin-amd64.out 2>&1
cd ../..

rm -f canu-$version/linux.sh

echo >> canu-$version/linux.sh  \#\!/bin/bash
echo >> canu-$version/linux.sh  yum install -y git
echo >> canu-$version/linux.sh  cd /build/canu-$version/src
echo >> canu-$version/linux.sh  gmake -j 12 \> ../Linux-amd64.out 2\>\&1
echo >> canu-$version/linux.sh  cd ../..
echo >> canu-$version/linux.sh  tar -cf canu-$version/README* canu-$version.Darwin-amd64.tar canu-$version/Darwin-amd64
echo >> canu-$version/linux.sh  tar -cf canu-$version/README* canu-$version.Linux-amd64.tar  canu-$version/Linux-amd64

chmod 755 canu-$version/linux.sh

docker run -v `pwd`:/build -t -i --rm phusion/holy-build-box-64:latest /hbb_exe/activate-exec bash /build/canu-$version/linux.sh

rm -rf canu-$version/*-amd64/obj

xz -9v canu-$version.Darwin-amd64.tar
xz -9v canu-$version.Linux-amd64.tar

exit
