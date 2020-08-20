#!/bin/sh
#
#  Before building a release:
#
#  Make a place to work, grab the bits you want to release:
#    git clone git@github.com:marbl/canu canu-release
#    cd canu-release
#
#  Commit to master:
#    Increase version in documentation/source/conf.py
#    Increase version in scripts/version_update.pl
#
#  Build.  This pulls in submodule code.  This build isn't used
#  for release.
#    cd src && gmake
#
#  Make a branch:
#    git checkout -b v2.1-maintenance
#
#  Commit to branch:
#    Change 'snapshot' to 'release' in scripts/version_update.pl
#
#  Run this script:
#    scripts/buildRelease.sh 2.1
#

version=$1

if [ x$version = x ] ; then
  echo usage: $0 numeric-version
  exit
fi

#
#  Cleanup any old build, make space for the new one, and initialize scripts.
#

echo Preparing build trees.

rm -rf build
rm -rf build-darwin build-darwin.out
rm -rf build-linux  build-linux.out
rm -rf build-src

rm  -f build-linux.sh

rm  -f canu-${version}.Darwin-amd64.tar canu-${version}.Darwin-amd64.tar.xz
rm  -f canu-${version}.Linux-amd64.tar  canu-${version}.Linux-amd64.tar.xz
rm  -f canu-${version}.tar  canu-${version}.tar.xz

mkdir -p build-src/scripts
mkdir -p build-darwin/scripts
mkdir -p build-linux/scripts

rsync -a src/ build-src/src
rsync -a src/ build-darwin/src
rsync -a src/ build-linux/src

cp -p README* build-src/
cp -p README* build-darwin/
cp -p README* build-linux/

cp -p scripts/version_update.pl build-src/scripts/
cp -p scripts/version_update.pl build-darwin/scripts/
cp -p scripts/version_update.pl build-linux/scripts/

echo >> build-linux.sh  "#!/bin/bash"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "rm -rf /dock/build"
echo >> build-linux.sh  "cd /dock/src"
echo >> build-linux.sh  "gmake -j 12 > ../build-linux.out 2>&1"
echo >> build-linux.sh  "cd .."
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build/* build-linux/"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "rm -rf build-darwin/obj"
echo >> build-linux.sh  "rm -rf build-linux/obj"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-darwin canu-$version"
echo >> build-linux.sh  "tar -cf canu-$version.Darwin-amd64.tar canu-$version/README* canu-$version/bin canu-$version/lib canu-$version/share"
echo >> build-linux.sh  "mv canu-$version build-darwin"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-linux canu-$version"
echo >> build-linux.sh  "tar -cf canu-$version.Linux-amd64.tar  canu-$version/README*  canu-$version/bin  canu-$version/lib  canu-$version/share"
echo >> build-linux.sh  "mv canu-$version build-linux"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-src canu-$version"
echo >> build-linux.sh  "tar -cf canu-$version.tar              canu-$version/README*  canu-$version/src  canu-$version/scripts"
echo >> build-linux.sh  "mv canu-$version build-src"
echo >> build-linux.sh  ""
echo >> build-linux.sh  ""

chmod 755 build-linux.sh

#
#
#

echo Build for MacOS.

cd src
gmake -j 12 > ../build-darwin.out 2>&1
cd ..

mv build/* build-darwin/

echo Make static binaries for MacOS.

cd build-darwin
python ../scripts/statifyOSX.py bin lib true true >> ../build-darwin.out 2>&1
python ../scripts/statifyOSX.py lib lib true true >> ../build-darwin.out 2>&1
cd ..

#
#
#

echo Build for Linux and make tarballs.

echo \
docker run -v `pwd`:/dock -t -i --rm phusion/holy-build-box-64:latest /hbb_exe/activate-exec bash /dock/build-linux.sh
docker run -v `pwd`:/dock -t -i --rm phusion/holy-build-box-64:latest /hbb_exe/activate-exec bash /dock/build-linux.sh

#  strip --only-keep-debug

echo Compress.

xz -9v canu-$version.Darwin-amd64.tar
xz -9v canu-$version.Linux-amd64.tar
xz -9v canu-$version.tar

exit
