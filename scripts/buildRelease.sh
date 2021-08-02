#!/bin/sh
#
#  Before building a release:
#
#  Make a place to work, grab the bits you want to release:
#    git clone git@github.com:marbl/PACKAGE PACKAGE-release
#    cd PACKAGE-release
#
#  Commit to master:
#    Increase version in documentation/source/conf.py  (if present)
#    Increase version in scripts/version_update.pl
#
#  Build.  This pulls in submodule code.  This build isn't used
#  for creating the release and can be deleted or aborted (once
#  submodules are populated):
#    cd src && gmake
#
#  Tag the next release development at the tip in master:
#    git tag -a v2.2-development -m "Development for v2.2."
#    git push --follow-tags
#
#  Make a branch:
#    git checkout -b v2.1-maintenance
#
#  Commit to branch:
#    Change 'snapshot' to 'release' in scripts/version_update.pl
#    git push --set-upstream origin v2.1-maintenance
#
#  Run this script:
#    scripts/buildRelease.sh PACKAGE 2.1
#
#  Some notes:
#   - The source code from the branch made above will be used
#     to build tarballs with source and binaries.
#   - These tarballs (and scratch space) are called 'build-*'.
#   - The very first step is to rename '.git' to 'dot-git-directory';
#     catastrophic failure (of this script) could leave it
#     renamed.
#

package=$1
version=$2

if [ x$version = x ] ; then
  echo "usage: $0 package-name numeric-version"
  exit
fi

#
#  Cleanup any old build, make space for the new one, and initialize scripts.
#

if [ -e .git ] ; then
    echo Moving .git directory out of the way.
    mv .git dot-git-directory
fi

echo Preparing build trees.

rm -rf build
rm -rf build-darwin build-darwin.out
rm -rf build-linux  build-linux.out
rm -rf build-src

rm  -f build-linux.sh

rm  -f ${package}-${version}.Darwin-amd64.tar ${package}-${version}.Darwin-amd64.tar.xz
rm  -f ${package}-${version}.Linux-amd64.tar  ${package}-${version}.Linux-amd64.tar.xz
rm  -f ${package}-${version}.tar  ${package}-${version}.tar.xz

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
echo >> build-linux.sh  "mv build-darwin ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.Darwin-amd64.tar ${package}-${version}/README* ${package}-${version}/bin ${package}-${version}/lib ${package}-${version}/share"
echo >> build-linux.sh  "mv ${package}-${version} build-darwin"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-linux ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.Linux-amd64.tar  ${package}-${version}/README*  ${package}-${version}/bin  ${package}-${version}/lib  ${package}-${version}/share"
echo >> build-linux.sh  "mv ${package}-${version} build-linux"
echo >> build-linux.sh  ""
echo >> build-linux.sh  "mv build-src ${package}-${version}"
echo >> build-linux.sh  "tar -cf ${package}-${version}.tar              ${package}-${version}/README*  ${package}-${version}/src  ${package}-${version}/scripts"
echo >> build-linux.sh  "mv ${package}-${version} build-src"
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

xz -9v ${package}-${version}.Darwin-amd64.tar
xz -9v ${package}-${version}.Linux-amd64.tar
xz -9v ${package}-${version}.tar

if [ -e dot-git-directory ] ; then
    echo Restoring .git directory.
    mv dot-git-directory .git
fi

exit
