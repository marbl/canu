#!/bin/sh

#  Build Canu using all available and supported compilers.
#   - does not disturb the usual build/ directory.
#
#  Only works on FreeBSD; compiler names (and LDFLAGS) would
#  need to be adjusted for other systems.

if [ ! -e ./Makefile ] ; then
  echo ERROR: Must be run from within the src/ directory.
  exit
fi

for v in 9 10 11 12 13 ; do
  echo gcc-$v
  st=$(date +%s)
  rm -rf ../build-gcc$v
  export LDFLAGS="-rpath /usr/local/lib/gcc$v"
  gmake -j 999 DESTDIR=../build-gcc$v CC=gcc$v CXX=g++$v > ../build-gcc$v.out 2> ../build-gcc$v.err || echo "  FAIL"
  unset LDFLAGS
  ed=$(date +%s)
  ss=$(expr $ed - $st)
  echo "  $ss seconds"
done

for v in 10 11 12 13 14 15 16 17 ; do
  echo llvm-$v
  st=$(date +%s)
  rm -rf ../build-llvm$v
  gmake -j 999 DESTDIR=../build-llvm$v CC=clang$v CXX=clang++$v > ../build-llvm$v.out 2> ../build-llvm$v.err || echo "  FAIL"
  ed=$(date +%s)
  ss=$(expr $ed - $st)
  echo "  $ss seconds"
done

exit
