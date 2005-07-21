#!/bin/sh

#  Builds a bunch of distributable .tar.bz2 files in ...../dist/

mkdir dist
mkdir dist/sim4db
mkdir dist/meryl
mkdir dist/leaff
mkdir dist/atac
mkdir dist/ESTmapper
mkdir dist/snapper

cp -pr ../build configure.sh Make.include libutil libbio libsim4 sim4db sim4dbutils \
  dist/sim4db

cp -pr ../build configure.sh Make.include libutil libbio libmeryl meryl \
  dist/meryl

cp -pr ../build configure.sh Make.include libutil libbio leaff \
  dist/leaff

#  fails in atac
cp -pr ../build configure.sh Make.include libutil libbio libmeryl libkmer seatac atac atac-driver \
  dist/atac

cp -pr ../build configure.sh Make.include libutil libbio ESTmapper leaff libkmer libmeryl libsim4 meryl seagen sim4db sim4dbutils \
  dist/ESTmapper

cp -pr ../build configure.sh Make.include libutil libbio libkmer libmeryl libsim4 snapper \
  dist/snapper


version=`date +%Y%m%d-%k%M`

cd dist
tar -cf - sim4db    | bzip2 -9vc > sim4db-$version.tar.bz2
tar -cf - meryl     | bzip2 -9vc > meryl-$version.tar.bz2
tar -cf - leaff     | bzip2 -9vc > leaff-$version.tar.bz2
tar -cf - atac      | bzip2 -9vc > atac-$version.tar.bz2
tar -cf - ESTmapper | bzip2 -9vc > ESTmapper-$version.tar.bz2
tar -cf - snapper   | bzip2 -9vc > snapper-$version.tar.bz2

