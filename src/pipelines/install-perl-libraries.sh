#!/bin/sh

installdir=$1

#  Should detect which modules we really need to install.  For now, just install both.


#  Remove old junk

rm -rf File-Path-2.09
rm -rf Filesys-Df-0.92

#  Clean our environment

unset CC
unset CXX

#  Decide which perl we have.  Great, thanks for being standard.
#
#  Linux reports:
#    This is perl, v5.8.8 ...
#
#  FreeBSD reports:
#    This is perl 5, version 18, subversion 4 (v5.18.4) built for amd64-freebsd-thread-multi
#
#  OSX reports:
#    This is perl 5, version 16, subversion 3 (v5.16.3) built for darwin-thread-multi-2level

maj=`perl -v | grep "This is" | awk '{ if ($5 = "version") { print $9 } else { print $4 } }' | tr -d 'v()' | tr . \  | awk '{ print $1 }'`
min=`perl -v | grep "This is" | awk '{ if ($5 = "version") { print $9 } else { print $4 } }' | tr -d 'v()' | tr . \  | awk '{ print $2 }'`
ver=`perl -v | grep "This is" | awk '{ if ($5 = "version") { print $9 } else { print $4 } }' | tr -d 'v()' | tr . \  | awk '{ print $3 }'`
pvv="$maj.$min.$ver"

#  Our older perl (5.8.8) uses an older MakeMaker, which doesn't know INSTALL_BASE, and installs the
#  Filesys::Df package in the wrong place.
#
#  Our newer perl (5.10.1) knows INSTALL_BASE, and has File::Path installed already.

if [ $min -lt 10 ] ; then

  echo "Installing File-Path-2.09"
  tar xzf File-Path-2.09.tar.gz
  cd File-Path-2.09
  perl Makefile.PL PREFIX=../$installdir/ca3g > make.maker.err 2>&1
  make install > make.install.err 2>&1
  cd ..

  echo "Installing Filesys-Df-0.92"
  tar xzf Filesys-Df-0.92.tar.gz
  cd Filesys-Df-0.92
  perl Makefile.PL PREFIX=../$installdir/ca3g > make.maker.err 2>&1
  make install > make.install.err 2>&1
  cd ..

  mv $installdir/ca3g/lib64/perl5/site_perl/${pvv}/x86_64-linux-thread-multi/Filesys      $installdir/ca3g/lib64/perl5/${pvv}/x86_64-linux-thread-multi/Filesys
  mv $installdir/ca3g/lib64/perl5/site_perl/${pvv}/x86_64-linux-thread-multi/auto/Filesys $installdir/ca3g/lib64/perl5/${pvv}/x86_64-linux-thread-multi/auto/Filesys

  rm -rf $installdir/ca3g/lib64/perl5/site_perl

else

  echo "Installing Filesys-Df-0.92"
  tar xvzf Filesys-Df-0.92.tar.gz
  cd Filesys-Df-0.92
  perl Makefile.PL INSTALL_BASE=../$installdir/ca3g
  make install
  cd ..

fi

#  Remove the junk

rm -rf File-Path-2.09
rm -rf Filesys-Df-0.92

