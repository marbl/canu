#!/bin/sh

installdir=$1

#  Should detect which modules we really need to install.  For now, just install both.


#  Remove old junk

rm -rf File-Path-2.09
rm -rf Filesys-Df-0.92

#  Clean our environment

unset CC
unset CXX

#  Decide which perl we have.

perlversion=`perl -e '$x = ($] < 5.009) ? "old" : "new"; print $x'`

#  Our older perl (5.8.8) uses an older MakeMaker, which doesn't know INSTALL_BASE, and installs the
#  Filesys::Df package in the wrong place.
#
#  Our newer perl (5.10.1) knows INSTALL_BASE, and has File::Path installed already.

if [ $perlversion = "old" ] ; then

  echo "Installing File-Path-2.09"

  echo \
  tar xzf File-Path-2.09.tar.gz
  tar xzf File-Path-2.09.tar.gz

  cd File-Path-2.09

  echo \
  perl Makefile.PL PREFIX=$installdir
  perl Makefile.PL PREFIX=$installdir > make.maker.err 2>&1

  echo \
  make install
  make install > make.install.err 2>&1
  cd ..

  echo "Installing Filesys-Df-0.92"

  echo \
  tar xzf Filesys-Df-0.92.tar.gz
  tar xzf Filesys-Df-0.92.tar.gz

  cd Filesys-Df-0.92

  echo \
  perl Makefile.PL PREFIX=$installdir
  perl Makefile.PL PREFIX=$installdir > make.maker.err 2>&1

  echo \
  make install
  make install > make.install.err 2>&1

  cd ..

  echo \
  mv $installdir/lib64/perl5/site_perl/5*/*/Filesys      $installdir/lib64/perl5/5*/*/
  mv $installdir/lib64/perl5/site_perl/5*/*/Filesys      $installdir/lib64/perl5/5*/*/

  echo \
  mv $installdir/lib64/perl5/site_perl/5*/*/auto/Filesys $installdir/lib64/perl5/5*/*/auto/
  mv $installdir/lib64/perl5/site_perl/5*/*/auto/Filesys $installdir/lib64/perl5/5*/*/auto/

  echo \
  rm -rf $installdir/lib64/perl5/site_perl
  rm -rf $installdir/lib64/perl5/site_perl

else

  echo "Installing Filesys-Df-0.92"

  echo \
  tar xzf Filesys-Df-0.92.tar.gz
  tar xzf Filesys-Df-0.92.tar.gz

  cd Filesys-Df-0.92

  echo \
  perl Makefile.PL INSTALL_BASE=$installdir
  perl Makefile.PL INSTALL_BASE=$installdir > make.maker.err 2>& 1

  #  The toplevel GNU make is setting MAKEFLAGS to be "-j --jobserver-fds=3,4".  The BSD make
  #  invoked below (on FreeBSD and probably OS X) requires a value for -j.  So we just remove it.
  #  In GNU make, '-j' says to run as many tasks in parallel as possible; this isn't heavy lifting,
  #  so won't matter if it's sequential.  The jobserver-fds baloney is to track jobs in parallel,
  #  again, we don't care.

  export MAKEFLAGS=

  echo \
  make install
  make install > make.install.err 2>&1

  cd ..

fi

#  Remove the junk

rm -rf File-Path-2.09
rm -rf Filesys-Df-0.92

