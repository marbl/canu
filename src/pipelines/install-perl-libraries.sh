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
  perl Makefile.PL PREFIX=../$installdir/ca3g
  perl Makefile.PL PREFIX=../$installdir/ca3g > make.maker.err 2>&1

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
  perl Makefile.PL PREFIX=../$installdir/ca3g
  perl Makefile.PL PREFIX=../$installdir/ca3g > make.maker.err 2>&1

  echo \
  make install
  make install > make.install.err 2>&1

  cd ..

  echo \
  mv $installdir/ca3g/lib64/perl5/site_perl/5*/*/Filesys      $installdir/ca3g/lib64/perl5/5*/*/
  mv $installdir/ca3g/lib64/perl5/site_perl/5*/*/Filesys      $installdir/ca3g/lib64/perl5/5*/*/

  echo \
  mv $installdir/ca3g/lib64/perl5/site_perl/5*/*/auto/Filesys $installdir/ca3g/lib64/perl5/5*/*/auto/
  mv $installdir/ca3g/lib64/perl5/site_perl/5*/*/auto/Filesys $installdir/ca3g/lib64/perl5/5*/*/auto/

  echo \
  rm -rf $installdir/ca3g/lib64/perl5/site_perl
  rm -rf $installdir/ca3g/lib64/perl5/site_perl

else

  echo "Installing Filesys-Df-0.92"

  echo \
  tar xzf Filesys-Df-0.92.tar.gz
  tar xzf Filesys-Df-0.92.tar.gz

  cd Filesys-Df-0.92

  echo \
  perl Makefile.PL INSTALL_BASE=../$installdir/ca3g
  perl Makefile.PL INSTALL_BASE=../$installdir/ca3g > make.maker.err 2>& 1

  echo \
  make install
  make install > make.install.err 2>&1

  cd ..

fi

#  Remove the junk

rm -rf File-Path-2.09
rm -rf Filesys-Df-0.92

