#!/bin/sh

#  Checks that all domap.sh runs finished.  Resubmits those that didn't.

mm=`echo 1`
while [ $mm -lt 2128 ] ; do
  nn=`expr $mm - 1`

  xx=`expr $nn % 8 + 1`
  piece=`echo f1 f2 f3 f4 f5 f6 f7 f8 | cut -d' ' -f $xx`
  nn=`expr $nn / 8`

  xx=`expr $nn % 14 + 1`
  ignore=`echo 0000 0001 0002 0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096  | cut -d' ' -f $xx`
  nn=`expr $nn / 14`

  xx=`expr $nn % 19 + 1`
  skip=`echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19                 | cut -d' ' -f $xx`
  nn=`expr $nn / 19`

  name=`echo $piece.skp$skip.ign$ignore`
  #echo $mm $nn $name

  if [ ! -e /project/huref0/assembly/chr14/out14/chr14.$name.stats ] ; then
      echo qsub -t $mm domap14.sh
  fi

  mm=`expr $mm + 1`
done
