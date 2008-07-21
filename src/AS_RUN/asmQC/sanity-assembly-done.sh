#!/bin/sh

bindir="$1"
wrkdir="$2"
prefix="$3"
thisdate="$4"
lastdate="$5"
gooddate="ref"

if [ ! -e $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm ] ; then
  echo "Assembly result for `pwd`: FAILURE (no asm file)"
  echo ""
  cat $wrkdir/$thisdate/$prefix/runCA.sge.out.[0-9][0-9] | tail
  exit
fi

resultlast="none"
resultgood="none"

if [ -e $wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm ] ; then
  if diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm \
             $wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm ; then
    resultlast="same"
  else
    resultlast="differs"
  fi
fi

if [ -e $wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm ] ; then
  if diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm \
             $wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm ; then
    resultgood="same"
  else
    resultgood="differs"
  fi
fi

echo "Assembly result for `pwd`: SUCCESS (last: $resultlast) (reference: $resultgood)"
echo ""

if [ $resultlast = "differs" -o $resultgood = "differs" ] ; then
  perl $wrkdir/mergeqc.pl $wrkdir/$gooddate/$prefix/$prefix.qc $wrkdir/$lastdate/$prefix/$prefix.qc $wrkdir/$thisdate/$prefix/$prefix.qc
fi
