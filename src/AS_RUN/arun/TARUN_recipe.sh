## Assembly Script for test against TARUN.

## Files Required
##REQUIRED PREFIX.seq PREFIX.qual
#> TA
echo $PREFIX
echo $WORKDIR
echo $OPTIONS
/usr/local/packages/TIGR_Assembler/bin/TIGR_Assembler $OPTIONS  -q $WORKDIR/$PREFIX.qual  -n $PREFIX -a $WORKDIR/$PREFIX.align -f $WORKDIR/$PREFIX.fasta $WORKDIR/$PREFIX.scratch  < $WORKDIR/$PREFIX.seq > $WORKDIR/$PREFIX.tasm 2>$WORKDIR/$PREFIX.error

#> SUBTASM
echo "OPTIONS: "
echo $OPTIONS
