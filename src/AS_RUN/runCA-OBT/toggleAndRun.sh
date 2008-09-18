#!/bin/bash

shopt -s extglob || exit 1 # needed for frg cp match +(0) below

usage="$0 asmBinPath specFile
Must be run from the toplevel asm directory with the stores, 5-consensus and 9-terminator."

if [[ $# != 2 ]]
then
    echo $usage
    exit
fi

asmBin=$1
specFile=$2

toggler='markUniqueUnique.rb'
gatekeeper='gatekeeper'
ecrEdits='frg.ECREdits.txt'

if [[ ! -d "9-terminator" ]]
then
    echo 9-terminator not found exiting
    exit 1
fi

prefix=`echo *.gkpStore | sed 's/.gkpStore//'`
if [[ $prefix == '' ]]
then
    echo 'prefix unset from echo *.gkpStore'
    exit 2
fi

newDir=toggledAsm
mkdir $newDir
cd $newDir    || exit

# if the specFile is local link it
if [[ -e ../$specFile ]]
then
    ln -s ../$specFile .
fi

# frg file needed by caqc.pl
ln ../$prefix.frg .

# link the stores for space savings
ln -s ../*.ovlStore .
gkp="$prefix.gkpStore"
mkdir $gkp
(cd $gkp && ln -s ../../$gkp/??? .)

# but the frg store is rewritten by cgw, so reset the ECR clear-ranges
rm -f $gkp/frg
cp ../$gkp/frg $gkp/frg || exit
$asmBin/gatekeeper -dumpfragments -tabular -allreads -clear OBT $gkp | grep -v "UID" | awk '{print "frg uid "$1" ECR1 ALL "$12" "$13}' > $gkp/$ecrEdits || exit
$asmBin/gatekeeper --edit $gkp/$ecrEdits $gkp > /dev/null || exit

# runCA looks for the 5-consensus *.err files at some point
conDir=5-consensus
mkdir $conDir
(cd $conDir && ln -s ../../$conDir/*.err .)

# create the toggled cgi file
cgiInput=$conDir/$prefix_*.cgi
cgi=$conDir/$prefix.cgi
if ($asmBin/$toggler ../*.asm ../$cgiInput > $cgi 2> toggle.err == 0)
then
    echo No toggling occured. Finished.
    exit 0
fi

# restart the asm at cgw, with the given spec file so options get passed
$asmBin/runCA -s $specFile -p $prefix -d . $cgi > runCA.out 2>&1 &
