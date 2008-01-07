#!/bin/bash

shopt -s extglob || exit 1 # needed for frg cp match below

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

ln -s ../*.ovlStore .
gkp="$prefix.gkpStore"
mkdir $gkp
(cd $gkp && ln -s ../../$gkp/??? .) 
rm -f $gkp/frg
cp ../$gkp/frg.before-7-1-ECR-scaffold.+(0) $gkp/frg || exit

$asmBin/$toggler ../*.asm ../5-consensus/$prefix.cgi > $prefix.cgi 2> toggle.err
if cmp -s ../5-consensus/$prefix.cgi $prefix.cgi
then
    echo No toggling occured. Finished.
    exit 0
fi
$asmBin/runCA -s $specFile -p $prefix -d . *.cgi > runCA.out 2>&1 &
