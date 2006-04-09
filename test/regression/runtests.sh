#!/bin/sh

###################################################################
# Wrapper to run the regression tests for the assembler
###################################################################

BINDIR=""
DIRBIN_LIST="AS_GKP:gatekeeper AS_CNS:consensus"

LOCAL_OSTYPE=$(uname)
MACHINE_TYPE=$(uname -m)

if [ "$MACHINE_TYPE" = "x86_64" ]; then
  LOCAL_OSTYPE="Linux64"
fi

set_BINDIR() {
  pushd ../../ > /dev/null 2>&1
  top=$(pwd)
  popd > /dev/null 2>&1

  BINDIR=${top}/${LOCAL_OSTYPE}/bin
}

check_BINS() {
  local dirbin=""
  local bin=""

  for dirbin in $DIRBIN_LIST
  do
    bin=${BINDIR}/$(echo $dirbin | awk -F: '{ print $2 }') 

    if [ ! -f $bin ]; then
      echo "check_BINS: file $bin does not exist"
      exit 1
    fi

    if [ ! -x $bin ]; then
      echo "check_BINS: file $bin is not executable"
      exit 1
    fi
  done
}

####
# MAIN:
####
set_BINDIR
check_BINS

for dirbin in $DIRBIN_LIST
do
  dir=$(echo $dirbin | awk -F: '{ print $1 }') 
  bin=${BINDIR}/$(echo $dirbin | awk -F: '{ print $2 }') 

  if [ -d $dir ]; then
    pushd $dir > /dev/null 2>&1
    for test in $(ls -1 [0-9][0-9][0-9].sh)
    do
      echo "[*] Running $dir/$test ..."
      cwd=$(pwd)
      ${cwd}/$test $bin
      test_status=$?
      if [ $test_status != "0" ]; then
        echo "Test $test failed"
      fi
    done
    popd > /dev/null 2>&1
  else
    echo "Error: directory $dir does not exist"
  fi
done

