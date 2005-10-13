#!/usr/local/bin/bash

cwd=`pwd`
baseTestDir=$cwd
testbenchDir="${baseTestDir}/testbench"
script="${baseTestDir}/apply_assembler.pl"
# ordered from small to large
testDirs=(defect small wolbachiaPipientis yersinaPestis aspergillusFlavus brugiaMalayi trichonomasVaginalis)

if [[ -z "$1" || -z "$2" ]]; then
  echo "Usage: $0  version  platform"
  exit
fi
version=$1
platform=$2

for testDir in "${testDirs[@]}"; do
  cd ${testbenchDir}/${testDir}
  echo "Running ${script} on ${testDir} genome(s)"
  perl ${script} ${version} ${platform} ${baseTestDir} > ${platform}.${version}.out 2> ${platform}.${version}.out
done
