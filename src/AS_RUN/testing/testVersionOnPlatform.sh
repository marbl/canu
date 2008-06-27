#!/usr/bin/env bash
#
###########################################################################
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
# $Id: testVersionOnPlatform.sh,v 1.5 2008-06-27 06:29:21 brianwalenz Exp $
#

cwd=`pwd`
baseTestDir=$cwd
testbenchDir="${baseTestDir}/testbench"
script="apply_assembler"
# ordered from small to large at TIGR
testDirs=(defect small wolbachiaPipientis trypanosomaBrucei yersinaPestis aspergillusFlavus brugiaMalayi trichonomasVaginalis)

if [[ -z "$1" || -z "$2" ]]; then
  echo "Usage: $0  version  platform"
  exit
fi
version=$1
platform=$2

# test for needed files
if [ ! -d ${testbenchDir} ]; then
  echo "Directory ${testbenchDir} must exist!"
  exit
fi

if [ ! -f ${baseTestDir}/run_CA.pl ]; then
  echo "File ${baseTestDir}/run_CA.pl must exit!"
  exit
fi


for testDir in "${testDirs[@]}"; do
  if [ -d ${testbenchDir}/${testDir} ]; then
    echo "Running ${script} on ${testDir} genome(s)"
    cd ${testbenchDir}/${testDir}
    ${script} ${version} ${platform} ${baseTestDir} > ${platform}.${version}.out 2> ${platform}.${version}.out
  else
    echo "Directory ${testbenchDir}/${testDir} doesn't exist!"
  fi
done
