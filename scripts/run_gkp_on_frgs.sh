#!/usr/local/bin/bash
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#


# Runtime settings
umask 2
ulimit -c 0



#############################################################################
# Function:   System
# Purpose:    executes a unix command & exits on failure
# Parameters: command string
#############################################################################
function System
{
  if [ -z "$1" ]; then
    echo "Warning: System function called without parameter"
  else
    echo "Executing ${1}"
    if ! ${1} 2>> ${prefix}.log; then
      echo "Failed to execute: $1" > ${prefix}.FAILURE
      exit 1
    fi
  fi
}


#############################################################################
# Function:   SafeSystem
# Purpose:    executes a unix command - & returns 0/1
# Parameters: command string
#############################################################################
function SafeSystem
{
  if [ -z "$1" ]; then
    echo "Warning: System function called without parameter"
  else
    echo "Executing ${1}"
    if ! ${1} 2>> ${prefix}.log; then
      echo "Failed to execute: $1" >> ${prefix}.log
      return 1
    fi
  fi
  return 0
}


#############################################################################
# Function:   CheckBinary
# Purpose:    checks existence & executability of a file
# Parameters: path/file
#############################################################################
function CheckBinary
{
  if [ -f ${1} ]; then
    if [ ! -x ${1} ]; then
      echo "${1} exists but cannot be executed!"
      exit 1
    fi
  else
    echo "${1} does not exist, but it is a required binary."
    exit 1
  fi
}

#############################################################################
# Function:   'main'
# Purpose:    runs the screener in 'parallel' using LSF
#             creates a temporary working directory
# Parameters: see Usage function, above
#############################################################################

###############################################
# Check input parameters
if [ -z ${1} ] || [ -z ${2} ] || [ -z ${3} ]; then
  echo "Usage: ${0} <component name>"
  echo "            <binaries directory>"
  echo "            <LSF queue for screener & overlapper>"
  echo "            <optional restart point>"
  exit 1
fi

prefix=${1}
export AS_BIN=${2}
queue_name=${3}

# Replace . characters with _ in prefix (cgw scrogs names with periods...)
new_prefix=`echo ${prefix} | sed "s/\./_/g"`
if [ "${new_prefix}" != "${prefix}" ]; then
  prefix=${new_prefix}
  echo "Script has replaced prefix periods with underscores." >> ${prefix}.log
fi

declare -i restart=0
# set the restart point, if specified
if [ ! -z ${4} ]; then
  restart=${4}
fi
###############################################


###############################################
# Check presence of required binaries
CheckBinary ${AS_BIN}/gatekeeper
CheckBinary ${AS_BIN}/PopulateFragStore
CheckBinary ${AS_BIN}/make_range_file
CheckBinary ${AS_BIN}/make-ovl-script
CheckBinary ${AS_BIN}/overlap
CheckBinary ${AS_BIN}/unitigger

CheckBinary ${AS_BIN}/consensus
CheckBinary ${AS_BIN}/cgw
CheckBinary ${AS_BIN}/terminator
###############################################


###############################################
# Feedback on startup & key variables
echo "Assembler started: `date`" >> ${prefix}.log
echo "Prefix=${prefix}" >> ${prefix}.log
echo "AS_BIN=${AS_BIN}" >> ${prefix}.log
echo "lsf queue=${queue_name}" >> ${prefix}.log
echo "LSF service & project = ${LSF_SERVICE}" >> ${prefix}.log
###############################################


###############################################
# Clean up any previous failure
if [ -f ${prefix}.FAILURE ]; then
  nc="rm -f ${prefix}.FAILURE "
  System "${nc}"
fi
###############################################


###############################################
# gatekeeper 
declare -i num_frg_files=0
if [ ${restart} -le 0 ]; then
  for frg_file in `ls *.frg`; do
    if [ $num_frg_files -eq 0 ]; then
      nc="time ${AS_BIN}/gatekeeper -n ${prefix}  -Q -e 4 -f -c -o ${prefix}.gkpStore ${frg_file} "
    else
      nc="time ${AS_BIN}/gatekeeper -n ${prefix} -Q  -e 4 -a -i ${prefix}.gkpStore ${frg_file} "
    fi
    System "${nc}"
    num_frg_files=${num_frg_files}+1
  done
fi
exit 0
#############################################################################
