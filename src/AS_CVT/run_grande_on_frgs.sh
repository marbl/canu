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
# $Id: run_grande_on_frgs.sh,v 1.2 2004-09-23 20:25:21 mcschatz Exp $
#

# Parameter globals
FGB_ALWAYS="-P -A 1 -f"
FGB_REAPER="-x 1 -z 10 -d 0 -M 0"

CGB_ALWAYS="-P -A 1 -j 5 -b 0"
CGB_BREAKERS="-s -U 0"
CGB_BUBBLES="-U 1"
CGB_FINAL="-U 0"

CNS_ALWAYS="-P -U"
#CNS_BUBBLES="-a L"
CNS_BUBBLES=" "

# Runtime settings
umask 2
ulimit -c 0

# LSF settings
LSF_SERVICE="IR:HUM_ASM:H"
if [ -z ${LSF_ENVDIR} ]; then
  . /usr/local/src/lsf/profile.lsf
fi
LSF_BIN="${LSF_ENVDIR}/../bin"


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
# Function:   RunUnitigger
# Purpose:    runs unitigger with various options
# Parameters: ${1} = assembly name prefix
#             ${2} = cgb custom parameters
#             ${3} = list of overlap files for fgb
#             ${4} = run count: 1 = first run - needed for output names
#             ${5} = last run count - needed for input names
#             ${6} = flag > 0 means fgb should be run
#############################################################################
function RunUnitigger
{
#  echo "Running unitigger on assembly ${1}"
#  echo "cgb basic parameters: ${CGB_ALWAYS}"
#  echo "      cgb parameters: ${2}"
#  echo "       overlap files: ${3}"
#  echo "          run number: ${4}"
#  echo "            last run: ${5}"
#  echo "            run fgb?: ${6}"

  # set up output store/file names  
  out_reaper_store="${1}.reaperStore${4}"
  out_fgb_store="${1}.fgbStore${4}"
  out_cgb_file="${1}.cgb${4}"

  # if parameter 5 is not zero, run fgb
  # won't rerun fgb if no spur overlaps
  if [ ${6} -ne 0 ]; then
    # define the I/O portion of the fgb reaper command line
    if [ ${4} -eq 1 ]; then
      io_string="-f -c -o ${out_reaper_store}"
    else
      io_string="-i ${1}.reaperStore${5} -o ${out_reaper_store}"
    fi

    # run fgb reaper
    nc="time ${AS_BIN}/fgb ${FGB_ALWAYS} ${FGB_REAPER} ${io_string} ${3} "
    System "$nc"

    # convert fgb reaper store to fgb store
    nc="time ${AS_BIN}/fgb ${FGB_ALWAYS} -i ${out_reaper_store} -o ${out_fgb_store} "
    System "$nc"
  else
    # set up links, so stores will be present for cgb & subsequent fgbs
    nc="ln -s ${1}.reaperStore${5} ${out_reaper_store}"
    System "$nc"
    nc="ln -s ${1}.fgbStore${5} ${out_fgb_store}"
    System "$nc"
  fi
  
  # run cgb with specified flags
  nc="time ${AS_BIN}/cgb ${CGB_ALWAYS} ${2} ${1}.frgStore ${out_fgb_store}"
  System "$nc"

  # mv the cgb output file
  nc="mv ${1}.cgb ${1}.cgb${4}"
  System "$nc"
}


#############################################################################
# Function:   DeleteUnitiggerFiles
# Purpose:    Deletes intermediate reaperStores, fgbStores, & cgb files
# Parameters: prefix & number of the one to keep  
#############################################################################
function DeleteUnitiggerFiles
{
  # number 4 of thes sets of stores/files may not be present, so do it by ls
  
  # delete the reaper stores
  for reapers in `ls |egrep "${prefix}.reaperStore[1-4]"`; do
    if [ -d ${reapers} ] && [ "${reapers}" != "${1}.reaperStore${2}" ]; then
      nc="rm -rf ${reapers}"
      SafeSystem "$nc"
    fi
    if [ -f ${reapers} ] && [ "${reapers}" != "${1}.reaperStore${2}.fge" ] && [ "${reapers}" != "${1}.reaperStore${2}.fgv" ]; then
      nc="rm -f ${reapers}"
      SafeSystem "$nc"
    fi
  done
  
  # delete the fgb stores
  for fgbs in `ls |grep "${prefix}.fgbStore[1-4]"`; do
    if [ -d ${fgbs} ] && [ "${fgbs}" != "${1}.fgbStore${2}" ]; then
      nc="rm -rf ${fgbs}"
      SafeSystem "$nc"
    fi
    if [ -f ${fgbs} ] && [ "${fgbs}" != "${1}.fgbStore${2}.fge" ] && [ "${fgbs}" != "${1}.fgbStore${2}.fgv" ]; then
      nc="rm -f ${fgbs}"
      SafeSystem "$nc"
    fi
  done

  # delete the cgb files
  for cgbs in `ls |egrep "${prefix}.cgb[1-4]"`; do
    if [ -f ${cgbs} ] && [ "${cgbs}" != "${1}.cgb${2}" ]; then
      nc="rm -f ${cgbs}"
      SafeSystem "$nc"
    fi
  done
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
CheckBinary ${AS_BIN}/lsf_screener.sh
CheckBinary ${AS_BIN}/PopulateFragStore
CheckBinary ${AS_BIN}/make_range_file
CheckBinary ${AS_BIN}/make-ovl-script
CheckBinary ${AS_BIN}/overlap
CheckBinary ${AS_BIN}/fgb
CheckBinary ${AS_BIN}/cgb
CheckBinary ${AS_BIN}/repair_breakers
CheckBinary ${AS_BIN}/consensus
CheckBinary ${AS_BIN}/cgw
CheckBinary ${AS_BIN}/terminator
CheckBinary ${LSF_BIN}/bsub
CheckBinary ${LSF_BIN}/bdel
CheckBinary ${LSF_BIN}/bkill
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
# gatekeeper - run with quality checking
declare -i num_frg_files=0
if [ ${restart} -le 0 ]; then
  for frg_file in `ls *.frg`; do
    if [ $num_frg_files -eq 0 ]; then
      nc="time ${AS_BIN}/gatekeeper -n ${prefix} -NP -e 10000 -f -c -o ${prefix}.gkpStore ${frg_file} "
    else
      nc="time ${AS_BIN}/gatekeeper -n ${prefix} -NP -e 10000 -a -i ${prefix}.gkpStore ${frg_file} "
    fi
    System "${nc}"
    num_frg_files=${num_frg_files}+1
  done
fi
###############################################

  
###############################################
# PopulateFragStore - relies on lexico ordering of ls
declare -i num_inp_files=0
if [ ${restart} -le 1 ]; then
  for inp_file in `ls *.inp`; do
    if [ $num_inp_files -eq 0 ]; then
      nc="time ${AS_BIN}/PopulateFragStore -f -c -o ${prefix}.frgStore ${inp_file} "
    else
      nc="time ${AS_BIN}/PopulateFragStore -a -i ${prefix}.frgStore ${inp_file} "
    fi
    System "$nc"
    num_inp_files=${num_inp_files}+1
  done
fi
###############################################


###############################################
# make range file that goes from 1 to n - relies on lexico ordering of ls
if [ ${restart} -le 2 ]; then
  for rng_file in `ls *_*.range`; do
    nc="time ${AS_BIN}/make_range_file ${rng_file} ${prefix}.range "
    System "$nc"
  done
fi
###############################################


###############################################
# make overlap script
if [ ${restart} -le 3 ]; then
  nc="time ${AS_BIN}/make-ovl-script -L ${LSF_SERVICE} -ePw -o ${prefix}_batch -q ${queue_name} -s ${prefix}.script ${prefix}.frgStore ${prefix}.range "
  System "$nc"
fi
###############################################

    
###############################################
# overlap
if [ ${restart} -le 4 ]; then
  nc="time ${prefix}.script "
  System "$nc"
fi
###############################################


###############################################
# Unitigger: run as follows
# 1. fgb reaper, fgb converter, cgb with -s -U 0 to generate cgb.crappies
# 1a. repair_breakers to generate breaker_overlaps.ovl
# 2. fgb reaper, fgb converter, cgb with -U 1 to generate ${prefix}.bubble_edges.ovl
# 3. fgb reaper, fgb converter, cgb with -U 0 for 'final' .cgb file
# Then try consensus with rollback
# consensus -aL on .cgb file #3
# if failure, consensus on .cgb file #2
#   if failure, run fgb with bubble_edges.ovl using initial fgbStore &
#      cgb -U 0 for .cgb file #4
#      run consensus -aL on .cgb file #4
#      if failure, run consensus on .cgb file #1
declare -i num_delta=0

# Run fgb/cgb 1st time for fragments in overlapper-determined overlaps
if [ ${restart} -le 5 ]; then
  RunUnitigger ${prefix} "${CGB_BREAKERS}" "${prefix}_0*.ovl ${prefix}_batch*.ovl" 1 0 1
fi

# generate spur smoothing overlaps, if necessary
if [ ${restart} -le 6 ]; then
  num_delta=`wc cgb.crappies | awk '{print $1}'`
  echo "$num_delta spur fragments" >> ${prefix}.log
  if [ ${num_delta} -gt 0 ]; then
    nc="time ${AS_BIN}/repair_breakers -S cgb.crappies -F ${prefix}.frgStore -O breaker_overlaps.ovl ${prefix}.cgb1 "
    if ! SafeSystem "$nc"; then
      echo "***** repair_breakers failed. Continuing..." >> ${prefix}.log
    fi
  fi
fi

# run unitigger second time - take in spur overlaps & generate bubble overlaps
if [ ${restart} -le 7 ]; then
  if [ ! -f breaker_overlaps.ovl ]; then
    num_delta=0
  else
    # update the number of spurs - now the number of overlaps found
    num_delta=`egrep -c "\{OVL" breaker_overlaps.ovl | awk '{print $1}'`
  fi

  # NOTE - if num_delta is 0, then breaker_overlaps.ovl will be ignored - so this is ok
  RunUnitigger ${prefix} "${CGB_BUBBLES}" breaker_overlaps.ovl 2 1 ${num_delta}
fi

# run unitigger third time - take in bubble overlaps & generate final cgb
if [ ${restart} -le 8 ]; then
  # see if bubble edge overlaps were generated
  if [ -f ${prefix}.bubble_edges.ovl ]; then
    num_delta=`egrep -c "\{OVL" ${prefix}.bubble_edges.ovl | awk '{print $1}'`
    RunUnitigger ${prefix} "${CGB_FINAL}" ${prefix}.bubble_edges.ovl 3 2 ${num_delta}
  fi
fi
###############################################


###############################################
# consensus - with rollback to previous unitigger runs
declare -i cns_count=4

if [ ${restart} -le 9 ]; then

  # need to figure out which is last cgb file
  while [ ${cns_count} -gt 0 ] && [ ! -f ${prefix}.cgb${cns_count} ]; do
    cns_count=${cns_count}-1
  done
  if [ ${cns_count} -eq 0 ]; then
    echo "Post-cgb consensus cannot run - there are no usable .cgb files"
    exit 1
  fi
  
  # run consensus until success or total failure
  declare -i cns_failure=1
  while [ ${cns_failure} -eq 1 ]; do

    echo "Running consensus on file cgb file ${cns_count}" >> ${prefix}.log
    # check which .cgb file & consensus options to invoke & how to handle failure
    case "$cns_count" in
    
      # 1 is the last resort
      1 ) nc="time ${AS_BIN}/consensus ${CNS_ALWAYS} ${prefix}.frgStore ${prefix}.cgb${cns_count} "
          System "$nc"
          cns_failure=0
          DeleteUnitiggerFiles ${prefix} 1
          ;;

      # 2 is spurs only - failure gets us to redoing unitigger with just bubbles (#4) or with nothing (#1)
      2 ) nc="time ${AS_BIN}/consensus ${CNS_ALWAYS} ${prefix}.frgStore ${prefix}.cgb${cns_count} "
          if ! SafeSystem "$nc"; then
            echo "Consensus with just spurs failed. Rerunning cgb to add bubbles." >> ${prefix}.log
            if [ -f ${prefix}.bubble_edges.ovl ]; then
              num_delta=`egrep -c "\{OVL" ${prefix}.bubble_edges.ovl | awk '{print $1}'`
              if [ ${num_delta} -gt 0 ]; then
                # run unitigger fourth time - take in bubble overlaps & generate final cgb
                RunUnitigger ${prefix} "${CGB_FINAL}" ${prefix}.bubble_edges.ovl 4 1 ${num_delta}
                cns_count=4
              else
                echo "No bubbles generated. Rolling back to basic cgb." >> ${prefix}.log
                cns_count=1
              fi
            else
              echo "Bubbles overlap file doesn't exist. Rolling back to basic cgb." >> ${prefix}.log
              cns_count=1
            fi
          else
            cns_failure=0
            DeleteUnitiggerFiles ${prefix} 2
          fi
          ;;

      # 3 is bubbles & spurs - failure gets us to spurs only (#2)
      3 ) nc="time ${AS_BIN}/consensus ${CNS_ALWAYS} ${CNS_BUBBLES} ${prefix}.frgStore ${prefix}.cgb${cns_count} "
          if ! SafeSystem "$nc"; then
            echo "Consensus with spurs and bubbles failed. Rolling back to just spurs." >> ${prefix}.log
            cns_count=2
          else
            cns_failure=0
            DeleteUnitiggerFiles ${prefix} 3
          fi
          ;;

      # 4 is bubbles only - failure gets us to #1
      4 ) nc="time ${AS_BIN}/consensus ${CNS_ALWAYS} ${CNS_BUBBLES} ${prefix}.frgStore ${prefix}.cgb${cns_count} "
          if ! SafeSystem "$nc"; then
            echo "Consensus with just bubbles failed. Rolling back to basic cgb." >> ${prefix}.log
            cns_count=1
          else
            cns_failure=0
            DeleteUnitiggerFiles ${prefix} 4
          fi
          ;;
    esac
  done
fi
###############################################


###############################################
# cgw
if [ ${restart} -le 10 ]; then
  nc="time ${AS_BIN}/cgw -m 100000000 -c -j 1 -k 5 -r 4 -s 2 -w 1 -T -P -f ${prefix}.frgStore -g ${prefix}.gkpStore -o ${prefix} ${prefix}.cgi "
  System "$nc"
fi
###############################################


###############################################
# post-cgw concatenation
if [ ${restart} -le 11 ]; then
  if ! cat ${prefix}.cgw ${prefix}.cgw_contigs ${prefix}.cgw_scaffolds > ${prefix}.cgw_total; then
    echo "Failed to cat cgw output files together!" >> ${prefix}.log
    exit 1
  fi
fi
###############################################


###############################################
# consensus
if [ ${restart} -le 12 ]; then
  nc="time ${AS_BIN}/consensus ${prefix}.frgStore ${prefix}.cgw_total "
  System "$nc"
fi
###############################################


###############################################
# terminator
if [ ${restart} -le 13 ]; then
  nc="time ${AS_BIN}/terminator -u -P -g ${prefix}.gkpStore -f ${prefix}.frgStore -i ${prefix}.cns -o ${prefix}.asm -m ${prefix}.map "
  System "$nc"
fi
###############################################


echo "Assembler completed: `date`" >> ${prefix}.log

echo "Removing several intermediate files"
nc="rm -rf *.inp *.urc *.rez.* *.crocks.* *.stone.* *.cstones.* *.ckp.* *.SeqStore *.cgw_total"
SafeSystem "$nc" 

exit 0
#############################################################################
