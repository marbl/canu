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
prefix=${1}
export AS_BIN=${2}
queue_name=${3}
declare -i restart=0

if [ ! -z ${4} ]; then
  restart=${4}
fi

echo "Assembler started: `date`" >> ${prefix}.log
echo "Prefix=${prefix}" >> ${prefix}.log
echo "AS_BIN=${AS_BIN}" >> ${prefix}.log
echo "lsf queue=${queue_name}" >> ${prefix}.log
service_project_priority="IR:HUM_ASM:H"
echo "LSF service & project = ${service_project_priority}" >> ${prefix}.log

#############################################################################
# Function:   System
# Purpose:    executes a unix command
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
# Function:   'main'
# Purpose:    runs the screener in 'parallel' using LSF
#             creates a temporary working directory
# Parameters: see Usage function, above
#############################################################################
umask 2
ulimit -c 0

if [ -f ${prefix}.FAILURE ]; then
  next_command="rm -f ${prefix}.FAILURE"
  System "${next_command}"
fi

###############################################
# gatekeeper
# repeats file
if [ ${restart} -le 1 ]; then
  next_command="time ${AS_BIN}/gatekeeper -n ${prefix} -f -e 1000000 -c -o ${prefix}.gkpStore dros1.frg"
  System "${next_command}"
fi
  
# Celera internal data file 
if [ ${restart} -le 2 ]; then
  next_command="time ${AS_BIN}/gatekeeper -n ${prefix} -NP -e 1000000 -a -i ${prefix}.gkpStore dros2.frg"
  System "${next_command}"
fi
#
# External data files
#if [ ${restart} -le 3 ]; then
#  next_command="time ${AS_BIN}/gatekeeper -n ${prefix} -NP -e 100000 -a -i ${prefix}.gkpStore new_vanilla.frg"
#  System "${next_command}"
#  next_command="time ${AS_BIN}/gatekeeper -n ${prefix} -NP -e 100000 -a -i ${prefix}.gkpStore new_deluxe.frg"
#  System "${next_command}"
#fi
###############################################

  
###############################################
# screener
# make a screener store
#if [ ${restart} -le 4 ]; then
#  if [ -d ${prefix}.urcStore ]; then
#    next_command="rm -rf ${prefix}.urcStore"
#    System "${next_command}"
#  fi
#  next_command="mkdir ${prefix}.urcStore"
#  System "${next_command}"
#fi
#
# use the update program for the repeats file
#if [ ${restart} -le 5 ]; then
#  next_command="time ${AS_BIN}/update_screener_lib -f -r ${prefix}.urcStore/screen_items.lib ${prefix}_00001.inp "
#  System "${next_command}"
#fi

# use the lsf screener on internal celera file
#if [ ${restart} -le 6 ]; then
#  next_command="time ${AS_BIN}/lsf_screener.sh -P ${service_project_priority} -a -i ${prefix}.urcStore -R -q ${queue_name} -n 5000 -d ${AS_BIN} ${prefix}_00002.inp"
#  System "$next_command"
#fi

# use the lsf screener on external data file
#if [ ${restart} -le 7 ]; then
#  next_command="time ${AS_BIN}/lsf_screener.sh -P ${service_project_priority} -a -i ${prefix}.urcStore -R -q ${queue_name} -n 5000 -d ${AS_BIN} ${prefix}_00003.inp"
#  System "$next_command"
#  next_command="time ${AS_BIN}/lsf_screener.sh -P ${service_project_priority} -a -i ${prefix}.urcStore -R -q ${queue_name} -n 5000 -d ${AS_BIN} ${prefix}_00004.inp"
#  System "$next_command"
#fi
###############################################


###############################################
# PopulateFragStore
if [ ${restart} -le 8 ]; then
  next_command="time ${AS_BIN}/PopulateFragStore -c -o ${prefix}.frgStore ${prefix}_00001.inp "
  System "$next_command"
fi

if [ ${restart} -le 9 ]; then
  next_command="time ${AS_BIN}/PopulateFragStore -a -i ${prefix}.frgStore ${prefix}_00002.inp "
  #System "$next_command"
  #next_command="time ${AS_BIN}/PopulateFragStore -a -i ${prefix}.frgStore ${prefix}_00004.inp "
  System "$next_command"
fi
###############################################


###############################################
# make range file that goes from 1 to n
#if [ ${restart} -le 10 ]; then
#  next_command="time ${AS_BIN}/make_range_file ${prefix}_00004.range ${prefix}.range "
#  System "$next_command"
#fi
###############################################


###############################################
# make overlap script
#if [ ${restart} -le 11 ]; then
#  next_command="time ${AS_BIN}/make-ovl-script -L ${service_project_priority} -ePw -o ${prefix}_batch -q ${queue_name} -s ${prefix}.script ${prefix}.frgStore ${prefix}.range "
#  System "$next_command"
#fi
###############################################

    
###############################################
# overlap
#if [ ${restart} -le 12 ]; then
#  next_command="time ${prefix}.script"
#  System "$next_command"
#fi
###############################################

###############################################
# fgb (The executable is kludged to be the Reaper. )
REAPER="-x 1 -z 10 -d 0 -M 0"
if [ ${restart} -le 15 ]; then
#  next_command="time ${AS_BIN}/fgb -P -A 1 ${REAPER} -f -c -o ${prefix}.fgbStore.reaper ${prefix}_00002.ovl ${prefix}_00003.ovl ${prefix}_00004.ovl ${prefix}_batch*.ovl "
  next_command="time ${AS_BIN}/fgb -P -A 1 ${REAPER} -f -c -o ${prefix}.fgbStore.reaper ${prefix}_00001.ovl ${prefix}_00002.ovl ${prefix}_batch*.ovl "
  System "$next_command"
fi
###############################################

###############################################
# fgb (first time: convert a Reaper store to a FGB store)
if [ ${restart} -le 16 ]; then
  next_command="time ${AS_BIN}/fgb -P -A 1 -i ${prefix}.fgbStore.reaper -o ${prefix}.fgbStore"
    System "$next_command"
fi
###############################################

###############################################
# cgb (1st time for spur fragment smoothing)
if [ ${restart} -le 17 ]; then
  next_command="time ${AS_BIN}/cgb -P -A 1 -j 5 -s -b 0 -U 0  ${prefix}.frgStore ${prefix}.fgbStore "
  System "$next_command"
fi
###############################################


if [ ${restart} -le 18 ]; then
  declare -i num_crappies=`wc cgb.crappies | awk '{print $1}'`
  echo "$num_crappies crappies"

  if [ -f cgb.crappies ] && [ ${num_crappies} -gt 0 ]; then
    ###############################################
    # repair breakers
    next_command="time ${AS_BIN}/repair_breakers -S cgb.crappies -F ${prefix}.frgStore -O breaker_overlaps.ovl ${prefix}.cgb "
    System "$next_command"
    ###############################################

    ###############################################
    # fgb (2nd time)
    next_command="time ${AS_BIN}/fgb -P -A 1 -a -i  ${prefix}.fgbStore breaker_overlaps.ovl "
    System "$next_command"
    ###############################################
  fi
fi


###############################################
# cgb (2nd time)
if [ ${restart} -le 19 ]; then
    next_command="time ${AS_BIN}/cgb -P -A 1 -j 5 -b 0 -U 1 ${prefix}.frgStore ${prefix}.fgbStore "
    # run_cgb_with_bubbles.pl -C "-A 1 -j 5 -b 0" ${prefix}.frgStore ${prefix}.fgbStore
  System "$next_command"
fi
###############################################

# Bubble smoothing goes here.


###############################################
# consensus - handle failure as if it's repair_breakers or bubble smoothing related

if ! [ -f ${AS_BIN}/consensus ]; then
   echo "consensus executable not present" >> ${prefix}.log ; exit 1
fi

if [ ${restart} -le 23 ]; then
    echo "Executing ${AS_BIN}/consensus -P -U -z ${prefix}.frgStore ${prefix}.cgb"
  if ! ${AS_BIN}/consensus -P -U -z ${prefix}.frgStore ${prefix}.cgb 2>> ${prefix}.log; then
    
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >> ${prefix}.log
    echo "Consensus failure attributed to spur or bubble smoothing. Rerunning from fgb." >> ${prefix}.log
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" >> ${prefix}.log
    
    next_command="time ${AS_BIN}/fgb -P -A 1 -i ${prefix}.fgbStore.reaper -f -o ${prefix}"
    System "$next_command"
    next_command="time ${AS_BIN}/cgb -P -A 1 -j 5 -b 0 -U 0 ${prefix}.frgStore ${prefix}.fgbStore "
    System "$next_command"
    next_command="time ${AS_BIN}/consensus -P -U -z ${prefix}.frgStore ${prefix}.cgb "
    System "$next_command"
  fi
fi
###############################################


# NOTE: uses .cgi as input, not .utg!!!!!!!!!!!
###############################################
# cgw 
if [ ${restart} -le 24 ]; then
  next_command="time ${AS_BIN}/cgw -c -j 1 -k 5 -r 4 -s 2 -w 1 -T -P -f ${prefix}.frgStore -g ${prefix}.gkpStore -o ${prefix} ${prefix}.cgi "
  System "$next_command"
fi
###############################################


###############################################
# consensus
if [ ${restart} -le 25 ]; then
  # cat the multiple cgw files together
  cat ${prefix}.cgw ${prefix}.cgw_contigs ${prefix}.cgw_scaffolds > ${prefix}.cgw_total

  next_command="time ${AS_BIN}/consensus ${prefix}.frgStore ${prefix}.cgw_total "
  System "$next_command"
fi
###############################################


###############################################
# terminator
if [ ${restart} -le 26 ]; then
  next_command="time ${AS_BIN}/terminator -u -P -g ${prefix}.gkpStore -f ${prefix}.frgStore -i ${prefix}.cns -o ${prefix}.asm -m ${prefix}.map "
  System "$next_command"
fi
###############################################


echo "Assembler completed: `date`" >> ${prefix}.log

echo "Removing .inp, .urc, and *.ckp.* files"
rm -f *.inp *.urc *.ckp.*

exit 0
#############################################################################
