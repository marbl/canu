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
# $Id: run_components_thru_grande.sh,v 1.3 2005-03-22 19:04:52 jason_miller Exp $
#
#
#set:
# ulimit -c 0
#
# AS_BIN
#
#parameters:
# queue_name
# AS_BIN
# regional_bactig_gatekeeper_store
# regional_fragment_gatekeeper_store
# grande_gatekeeper_store
# grande_fragment_store
# set of comps-input/.icns files
# repeats filename
#
#binaries needed:
# extract_component
# make_cns_list
# bsub
# all assembler programs
# assembler script
# make_range_file
#
#processing:
# check command line parameters
# check for presence of binaries & scripts
# for each comps-input/.icns file on command line:
#   make prefix directory
#   cd prefix directory
#   make_cns_list input.icns $prefix_cns_list.txt
#   extract_component
#   copy repeats filename
#   bsub assembler without blocking script specifying:
#     AS_BIN
#     data directory
#
#assembler script:
# export AS_BIN=command-line parameter
# prefix=command-line parameter
# $AS_BIN/gatekeeper -c -o prefix.gkpStore with repeats.frg 2>> prefix.log
# $AS_BIN/gatekeeper -NQP -a -i prefix.gkpStore with INT*.frg 2>> prefix.log
# $AS_BIN/gatekeeper -NQP -a -i prefix.gkpStore with EXT*.frg 2>> prefix.log
# mkdir $prefix.urcStore
# $AS_BIN/update_screener_lib -f -r $prefix.urcStore/screen_items.lib $prefix_00001.inp
# $AS_BIN/lsf_screener.sh -a -i $prefix.urcStore R -q $queue_name -n 10000 -d $AS_BIN $prefix_00002.inp
# $AS_BIN/lsf_screener.sh -a -i $prefix.urcStore R -q $queue_name -n 10000 -d $AS_BIN $prefix_00003.inp
# $AS_BIN/PopulateFragStore -c $prefix.frgStore $prefix_00002.urc
# $AS_BIN/PopulateFragStore -a -i $prefix.frgStore $prefix_00003.urc
# $AS_BIN/make_range_file $prefix_00003.range $prefix.range
# $AS_BIN/make-ovl-script -ePw -o $prefix_batch -q $queue_name -s $prefix.script $prefix.frgStore $prefix.range
# // cp $AS_BIN/overlap OlapDir*
# $prefix.script
# $AS_BIN/fgb -P -A 1 -f -c -o $prefix.fgbStore $prefix_00002.ovl $prefix_00003.ovl $prefix_batch*.ovl
# $AS_BIN/cgb -A 1 -S
# $AS_BIN/consensus -U
# $AS_BIN/AS_CGB_fom2uom
# $AS_BIN/cgw -c -j 1 -k 10 -r 4


# Globals - defined here for convenience
CVS_REVISION="\$Revision: 1.3 $"

library_filename_stub="screen_items.lib"
lsf_job_group="/Assembly"

lsf_dir="/usr/local/lsf/bin"
bsub_program=${lsf_dir}"/bsub -P IR:HUM_ASM"
bdel_program=${lsf_dir}"/bdel"

ls_pipe="ls *${urc_suffix} | sort -n"
cns_list="cns_list.txt"

curr_dir=`pwd`

ulimit -c 0

#############################################################################
# Function:   Usage
# Purpose:    prints usage instructions & exits
# Parameters: pass a string to print a message - not required
#############################################################################
function Usage
{
  if [ $# -ne 0 ]; then
    echo -e "\nUsage error!"
    echo -n $1
    echo -e "\n"
  fi
  echo -e "\n${0##*/} [-h] [-d AS_BIN] [-q LSF queue1] -[-Q LSF queue2]"
  echo -e "\t\t[-b Regional bactig & gatekeeper store]"
  echo -e "\t\t[-f Regional fragment & gatekeeper store]"
  echo -e "\t\t[-G Grande gatekeeper store]"
  echo -e "\t\t[-F Grande fragment store]"
  echo -e "\t\t[-r Repeat screen items .frg file]"
  echo -e "\t\t[-a | -e]"
  echo -e "\t\t<input (comps-input)/.icns filenames>\n"
  echo "-h              print this message"
  echo "-d              directory where screener binaries are located"
  echo "                  use only to override AS_BIN environment variable"
  echo "-q LSF queue    name of LSF queue for running assembler"
  echo "-Q LSF queue    name of LSF queue for running screener & overlapper"
  echo "-b store        name of regional bactig/gatekeeper store directory"
  echo "-f store        name of regional fragment/gatekeeper store directory"
  echo "-G store        name of grande gatekeeper store directory"
  echo "-F store        name of grande fragment store directory"
  echo "-r filename     name of repeat screen items file"
  echo "-a              assemble only (extraction already done)"
  echo "-e              extract only (assembly to be done separately"
  echo "input files     arbitrary number of .icns files"
  exit 1
}


#############################################################################
# Function:   CheckParameters
# Purpose:    checks command-line parameters. If there's a problem
#             Usage is called with a message to print
# Parameters: none
#############################################################################
function CheckParameters
{
  ###################################################
  # don't extract & don't assemble doesn't make sense
  if [ ${do_extraction} -eq 0 ] && [ ${do_assembly} -eq 0 ]; then
    Usage "Please either extract or assemble or both."
  fi
  ###################################################
  

  ###################################################
  # LSF-related checks  

  # queue name1
  if [ -z "$queue_name1" ]; then
    Usage "Please specify an LSF queue name for the assembler."
  fi
  
  # queue name2
  if [ -z "$queue_name2" ]; then
    Usage "Please specify an LSF queue name for the screener & overlapper."
  fi

  # must not be the same queue!
  if [ "$queue_name1" = "$queue_name2" ]; then
    Usage "Please specify two different queues."
  fi
  ###################################################


  ###################################################
  # file & store existence checks
  
  # regional bactig/gatekeeper store
  if [ -z $r_btg_store ]; then
    Usage "Please specify a regional bactig/gatekeeper store."
  fi
  if [ ! -d $r_btg_store ]; then
    Usage "The regional bactig/gatekeeper store ${r_btg_store} does not exist."
  fi
  
  # regional fragment/gatekeeper store
  if [ -z $r_frg_store ]; then
    Usage "Please specify a regional fragment/gatekeeper store."
  fi
  if [ ! -d $r_frg_store ]; then
    Usage "The regional fragment/gatekeeper store ${r_frg_store} does not exist."
  fi
  
  # grande gatekeeper store
  if [ -z $g_gkp_store ]; then
    Usage "Please specify a grande gatekeeper store."
  fi
  if [ ! -d $g_gkp_store ]; then
    Usage "The grande gatekeeper store ${g_gkp_store} does not exist."
  fi
  
  # grande fragment store
  if [ -z $g_frg_store ]; then
    Usage "Please specify a grande fragment store."
  fi
  if [ ! -d $g_frg_store ]; then
    Usage "The grande fragment store ${g_frg_store} does not exist."
  fi
  
  # repeats file
  if [ -z $repeats_file ]; then
    Usage "Please specify a repeats file."
  fi
  if [ ! -f $repeats_file ]; then
    Usage "The repeats file ${repeats_file} does not exist."
  fi
  ###################################################
}


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
    if ! ${1}; then
      echo "Failed to execute: $1"
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
# initialize variables as integers or strings

r_btg_store=""
r_frg_store=""
g_gkp_store=""
g_frg_store=""
repeats_file=""
queue_name1=""
queue_name2=""

declare -i do_assembly=1
declare -i do_extraction=1
declare -i all_is_well=0

# parse command line - collect parameters
while getopts ":hd:q:Q:b:f:G:F:r:ae" opt; do
  case $opt in
    h ) Usage;;
    q ) queue_name1=$OPTARG;;
    Q ) queue_name2=$OPTARG;;
    d ) AS_BIN=$OPTARG;;
    b ) r_btg_store=$OPTARG;;
    f ) r_frg_store=$OPTARG;;
    G ) g_gkp_store=$OPTARG;;
    F ) g_frg_store=$OPTARG;;
    r ) repeats_file=$OPTARG;;
    a ) do_extraction=0;;
    e ) do_assembly=0;;
   \? ) Usage "Bad parameter(s)"
  esac
done

# make the paths absolute
if [ "${AS_BIN##/*}" = "${AS_BIN}" ]; then
  AS_BIN=${curr_dir}/${AS_BIN}
fi
if [ "${r_btg_store##/*}" = "${r_btg_store}" ]; then
  r_btg_store=${curr_dir}/${r_btg_store}
fi
if [ "${r_frg_store##/*}" = "${r_frg_store}" ]; then
  r_frg_store=${curr_dir}/${r_frg_store}
fi
if [ "${g_gkp_store##/*}" = "${g_gkp_store}" ]; then
  g_gkp_store=${curr_dir}/${g_gkp_store}
fi
if [ "${g_frg_store##/*}" = "${g_frg_store}" ]; then
  g_frg_store=${curr_dir}/${g_frg_store}
fi
if [ "${repeats_file##/*}" = "${repeats_file}" ]; then
  repeats_file=${curr_dir}/${repeats_file}
fi

echo -e "Command line: $@\n"
echo -e "Using AS_BIN=${AS_BIN}\n"

# Check parameters after parsing
CheckParameters

# the final, unflagged parameters are the input filename
shift $(($OPTIND - 1))

# make sure input files were specified
if [ -z "$@" ]; then
  Usage "Need one or more input .icns files"
  exit 1
fi

# loop over the input files
for filename in "$@"; do

  # determine .icns filename prefix - first make name absolute
  if [ "${filename##/*}" = "${filename}" ]; then
    filename=${curr_dir}/${filename}
  fi
  prefix=${filename##/*/}
  prefix=${prefix%%.*}
  echo "Processing ${prefix}"

  # if extracting, need to set up directory & extract in it
  if [ ${do_extraction} -eq 1 ]; then
  
    # if directory is there, assume it's already been extracted
    if [ -d ${prefix} ]; then
      echo "Directory ${prefix} already exists. Assuming already extracted."
      
      # move to the existing subdirectory in case we're also assembling
      next_command="cd ${prefix}"
      System "${next_command}"
      
      all_is_well=1
    else
      # if here, do the extraction
      
      # make prefix directory
      next_command="mkdir ${prefix}"
      System "${next_command}"
  
      # move to new directory
      next_command="cd ${prefix}"
      System "${next_command}"

      # make a file listing the component .cns files of this .icns file
      next_command="${AS_BIN}/make_cns_list ${filename} ${cns_list}"
      if ! ${next_command}; then
        echo "Failed to read .icns file: ${filename}."
        echo "Removing directory ${prefix} and continuing."
    
        # move to previous directory
        next_command="cd ${curr_dir}"
        System "${next_command}"

        # remove the new directory
        next_command="rm -rf ${prefix}"
        System "${next_command}"
      
        all_is_well=0
      else
        # if here, cns_list file was created ok
        
        # run extraction program
        # don't use System, since this stupid program core dumps
        # but it produces good output - so far as we know
        echo "Executing ${AS_BIN}/extract_component"
        time ${AS_BIN}/extract_component -b ${r_btg_store} -f ${r_frg_store} -G ${g_gkp_store} -F ${g_frg_store} -l ${cns_list} -m >> data_extraction.log

        # copy the repeats file to the local directory
        next_command="cp ${repeats_file} repeats.frg"
        System "${next_command}"
      
        all_is_well=1
      fi
    fi
  else
    # if here, just assembling
    
    # move to new directory
    next_command="cd ${prefix}"
    System "${next_command}"
  
    all_is_well=1
  fi

  # launch the assembler if all is well and do_assembly
  if [ ${all_is_well} -eq 1 ] && [ ${do_assembly} -eq 1 ]; then  
    # run the assembler without blocking
    # need to keep in mind memory requirements for fgb, cgb
    # so far, hasn't exceeded 600MB
    next_command="bsub -q ${queue_name1} -o /dev/null -R \"select[physmem>800]rusage[physmem=800]\" ${AS_BIN}/run_grande_on_component.sh ${prefix} ${AS_BIN} ${queue_name2}"
    System "$next_command"
  fi

  # move to initial directory in any event
  next_command="cd ${curr_dir}"
  System "${next_command}"
  
done
# done looping over input files


exit 0
################################# THE END #################################
