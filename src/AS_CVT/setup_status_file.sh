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
# $Id: setup_status_file.sh,v 1.2 2004-09-23 20:25:21 mcschatz Exp $
#
umask 2

curr_dir=`pwd`

# for each .asm & .FAILURE file in subdirectories
# get it's prefix, filename, time stamp & print to tab-delimited file


delivery=${curr_dir##/*/}
#delivery=${curr_dir##*_}
#delivery=`echo ${delivery} | cut -b1-8`
delivery="2000-${delivery}"

tag=""
for note in "$@"; do
  tag="${tag}\t${note}"
done

ls * > all_files.txt

for filename in `cat all_files.txt`; do

  prefix=${filename%%.*}
  suffix=${filename##*.}

  filename=${prefix}/${filename}

  if [ "$suffix" = "asm" ] || [ "$suffix" = "FAILURE" ]; then
  
    # determine .asm filename prefix - first make name absolute
    if [ "${filename##/*}" = "${filename}" ]; then
      filename=${curr_dir}/${filename}
    fi
    prefix=${filename##/*/}
    prefix=${prefix%%.*}

    file_month=`ls -l ${filename} | awk '{print $6}'`
    case $file_month in
      Jan ) file_month="01";;
      Feb ) file_month="02";;
      Mar ) file_month="03";;
      Apr ) file_month="04";;
      May ) file_month="05";;
      Jun ) file_month="06";;
      Jul ) file_month="07";;
      Aug ) file_month="08";;
      Sep ) file_month="09";;
      Oct ) file_month="10";;
      Nov ) file_month="11";;
      Dec ) file_month="12";;
    esac

    file_day=`ls -l ${filename} | awk '{print $7}'`
    file_hour=`ls -l ${filename} | awk '{print $8}' | cut -f 1 -d ':'`
  
    if [ "$suffix" = "FAILURE" ]; then
      echo -e "${delivery}\t${prefix}\t${filename%%.*}\tfailed   \t2000-${file_month}-${file_day}-${file_hour}${tag}"
    else
      echo -e "${delivery}\t${prefix}\t${filename}\tsucceeded\t2000-${file_month}-${file_day}-${file_hour}${tag}"
    fi
  fi
  
done

exit 0
#############################################################################
