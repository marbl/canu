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
# $Id: getIntraResults.sh,v 1.3 2005-03-22 19:05:53 jason_miller Exp $
#

AS=${1}
Type=(stretched compressed outtie normal antinormal inversion transposition)
Status=(raw confirmed fewer)

# process files in the directory
# create files for each unsatisfied mate pair type
#  1. ${AS}_${t}.csv - comma-delimited file listing raw, confirmed, & fewer
#        for import into excel spreadsheet
#    header line: chromosome, raw, confirmed, fewer
#    other lines:   # , # , # , #
#  2. ${AS}_${t}_lengths_weights.txt file summarizing lengths & weights of confirmed bad intervals
# Also create a ${AS}_satisfied.csv file

function ProcessDirFiles
{

  # count the number of chromosome input files - if it's only 1, treat differently
  declare -i numFiles=`ls [0-9][0-9][0-9].txt | gawk 'END{print NR}'`
  
  # process each type of unsatisfied mate pair/bad interval
  for t in "${Type[@]}"; do
    string=""
    echo "Scanning ${t}"
    
    # do each category of mate pairs within this type
    for s in "${Status[@]}"; do
      echo "  scanning ${s}"

      # delete temporary files to be written to
      if [ -f ${AS}_${t}_${s}.txtA ]; then
        rm ${AS}_${t}_${s}.txtA
      fi
      if [ -f ${AS}_${t}_${s}.txt ]; then
        rm ${AS}_${t}_${s}.txt
      fi
      echo "${s}" > ${AS}_${t}_${s}.txt

      # do raw numbers separately
      if [ ${s} == "raw" ]; then

        # inversions & transpositions don't have raw mate pairs
        if [ ${t} != "inversion" ] && [ ${t} != "transposition" ]; then

          # create line with two lines: chromosome# and count of raw mate pairs
          wc -l ${AS}.[0-9][0-9][0-9].${t}.raw | sed 's/[.]/ /g' | gawk '{print $3+1,$1}' >> ${AS}_${t}_${s}.txtA

          # count number of chromosomes/lines in the temporary file
          # if there are more than 1 .raw files, wc will add a total line
          # declare -i a=`wc -l ${AS}_${t}_${s}.txtA | gawk '{print $1}'`
          if [ ${numFiles} -gt 1 ]; then
            head -n ${numFiles} ${AS}_${t}_${s}.txtA >> ${AS}_${t}_${s}.txt
          else
            cat ${AS}_${t}_${s}.txtA >> ${AS}_${t}_${s}.txt
          fi
          rm ${AS}_${t}_${s}.txtA
        else

          # for inversions & transpositions, create dummy lines
          # use gawk because numbers may confuse bash with leading 0s
          for file in `ls [0-9][0-9][0-9].txt`; do
            chr=${file%%.*}
            gawk -v c=${chr} 'BEGIN{print c+1, 0}' >> ${AS}_${t}_${s}.txt
          done
        fi
        
      else

        # for non-raw types, scan .err files
        # if there multiple files, egrep will prefix the matching line with the filename
        if [ ${numFiles} -gt 1 ]; then
          egrep " ${t}" [0-9][0-9][0-9].err |egrep ${s} |sed 's/[:.]/ /g' |gawk '{print $1+1, $3}' >> ${AS}_${t}_${s}.txt
        else
          filename=`ls [0-9][0-9][0-9].txt`
          chr=${filename%%.*}
          egrep " ${t}" [0-9][0-9][0-9].err |egrep ${s} |gawk -v c=${chr} '{print c+1, $1}' >> ${AS}_${t}_${s}.txt
        fi
      fi
      string="${string} ${AS}_${t}_${s}.txt"
    done

    # paste together the raw, confirmed, & fewer temporary files
    # and strip out extra columns identifying the chromosome number
    paste -d ' ' ${string} |gawk '{if(NR==1){print "chromosome ,", $1, "," $2, "," $3}else{print $1, ",", $2, ",", $4, ",", $6}}' > ${AS}_${t}.csv
    # rm ${string}
  
    echo "  getting lengths & weights"
    if [ ${numFiles} -gt 1 ]; then
      egrep weight ${AS}.[0-9][0-9][0-9].${t}.ata |sed 's/[=.]/ /g' |gawk '{print $2+1, $9, $14}' > ${AS}_${t}_lengths_weights.txt
      egrep "raw satisfied" [0-9][0-9][0-9].err |sed 's/[:.]/ /g' |gawk '{print $1+1, ",", $3}' > ${AS}_satisfied.csv
    else
      filename=`ls [0-9][0-9][0-9].txt`
      chr=${filename%%.*}
      egrep weight ${AS}.[0-9][0-9][0-9].${t}.ata |sed 's/[=.]/ /g' |gawk -v c=${chr} '{print c, $6, $11}' > ${AS}_${t}_lengths_weights.txt
      egrep "raw satisfied" [0-9][0-9][0-9].err |sed 's/[:.]/ /g' |gawk -v c=${chr} '{print c, ",", $1}' > ${AS}_satisfied.csv

    fi
  done
}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ]; then
  echo "Please identify assembly name"
  return
fi

ProcessDirFiles

getLibSpecifics.sh ${AS}