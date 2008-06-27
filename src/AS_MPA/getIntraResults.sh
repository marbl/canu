#!/usr/bin/env bash
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
# $Id: getIntraResults.sh,v 1.7 2008-06-27 06:29:17 brianwalenz Exp $
#

AS=${1}
LIBFILE=${2}

Type=(stretched compressed outtie normal antinormal inversion transposition)
Status=(raw confirmed fewer)

# process files in the directory
# create files for each unsatisfied mate pair type
#  1. ${AS}_${t}.csv - comma-delimited file listing raw, confirmed, & fewer
#        for import into excel spreadsheet
#    header line: chromosome, raw, confirmed, fewer
#    other lines:   # , # , # , #
#  2. ${AS}.${t}.lengthsWeights.txt file summarizing lengths & weights of confirmed bad intervals
# Also create a ${AS}_satisfied.csv file

function ProcessDirFiles
{
  # count the number of chromosome input files - if it's only 1, treat differently
  declare -i numFiles=`ls | egrep "${AS}_(.+)_intra.txt" | gawk 'END{print NR}'`

  # process each type of unsatisfied mate pair/bad interval
  for t in "${Type[@]}"; do
    string=""
    echo "Scanning ${t}"

    # do each category of mate pairs within this type
    for s in "${Status[@]}"; do
      echo "  scanning ${s}"

      # a couple filenames
      fn="${AS}.${t}.${s}.txt"
      fnA="${fn}A"

      # delete temporary files to be written to
      if [ -f ${fnA} ]; then
        rm ${fnA}
      fi
      if [ -f ${fn} ]; then
        rm ${fn}
      fi
      echo "${s}" > ${fn}

      # do raw numbers separately
      if [ ${s} == "raw" ]; then

        # inversions & transpositions don't have raw mate pairs
        if [ ${t} != "inversion" ] && [ ${t} != "transposition" ]; then

          # create line with two lines: chromosome# and count of raw mate pairs
          wc -l `ls | egrep "${AS}\.(.+)\.${t}\.raw"` | sed 's/[.]/ /g' | gawk '{print $3,$1}'  >> ${fnA}

          # count number of chromosomes/lines in the temporary file
          # if there are more than 1 .raw files, wc will add a total line
          # declare -i a=`wc -l ${fnA} | gawk '{print $1}'`
          if [ ${numFiles} -gt 1 ]; then
            head -n ${numFiles} ${fnA} |sort -n >> ${fn}
          else
            cat ${fnA} >> ${fn}
          fi
          # rm ${fnA}
        else

          # for inversions & transpositions, create dummy lines
          # use gawk because numbers may confuse bash with leading 0s
          for file in `ls | egrep "${AS}_(.+)_intra.txt"`; do
            chr=${file%_*}
            chr=${chr##*_}
            gawk -v c=${chr} 'BEGIN{print c, 0}' >> ${fn}
          done
        fi

      else

        # for non-raw types, scan .err files
        # if there multiple files, egrep will prefix the matching line with the filename
        if [ ${numFiles} -gt 1 ]; then
          egrep " ${t}" ${AS}.*.intra.err |egrep ${s} |sed 's/[:.]/ /g' |gawk '{print $2, $5}' |sort -n >> ${fn}
        else
          filename=`ls |egrep "${AS}_(.+)_intra.txt"`
          chr=${filename%_*}
          chr=${chr##*_}
          egrep " ${t}" ${AS}.*.intra.err |egrep ${s} |gawk -v c=${chr} '{print c, $1}' >> ${fn}
        fi
      fi
      string="${string} ${fn}"
    done

    # paste together the raw, confirmed, & fewer temporary files
    # and strip out extra columns identifying the chromosome number
    paste -d ' ' ${string} |gawk '{if(NR==1){print "chromosome ,", $1, "," $2, "," $3}else{print $1, ",", $2, ",", $4, ",", $6}}' > ${AS}_${t}.csv
    # rm ${string}

    echo "  getting lengths & weights"
    if [ ${numFiles} -gt 1 ]; then
      egrep weight ${AS}.*.${t}.ata |sed 's/[=.]/ /g' |gawk '{print $2, $9, $14}' > ${AS}.${t}.lengthsWeights.txt
      egrep "raw satisfied" ${AS}.*.intra.err |sed 's/[:.]/ /g' |gawk '{print $1, ",", $3}' > ${AS}_satisfied.csv
    else
      filename=`ls | egrep "${AS}_(.+)_intra.txt"`
      chr=${filename%_*}
      chr=${chr##*_}
      egrep weight ${AS}.*.${t}.ata |sed 's/[=.]/ /g' |gawk -v c=${chr} '{print c, $6, $11}' > ${AS}.${t}.lengthsWeights.txt
      egrep "raw satisfied" ${AS}.*.intra.err |sed 's/[:.]/ /g' |gawk -v c=${chr} '{print c, ",", $1}' > ${AS}_satisfied.csv

    fi
  done
}

ProcessDirFiles

getLibSpecifics.sh ${AS} ${LIBFILE}