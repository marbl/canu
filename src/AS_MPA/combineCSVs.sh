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
# $Id: combineCSVs.sh,v 1.1.1.1 2004-04-14 13:51:58 catmandew Exp $
#


######################################################################
function MakeInterSummaryFile
{
  cp ${ttFile} ${1}
  gawk -v s=${2} -v l=${3} '{if(NR>1){printf("%d ", $1);for(i=s;i<=NF;i+=l)printf(", %d ", $i);printf("\n")}}' tempWorking.txt >> ${1}
}


######################################################################
function PopulateWorkingFile
{
  if [ -f tempWorking.txt ]; then
    rm tempWorking.txt
  fi
  touch tempWorking.txt
  for assembly in "${AllAssemblies[@]}"; do
    paste -d ' ' tempWorking.txt ${assembly}/${2}/${assembly}_${1}.csv > tempWorking2.txt
    mv tempWorking2.txt tempWorking.txt
  done
}


######################################################################
function DoRaw
{
  PopulateWorkingFile ${1} intraChromosome
  
  gawk '{printf("%d ", $1);for(i=3;i<=NF;i+=7){printf(", %d ", $i)}printf("\n")}' tempWorking.txt > ALL_${1}Raw.csv
  
  # do unmapped ${1} Raw for Celera assemblies
  rm tempWorking.txt
  for assembly in "${CeleraAssemblies[@]}"; do
    gawk 'BEGIN{a=0}{a+=$3}END{print a}' ${assembly}/intraChromosome/unmapped/${assembly}_${1}.csv >> tempWorking.txt
  done
  gawk '{a[NR]=$1}END{printf("25 "); for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking2.txt
  cat ALL_${1}Raw.csv tempWorking2.txt > tempWorking.txt
  tail +2 tempWorking.txt > tempWorking2.txt
  cat ${ttFile} tempWorking2.txt > ALL_${1}Raw.csv
}


######################################################################
function DoConfirmed
{
  PopulateWorkingFile ${1} intraChromosome
  
  gawk '{printf("%d ", $1);for(i=5;i<=NF;i+=7){printf(", %d ", $i)}printf("\n")}' tempWorking.txt > ALL_${1}Confirmed.csv
  
  # do unmapped ${1} Confirmed for Celera assemblies
  rm tempWorking.txt
  for assembly in "${CeleraAssemblies[@]}"; do
    gawk 'BEGIN{a=0}{a+=$5}END{print a}' ${assembly}/intraChromosome/unmapped/${assembly}_${1}.csv >> tempWorking.txt
  done
  gawk '{a[NR]=$1}END{printf("25 "); for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking2.txt
  cat ALL_${1}Confirmed.csv tempWorking2.txt > tempWorking.txt
  tail +2 tempWorking.txt > tempWorking2.txt
  cat ${ttFile} tempWorking2.txt > ALL_${1}Confirmed.csv
}


######################################################################
function DoConfirmedTranspositions
{
  PopulateWorkingFile transposition intraChromosome
  
  gawk '{printf("%d ", $1);for(i=5;i<=NF;i+=6){printf(", %d ", $i)}printf("\n")}' tempWorking.txt > ALL_${1}Confirmed.csv
  
  # do unmapped ${1} Confirmed for Celera assemblies
  rm tempWorking.txt
  for assembly in "${CeleraAssemblies[@]}"; do
    gawk 'BEGIN{a=0}{a+=$5}END{print a}' ${assembly}/intraChromosome/unmapped/${assembly}_${1}.csv >> tempWorking.txt
  done
  gawk '{a[NR]=$1}END{printf("25 "); for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking2.txt
  cat ALL_${1}Confirmed.csv tempWorking2.txt > tempWorking.txt
  tail +2 tempWorking.txt > tempWorking2.txt
  cat ${ttFile} tempWorking2.txt > ALL_${1}Confirmed.csv
}


######################################################################
function DoStretchedCompressed
{
  DoRaw ${1}

  DoConfirmed ${1}

  # count double-counted in transpositions
  # do celera separately from public
  rm tempWorking.txt
  for assembly in "${CeleraAssemblies[@]}"; do
    gawk -v f=${1} 'BEGIN{a=0}{if($2==f)a+=$1}END{print a}' ${assembly}/intraChromosome/${assembly}_doubleCountedInTranspositions.txt ${assembly}/intraChromosome/unmapped/${assembly}_doubleCountedInTranspositions.txt >> tempWorking.txt
  done
  gawk '{a[NR]=$1}END{printf("double counted in transpositions "); for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking2.txt

  # do public double-counted in transpositions
  rm tempWorking.txt
  for assembly in "${PublicAssemblies[@]}"; do
    gawk -v f=${1} 'BEGIN{a=0}{if($2==f)a+=$1}END{print a}' ${assembly}/intraChromosome/${assembly}_doubleCountedInTranspositions.txt >> tempWorking.txt
  done
  gawk '{a[NR]=$1}END{for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking3.txt
  paste -d ' ' tempWorking2.txt tempWorking3.txt > tempWorking.txt
  cat ALL_${1}Confirmed.csv tempWorking.txt > tempWorking2.txt
  cat ${ttFile} tempWorking2.txt > ALL_${1}Confirmed.csv
  rm tempWorking3.txt

  # get lib-specific numbers for all - celera separate from public
  rm tempWorking.txt
  touch tempWorking.txt
  for assembly in "${CeleraAssemblies[@]}"; do
    paste -d ' ' tempWorking.txt ${assembly}/intraChromosome/${assembly}_libSummary.csv > tempWorking2.txt
    paste -d ' ' tempWorking2.txt ${assembly}/intraChromosome/unmapped/${assembly}_libSummary.csv > tempWorking.txt
  done
  if [ ${1} == "stretched" ]; then
    gawk '{if(NR>1){printf("%s , %s ", $1, $3);for(i=9;i<=NF;i+=22)printf(", %d ", $i+$(i+11));printf("\n")}}' tempWorking.txt > ALL_${1}ByLib.csv
  else
    gawk '{if(NR>1){printf("%s , %s ", $1, $3);for(i=5;i<=NF;i+=22)printf(", %d ", $i+$(i+11));printf("\n")}}' tempWorking.txt > ALL_${1}ByLib.csv
  fi

  # do public
  rm tempWorking.txt
  touch tempWorking.txt
  for assembly in "${PublicAssemblies[@]}"; do
    paste -d ' ' tempWorking.txt ${assembly}/intraChromosome/${assembly}_libSummary.csv > tempWorking2.txt
    mv tempWorking2.txt tempWorking.txt
  done
  if [ ${1} == "stretched" ]; then
    gawk '{if(NR>1){for(i=9;i<=NF;i+=11)printf(", %d ", $i);printf("\n")}}' tempWorking.txt > tempWorking3.txt
  else
    gawk '{if(NR>1){for(i=5;i<=NF;i+=11)printf(", %d ", $i);printf("\n")}}' tempWorking.txt > tempWorking3.txt
  fi
  paste -d ' ' ALL_${1}ByLib.csv tempWorking3.txt > tempWorking.txt
  cat ${ttFile} tempWorking.txt > ALL_${1}ByLib.csv
}


if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

# mapped & unmapped for Celera
CeleraAssemblies=(CSAB WGAB VAN VAN_asm R26 R27)

# just mapped for Public
PublicAssemblies=(HG05 HG06 B28 B33A B34)

AllAssemblies=(CSAB WGAB VAN VAN_asm R26 R27 HG05 HG06 B28 B33A B34)

# satisfied raw
# stretched raw & confirmed & by library & double-counted
# compressed raw & confirmed & by library & double-counted
# normals & antinormal raw & confirmed inversion
# outtie raw & confirmed transposition

# inter-chromosome intervals & bps
# mappable intervals & bps

cd ${DATA_DIR}
ttFile=assemblyTitlesFile.txt

######################################################################
# make a titles file for concatenation
if [ -f ${ttFile} ]; then
  rm ${ttFile}
fi
echo -n "Chromosome " >> ${ttFile}
for assembly in "${AllAssemblies[@]}"; do
  echo -n ", ${assembly} " >> ${ttFile}
done
echo "" >> ${ttFile}
  

######################################################################
# satisfied: two columns - chrom & count
PopulateWorkingFile satisfied intraChromosome

gawk '{printf("%d ", $1);for(i=3;i<=NF;i+=3){printf(", %d ", $i)}printf("\n")}' tempWorking.txt > tempWorking2.txt
cat ${ttFile} tempWorking2.txt > ALL_satisfied.csv

# do unmapped satisfieds for Celera assemblies
rm tempWorking.txt
for assembly in "${CeleraAssemblies[@]}"; do
  gawk 'BEGIN{a=0}{a+=$3}END{print a}' ${assembly}/intraChromosome/unmapped/${assembly}_satisfied.csv >> tempWorking.txt
done

gawk '{a[NR]=$1}END{printf("25 "); for(i=1;i<=NR;i++){printf(", %d ", a[i])}printf("\n")}' tempWorking.txt > tempWorking2.txt
cat ${ttFile} ALL_satisfied.csv tempWorking2.txt > tempWorking.txt
mv tempWorking.txt ALL_satisfied.csv

######################################################################
# stretched: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoStretchedCompressed stretched

######################################################################
# compressed: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoStretchedCompressed compressed

######################################################################
# antinormal: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoRaw antinormal

######################################################################
# normal: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoRaw normal

######################################################################
# inversion: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoConfirmed inversion

######################################################################
# outtie: 4 columns (with header line) - chrom, raw, confirmed, fewer
DoRaw outtie

######################################################################
# transposition: 3 columns (with header line) - chrom, raw, confirmed
DoConfirmedTranspositions transposition

######################################################################
# inter-chromosome intervals & bps
PopulateWorkingFile mappedInterSummary interChromosome
MakeInterSummaryFile ALL_mappedInterChromosomeToIntervals.csv 3 17
MakeInterSummaryFile ALL_mappedInterChromosomeToBasepairs.csv 5 17

MakeInterSummaryFile ALL_mappedInterChromosomeFromIntervals.csv 7 17
MakeInterSummaryFile ALL_mappedInterChromosomeFromBasepairs.csv 9 17

MakeInterSummaryFile ALL_mappedInterChromosomeRepeatIntervals.csv 11 17
MakeInterSummaryFile ALL_mappedInterChromosomeRepeatBasepairs.csv 13 17

MakeInterSummaryFile ALL_mappedInterChromosomeUnknownIntervals.csv 15 17
MakeInterSummaryFile ALL_mappedInterChromosomeUnknownBasepairs.csv 17 17

######################################################################
# mappable intervals & bps
PopulateWorkingFile unmappedInterSummary interChromosome/unmapped
MakeInterSummaryFile ALL_unmappedInterChromosomeToIntervals.csv 3 17
MakeInterSummaryFile ALL_unmappedInterChromosomeToBasepairs.csv 5 17

MakeInterSummaryFile ALL_unmappedInterChromosomeFromIntervals.csv 7 17
MakeInterSummaryFile ALL_unmappedInterChromosomeFromBasepairs.csv 9 17

MakeInterSummaryFile ALL_unmappedInterChromosomeRepeatIntervals.csv 11 17
MakeInterSummaryFile ALL_unmappedInterChromosomeRepeatBasepairs.csv 13 17

MakeInterSummaryFile ALL_unmappedInterChromosomeUnknownIntervals.csv 15 17
MakeInterSummaryFile ALL_unmappedInterChromosomeUnknownBasepairs.csv 17 17
