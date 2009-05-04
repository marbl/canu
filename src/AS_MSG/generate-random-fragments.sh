#!/bin/sh

#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2009, J. Craig Venter Instititue.
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
#  Generates a 'random' genomic sequence, and samples fragments of
#  some length at some error (mismatch and indel) from it.
#
#  User unfriendly.
#
#  Parameters, and suggestions:
#
#    genome size 1000000
#    read length 400
#    coverage    10
#    error rate  0.020

if [ $# != 4 ] ; then
  echo "usage $0 <genomesize> <readlength> <coverage> <errorrate>"
  exit
fi

genomesize=$1
readlength=$2
coverage=$3
errorrate=$4

leaff="/work/wgs/kmer/leaff/leaff"
cvtto="perl /work/wgs/src/AS_MSG/convert-fasta-to-v2.pl"

numreads=`expr $genomesize \* $coverage / $readlength`

if [ ! -e genome-${genomesize}bp.fasta ] ; then
  $leaff -G 1 $genomesize $genomesize > genome-${genomesize}bp.fasta
fi

if [ ! -e genome-${genomesize}bp-${coverage}x-${errorrate}percent.frg ] ; then
  $leaff --errors $readlength $numreads 1 $errorrate genome-${genomesize}bp.fasta > genome-${genomesize}bp-${coverage}x-${errorrate}percent.fasta

  tr 'ACGT' 'NNNN' < genome-${genomesize}bp-${coverage}x-${errorrate}percent.fasta > genome-${genomesize}bp-${coverage}x-${errorrate}percent.qlt

  $cvtto -l RAND -mean 0 -stddev 0 \
    -s genome-${genomesize}bp-${coverage}x-${errorrate}percent.fasta \
    -q genome-${genomesize}bp-${coverage}x-${errorrate}percent.qlt \
    >  genome-${genomesize}bp-${coverage}x-${errorrate}percent.frg

  #rm genome-${genomesize}bp-${coverage}x-${errorrate}percent.qlt
  #rm genome-${genomesize}bp-${coverage}x-${errorrate}percent.fasta
fi

