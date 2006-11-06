#!/bin/sh
#$ -p -333
#$ -j y
#$ -o crap.$TASK_ID
#$ -cwd
#$ -N snapper2test
#$ -A snapper2test

#  mersize -- fixed at 28
#  piece   -- pick one of 8 different input sets
#  ignore  -- pick one of 14 different mer thresholds
#  skip    -- pick one of 19 different mer skips

nn=`expr $SGE_TASK_ID - 1`

xx=`expr $nn % 8 + 1`
piece=`echo f1 f2 f3 f4 f5 f6 f7 f8 | cut -d' ' -f $xx`
nn=`expr $nn / 8`

xx=`expr $nn % 14 + 1`
ignore=`echo 0000 0001 0002 0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096  | cut -d' ' -f $xx`
nn=`expr $nn / 14`

xx=`expr $nn % 19 + 1`
skip=`echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19                 | cut -d' ' -f $xx`
nn=`expr $nn / 19`

name=`echo $piece.skp$skip.ign$ignore`
echo $name

if [ ! -e /project/huref0/assembly/chr14/out14/chr14.$name.stats ] ; then
  /project/huref0/assembly/chr14/src/genomics/snapper/snapper2 \
    -verbose \
    -mersize 28 -merskip $skip \
    -ignore $ignore \
    -minmatchidentity 95 \
    -minmatchcoverage 50 \
    -queries /project/huref0/assembly/chr14/$piece.fasta \
    -genomic /project/huref0/assembly/chr14/chr14.fasta \
    -numthreads 2 \
    -noaligns \
    -output /scratch/$name.sim4db \
    -stats  /scratch/$name.stats \
  && \
  bzip2 -9v /scratch/$name.sim4db \
  && \
  mv /scratch/$name.sim4db.bz2 /project/huref0/assembly/chr14/out14/chr14.$name.sim4db.bz2 \
  && \
  mv /scratch/$name.stats      /project/huref0/assembly/chr14/out14/chr14.$name.stats
fi
