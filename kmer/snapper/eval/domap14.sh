#!/bin/sh
#$ -p -666
#$ -j y
#$ -o domap14.$TASK_ID.out
#$ -cwd
#$ -N snapper2test
#$ -A snapper2test

#  gn      -- pick one of 1 different genome sequences
#  fg      -- pick one of 8 different fragment sets
#  ms      -- fixed at 28
#  mk      -- pick one of 19 different mer skips
#  ig      -- pick one of 14 different mer thresholds

workd=/project/huref0/assembly/chr14

decodeJob() {
    nn=$1

    xx=`expr $nn % 1 + 1`
    gn=`echo chr14                   | cut -d' ' -f $xx`
    nn=`expr $nn / 1`

    xx=`expr $nn % 8 + 1`
    fg=`echo f1 f2 f3 f4 f5 f6 f7 f8 | cut -d' ' -f $xx`
    nn=`expr $nn / 8`

    xx=`expr $nn % 1 + 1`
    ms=`echo 28                                                                     | cut -d' ' -f $xx`
    nn=`expr $nn / 1`

    xx=`expr $nn % 19 + 1`
    mk=`echo 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19               | cut -d' ' -f $xx`
    nn=`expr $nn / 19`

    xx=`expr $nn % 14 + 1`
    ig=`echo 0000 0001 0002 0004 0008 0016 0032 0064 0128 0256 0512 1024 2048 4096  | cut -d' ' -f $xx`
    nn=`expr $nn / 14`

    name=`echo $gn.$fg.ms$ms.mk$mk.ig$ig`
}



if [ x$1 = "xcheck" ] ; then
    shift
    mm=1
    jj=`expr 1 \* 8 \* 1 \* 19 \* 14`
    echo "Checking $jj jobs."
    while [ $mm -lt $jj ] ; do
        decodeJob $mm
        if [ ! -e $workd/$name.stats ] ; then
            echo qsub -t $mm domap14.sh
        fi
        mm=`expr $mm + 1`
    done
    exit
fi



nn=`expr $SGE_TASK_ID - 1`
decodeJob $nn

if [ ! -e $workd/$name.stats ] ; then
  $workd/src/genomics/snapper/snapper2 \
    -verbose \
    -queries $workd/$fg.fasta \
    -genomic $workd/$gn.fasta \
    -mersize $ms \
    -merskip $mk \
    -ignore  $ig \
    -minmatchidentity 95 \
    -minmatchcoverage 50 \
    -numthreads 2 \
    -noaligns \
    -output   /scratch/$name.sim4db \
    -validate /scratch/$name.validate \
    -stats    /scratch/$name.stats \
  && \
  bzip2 -9v /scratch/$name.sim4db \
  && \
  mv /scratch/$name.sim4db.bz2 $workd/$name.sim4db.bz2 \
  && \
  mv /scratch/$name.validate   $workd/$name.validate \
  && \
  mv /scratch/$name.stats      $workd/$name.stats
fi
