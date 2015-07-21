#!/bin/sh

bgn=`expr $SGE_TASK_ID + 0`
end=`expr $SGE_TASK_ID + 1`

./bin-bound-test $bgn $end


for er in 0100 0200 0300 0400 0500 0600 0700 0800 0900 1000 \
          1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 \
          2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 \
          3100 3200 3300 3400 3500 3600 3700 3800 3900 4000 \
          4100 4200 4300 4400 4500 4600 4700 4800 4900 5000 ; do
  qsub -cwd -j y -o /dev/null -b y ./prefixEditDistance-matchLimit $er
done
