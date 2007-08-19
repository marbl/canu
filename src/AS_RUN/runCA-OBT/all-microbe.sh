#!/bin/sh

# qsub -cwd -j y -o /dev/null -A assembly -t 1-108 asmmicrobe.sh

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined ]; then
  jobid=$1
fi
if [ x$jobid = x ]; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line.
  exit 1
fi

if [ -d /home/work/FRAG/$jobid ]; then
  name=$jobid
else
name=`echo M001 M002 M003 M005 M006 M007 M008 M010 M011 M012 M013 M014 M015 M016 M017 M018 M019 M020 M022 M023 M024 M025 M026 M027 M028 M029 M030 M032 M033 M035 M036 M037 M038 M039 M040 M043 M046 M047 M048 M049 M051 M053 M054 M055 M056 M057 M058 M059 M061 M062 M063 M065 M066 M068 M069 M070 M071 M072 M073 M074 M075 M076 M078 M081 M082 M083 M084 M085 M086 M087 M088 M092 M093 M095 M096 M098 M102 M103 M104 M105 M106 M108 M109 M110 M111 M112 M113 M114 M115 M116 M117 M119 M121 M122 M123 M127 M128 M129 M133 M134 M138 M139 M140 M141 M143 M144 M145 M159 \
  | \
  cut -d \  -f $jobid`
fi

#  CHANGE BELOW

outdir=/home/work/moore
bindir=/home/work/wgs/FreeBSD-amd64

mkdir -p $outdir/$name

SGE_TASK_ID=undefined

perl $bindir/bin/runCA-OBT.pl \
  -p $name \
  -d $outdir/$name \
  fakeUIDs=1 \
  vectorIntersect=/home/work/FRAG/moore/$name/$name.vector \
  /home/work/FRAG/moore/$name/$name*.frg.bz2 \
> $outdir/$name/runCA.out 2>&1

