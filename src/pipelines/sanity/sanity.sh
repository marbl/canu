#!/bin/sh

#  (re)Load the sge config.
. $SGE_ROOT/$SGE_CELL/common/settings.sh

#
#  Master controller of nightly sanity checks.  Optional date on command line.
#
#  To start the next nightly, launch with something like:
#    qsub -cwd -j y -o init.err -A assembly-nightly -l fast -a 200907180001.00 -b n sanity.sh 2009-07-18-0001 grid
#
#  To start a previous nightly, the sge hold isn't needed, well, neither is sge for that matter.
#  Just run sanity.sh with the date/time you want to start at.
#
#  CAREFUL!  Submitting something in the past will automagically submit every nightly since that date!
#

date=$1
bins=/work/canu/src/pipelines/sanity
spec=/work/canu/src/pipelines/sanity

if [ x$date = x ] ; then
  date=`date +%Y-%m-%d-%H%M`
fi

echo "SANITY BEGINS for $date at `date`"

#  Remove old versions
#perl $bins/sanity-purge-old.pl purge
#rm -rf DEL

perl $bins/sanity.pl fetch                #  Update the local repo.
perl $bins/sanity.pl checkout $date       #  Checkout from that repo.
perl $bins/sanity.pl build    $date       #  Compile.
perl $bins/sanity.pl submit   $date       #  Submit the next iteration.

#  Run small stuff daily.

for ss in small.bibersteinia_trehalosi.pacbio.h5-1000 \
          small.bibersteinia_trehalosi.pacbio.h5-5000 \
          small.bibersteinia_trehalosi.pacbio.sra-3000 \
          small.escherichia_coli_k12.pacbio.p6 \
          small.escherichia_coli_ne92.pacbio.p4 \
          small.escherichia_coli_ne92.pacbio.p5 \
          small.escherichia_coli_o157_h7_str_f8092b.pacbio.p4c2.average \
          small.escherichia_coli_o157_h7_str_f8092b.pacbio.p4c2.long \
          small.francisella_tularensis.pacbio \
          small.saccharomyces_cerevisiae_glbrcy22-3.pacbio \
          small.saccharomyces_cerevisiae_glbrcy22-3.pacbio.sra \
          small.saccharomyces_cerevisiae_s288c.pacbio
  do
  perl $bins/sanity.pl assemble $date $spec/$ss.spec
done

for ss in small.escherichia_coli_k12.nanopore.all.2d \
          small.escherichia_coli_k12.nanopore.map006-1.2d \
          small.escherichia_coli_k12.nanopore.map006-2.2d \
          small.escherichia_coli_k12.nanopore.map006-pcr-1.2d \
          small.escherichia_coli_k12.nanopore.map006-pcr-2.2d \
          small.escherichia_coli_k12.nanopore.r9.SpotOn.1d
  do
  perl $bins/sanity.pl assemble $date $spec/$ss.spec
done

#  Run big stuff weekly.

if [ `date +%u` = 6 ] ; then
  for ss in medium.arabidopsis_thaliana.pacbio.p4c2 \
            medium.arabidopsis_thaliana.pacbio.p5c3 \
            medium.drosophila_melanogaster.pacbio.p5c3
    do
    perl $bins/sanity.pl assemble $date $spec/$ss.spec
  done
fi

