#!/bin/sh
#
#  Given a set of trimmed reads, map them to reference, filter out any that do not
#  map completely, and generate an output set of reads that is ordered and oriented.
#

if [ -d "/work/scripts" ] ; then
  bin=/work/wgspb/FreeBSD-amd64/bin
  scp=/work/pacbio-overlapper/scripts
  export PATH=${PATH}:/work/scripts
fi

if [ -d "/usr/local/projects/CELERA/bri/scripts" ] ; then
  bin=/usr/local/projects/CELERA/bri/wgspb/Linux-amd64/bin
  scp=/usr/local/projects/CELERA/bri/pacbio-assembler/scripts
  export PATH=${PATH}:/usr/local/projects/CELERA/bri/scripts
fi

if [ -z $bin ] ; then
  echo scripts not found.
  exit
fi

if [ -z $1 ] ; then
  echo usage: $0 work-directory reads.fastq reference.fasta
  exit
fi

wrk=$1      #  Working directory name
pfx="test"  #  Prefix of the files
inp=$2      #  Path to input FASTQ reads
ref=$3      #  Path to input FASTA reference

if [ ! -d $wrk ] ; then
  mkdir $wrk
fi

if [ ! -e $inp ] ; then
  echo "Failed to find input FASTQ reads in '$inp'"
  exit
fi

if [ ! -e $ref ] ; then
  echo ""
  echo "WARNING: no reference, will not build ordered assembly."
  echo
fi


if [ ! -e "$wrk/build-ordered-reads.spec" ] ; then
  echo  $wrk/build-ordered-reads.spec
  cat > $wrk/build-ordered-reads.spec \
<<EOF
fakeUIDs      = 1

useGrid       = 0
scriptOnGrid  = 0

sge                   = -l medium
sgeScript             = -pe threaded 1 -l memory=8g
sgeOverlap            = -pe threaded 4 -l memory=2g
sgeMerOverlapSeed     = DIE
sgeMerOverlapExtend   = DIE
sgeConsensus          = -pe threaded 1 -l memory=6g
sgeFragmentCorrection = DIE
sgeOverlapCorrection  = DIE

ovlConcurrency        = 3
cnsConcurrency        = 12


ovlErrorRate          = 0.30

utgGraphErrorRate     = 0.30
utgMergeErrorRate     = 0.30

cgwErrorRate          = 0.30
cnsErrorRate          = 0.30


merylMemory           = DIE
merylThreads          = 4

merThreshold          = 0

merSize               = 17

doOBT                 = 0
doFragmentCorrection  = 0

ovlThreads            = 4
ovlHashBits           = 23
ovlHashBlockLength    = 10000000
ovlRefBlockSize       = 1000

unitigger             = bogart

cnsMinFrags           = 10
cnsPartitions         = 256
EOF
fi


########################################
#
#  Build a gkpStore for the reads, then dump for mapping
#
#  Map reads to reference, build overlaps, filter

if [ ! -d $wrk/BO ] ; then
  mkdir $wrk/BO
fi

if [ ! -d "$wrk/BO/$pfx.gkpStore" ] ; then
  echo fastqToCA
  $bin/fastqToCA \
    -libraryname L \
    -technology sanger \
    -type sanger \
    -reads $inp \
  > $wrk/BO/$pfx.frg

  echo gatekeeper -T -F -o $wrk/BO/$pfx.gkpStore $wrk/BO/$pfx.frg
  $bin/gatekeeper -T -F -o $wrk/BO/$pfx.gkpStore $wrk/BO/$pfx.frg
fi

if [ ! -e "$wrk/BO/$pfx.fastq" ] ; then
  echo gatekeeper -dumpfastq $wrk/BO/$pfx $wrk/BO/$pfx.gkpStore
  $bin/gatekeeper -dumpfastq $wrk/BO/$pfx $wrk/BO/$pfx.gkpStore

  awk '{ print $1 }' < $wrk/BO/$pfx.unmated.fastq > $wrk/BO/$pfx.fastq

  rm -f $wrk/BO/$pfx.unmated.fastq
  rm -f $wrk/BO/$pfx.paired.fastq
  rm -f $wrk/BO/$pfx.1.fastq
  rm -f $wrk/BO/$pfx.2.fastq

  echo replaceUIDwithName $wrk/BO/$pfx.gkpStore.fastqUIDmap $wrk/BO/$pfx.fastq
  $bin/replaceUIDwithName $wrk/BO/$pfx.gkpStore.fastqUIDmap $wrk/BO/$pfx.fastq
fi

#
#  OLD blasr options
#
#    -maxLCPLength    15
#    -nCandidates     25
#    -maxScore       -500
#

if [ -e $ref ] ; then
  if [ ! -e "$wrk/BO/$pfx.blasr.badnm.sam" ] ; then
    echo \
    blasr \
      -noSplitSubreads \
      -nproc           4 \
      -minMatch        12 \
      -bestn           10 \
      -nCandidates     25 \
      -minPctIdentity  65.0 \
      -sam -clipping soft \
      -out             $wrk/BO/$pfx.blasr.badnm.sam \
      $wrk/BO/$pfx.fastq \
      $ref

    blasr \
      -noSplitSubreads \
      -nproc           4 \
      -minMatch        12 \
      -bestn           10 \
      -minPctIdentity  65.0 \
      -sam -clipping soft \
      -out             $wrk/BO/$pfx.blasr.badnm.sam \
      $wrk/BO/$pfx.fastq \
      $ref
  fi


  if [ ! -e "$wrk/BO/$pfx.blasr.sam" ] ; then
    echo a
    echo \
    samtools calmd \
      -S $wrk/BO/$pfx.blasr.badnm.sam \
      $ref \
    \>  $wrk/BO/$pfx.blasr.sam

    samtools calmd \
      -S $wrk/BO/$pfx.blasr.badnm.sam \
      $ref \
    >  $wrk/BO/$pfx.blasr.sam
  fi


  if [ ! -e "$wrk/BO/$pfx.blasr.sam.coords" ] ; then
   echo bowtie2-to-nucmercoords.pl $ref \< $wrk/BO/$pfx.blasr.sam \> $wrk/BO/$pfx.blasr.sam.coords
   bowtie2-to-nucmercoords.pl $ref < $wrk/BO/$pfx.blasr.sam > $wrk/BO/$pfx.blasr.sam.coords
  fi

  if [ ! -e "$wrk/BO/$pfx.blasr.sam.coords.ova" ] ; then
    echo infer-olaps-from-coords
    perl $scp/infer-olaps-from-genomic-coords.pl \
      $wrk/BO/$pfx.blasr.sam.coords \
      $wrk/BO/$pfx.fastq \
      $wrk/BO/$pfx.blasr.sam.coords \
      $wrk/BO/$pfx.gkpStore.fastqUIDmap
  fi
fi



########################################
#
#  Start an assembly using the blasr overlaps for the perfectly mapping reads.
#
#  We need to rebuild the overlaps to get the IIDs correct for this assembly.
#  If we were to just use the BO/gkpStore and BO/ovlStore, all the reads that
#  don't map perfectly end up as singletons.
#
#  For this to work, it is critical that the reads have their original names, NOT UIDs.
#
#  NOTE!  Reading coords and lengths from the BL work above, and using reads from the BO work here.
#

if [ -e $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq ] ; then
  if [ ! -d $wrk/BL ] ; then
    mkdir $wrk/BL
  fi

  if [ ! -d "$wrk/BL/$pfx.gkpStore" ] ; then
    $bin/fastqToCA \
      -libraryname L \
      -technology sanger \
      -type sanger \
      -reads $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq \
    > $wrk/BL/$pfx.frg

    $bin/gatekeeper -T -F \
      -o $wrk/BL/$pfx.gkpStore \
      $wrk/BL/$pfx.frg
  fi

  if [ ! -e "$wrk/BL/$pfx.blasr.sam.coords.ova" ] ; then
    perl $scp/infer-olaps-from-genomic-coords.pl \
      $wrk/BL/$pfx.blasr.sam.coords \
      $wrk/BO/$pfx.fastq \
      $wrk/BO/$pfx.blasr.sam.coords \
      $wrk/BL/$pfx.gkpStore.fastqUIDmap
  fi

  if [ ! -d "$wrk/BL/$pfx.ovlStore" ] ; then
    $bin/convertOverlap -ovl \
      < $wrk/BL/$pfx.blasr.sam.coords.ova \
      > $wrk/BL/$pfx.blasr.sam.coords.ovb

    $bin/overlapStoreBuild \
      -o $wrk/BL/$pfx.ovlStore \
      -g $wrk/BL/$pfx.gkpStore \
      -F 1 \
      $wrk/BL/$pfx.blasr.sam.coords.ovb
  fi

  if [ ! -d "$wrk/BL/$pfx.tigStore" ] ; then
    perl $bin/runCA -p $pfx -d $wrk/BL -s $wrk/build-ordered-reads.spec \
      useGrid=0 scriptOnGrid=0 \
      stopAfter=unitigger \
      $wrk/BL/$pfx.frg
  else
    if [ ! -e "$wrk/BL/$pfx.qc" ] ; then
      perl $bin/runCA -p $pfx -d $wrk/BL -s $wrk/build-ordered-reads.spec \
        useGrid=0 scriptOnGrid=0 \
        $wrk/BL/$pfx.frg
    fi
  fi
fi



########################################
#
#  Assembly with all reads.
#
if [ ! -d $wrk/CA ] ; then
  mkdir $wrk/CA
fi

if [ ! -d "$wrk/CA/$pfx.frg" ] ; then
  $bin/fastqToCA \
    -libraryname L \
    -technology sanger \
    -type sanger \
    -reads $inp \
  > $wrk/CA/$pfx.frg
fi

if [ ! -d "$wrk/CA/$pfx.tigStore" ] ; then
  perl $bin/runCA -p $pfx -d $wrk/CA -s $wrk/build-ordered-reads.spec \
    useGrid=0 scriptOnGrid=0 \
    stopAfter=unitigger \
    $wrk/CA/$pfx.frg
else
  if [ ! -e "$wrk/CA/$pfx.qc" ] ; then
    perl $bin/runCA -p $pfx -d $wrk/CA -s $wrk/build-ordered-reads.spec \
      useGrid=0 scriptOnGrid=0 \
      $wrk/CA/$pfx.frg
  fi
fi



########################################
#
#  Assembly with all reads.
#
if [ ! -d $wrk/CAcorrected ] ; then
  mkdir $wrk/CAcorrected
fi

if [ ! -d "$wrk/CAcorrected/$pfx.frg" ] ; then
  $bin/fastqToCA \
    -libraryname L \
    -technology sanger \
    -type sanger \
    -reads $inp \
  > $wrk/CAcorrected/$pfx.frg
fi

if [ ! -d "$wrk/CAcorrected/$pfx.tigStore" ] ; then
  perl $bin/runCA -p $pfx -d $wrk/CAcorrected -s $wrk/build-ordered-reads.spec \
    useGrid=0 scriptOnGrid=0 \
    doFragmentCorrection=1 \
    stopAfter=unitigger \
    $wrk/CAcorrected/$pfx.frg
else
  if [ ! -e "$wrk/CAcorrected/$pfx.qc" ] ; then
    perl $bin/runCA -p $pfx -d $wrk/CAcorrected -s $wrk/build-ordered-reads.spec \
      useGrid=0 scriptOnGrid=0 \
      doFragmentCorrection=1 \
      $wrk/CAcorrected/$pfx.frg
  fi
fi





########################################
#
#  Assembly with the perfectly mapping reads.
#

if [ -e $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq ] ; then
  if [ ! -d $wrk/CAordered ] ; then
    mkdir $wrk/CAordered
  fi

  if [ ! -d "$wrk/CAordered/$pfx.frg" ] ; then
    $bin/fastqToCA \
      -libraryname L \
      -technology sanger \
      -type sanger \
      -reads $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq \
    > $wrk/CAordered/$pfx.frg
  fi

  if [ ! -d "$wrk/CAordered/$pfx.tigStore" ] ; then
    perl $bin/runCA -p $pfx -d $wrk/CAordered -s $wrk/build-ordered-reads.spec \
      useGrid=0 scriptOnGrid=0 \
      stopAfter=unitigger \
      $wrk/CAordered/$pfx.frg
  else
    if [ ! -e "$wrk/CAordered/$pfx.qc" ] ; then
      perl $bin/runCA -p $pfx -d $wrk/CAordered -s $wrk/build-ordered-reads.spec \
        useGrid=0 scriptOnGrid=0 \
        $wrk/CAordered/$pfx.frg
    fi
  fi
fi



########################################
#
#  Assembly with the perfectly mapping reads.
#

if [ -e $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq ] ; then
  if [ ! -d $wrk/CAorderedCorrected ] ; then
    mkdir $wrk/CAorderedCorrected
  fi

   if [ ! -d "$wrk/CAorderedCorrected/$pfx.frg" ] ; then
    $bin/fastqToCA \
      -libraryname L \
      -technology sanger \
      -type sanger \
      -reads $wrk/BO/$pfx.blasr.sam.coords.mapped.ordered.fastq \
    > $wrk/CAorderedCorrected/$pfx.frg
  fi

  if [ ! -d "$wrk/CAorderedCorrected/$pfx.tigStore" ] ; then
    perl $bin/runCA -p $pfx -d $wrk/CAorderedCorrected -s $wrk/build-ordered-reads.spec \
      useGrid=0 scriptOnGrid=0 \
      doFragmentCorrection=1 \
      stopAfter=unitigger \
      $wrk/CAorderedCorrected/$pfx.frg
  else
    if [ ! -e "$wrk/CAorderedCorrected/$pfx.qc" ] ; then
      perl $bin/runCA -p $pfx -d $wrk/CAorderedCorrected -s $wrk/build-ordered-reads.spec \
        useGrid=0 scriptOnGrid=0 \
        doFragmentCorrection=1 \
        $wrk/CAorderedCorrected/$pfx.frg
    fi
  fi
fi
