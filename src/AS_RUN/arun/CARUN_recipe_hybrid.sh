## Assembly Script.
## Hybrid script (Using Sanger and 454 reads) ver. 1.0.1

## Files Required
##REQUIRED PREFIX.seq.features PREFIX.catmap PREFIX.frg PREFIX.immutable

#> Gatekeeper

mkdir -p $WORKDIR/0-preoverlap && cd $WORKDIR/0-preoverlap
$CA_BIN/gatekeeper -X -C -P -e 10000000 -Q -T -N \
   -f $WORKDIR/$PREFIX.gkpStore \
   $WORKDIR/$PREFIX.frg \
   > $WORKDIR/0-preoverlap/gatekeeper.out \
   2> $WORKDIR/0-preoverlap/gatekeeper.err

#> Populate Frag Store

cd $WORKDIR/0-preoverlap
$CA_BIN/PopulateFragStore -P -c -f \
   -o $WORKDIR/$PREFIX.frgStore \
   -V $WORKDIR/0-preoverlap/$PREFIX.ofg \
   $WORKDIR/$PREFIX.inp \
   > $WORKDIR/0-preoverlap/populatefragstore.out \
   2> $WORKDIR/0-preoverlap/populatefragstore.err

#> Initial Trim

awk '{print $1,$5,$6}' $WORKDIR/$PREFIX.seq.features > $WORKDIR/$PREFIX.clv

cd $WORKDIR/0-preoverlap
$CA_BIN/initialTrim -update -q 12  \
   -immutable $WORKDIR/$PREFIX.immutable \
   -vector $WORKDIR/$PREFIX.clv  \
   -log $PREFIX.initialTrimLog  \
   -frg $WORKDIR/$PREFIX.frgStore  \
   > initialTrim.err \
   2>&1

#> Meryl

cd $WORKDIR/0-preoverlap
$CA_BIN/meryl -m $KMER -n $NMER -K 10 \
   -s $WORKDIR/$PREFIX.frgStore \
   -o $PREFIX.${KMER}mers${NMER}.fasta \
   > meryl.out \
   2> meryl.err

#> Overlapper, first run  

mkdir -p $WORKDIR/0-overlap && cd $WORKDIR/0-overlap
$CA_BIN/lastfraginstore $WORKDIR/$PREFIX.frgStore  >lastfraginstore.out
HIFRAGID=`sed 's/Last frag in store is iid = //'  < lastfraginstore.out `
# Required options: -M, -h, -r, -o. Do partial overlaps (-G).
$CA_BIN/overlap -P -G -M $MEMORY -t 2 -h 1-$HIFRAGID -r 0-$HIFRAGID \
  -k $WORKDIR/0-preoverlap/$PREFIX.${KMER}mers${NMER}.fasta \
  -o $PREFIX.ovl \
  $WORKDIR/$PREFIX.frgStore  

#> Convert overlaps from text to binary

cd $WORKDIR/0-overlap
$CA_BIN/overlap-convert -b $PREFIX.ovl $PREFIX.ovb 

#> Sort overlaps.

# Also list AvsB and BvsA for conslidate and chimera.
mkdir -p $WORKDIR/0-trim && cd $WORKDIR/0-trim
$CA_BIN/sort-overlaps -memory 1024  -maxiid $HIFRAGID  \
  $WORKDIR/0-overlap/$PREFIX.ovb \
  > $PREFIX.ovl.sorted \
  2> $PREFIX.ovl.sorted.err

#> Consolidate overlaps.

cd $WORKDIR/0-trim
$CA_BIN/consolidate \
  < $PREFIX.ovl.sorted > $PREFIX.ovl.consolidated

#> Overlap-based trimming

cd $WORKDIR/0-trim
$CA_BIN/merge-trimming \
   -immutable $WORKDIR/$PREFIX.immutable \
  -log $PREFIX.mergeLog \
  -frg $WORKDIR/$PREFIX.frgStore \
  -ovl $PREFIX.ovl.consolidated \
  > $PREFIX.ovl.consolidated.err 2>&1

#> Overlap-based chimer detection

cd $WORKDIR/0-trim
$CA_BIN/chimera  \
  -frg $WORKDIR/$PREFIX.frgStore  \
  -immutable $WORKDIR/$PREFIX.immutable \
  -summary $PREFIX.chimera.summary  \
  -report $PREFIX.chimera.report  \
  < $PREFIX.ovl.sorted  2> $PREFIX.chimera.err 

#> Delete links for deleted frags

cd $WORKDIR/0-trim
$CA_BIN/deleteLinks  \
  -f $WORKDIR/$PREFIX.frgStore \ 
  -g $WORKDIR/$PREFIX.gkpStore  \
  > $PREFIX.deletelinks.out 2>&1


#> Overlapper, second run  

mkdir -p $WORKDIR/1-overlap && cd $WORKDIR/1-overlap
$CA_BIN/lastfraginstore $WORKDIR/$PREFIX.frgStore  >lastfraginstore.out
HIFRAGID=`sed 's/Last frag in store is iid = //'  < lastfraginstore.out `
# Required options: -M, -h, -r, -o. No partial overlaps (-G).
$CA_BIN/overlap -P -M $MEMORY -t 2 -h 1-$HIFRAGID -r 0-$HIFRAGID \
  -k $WORKDIR/0-preoverlap/$PREFIX.${KMER}mers${NMER}.fasta \
  -o $PREFIX.ovl \
  $WORKDIR/$PREFIX.frgStore  

#> Create overlap store

cd $WORKDIR/1-overlap
$CA_BIN/grow-olap-store -cfS -M 1024 \
  -o $WORKDIR/$PREFIX.ovlStore   $PREFIX.ovl  \
  > grow-olap-store.out 2> grow-olap-store.err 

#> Fragment correction  

mkdir -p $WORKDIR/2-frgcorr && cd $WORKDIR/2-frgcorr
$CA_BIN/correct-frags -t 2 \
  -S $WORKDIR/$PREFIX.ovlStore -o $PREFIX.frgcorr  \
  $WORKDIR/$PREFIX.frgStore  1  $HIFRAGID  > correct-frags.err 2>&1 
$CA_BIN/cat-corrects  -o $PREFIX.corr  \
  $PREFIX.frgcorr > cat-corrects.out 2> cat-corrects.err

#> Overlap error correction

mkdir -p $WORKDIR/3-ovlcorr && cd $WORKDIR/3-ovlcorr 
$CA_BIN/correct-olaps \
  -S $WORKDIR/$PREFIX.ovlStore   -e $PREFIX.erate \
  $WORKDIR/$PREFIX.frgStore  $WORKDIR/2-frgcorr/$PREFIX.corr  \
  1  $HIFRAGID  > correct-olaps.err 2>&1 
# In single-partition mode, move file to simulate cat-erates
cd $WORKDIR/3-ovlcorr && mv $PREFIX.erate $PREFIX.erates

#> Update overlap error rates

cd $WORKDIR/3-ovlcorr 
$CA_BIN/update-erates $WORKDIR/$PREFIX.ovlStore $PREFIX.erates \
  > update-erates.err 2>&1

#> Prepare files for unitigger

mkdir -p $WORKDIR/4-unitigger && cd $WORKDIR/4-unitigger
$CA_BIN/make_OFG_from_FragStore $WORKDIR/$PREFIX.frgStore  \
  > $PREFIX.ofg  2> make_OFG.err
# Even in single-partition mode, we need a list of files.
ls -1 *.ofg > $PREFIX.ofgList

#> Unitigger

# Note Brian uses -B 75000. Code supports B=0 => single partition.
# The -e 15 parameter equates to 1.5% assumed sequencing error.
# This generates $PREFIX.cgb file.
cd $WORKDIR/4-unitigger 
$CA_BIN/unitigger -c -P -d 1 -x 1 -z 10 -j 5 $GENOMELENGTH -U $BUBBLE -e $ERATE -f \
  -F $WORKDIR/$PREFIX.frgStore  -o $PREFIX.fgbStore  \
  -L $PREFIX.ofgList  -I $WORKDIR/$PREFIX.ovlStore  \
  > unitigger.out  2> unitigger.err 
cp $WORKDIR/4-unitigger/$PREFIX.cga.0 $WORKDIR/.

#> Consensus on unitigs

# For single-partition case, we do not call partitionFragStore.
# Note Brian uses -S for frag store partition.
mkdir -p $WORKDIR/5-consensus && cd $WORKDIR/5-consensus
$CA_BIN/consensus -P -m -U   -o $PREFIX.cgi \
  $WORKDIR/$PREFIX.frgStore  $WORKDIR/4-unitigger/$PREFIX.cgb  \
  > consensus.err 2>&1
mv $PREFIX.cgi $WORKDIR
ln -s $WORKDIR/$PREFIX.cgi .

#> Preliminary scaffolds, estimate mate distances, build SeqStore

mkdir -p $WORKDIR/6-distances && cd $WORKDIR/6-distances
# cgw options: -c for checkpoints, -j and -k for Astat cutoffs,
# -r for rocks, -s for stones, -w for walks,
# -T to ignore transchunks, -P for text output.
# The -o creates $PROJECT.SeqStore directory.
$CA_BIN/cgw  -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 2 -w 0 -T -P    \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.out 2>&1 
# SeqStore needed by subsequent cgw calls. Move it up one level.
mv $PREFIX.SeqStore $WORKDIR
ln -s $WORKDIR/$PREFIX.SeqStore .

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Dump distance estimates

cd $WORKDIR/6-distances
# Next command reads ckp files and SeqStore directory.
$CA_BIN/dumpDistanceEstimates  \
  -u  -n $MAXCHKPT  -p $PREFIX   \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
   > update.dst  2> dumpDistanceEstimates.err

#> Store distance estimates with gatekeeper

cd $WORKDIR/6-distances
$CA_BIN/gatekeeper  -X -Q -C -P \
  -a $WORKDIR/$PREFIX.gkpStore  update.dst  \
  > gatekeeper.err 2>&1

#> Initial scaffolds

mkdir -p $WORKDIR/7-0-scaffold && cd $WORKDIR/7-0-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
$CA_BIN/cgw  -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 0 -w 0 -T -P  \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX    $WORKDIR/$PREFIX.cgi  \
  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Initial extend clear ranges

mkdir -p $WORKDIR/7-1-extend && cd $WORKDIR/7-1-extend
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-0-scaffold/$PREFIX.ckp.$MAXCHKPT .
# Brian specifies a range of scaffolds (-b and -e) but we do not.
# Next command reads ckp files and SeqStore directory.
$CA_BIN/extendClearRanges  -B  -c $PREFIX  -n $MAXCHKPT    \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  > extendClearRanges.err 2>&1
assert_exists $PREFIX.ckp.* > /dev/null
MAXCHKPT=`ls -1 $PREFIX.ckp.* | sed 's/.ckp./ /' | awk '{print $2 ;}' | sort -nr | head -n 1`

#> Secondary scaffolds

mkdir -p $WORKDIR/7-2-scaffold && cd $WORKDIR/7-2-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-1-extend/$PREFIX.ckp.$MAXCHKPT .
# The -R option says to start from previous checkpoint.
# We use -N 3 for all cgw runs except first and last.
$CA_BIN/cgw -y -N 3 -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 0 -w 0 -T -P   \
  -R $MAXCHKPT   \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore   \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Secondary extend clear ranges

mkdir -p $WORKDIR/7-3-extend && cd $WORKDIR/7-3-extend
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-2-scaffold/$PREFIX.ckp.$MAXCHKPT .
# Brian specifies a range of scaffolds (-b and -e) but we do not.
# Next command reads ckp files and SeqStore directory.
# Next command generates another checkpoint file, not used.
$CA_BIN/extendClearRanges  -B  -c $PREFIX  -n $MAXCHKPT     \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore   \
  > extendClearRanges.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Tertiary scaffolds

mkdir -p $WORKDIR/7-4-scaffold && cd $WORKDIR/7-4-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-3-extend/$PREFIX.ckp.$MAXCHKPT .
# The -R option says to start from previous checkpoint.
# We use -N 3 for all cgw runs except first and last.
$CA_BIN/cgw -y -N 3 -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 2 -w 0 -T -P   \
  -R $MAXCHKPT   \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX  $WORKDIR/$PREFIX.cgi  > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Resolve surrogates

mkdir -p $WORKDIR/7-5-resolve && cd $WORKDIR/7-5-resolve 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-4-scaffold/$PREFIX.ckp.$MAXCHKPT .
$CA_BIN/resolveSurrogates   -c $PREFIX  -n $MAXCHKPT  -1   \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  > resolveSurrogates.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Final scaffolds

mkdir -p $WORKDIR/7-6-scaffold && cd $WORKDIR/7-6-scaffold 
ln -s $WORKDIR/$PREFIX.SeqStore .
ln -s $WORKDIR/7-5-resolve/$PREFIX.ckp.$MAXCHKPT .
# We use -N 14 for the very last cgw run.
$CA_BIN/cgw -y -N 14 -c -j $ASTATLOW -k $ASTATHIGH -r 5 -s 2 -w 0 -T -P   \
  -R $MAXCHKPT   \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore   \
  -o $PREFIX   $WORKDIR/$PREFIX.cgi   > cgw.err 2>&1

# Begin : determine max checkpoint number
ls -1 $PREFIX.ckp.*   > pipe01
sed 's/.ckp./ /'      < pipe01 > pipe02 
awk '{print $2 ;}'    < pipe02 > pipe03
sort -nr              < pipe03 > pipe04
MAXCHKPT=`head -n 1 pipe04`
rm     pipe01 pipe02 pipe03 pipe04
# End : determine max checkpoint number

#> Consensus on scaffolds

mkdir -p $WORKDIR/8-consensus && cd $WORKDIR/8-consensus 
ln -s $WORKDIR/7-6-scaffold/$PREFIX.cgw_contigs .
# Brian uses -p and -S and -m to specify partition, but we do not.
$CA_BIN/consensus \
  -P -s $WORKDIR/$PREFIX.SeqStore  -V $MAXCHKPT    \
  -o $PREFIX.cns_contigs  $WORKDIR/$PREFIX.frgStore   \
  $PREFIX.cgw_contigs > consensus.err 2>&1 

#> Terminator

mkdir -p $WORKDIR/9-terminator && cd $WORKDIR/9-terminator 
assert_exists $WORKDIR/7-6-scaffold/$PREFIX.cgw  \
    $WORKDIR/8-consensus/$PREFIX.cns_contigs   \
    $WORKDIR/7-6-scaffold/$PREFIX.cgw_scaffolds > /dev/null
cat $WORKDIR/7-6-scaffold/$PREFIX.cgw  \
    $WORKDIR/8-consensus/$PREFIX.cns_contigs   \
    $WORKDIR/7-6-scaffold/$PREFIX.cgw_scaffolds  \
  | $CA_BIN/terminator -P $EUIDSERVICE    \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore  \
  -o $PREFIX.asm  -m $PREFIX.map  > terminator.err 2>&1 
cd $WORKDIR
ln -s 9-terminator/$PREFIX.asm .

#> Generate QC report
cd $WORKDIR
$CA_BIN/caqc.pl $GENOMELENGTH -metrics $PREFIX.asm > caqc.out 2>&1

#> Generate scaffold FASTA files
cd $WORKDIR
$CA_BIN/asmProcessScaffolds_TER -f   \
  $PREFIX.scaffold.fasta   < $PREFIX.asm
ln -s $PREFIX.scaffold.fasta $PREFIX.scaffolds.fasta

#> Generate singleton FASTA files
cd $WORKDIR
$CA_BIN/dumpSingletons  \
  -f $WORKDIR/$PREFIX.frgStore  -g $WORKDIR/$PREFIX.gkpStore   \
  -c $WORKDIR/7-6-scaffold/$PREFIX -n $MAXCHKPT -U -S  \
  > $PREFIX.singleton.fasta   2> dumpSingletons.err 

#> Generate contig FASTA files (includes degenerates)
cd $WORKDIR
$CA_BIN/asmOutputContigsFasta -d  < $PREFIX.asm > $PREFIX.fasta
