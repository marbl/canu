## Preloader Script ver. 1.2
## By Marwan Oweis, Feb 2007
##
## Files Required
##REQUIRED PREFIX.asm PREFIX.seq.features PREFIX.frg PREFIX.insert

cd $WORKDIR

#Produce a seq & qual file
$COMMON_BIN/subfrag $WORKDIR/$PREFIX.frg -o $PREFIX.seq,$PREFIX.qual

#> Generate scaff file
$COMMON_BIN/ca2scaff -i $PREFIX.asm -s -o $PREFIX

#> Produce BCP file (3.1)
$COMMON_BIN/reclip $WORKDIR/$PREFIX.seq.features $WORKDIR/$PREFIX.asm -o $PREFIX --ref rid --metrics

#> Obtain .tasm file (3.2.1)
$COMMON_BIN/ca2ta -debug 0 -report -coveragestats -o $PREFIX $WORKDIR/$PREFIX.asm
cp $PREFIX.tasm $PREFIX.tasm.orig

ln -s $PREFIX.tasm tmpA.tasm
ln -s $PREFIX.surrogates.fasta tmpA.surrogates.fasta
ln -s $PREFIX.fasta tmpA.fasta
ln -s $PREFIX.surrogates.contig tmpA.surrogates.contig
ln -s $PREFIX.contig tmpA.contig
ln -s $PREFIX.nuu.tasm tmpA.nuu.tasm
#Generate a list of surrogate contig ids
grep '>' tmpA.surrogates.fasta | sed 's/>//g' | awk '{print $1}' | sort > tmpA.surrogates.contigs
#Generate a list of main contig ids
grep '>' tmpA.fasta | sed 's/>//g' | awk '{print $1}' | sort > tmpA.contigs
#Generate a list of reads placed in surrogate contigs
grep '#' tmpA.surrogates.contig | grep -v '##' | sed 's/#//g' | cut -d '(' -f 1 | sort > tmpA.surrogates.seqs
#enerate a list of reads placed in contigs
grep '#' tmpA.contig | grep -v '##' | sed 's/#//g' | cut -d '(' -f 1 | sort > tmpA.seqs
#Find the intersection of reads between the contigs and the surrogates
comm -12 tmpA.seqs tmpA.surrogates.seqs > tmpA.resolved.seqs
#Check if there are resolved reads
if [ -s tmpA.resolved.seqs ] ; \
then echo 'There are redunantly place (resolved) reads. Cleaning the contigs .tasm file'; \
grep asmbl_id tmpA.nuu.tasm | awk '{ print $2 }' | sort > tmpA.nuu.contigs ;\
$COMMON_BIN/subtasm -X tmpA.nuu.contigs tmpA.tasm -o tmpB ;\
else \
ln -s tmpA.tasm tmpB.tasm; echo 'There are no redundantly placed (resolved) reads'; \
fi

#> Produce the surrogates .tasm file (3.2.3)
ln -s $PREFIX.seq tmpA.seq
#Generate a .tasm file from the surrogates
$COMMON_BIN/convertContig --comment "CA_FREE"  tmpA.surrogates.contig tmpA.surrogates.tasm > convertContig.log
#Strip out the resolved reads from the surrogates
$COMMON_BIN/subtasm tmpA.surrogates.tasm -x tmpA.resolved.seqs -o tmpB.surrogates

#> Combine the .tasm files for uploading (3.2.4)
$COMMON_BIN/subtasm tmpB.tasm tmpB.surrogates.tasm -o tmpC
mv tmpC.tasm $PREFIX.tasm

#> Produce the .library file (3.3)
$COMMON_BIN/restat $PREFIX.asm $PREFIX.insert -o $PREFIX.insert,$PREFIX.restat.metrics,$PREFIX.library

#> Produce .sq.contigs file (3.6)
#Create the scaffold layout file ($PREFIX.sq.contigs)
#casper doesn't accept absolute/relative path for asm file
#requires $PREFIX.asm and $PREFIX.fasta
$COMMON_BIN/casper -l $PREFIX
#Build the insert-contig layout files
export ACTMPDIR=/local/aserver_new/var/AC_WORK
$COMMON_BIN/buildTiling -l $PREFIX > buildTiling.log

#> Remove temporary files
rm -f tmpA.* tmpB.* *.frg *.asm *.insert *.seq.features
mv $PREFIX.nuu.tasm $PREFIX.tasm.nuu
mv $PREFIX.orig.contig $PREFIX.contig.orig
mv $PREFIX.sq.fasta $PREFIX.fasta.sq