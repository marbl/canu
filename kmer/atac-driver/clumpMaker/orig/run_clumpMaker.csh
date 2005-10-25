#!/usr/local/bin/tcsh -efx

set cmddir=`echo $0|awk -F/ '{for(i=1;i<NF;i++)printf("%s/",$i)}'`

grep '^M [ug]' $1 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > tmp.clumpMaker.$$
set numhits=`cat tmp.clumpMaker.$$ | wc -l`
@ numhits += 1000000
cat tmp.clumpMaker.$$ | $cmddir/clumpHits_${OSTYPE} -n $numhits > $1.clumps
rm tmp.clumpMaker.$$
awk -f  $cmddir/clump_diags.awk $1.clumps > $1.diags
cat $1.diags $1.clumps | awk -v fudge=200000 -f  $cmddir/near_clumps.awk > $1.near
awk '$14!=-1{inclump+=$7}$14==-1{if($15~/near/){if($12==-1){revnear+=$7}else{fwdnear+=$7}}else{far+=$7}}END{print inclump+revnear+fwdnear+far,inclump,fwdnear,revnear,far}' $1.near
