#!/bin/sh

#  1) run clumpMaker
#  2) I expect the prefix of the *.clumps file

awk -f ../clump_diags.awk $1.clumps > $1.diags
awk -v fudge=200000 -f ../near_clumps.awk < $1.clumps > $1.near
awk '$14!=-1{inclump+=$7}$14==-1{if($15~/near/){if($12==-1){revnear+=$7}else{fwdnear+=$7}}else{far+=$7}}END{print inclump+revnear+fwdnear+far,inclump,fwdnear,revnear,far}' $1.near
