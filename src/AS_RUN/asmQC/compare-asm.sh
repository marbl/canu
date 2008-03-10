#!/bin/sh

#  Hack to diff two assemblies.

asm1=$1
asm2=$2
pfx=$3

diff="diff -b --side-by-side --suppress-common-lines"

echo "--------------------------------------------------------------------------------"
echo "initialTrimLog"
$diff $asm1/0-overlaptrim/$pfx.initialTrimLog   $asm2/0-overlaptrim/$pfx.initialTrimLog
echo "--------------------------------------------------------------------------------"
echo "consolidated"
$diff $asm1/0-overlaptrim/$pfx.ovl.consolidated $asm2/0-overlaptrim/$pfx.ovl.consolidated
echo "--------------------------------------------------------------------------------"
echo "mergeLog"
$diff $asm1/0-overlaptrim/$pfx.mergeLog         $asm2/0-overlaptrim/$pfx.mergeLog
echo "--------------------------------------------------------------------------------"
echo "chimera.summary"
$diff $asm1/0-overlaptrim/$pfx.chimera.summary  $asm2/0-overlaptrim/$pfx.chimera.summary
echo "--------------------------------------------------------------------------------"
echo "chimera.report"
cat $asm1/0-overlaptrim/$pfx.chimera.report | sed s/gatekeeper\ store// | sed s/fragStore// > $asm1/0-overlaptrim/$pfx.chimera.report.fixed
cat $asm2/0-overlaptrim/$pfx.chimera.report | sed s/gatekeeper\ store// | sed s/fragStore// > $asm2/0-overlaptrim/$pfx.chimera.report.fixed
$diff $asm1/0-overlaptrim/$pfx.chimera.report.fixed   $asm2/0-overlaptrim/$pfx.chimera.report.fixed

echo "--------------------------------------------------------------------------------"
echo "unitigger"
$diff $asm1/4-unitigger/$pfx.cga.0              $asm2/4-unitigger/$pfx.cga.0

echo "--------------------------------------------------------------------------------"
echo "consensus stats"
$diff $asm1/5-consensus/consensus.stats.summary $asm2/5-consensus/consensus.stats.summary
echo "--------------------------------------------------------------------------------"
echo "consensus"
$diff $asm1/5-consensus/$pfx.cgi                $asm2/5-consensus/$pfx.cgi


echo "--------------------------------------------------------------------------------"
echo "cgw final (cgw)"
$diff $asm1/7-CGW/$pfx.cgw                      $asm2/7-CGW/$pfx.cgw
echo "--------------------------------------------------------------------------------"
echo "cgw final (cgw_contigs)"
$diff $asm1/7-CGW/$pfx.cgw_contigs              $asm2/7-CGW/$pfx.cgw_contigs
echo "--------------------------------------------------------------------------------"
echo "cgw final (cgw_scaffolds)"
$diff $asm1/7-CGW/$pfx.cgw_scaffolds            $asm2/7-CGW/$pfx.cgw_scaffolds


echo "--------------------------------------------------------------------------------"
echo "consensus stats"
$diff $asm1/8-consensus/consensus.stats.summary $asm2/8-consensus/consensus.stats.summary

echo "--------------------------------------------------------------------------------"
echo "posmap.frgctg"
$diff $asm1/9-terminator/$pfx.posmap.frgctg     $asm2/9-terminator/$pfx.posmap.frgctg
echo "--------------------------------------------------------------------------------"
echo "posmap.frgscf"
$diff $asm1/9-terminator/$pfx.posmap.frgscf     $asm2/9-terminator/$pfx.posmap.frgscf
echo "--------------------------------------------------------------------------------"
echo "posmap.varscf"
$diff $asm1/9-terminator/$pfx.posmap.varscf     $asm2/9-terminator/$pfx.posmap.varscf

echo "--------------------------------------------------------------------------------"
echo "qc"
$diff $asm1/9-terminator/$pfx.qc                $asm2/9-terminator/$pfx.qc

echo "--------------------------------------------------------------------------------"
echo "asm"
$diff $asm1/9-terminator/$pfx.asm               $asm2/9-terminator/$pfx.asm

