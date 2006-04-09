#!/bin/sh

######################################################################
# Verify that gatekeeper will abort when fed 
# a fragment with an invalid clear range
######################################################################

###
# Name by which this script was invoked.
###
progname=`echo "$0" | sed -e 's/[^\/]*\///g'`

usage="Usage: $progname <gkpbin>

<gkpbin>            Path to gatekeeper binary to test
"

PROGFILE=$(basename $progname)
PROGPRFX=${PROGFILE%%.sh}
RUN_DIR=${PROGPRFX}_run

FRAG_PREF="inv_clr"
FRAG_FILE=${FRAG_PREF}.frg
ERR_FILE=${FRAG_PREF}.err

GKP_BIN=""

PASSED() {
  echo "-----------------------------"; 
  echo "*** SUCCESS - Passed Test ***"; 
  echo "-----------------------------"; echo ""
  exit 0
}

FAILED() {
  echo "-----------------------------"; 
  echo "*** FAILURE - Failed Test ***"; 
  echo "-----------------------------"; echo ""
  exit 1
}

CLEANUP() {
  rm -f ${RUN_DIR}/$FRAG_FILE
  rm -f ${RUN_DIR}/$ERR_FILE
  rmdir ${RUN_DIR}
}

INIT() {
  rm -rf $RUN_DIR
  mkdir $RUN_DIR

  create_frag_file
}

#####
# Function
#   RUN_TEST
#
# Purpose
#   Run the given gatekeeper binary against
#   the frag file within this script. Verify
#   that gatekeeper reports an error. It should
#   since the frag file has an invalid clear range.
#
# Arguments
#   NONE
#
# Return
#   success - NONE
#   failure - exit 1
#####
RUN_TEST() {
  if [ ! -f $GKP_BIN ]; then
    echo "RUN_TEST: $(pwd)/gatekeeper executable not found"
    FAILED
  fi

  if [ ! -f ${RUN_DIR}/$FRAG_FILE ]; then
    echo "RUN_TEST: $(pwd)/${RUN_DIR}/$FRAG_FILE not found"
    FAILED
  fi

  pushd ${RUN_DIR} > /dev/null 2>&1

  # clean up from any previous runs
  rm -rf ${FRAG_PREF}.gkpStore

  $GKP_BIN -X -C -N -Q -P -f \
    ${FRAG_PREF}.gkpStore ${FRAG_PREF}.frg > /dev/null 2>&1

  ret=$?

  if [ "$ret" -eq "0" ]; then
    echo "RUN_TEST: $(pwd)/gatekeeper should have failed on this frag file!"
    exit 1
  fi

  if [ ! -f ${FRAG_PREF}.err ]; then
    echo "RUN_TEST: $(pwd)/gatekeeper should have created a *.err file!"
    FAILED
  fi

  popd > /dev/null 2>&1
}

#####
# Function
#   create_frag_file
#
# Purpose
#   Create the frag file to feed to gatekeeper. This
#   frag file's clear range (clr:1935,2254) is outside
#   the 2k limit and should generate an error from
#   gatekeeper
#
# Arguments
#   NONE
#
# Return
#   success - NONE
#   failure - exit 1
#####
create_frag_file() {
cat > ${RUN_DIR}/$FRAG_FILE <<EOF
{BAT
bna:Celera Assembler
crt:1140792660
acc:1100042561717
com:
scratch assembly
.
}
{FRG
act:A
acc:1100048178561
typ:R
src:
TYABL06TF
.
etm:1057604340
seq:
acagtccgactcagtcatgaacttgaacagagtcataactgtcactagctgaatgcactgaacgcagatc
tgacgcgcatctaagacgatcagacatctgcagtagtttgtatgacctgtagatcagtctaaataggacg
tagctcttgaacttacgcttatgacatcggccctctagcgctcgatattgtcagctctctctctgaatct
cgctgcgtagctcattgccctagcgtatgtgatttatcgtatcattacactggcttgctgtagtgttcct
tagccgttttatcgctgtgtcattctatgacttatatcccgcagtggctgtgcgctcgtcgtttccctac
gggtagctaccacttcgtagttgatatcgatcagactacattaacgcagaggcgcgtgtttcaacagtcg
atgtcagctgatattcagacttcaccattttgctatgtactgtactcgaacacccctcctcgattgccgc
tgctctactacgttcgttattcgcccccaccaccgggtgattgcatatgactgtgcgcacctggcgcatg
agtcagtactagtactaatgaacgttctgttgcatctatatatgtataacactagcggtgctgaagtatc
atctatatgatcgttgtaatttatgattgctactctttttactgctgttctacgtacatgtgagtttgcg
ttgtagatgacttgcgttttgactgactgctatcatcgacgctgttgactcacatacatgtcgctagctt
cgacgtattatgattacgcatatgtgcacttatacgatgcttaattcttgtgctatagagctttcctgct
attgcgatgtatgttctgatcatatatgctttttacgctttttattgtatatgtgttattctttttctct
atctctgctactattctgtctatcctgatgacatgccattgtattttctactcgtttcggcataatgtac
tgtgtgcgtgtttgttctactttgttagctattgcttcttgccaggcgctcgcgtgatcgacgtagctgt
cactctcgctcgcgcgccccctcggtgcatcctgtactcaagcattgtgtcagaggtctttaaggtgtat
gtgtggttttatcgggctgttctttgtggcggtggcttgtttatctacttgtacgctatcggggggcgca
tggataatatacctatactattcctgcgatgtaaacatgcttcttgatcctccacggtattcctctatgc
gcgaattttgtcgtgtaagcgattgatggagtacatctttataattatgtgcgggctgagcggcgatcct
tttagccaccttttctttttgtttgcgtgatgatcactacatatgtttcttgttcgtgcacattcggtgc
ttgccgtgtcgctatttcaattcatgattgtattctgtcttgtctgcttcgtgtgatatggtagtatcct
gcgtctgcttttgcatgcattcggttgcgtgttctgtgatagaactgtgccatgctacagatgtgtacga
tgtcacgtgttacttgattgttatacttcagtcatcgacatatattacttgctataatgactgtcttatc
aatatagaatatctacgtatgatgcacgcgttgatttgtcgtatccgtagatttttattgctacgaatac
tcaagtattatgtgggtttgctaaagctgctctttatagatctcatagctgtatgttatatgtacttgcg
tgttgttcagttctctgctataaatatcacctccctatcttacttactacgtatttctccatctgacttt
agaggaatcacgcgacgtatccgagaacgaccatatcttctacctagacttttatggtattatgctgcgt
cgtatgttgtgtactcatcatgttgatatgcactcctagggtaatatctttgattgctcttgtttctcgt
tagaatcaacacggtcgctgctgttcgtcactgctagtctatcgtccatctaattatgaatcgtggcgag
attcaatgcatctcatgt
.
qlt:
:9;;;;98996996::<::7999:6<==>=8A979;799:::66:9999:96<:::;69967999?@;:9
9<:99998;:;;>=::699969896:99:899999::9997:<:7:7<99::;;9:;9<:9:9:798696
8997<:<:9:9699::99:A9;99699696999<9:::99;996969AA<:9:999668>:78?:68:<8
7:;:9::<<<6;<<:=999::9;99<9::99;999::;9999986999999;99<9;::<::;<9996=9
9::796>;998:;;:869999898898;;<><;:998899:99:::9<;::969996:<9:998669988
;99::9<:69999:6:;;8<::;<9:7;7<:9<78:<::9:::::8699999;996>:9:;;99::::8;
;968769<<8;:9:6::::99998:<:9;9998<:><9:99999998<999:7696989:9::9999999
9968:99699;999899999999:87:9:7:96966999699;;99998;99::69::89:<<::::999
999:<<=::<;99;69:6966669999;9::996:98:6:9;9<@=99:998999:99:::;99999699
999:99:;<9999<:69:::>:999=9:99:<;;:9999:88;69:996996:9<:7:::9:8:999:97
99999:968;;999699999998<;88::66;:9::9;:<:;99:999996:9;::6969;69::99996
:9999899:99699;:699;7:<=>79999;:88::@<=<<=:<999:::9=;<::<9789::9;78;89
9:::9<<>AB:<:::699:979999<:6;6999999:;;:<99<@>9<?::<699998:89;:9::9969
96696:999999977:99;9::<:8=7;;6999<;;79<?6:A:9::9999=9<:9::9996;87:::;9
;999:=:6:9:6::<<8<9999:69:9;7:9:999996==9999::<996;9998<;<7899988;:9:9
99:888878;988;>@@@@KK@>>>>>CKKKKK<KK<@<><>><><<CCK<KC<<<<<<@<<<9<<<K@:
::::::<><>><<@@<@<@@@<>><<<><<>>>>><>C@@<@@@<KK<CCCCCKK<<KKKKKKK<KKK<@
>>>><>><>>>>>>>C<K@<<<<<<>CKTT@>>><><T@>>>>>C::::::<<<::::::::::<KKKKK
K<@@>><<K<><<<<<<<<:777777:<<<KKT<T@<<<<<<>>>>KKKC>>>>><<K<@<>><><K<<<
::::::::<<::::::<>>>>>CC<<K<KKKK<KKKKK@<<<<<@C<<<<::::::@<<<<<<<<<<<@@
C<>>>>>CKKKKK@>>>>>>>><>K<<<<<<<<<<<::::::><C>>>><>KTT<<K@<<<<<@K<<<<<
<<KKK<KC@<@@@G<KK<CCCC>><><<CK<KK<K<KKKKK<<KKKKK<<<KK<KTK<G<G><>><>G<T
TTTTT<TT<T<T<GGG<GG<TTTTT<TT<TTTT<T<<T<TTTTTTG<>>>><<<T<T<TTT<<<G@<<@<
@<<T<<<TTTTTTT<<TTTTTTTTTT<TTTTT<T<TT<TTTTTTT<TTTCC::::<<<T<TTTTTTTTT<
TTTT<TTTTTT<TTTTT<TK<KCC<K<KTTTTTTTTTTT<T<TT<TTTTTTTTTTTTTTTTKK<T<KKK<
KKKKKKTTC@@@<@GT<TTTTT<CCCCCTTTT<TT<<TTTTTTTT<TT<TK<<<<<<CTTKK<KK><777
77:KK<TT<TT<T<<TTTTT<TTT<T<T<<TTTT<TTTTTT<TTTTK>>>>>>KKTTTTT<<TTTTTTTT
TTTTTT<TTTTTTTTT<<TT<TTTT<TTTTTTTTTT<@77777:<TT<<K<K<KKTTTTTT<TTTTTT<T
TTTTTTTTTTT<T<T<<TT<T<TTT<TT<TTT<<TT<T<TT<<TT<TTTT<TTTTTT<TTTTTTTTT<<T
TTTTTTTTTTT<TTTTTT
.
clr:1935,2254
}
EOF

}

#####
# MAIN:
#####
if [ "$#" -ne 1 ]; then
  echo "$usage" 1>&2
  exit 1
fi

GKP_BIN=$1

###
# INIT
###
INIT

###
# BODY
###
RUN_TEST

###
# DONE
###
CLEANUP
PASSED
