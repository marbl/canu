#!/bin/sh

leaff="./leaff"
source="leaff-test.tmp.fasta"
tmpX="leaff-test.xxx.fasta"

#  Make some sequences.  'U' is missing because it isn't a valid DNA
#  ambiguitity code.
#
cat > ${source} <<EOF
>first sequence seq1 big-unique-thing-on-seq-1
THISISTHEFIRSTSEQENCE
>second sequence seq2 with junk at the end of the defline seq-2-gets-another-unique-thing
THISISTHESECONDSEQENCE
>third sequence seq3 finally-seq-3-has-something
THISISTHETHIRDSEQENCE
>fourth sequence seq4 cant-forget-sequence-number-four
THISISHTEFORTHSEQENCE
>fifth sequence seq5 with junk at the end of the defline but-what-will-we-put-for-sequence-five
THISISTHEFIFTHSEQENCE
EOF


#  Things without an index
#

#  Print ten short random sequences.  The first and third should be
#  errors, but leaff accepts it, silently increasing s to be 1, and
#  rearranging so that s <= l
#
echo ""
echo "Print three sets of ten short random sequences"
${leaff} -G 10  0 12
${leaff} -G 10 12 12
${leaff} -G 10 12 0

echo ""
echo "Dump the test file"
${leaff} -f ${source} -W

#  Should print only seq1, seq3 and seq5
#
echo ""
echo "Printing sequences < 22bp (1, 3, 4 and 5)"
${leaff} -f ${source} -L 0 22

#  Should print only seq2, seq4
#
echo ""
echo "Printing sequences >= 22bp (seq2)"
${leaff} -f ${source} -L 22 100


#  -6 is untested
#  -u is untested

echo ""
echo "Original file checksum"
cat ${source} | md5

echo ""
echo "Twice reversed checksum"
${leaff} -f ${source} -R -W > ${tmpX}
${leaff} -f ${tmpX} -R -W | md5

echo ""
echo "Twice complemented checksum"
${leaff} -f ${source} -C -W > ${tmpX}
${leaff} -f ${tmpX} -C -W | md5

echo ""
echo "Twice reverse-complemented checksum"
${leaff} -f ${source} -R -C -W > ${tmpX}
${leaff} -f ${tmpX} -R -C -W | md5

echo ""
echo "Flipping defline flags on/off/on"
${leaff} -f ${source} -h tmp -H -h tmp2 -H -H -W | md5

#
#  internal id index tests
#
echo ""
echo "Number of sequences in the file"
${leaff} -F ${source} -d

echo ""
echo "Printing all sequences in various orders; all checksums should agree"
sort ${source} | md5
${leaff} -F ${source} -s 0 -s 1 -s 2 -s 3 -s 4 | sort | md5
${leaff} -F ${source} -s 4 -s 3 -s 2 -s 1 -s 0 | sort | md5
${leaff} -F ${source} -s 3 -s 4 -s 1 -s 0 -s 2 | sort | md5
${leaff} -F ${source} -s 0 -s 4 -s 1 -s 3 -s 2 | sort | md5
${leaff} -F ${source} -s 0 -s 2 -s 4 -s 3 -s 1 | sort | md5

echo ""
echo "Printing all sequences using a range"
${leaff} -F ${source} -S 0 4 | sort | md5
${leaff} -F ${source} -S 4 0 | sort | md5

echo ""
echo "Printing all sequences using randomness"
${leaff} -F ${source} -r 5

echo ""
echo "Printing three sequences at random"
${leaff} -F ${source} -r 3

echo ""
echo "Printing three sequences using id's from a file (seq5, seq1 and seq2)"
cat > leaff-test.tmp <<EOF
4 0
1
EOF
${leaff} -F ${source} -q leaff-test.tmp
rm leaff-test.tmp



#
#  external id index tests
#

echo ""
echo "Printing three sequences using NAME id's from a file (seq5, seq1 and seq2)"
cat > leaff-test.tmp <<EOF
fifth
first second
EOF
${leaff} -Fn ${source} -q leaff-test.tmp
rm leaff-test.tmp

echo ""
echo "Printing three sequences using DEFLINE id's from a file (seq5, seq1 and seq2)"
cat > leaff-test.tmp <<EOF
sequence-five
unique-thing-on-seq-1 another
EOF
${leaff} -Fd ${source} -q leaff-test.tmp
rm leaff-test.tmp


