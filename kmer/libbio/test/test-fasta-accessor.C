#include <stdio.h>
#include <stdlib.h>
#include "fasta-accessor.H"



//  rev is the reverse complement of fwd, so if we make an accessor
//  from it, and set the reverse-complement flag we should get back
//  exactly fwd.
//
void
simpleTest(void) {

  //  Yup, real sequence.  I don't remember what it belongs to.

  char *fwd =
    "CGCTTTAATGGGCGAACAGCCCAACCCTTGGGAGAGACTCCACCCCCAGGATGCGACGAGCCGACATCGAGGTGCCAAAC"
    "CATGCCGTCGATATGGACTCTGAGCGACGCCGCTTCCACAAGCCAGCGCCGGGTCACTAGTTCCGACTTTCGTCCCTGCT"
    "CGACCTGCCGGTCTCGCAGTCAAGCTCCCTTGTGCACTTAGCCTCCGTTACCTTTTAGGAGGCAACCGCCCCAGTTAAAC"
    "TACCCACCAGGCAATGTCCCTGATCCGGATCACGGACCTAGGTTAGATATCCAGAACGACGCTTCACAGTCTCCCACCTA"
    "TCCTACACAAGCCGTACCGAACACCAATACCAAGCTATAGTAAAGGTCCCGGGGTCTTTCCGTCCTGCTGCGCGTAACGA"
    "CAGTGGAGAAGTCGTTACGCCATTCGTGCAGGTCGGAACTTACCCGACAAGGAATTTCGCTACCTTAGGATGGTTATAGT"
    "TACCACCGCCGTTTACTGGGTTAACCTTCCAGCACCGGGCAGGCGTCAGTCCGTATACATCGTCTTGCGACTTAGCACGG"
    "ACCTGTGTTTTTAGTAAACAGTCGCTTCTCCCTGGTCTCTCCCCTTCTCCCGAAGTTACGGGGGTATTTTGCCGAGTTCC"
    "TTAACCATGATTCACTCGATCGCCTTGGTATTCTCTACCTAACCACCTGAGTCGGTTTGGTAGGATCACCCTGCTTCCCG"
    "CATTCGCGGTCACTATCAGGTCTCAGGATATGTGTGAGACGGATTTGCCTATCTCACTCCCTACACCCTTGGACGTGGAC"
    "TTGACTACTACCAAATCGGGTCACGCGCTCCGCTCAACATTCCATCACCCGAAGGTGACAGAAAAAAGAGTTTTAGGCGT"
    "TTAGCATCAAAAGGTTCATCTCGACTACGCCTGTCGGCCTCGCCTTAGGTCCCGACTTACCCAGGGCAGATTAGCTTGAC"
    "CCTGGAACCCTTGGTTATTCGGCGGACGGGTTTCTCGCCC";

  char *rev =
    "GGGCGAGAAACCCGTCCGCCGAATAACCAAGGGTTCCAGGGTCAAGCTAATCTGCCCTGGGTAAGTCGGGACCTAAGGCG"
    "AGGCCGACAGGCGTAGTCGAGATGAACCTTTTGATGCTAAACGCCTAAAACTCTTTTTTCTGTCACCTTCGGGTGATGGA"
    "ATGTTGAGCGGAGCGCGTGACCCGATTTGGTAGTAGTCAAGTCCACGTCCAAGGGTGTAGGGAGTGAGATAGGCAAATCC"
    "GTCTCACACATATCCTGAGACCTGATAGTGACCGCGAATGCGGGAAGCAGGGTGATCCTACCAAACCGACTCAGGTGGTT"
    "AGGTAGAGAATACCAAGGCGATCGAGTGAATCATGGTTAAGGAACTCGGCAAAATACCCCCGTAACTTCGGGAGAAGGGG"
    "AGAGACCAGGGAGAAGCGACTGTTTACTAAAAACACAGGTCCGTGCTAAGTCGCAAGACGATGTATACGGACTGACGCCT"
    "GCCCGGTGCTGGAAGGTTAACCCAGTAAACGGCGGTGGTAACTATAACCATCCTAAGGTAGCGAAATTCCTTGTCGGGTA"
    "AGTTCCGACCTGCACGAATGGCGTAACGACTTCTCCACTGTCGTTACGCGCAGCAGGACGGAAAGACCCCGGGACCTTTA"
    "CTATAGCTTGGTATTGGTGTTCGGTACGGCTTGTGTAGGATAGGTGGGAGACTGTGAAGCGTCGTTCTGGATATCTAACC"
    "TAGGTCCGTGATCCGGATCAGGGACATTGCCTGGTGGGTAGTTTAACTGGGGCGGTTGCCTCCTAAAAGGTAACGGAGGC"
    "TAAGTGCACAAGGGAGCTTGACTGCGAGACCGGCAGGTCGAGCAGGGACGAAAGTCGGAACTAGTGACCCGGCGCTGGCT"
    "TGTGGAAGCGGCGTCGCTCAGAGTCCATATCGACGGCATGGTTTGGCACCTCGATGTCGGCTCGTCGCATCCTGGGGGTG"
    "GAGTCTCTCCCAAGGGTTGGGCTGTTCGCCCATTAAAGCG";

  FastAAccessor F(fwd, 1000, false);
  FastAAccessor R(rev, 1000, true);
  u32bit        C;

  for (u32bit i=0; i<1000; i++)
    if (F[i] != R[i])
      exit(1);

  F.setPosition(0);
  R.setPosition(0);
  for (C=0; F.isValid() && R.isValid(); ++C, ++F, ++R)
    if (*F != *R)
      exit(2);
  if (C != 1000)
    exit(3);

  F.setPosition(999);
  R.setPosition(999);
  for (C=0; F.isValid() && R.isValid(); ++C, --F, --R)
    if (*F != *R)
      exit(4);
  if (C != 1000)
    exit(5);
}



//  Test pulling out a subsequence.  We're given coordinates in the
//  forward direction, but want to pull out the reverse complement
//  sequence.
//
void
easierTest(void) {
  char sub[1000];
  int  i;

  //  100A 200N 100T 400N 200G
  //
  for (i=0; i<1000; i++)
    sub[i] = 'N';
  for (i=0; i<100; i++)
    sub[i] = 'A';
  for (i=300; i<400; i++)
    sub[i] = 'T';
  for (i=600; i<700; i++)
    sub[i] = 'R';
  for (i=800; i<1000; i++)
    sub[i] = 'G';


  //  Pull out the reverse-complement sequence from 300-400
  //
  //  Asking for sequence from 300 to 400 should give up exactly to
  //  'A' (reverse-complelent of T) block.
  //
  //  Without setting the range, we'd get back the sequence at
  //  700-600, the location when globally reverse-complemented.
  //
  FastAAccessor S(sub, 1000, true);

  S.setRange(300, 100);
  for (i=300; i<400; i++)
    if (S[i] != 'A')
      fprintf(stderr, "FAILED: got %c at pos %d\n", S[i], i), exit(5);

  S.setRange(0, 0);
  for (i=300; i<400; i++)
    if (S[i] != 'Y')
      fprintf(stderr, "FAILED: got %c at pos %d\n", S[i], i), exit(6);
}



//  A harder test: build an accessor to access sequence from 100 to
//  300, and grow/shrink the region.
//  
void
harderTest(void) {
  char sub[1000];
  int  i;
  int  e;

  //  100A 200C 300G 400N
  //
  for (i=0; i<1000; i++)
    sub[i] = 'N';
  for (i=0; i<100; i++)
    sub[i] = 'A';
  for (i=100; i<300; i++)
    sub[i] = 'C';
  for (i=300; i<600; i++)
    sub[i] = 'G';

  //  Try forward.

  {
    fprintf(stderr, "Forward setRange/setPosition\n");
    FastAAccessor  A(sub, 1000, false);

    A.setRange(100, 200);
    A.setPosition(100);

    fprintf(stderr, "Range: "u32bitFMT"-"u32bitFMT" len="u32bitFMT"\n",
            A.getRangeBegin(), A.getRangeEnd(), A.getRangeLength());
    if ((A.getRangeBegin() != 100) || (A.getRangeEnd() != 300) || (A.getRangeLength() != 200))
      fprintf(stderr, "FAILED.\n"), exit(1);

    e = 0;
    for (int j=0; j<200; j++) {
      fprintf(stderr, "%c", *A);
      if (*A != 'C')
        e++;
      ++A;
    }

    fprintf(stderr, "\n");

    if (e)
      fprintf(stderr, "FAILED forward setRange/setPosition test: %d errors\n", e), exit(1);

    //  Decrease the size of our region, using the extend operators, then shift it
    //  to the right/left.
    //
    for (int i=0; i<190; i++)
      A.extendLeft(-1);
    for (int i=0; i<10; i++)
      A.extendRight(1);
    for (int i=0; i<10; i++)
      A.extendLeft(-1);

    fprintf(stderr, "Range: "u32bitFMT"-"u32bitFMT" len="u32bitFMT"\n",
            A.getRangeBegin(), A.getRangeEnd(), A.getRangeLength());
    if ((A.getRangeBegin() != 300) || (A.getRangeEnd() != 310) || (A.getRangeLength() != 10))
      fprintf(stderr, "FAILED.\n"), exit(1);

    e = 0;
    for (int j=0; j<20; j++) {
      fprintf(stderr, "%c", *A);
      if (*A != 'G')
        e++;
      ++A;
    }
    fprintf(stderr, "\n");

    if (e)
      fprintf(stderr, "FAILED reverse extendRange test: %d errors\n", e), exit(1);
  }

  {
    fprintf(stderr, "Reverse setRange/setPosition\n");
    FastAAccessor  A(sub, 1000, true);
    A.setRange(100, 200);
    A.setPosition(100);

    fprintf(stderr, "Range: "u32bitFMT"-"u32bitFMT" len="u32bitFMT"\n",
            A.getRangeBegin(), A.getRangeEnd(), A.getRangeLength());
    if ((A.getRangeBegin() != 100) || (A.getRangeEnd() != 300) || (A.getRangeLength() != 200))
      fprintf(stderr, "FAILED.\n"), exit(1);

    e = 0;
    for (int j=0; j<200; j++) {
      fprintf(stderr, "%c", *A);
      if (*A != 'G')
        e++;
      ++A;
    }
    fprintf(stderr, "\n");

    if (e)
      fprintf(stderr, "FAILED reverse setRange/setPosition test: %d errors\n", e), exit(1);

    //  Decrease the size of our region, using the extend operators, then shift it
    //  to the right/left.
    //
    for (int i=0; i<190; i++)
      A.extendLeft(-1);
    for (int i=0; i<10; i++)
      A.extendRight(1);
    for (int i=0; i<10; i++)
      A.extendLeft(-1);

    fprintf(stderr, "Range: "u32bitFMT"-"u32bitFMT" len="u32bitFMT"\n",
            A.getRangeBegin(), A.getRangeEnd(), A.getRangeLength());
    if ((A.getRangeBegin() != 90) || (A.getRangeEnd() != 100) || (A.getRangeLength() != 10))
      fprintf(stderr, "FAILED.\n"), exit(1);

    e = 0;
    for (int j=0; j<20; j++) {
      fprintf(stderr, "%c", *A);
      if (*A != 'T')
        e++;
      ++A;
    }
    fprintf(stderr, "\n");

    if (e)
      fprintf(stderr, "FAILED reverse extendRange test: %d errors\n", e), exit(1);
  }
}







int
main(int argc, char **argv) {
  simpleTest();
  easierTest();
  harderTest();

  fprintf(stderr, "All tests OK!\n");
  exit(0);
}
