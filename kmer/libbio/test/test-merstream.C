#include <stdio.h>
#include <stdlib.h>

#include "bri++.H"

//  This is the input file.  Files are tested in test-merstream-speed.C
//
const char *sequence =
">(one)\n"
"TTTTTTTTTT\n"
"AAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT\n"
">(two)\n"
"TTTTTTTTAA\n"
">(thr)\n"
"TTTTTTTTAC\n"
">(fou)\n"
"T       T       T       T       T\n"
"T       T       T       A       G\n"
">(fiv)\n"
"T T T T T T T T A T C\n"
"\n"
">(six)\n"
"\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n"
"T\n";

//  This is the correct output
//
struct answers {
  u64bit  mer;
  u64bit  pos;
  u64bit  str;
};

//  These numbers seem to be correct, but I haven't rigorously
//  verified each one.
//
answers correct[22] = {
  { 0x00000000000fffff, 0,  0   },
  { 0x00000000000ffffc, 1,  1   },
  { 0x00000000000ffff0, 2,  2   },
  { 0x00000000000fffc0, 3,  3   },
  { 0x00000000000fff00, 4,  4   },
  { 0x00000000000ffc00, 5,  5   },
  { 0x00000000000ff000, 6,  6   },
  { 0x00000000000fc000, 7,  7   },
  { 0x00000000000f0000, 8,  8   },
  { 0x00000000000c0000, 9,  9   },
  { 0x0000000000000000, 10, 10  },
  { 0x0000000000000001, 11, 11  },
  { 0x0000000000000006, 12, 12  },
  { 0x00000000000fffff, 23, 23  },
  { 0x00000000000aaaaa, 34, 34  },
  { 0x00000000000fffff, 55, 55  },
  { 0x00000000000ffff0, 0,  66  },
  { 0x00000000000ffff1, 0,  77  },
  { 0x00000000000ffff2, 0,  88  },
  { 0x00000000000ffff3, 0,  99  },
  { 0x00000000000fffcd, 1,  100 },
  { 0x00000000000fffff, 0,  111 }
};


int
test1(merStream *S, char *id) {
  int c = 0;
  int e = 0;

  while (S->nextMer()) {
    if (c > 21)
      fprintf(stderr, "merStream(char): Too many mers in stream.\n"), e++;
    else if ((S->theFMer()               != correct[c].mer) ||
             (S->thePositionInSequence() != correct[c].pos) ||
             (S->thePositionInStream()   != correct[c].str))
      fprintf(stderr, "merStream(%s): "u64bitHEX"/"u64bitFMT"/"u64bitFMT" != correct: "u64bitHEX"/"u64bitFMT"/"u64bitFMT".\n",
              id,
              S->theFMer(), S->thePositionInSequence(), S->thePositionInStream(),
              correct[c].mer, correct[c].pos, correct[c].str), e++;
    c++;
  }

  return(e);
}


int
test2(merStream *MS1,
      merStream *MS2,
      merStream *MS3) {
  int e = 0;

  speedCounter C(" %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

  //  Run through all the mers, making sure they are all the same
  //
  while (MS1->nextMer() && MS2->nextMer() && MS3->nextMer()) {

#if 0
    //  If you're curious that things are actually non-zero...
    if (MS1->thePositionInSequence() == 666) {
      fprintf(stderr, "STAT: MS1: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS2: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS3: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
              MS1->theFMer(), MS1->theRMer(), MS1->thePositionInSequence(), MS1->thePositionInStream(), MS1->theSequenceNumber(),
              MS2->theFMer(), MS2->theRMer(), MS2->thePositionInSequence(), MS2->thePositionInStream(), MS2->theSequenceNumber(),
              MS3->theFMer(), MS3->theRMer(), MS3->thePositionInSequence(), MS3->thePositionInStream(), MS3->theSequenceNumber());
    }
#endif

    if ((MS1->theFMer()               != MS2->theFMer())               || (MS2->theFMer()               != MS3->theFMer()) ||
        (MS1->theRMer()               != MS2->theRMer())               || (MS2->theRMer()               != MS3->theRMer()) ||
        (MS1->thePositionInSequence() != MS2->thePositionInSequence()) || (MS2->thePositionInSequence() != MS3->thePositionInSequence()) ||
        (MS1->thePositionInStream()   != MS2->thePositionInStream())   || (MS2->thePositionInStream()   != MS3->thePositionInStream()) ||
        (MS1->theSequenceNumber()     != MS2->theSequenceNumber())     || (MS2->theSequenceNumber()     != MS3->theSequenceNumber())) {
      fprintf(stderr, "OOPS: MS1: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS2: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS3: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
              MS1->theFMer(), MS1->theRMer(), MS1->thePositionInSequence(), MS1->thePositionInStream(), MS1->theSequenceNumber(),
              MS2->theFMer(), MS2->theRMer(), MS2->thePositionInSequence(), MS2->thePositionInStream(), MS2->theSequenceNumber(),
              MS3->theFMer(), MS3->theRMer(), MS3->thePositionInSequence(), MS3->thePositionInStream(), MS3->theSequenceNumber());
    }
    C.tick();
  }
  C.finish();

  if (MS1->nextMer())
    fprintf(stderr, "merStream(msfr): Still has mers!\n"), e++;
  if (MS2->nextMer())
    fprintf(stderr, "merStream(file): Still has mers!\n"), e++;
  if (MS3->nextMer())
    fprintf(stderr, "merStream(char): Still has mers!\n"), e++;

  return(e);
}


int
main(int argc, char **argv) {
  merStream  *S = 0L;
  FILE       *F = 0L;
  u32bit      c = 0;
  u32bit      e = 0;


  //  This code does two tests.
  //
  //  The first, compares all three merStreams against a fixed and
  //  known input/output.  It's easy to debug.
  //
  //  The second compares all three merStreams against themselves,
  //  using whatever was passed in on the command line.  It's a
  //  heavier test, and probably not so easy to debug.  But I haven't
  //  had to try it.  :-)

  if (argc == 1) {

    ////////////////////////////////////////
    //
    //  Test the merstream using a character string as input
    //
    S = new merStream(10, sequence, strlen(sequence));
    e += test1(S, "char1");
    S->rewind();
    e += test1(S, "char2");
    delete S;

    if (e)
      return(e);

    ////////////////////////////////////////
    //
    //  Test the merstream using a fasta file as input
    //
    F = fopen("test-merstream.junk", "w");
    fwrite(sequence, sizeof(char), strlen(sequence) + 1, F);
    fclose(F);

    S = new merStream(10, "test-merstream.junk");
    c += test1(S, "file1");
    S->rewind();
    c += test1(S, "file2");
    delete S;

    if (e)
      return(e);

    ////////////////////////////////////////
    //
    //  Test the merstream using a merStreamFile as input
    //
    merStreamFileBuilder *B = new merStreamFileBuilder(10, "test-merstream.junk", "test-merstream.merstreamfile.junk");
    B->build();
    delete B;

    merStreamFileReader  *R = new merStreamFileReader("test-merstream.merstreamfile.junk");
    S = new merStream(R);
    c += test1(S, "msfr1");
    S->rewind();
    c += test1(S, "msfr1");
    delete S;
    delete R;

    unlink("test-merstream.junk");

    if (e)
      return(e);
  } else {

    //  We build three merStreams all from the same input.
    //    MS1 is built using a merStreamFile.
    //    MS2 is built using the FastA file directly.
    //    MS3 is built using the mmap()d file as a character string.
    //

    //  First, mmap() the file.  If it fails, we didn't waste time
    //  building the merStreamFile.
    //
    fprintf(stderr, "Mapping the input file to a char*\n");
    size_t  sequenceLen = 0;
    char   *sequence    = (char *)mapFile(argv[1], &sequenceLen, 'r');

    //  Then build the merStreamFile.
    //
    fprintf(stderr, "Building merStreamFile from file\n");
    merStreamFileBuilder  *MS1b = new merStreamFileBuilder(27, argv[1], "test-merstream.merstreamfile1.junk");
    MS1b->build();
    delete MS1b;

    fprintf(stderr, "Creating merStreamFileReader\n");
    merStreamFileReader   *MS1r = new merStreamFileReader("test-merstream.merstreamfile1.junk");


    //  Create some merStreams
    //
    fprintf(stderr, "Creating merStreams\n");
    merStream  *MS1 = new merStream(MS1r);
    merStream  *MS2 = new merStream(27, argv[1]);
    merStream  *MS3 = new merStream(27, sequence, sequenceLen);

    //  And another one, to test out merStreamFile buliding from a merStream
    //
#if 0
    fprintf(stderr, "Building merStreamFile from merStream of char*\n");
    merStreamFileBuilder  *MS4b = new merStreamFileBuilder(MS3, "test-merstream.merstreamfile2.junk");
    MS4b->build();
    delete MS4b;

    fprintf(stderr, "Creating merStreamFileReader\n");
    merStreamFileReader   *MS4r = new merStreamFileReader("test-merstream.merstreamfile2.junk");
    merStream             *MS4  = new merStream(MS4r);
#endif

    MS1->rewind();
    MS2->rewind();
    MS3->rewind();
#if 0
    MS4->rewind();
#endif

    e += test2(MS1, MS2, MS3);

    MS1->rewind();
    MS2->rewind();
    MS3->rewind();
#if 0
    MS4->rewind();
#endif

    e += test2(MS1, MS2, MS3);

    delete MS1;
    delete MS2;
    delete MS3;
#if 0
    delete MS4;
#endif
  }

  return(e);
}

