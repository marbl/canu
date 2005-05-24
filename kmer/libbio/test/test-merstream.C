#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"

//  XXX I have malloc errors.

//  XXX I should have (u64bit) cases removed before FMer() -- output
//  the mer as bases.

//  This is the input file.  Deflines are exactly 10 characters long
//  because the chained sequence uses 10 characters as a separator.
//  This allows us to use the same test harness for everything.
//
//  BUT, total space between sequences is 11 characters -- the newline
//  at the end of the sequence + defline + newline at the end of the
//  defline.
//
//  sequenceU is unsqueezed -- use it with the MSFR and the chainedSeq
//  sequenceS is squeezed   -- use it with the FastAstream and char*
//
const char *sequenceU =
">(one...)\n"
"TTTTTTTTTT\n"
"AAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT\n"
">(two...)\n"
"TTTTTTTTAA\n"
">(thr...)\n"
"TTTTTTTTAC\n"
">(fou...)\n"
"T       T       T       T       T\n"
"T       T       T       A       G\n"
">(fiv...)\n"
"T T T T T T T T A T C\n"
"\n"
">(six...)\n"
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

//1234567890123456789012345678901234567890123456789012345678901234567890
const char *sequenceS =
">(one...)\n"   // 10 long, including newline
"TTTTTTTTTTAAAAAAAAAACGNTTTTTTTTTTNGGGGGGGGGGNAAAAAAAAANTTTTTTTTTT\n"  //  76 long
">(two...)\n"
"TTTTTTTTAA\n"
">(thr...)\n"
"TTTTTTTTAC\n"
">(fou...)\n"
"TTTTTTTTAG\n"
">(fiv...)\n"
"TTTTTTTTATC\n"
">(six...)\n"
"TTTTTTTTTT\n";

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
  { 0x00000000000fffff, 0,  10  },  //  Yes, this entry is correct.
  { 0x00000000000ffffc, 1,  11  },
  { 0x00000000000ffff0, 2,  12  },
  { 0x00000000000fffc0, 3,  13  },
  { 0x00000000000fff00, 4,  14  },
  { 0x00000000000ffc00, 5,  15  },
  { 0x00000000000ff000, 6,  16  },
  { 0x00000000000fc000, 7,  17  },
  { 0x00000000000f0000, 8,  18  },
  { 0x00000000000c0000, 9,  19  },
  { 0x0000000000000000, 10, 20  },
  { 0x0000000000000001, 11, 21  },
  { 0x0000000000000006, 12, 22  },
  { 0x00000000000fffff, 23, 33  },
  { 0x00000000000aaaaa, 34, 44  },
  { 0x00000000000fffff, 55, 65  },
  { 0x00000000000ffff0, 0,  86  },
  { 0x00000000000ffff1, 0,  107 },
  { 0x00000000000ffff2, 0,  128 },
  { 0x00000000000ffff3, 0,  149 },
  { 0x00000000000fffcd, 1,  150 },
  { 0x00000000000fffff, 0,  171 }
};


//  The positionInStream() that a correct FastAstream should return.
//  Initialized in main()
//



int
test1(merStream *S, char *id, int offset) {
  int c = 0;
  int e = 0;

  while (S->nextMer()) {
    if (c > 21)
      fprintf(stderr, "merStream(char): Too many mers in stream.\n"), e++;
    else if (((u64bit)S->theFMer()       != correct[c].mer) ||
             (S->thePositionInSequence() != correct[c].pos) ||
             (S->thePositionInStream()   != correct[c].str - offset))
      fprintf(stderr, "merStream(%s): "u64bitHEX"/"u64bitFMT"/"u64bitFMT" != correct: "u64bitHEX"/"u64bitFMT"/"u64bitFMT".\n",
              id,
              (u64bit)S->theFMer(), S->thePositionInSequence(), S->thePositionInStream(),
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
              (u64bit)MS1->theFMer(), (u64bit)MS1->theRMer(), MS1->thePositionInSequence(), MS1->thePositionInStream(), MS1->theSequenceNumber(),
              (u64bit)MS2->theFMer(), (u64bit)MS2->theRMer(), MS2->thePositionInSequence(), MS2->thePositionInStream(), MS2->theSequenceNumber(),
              (u64bit)MS3->theFMer(), (u64bit)MS3->theRMer(), MS3->thePositionInSequence(), MS3->thePositionInStream(), MS3->theSequenceNumber());
    }
#endif

    if ((MS1->theFMer()               != MS2->theFMer())               || (MS2->theFMer()               != MS3->theFMer()) ||
        (MS1->theRMer()               != MS2->theRMer())               || (MS2->theRMer()               != MS3->theRMer()) ||
        (MS1->thePositionInSequence() != MS2->thePositionInSequence()) || (MS2->thePositionInSequence() != MS3->thePositionInSequence()) ||
        (MS1->thePositionInStream()   != MS2->thePositionInStream())   || (MS2->thePositionInStream()   != MS3->thePositionInStream()) ||
        (MS1->theSequenceNumber()     != MS2->theSequenceNumber())     || (MS2->theSequenceNumber()     != MS3->theSequenceNumber())) {
      fprintf(stderr, "OOPS: MS1: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS2: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MS3: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
              (u64bit)MS1->theFMer(), (u64bit)MS1->theRMer(), MS1->thePositionInSequence(), MS1->thePositionInStream(), MS1->theSequenceNumber(),
              (u64bit)MS2->theFMer(), (u64bit)MS2->theRMer(), MS2->thePositionInSequence(), MS2->thePositionInStream(), MS2->theSequenceNumber(),
              (u64bit)MS3->theFMer(), (u64bit)MS3->theRMer(), MS3->thePositionInSequence(), MS3->thePositionInStream(), MS3->theSequenceNumber());
      fprintf(stderr, "MS1:\n"); MS1->theFMer().dump(stderr);
      fprintf(stderr, "MS2:\n"); MS2->theFMer().dump(stderr);
      fprintf(stderr, "MS3:\n"); MS3->theFMer().dump(stderr);
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
  u32bit       e = 0;



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

    {
      FILE *F = fopen("junkU", "w");
      fwrite(sequenceU, sizeof(char), strlen(sequenceU), F);
      fclose(F);
    }

    {
      FILE *F = fopen("junkS", "w");
      fwrite(sequenceS, sizeof(char), strlen(sequenceS), F);
      fclose(F);
    }


    ////////////////////////////////////////
    //
    //  Test that the FastAstream works
    //
    {
      //  Initialize correct answers -- FastAstream only returns ACGT.
      //
      u32bit fastastreamcorrect[200];

      for (u32bit i=0, j=0; i<strlen(sequenceS); i++) {
        if ((sequenceS[i] == 'A') ||
            (sequenceS[i] == 'C') ||
            (sequenceS[i] == 'G') ||
            (sequenceS[i] == 'T')) {
          fastastreamcorrect[j++] = i;
        }
      }

      FastAstream *F = new FastAstream("junkS");
      unsigned char  r;

      e = 0;

      for (u32bit j=0; (r = F->nextSymbol()) != 0; ) {

        //  We need to skip the FastAstream new sequence and gap
        //  status codes (253 and 254).

        if (r < 128) {
          if (fastastreamcorrect[j] != F->thePositionInStream()) {
            fprintf(stderr, "FastAstream positions wrong: got "u64bitFMT" expected "u32bitFMT"\n",
                    F->thePositionInStream(),
                    fastastreamcorrect[j]);
            e++;
          }
          j++;
        }
      }
      delete F;

      if (e)
        return(e);
    }


    ////////////////////////////////////////
    //
    //  Test that the chained sequence
    //
    {
      //  Initialize correct answers -- chainedSequence returns all letters, but we only test NACGT.
      //
      u32bit fastastreamcorrect[200];

      for (u32bit i=0, j=0; i<strlen(sequenceS); i++) {
        if ((sequenceS[i] == 'N') ||
            (sequenceS[i] == 'A') ||
            (sequenceS[i] == 'C') ||
            (sequenceS[i] == 'G') ||
            (sequenceS[i] == 'T')) {
          fastastreamcorrect[j++] = i;
        }
      }

      chainedSequence *S = new chainedSequence();
      S->setSource("junkU");
      S->setSeparator('.');
      S->setSeparatorLength(11);
      S->finish();
      unsigned char  r;

      e = 0;

      for (u32bit j=0; (r = S->get()) != 0; ) {

        //  We need to skip the separator

        //fprintf(stderr, u64bitFMTW(3)" - "u64bitFMTW(3)" - %c\n", S->thePositionInStream(), S->thePositionInSequence(), r);

        if (r != '.') {
          if (fastastreamcorrect[j] - 10 != S->thePositionInStream()) {
            fprintf(stderr, "chainedSequence positions wrong:   got "u64bitFMT" expected "u32bitFMT"\n",
                    S->thePositionInStream(),
                    fastastreamcorrect[j] - 10);
            e++;
#if 0
          } else {
            fprintf(stderr, "chainedSequence positions correct: got "u64bitFMT" expected "u32bitFMT"\n",
                    S->thePositionInStream(),
                    fastastreamcorrect[j] - 10);
#endif
          }
          j++;
        }
      }
      delete S;

      if (e)
        return(e);
    }










    ////////////////////////////////////////
    //
    //  Test the merstream using a character string as input
    //
    {
      merStream *S = new merStream(10, sequenceS, strlen(sequenceS));
      e += test1(S, "char1s", 0);
      S->rewind();
      e += test1(S, "char2s", 0);
      delete S;

      if (e)
        return(e);
    }


    ////////////////////////////////////////
    //
    //  Test the merstream using a FastAstream as input
    //
    {
      FastAstream *F = new FastAstream("junkS");
      merStream   *S = new merStream(10, F);
      e += test1(S, "file1s", 0);
      S->rewind();
      e += test1(S, "file2s", 0);
      delete S;

      if (e)
        return(e);
    }


    ////////////////////////////////////////
    //
    //  Test the merstream using a chainedSequence as input
    //
    {
      chainedSequence *C = new chainedSequence;
      C->setSource("junkS");
      C->setSeparator('.');
      C->setSeparatorLength(11);
      C->finish();
      merStream       *S = new merStream(10, C);
      e += test1(S, "chain1s", 10);
      S->rewind();
      e += test1(S, "chain2s", 10);
      delete S;

      if (e)
        return(e);
    }

    {
      chainedSequence *C = new chainedSequence;
      C->setSource("junkU");
      C->setSeparator('.');
      C->setSeparatorLength(11);
      C->finish();
      merStream       *S = new merStream(10, C);
      e += test1(S, "chain1u", 10);
      S->rewind();
      e += test1(S, "chain2u", 10);
      delete S;

      if (e)
        return(e);
    }



    ////////////////////////////////////////
    //
    //  Test the merstream using a merStreamFile as input
    //
    {
      merStreamFileBuilder *B = new merStreamFileBuilder(10, "junkS", "junk.msf");
      B->build();
      delete B;

      merStreamFileReader  *R = new merStreamFileReader("junk.msf");
      merStream            *S = new merStream(R);
      e += test1(S, "msfr1s", 0);
      S->rewind();
      e += test1(S, "msfr1s", 0);
      delete S;
      delete R;

      if (e)
        return(e);
    }


    unlink("junkS");
    unlink("junkU");
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
    merStream   *MS1 = new merStream(MS1r);

    FastAstream *FS2 = new FastAstream(argv[1]);
    merStream   *MS2 = new merStream(27, FS2);

    merStream   *MS3 = new merStream(27, sequence, sequenceLen);

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
    delete FS2;
    delete MS3;
#if 0
    delete MS4;
#endif
  }

  return(e);
}

