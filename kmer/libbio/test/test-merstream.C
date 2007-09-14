#include <stdio.h>
#include <stdlib.h>

#include "bio++.H"



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
test2(merStream *MSa,
      merStream *MSb) {
  static char  stra[256], strb[256], strc[256], strd[256];
  int e = 0;

  speedCounter C(" %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, true);

  //  Run through all the mers, making sure they are all the same
  //
  while (MSa->nextMer() && MSb->nextMer()) {

#if 1
    //  If you're curious that things are actually non-zero...
    fprintf(stderr, "STAT: MSa: %s/%s @ "u64bitFMT"/"u64bitFMT","u64bitFMT"  MSb: %s/%s @ "u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
            MSa->theFMer().merToString(stra), MSa->theRMer().merToString(strb), MSa->thePositionInSequence(), MSa->thePositionInStream(), MSa->theSequenceNumber(),
            MSb->theFMer().merToString(strc), MSb->theRMer().merToString(strd), MSb->thePositionInSequence(), MSb->thePositionInStream(), MSb->theSequenceNumber());
    //fprintf(stderr, "STAT: MSa: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"  MSb: "u64bitHEX","u64bitHEX"@"u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
    //        (u64bit)MSa->theFMer(), (u64bit)MSa->theRMer(), MSa->thePositionInSequence(), MSa->thePositionInStream(), MSa->theSequenceNumber(),
    //        (u64bit)MSb->theFMer(), (u64bit)MSb->theRMer(), MSb->thePositionInSequence(), MSb->thePositionInStream(), MSb->theSequenceNumber());
#endif

    if ((MSa->theFMer()               != MSb->theFMer())               ||
        (MSa->theRMer()               != MSb->theRMer())               ||
        (MSa->thePositionInSequence() != MSb->thePositionInSequence()) ||
        (MSa->thePositionInStream()   != MSb->thePositionInStream())   ||
        (MSa->theSequenceNumber()     != MSb->theSequenceNumber())) {
      fprintf(stderr, "OOPS: MSa: %s/%s @ "u64bitFMT"/"u64bitFMT","u64bitFMT"  MSb: %s/%s @ "u64bitFMT"/"u64bitFMT","u64bitFMT"\n",
              MSa->theFMer().merToString(stra), MSa->theRMer().merToString(strb), MSa->thePositionInSequence(), MSa->thePositionInStream(), MSa->theSequenceNumber(),
              MSb->theFMer().merToString(strc), MSb->theRMer().merToString(strd), MSb->thePositionInSequence(), MSb->thePositionInStream(), MSb->theSequenceNumber());
      //fprintf(stderr, "MSa:\n"); MSa->theFMer().dump(stderr);
      //fprintf(stderr, "MSb:\n"); MSb->theFMer().dump(stderr);
    }
    C.tick();
  }
  C.finish();

  if (MSa->nextMer())
    fprintf(stderr, "merStream MSa still has mers!\n"), e++;
  if (MSb->nextMer())
    fprintf(stderr, "merStream MSb still has mers!\n"), e++;

  fprintf(stderr, "test2() finished.\n");

  return(e);
}


int
main(int argc, char **argv) {
  u32bit       e = 0;

  if (argc != 2) {
    fprintf(stderr, "usage: %s some.fasta\n", argv[0]);
    fprintf(stderr, "  if some.fasta has one sequence, the merStreamFile, chainedSeq and seqInCore ifaces are tested\n");
    fprintf(stderr, "  otherwise, only merStreamFile and chainedSeq are tested.\n");
    exit(1);
  }

  seqStore   *SS = new seqStore(argv[1], new seqStream(argv[1], true));
  seqStream  *CS = new seqStream(argv[1], true);

  kMerBuilder  KB1(27);
  merStream   *MS1 = new merStream(&KB1, SS);

  kMerBuilder  KB2(27);
  merStream   *MS2 = new merStream(&KB2, CS);

  fprintf(stderr, "test2()-- pass 1\n");
  e += test2(MS1, MS2);

  MS1->rewind();
  MS2->rewind();

  fprintf(stderr, "test2()-- pass 2\n");
  e += test2(MS1, MS2);

  seqFile  *SF = openSeqFile(argv[1]);
  SF->openIndex();

  if (SF->getNumberOfSequences() == 1) {
    seqInCore  *SC  = SF->getSequenceInCore();

    kMerBuilder  KB3(27);
    merStream  *MS3 = new merStream(&KB3, SC);

    MS1->rewind();
    MS2->rewind();
    MS3->rewind();

    fprintf(stderr, "test2()-- pass 3\n");
    e += test2(MS1, MS3);

    MS1->rewind();
    MS2->rewind();
    MS3->rewind();

    fprintf(stderr, "test2()-- pass 4\n");
    e += test2(MS2, MS3);

    delete MS3;
    delete SC;
  }

  delete MS2;
  delete MS1;
  delete SF;
  delete CS;
  delete SS;

  return(e);
}

