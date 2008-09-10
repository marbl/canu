#include "util.h"

#include "seqFile.H"
#include "sffFile.H"
#include "fastaFile.H"
#include "seqStore.H"
//#include "gkpStore.H"
#include "seqFactory.H"
#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

//  Include all header files directly; this is test code, right?

int
main(int argc, char **argv) {
  mt_s *mtctx = mtInit(451677);

  if (argc == 0) {
    fprintf(stderr, "usage: %s infile\n", argv[0]);
    exit(1);
  }

  seqFile  *SF = openSeqFile(argv[1]);


  fprintf(stderr, "source '%s' of type '%s' has "u32bitFMT" sequences.\n",
          SF->getSourceName(), SF->getFileTypeName(), SF->getNumberOfSequences());


  fprintf(stderr, "getSequenceLength() vs getSequence(full)\n");
  {
    char  *h = 0L;
    char  *s = 0L;
    u32bit hLen=0, hMax=0;
    u32bit sLen=0, sMax=0;

    for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
      SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

      if ((strlen(s) != SF->getSequenceLength(sid)) ||
          (strlen(s) != sLen) ||
          (SF->getSequenceLength(sid) != sLen)) {
        fprintf(stderr, "length differ for sid="u32bitFMT" h='%s' strlen(s)=%d sLen="u32bitFMT" getSequenceLength()="u32bitFMT"\n",
                sid, h, strlen(s), sLen, SF->getSequenceLength(sid));
      }
    }

    delete [] h;
    delete [] s;
  }


  fprintf(stderr, "getSequenceLength() vs getSequence(part)\n");
  {
    char  *p = new char [128 * 1024 * 1024];

    for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
      SF->getSequence(sid, 0, SF->getSequenceLength(sid), p);

      if (strlen(p) != SF->getSequenceLength(sid)) {
        fprintf(stderr, "length differ for sid="u32bitFMT" strlen(s)=%d getSequenceLength()="u32bitFMT"\n",
                sid, strlen(p), SF->getSequenceLength(sid));
      }
    }

    delete [] p;
  }



  fprintf(stderr, "getSequence(full) vs getSequence(part) & random access\n");
  {
    char  *h = 0L;
    char  *s = 0L;
    u32bit hLen=0, hMax=0;
    u32bit sLen=0, sMax=0;

    char **alls = new char * [SF->getNumberOfSequences()];

    for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
      SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

      alls[sid] = new char [sLen + 1];
      strcpy(alls[sid], s);
    }

    delete [] s;
    delete [] h;

    char  *p = new char [128 * 1024 * 1024];

    for (u32bit tl=0; tl<200 * SF->getNumberOfSequences(); tl++) {
      u32bit sid = mtRandom32(mtctx) % SF->getNumberOfSequences();

      u32bit b   = mtRandom32(mtctx) % SF->getSequenceLength(sid);
      u32bit e   = mtRandom32(mtctx) % SF->getSequenceLength(sid);
      if (b > e) {
        u32bit t=b;
        b = e;
        e = t;
      }

      SF->getSequence(sid, b, e, p);

      if (strlen(p) != e-b) {
        fprintf(stderr, "length differs; strlen(p)=%s requested="u32bitFMT"\n",
                strlen(p), e-b);
      }

      if (strncmp(alls[sid]+b, p, e-b)) {
        fprintf(stderr, "sequence differs.\n");
      }
    }

    delete [] p;

    for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++)
      delete [] alls[sid];

    delete [] alls;
  }




  //  Random get sequence - read all sequences using the stream, then
  //  read all backwards and randomly



#if 0

  fprintf(stderr, "seqInCore stream.\n");

  {
    seqInCore *SIC = SF->getSequenceInCore();
    while (SIC) {
      delete SIC;
      SIC = SF->getSequenceInCore();
    }
  }

  fprintf(stderr, "seqInCore iteration.\n");

  {
    for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++)
      seqInCore *SIC = SF->getSequenceInCore(sid);
      delete SIC;
    }
  }

#endif


  delete SF;

  return(0);
}

