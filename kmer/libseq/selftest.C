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

  u32bit nums = 0;
  char **allh = 0L;
  char **alls = 0L;

  if (argc == 0) {
    fprintf(stderr, "usage: %s infile\n", argv[0]);
    exit(1);
  }

  {
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


    fprintf(stderr, "loading allh/alls.\n");
    {
      char  *h = 0L;
      char  *s = 0L;
      u32bit hLen=0, hMax=0;
      u32bit sLen=0, sMax=0;

      nums = SF->getNumberOfSequences();
      allh = new char * [SF->getNumberOfSequences()];
      alls = new char * [SF->getNumberOfSequences()];

      for (u32bit sid=0; sid<SF->getNumberOfSequences(); sid++) {
        SF->getSequence(sid, h, hLen, hMax, s, sLen, sMax);

        allh[sid] = new char [hLen + 1];
        strcpy(allh[sid], h);

        alls[sid] = new char [sLen + 1];
        strcpy(alls[sid], s);
      }

      delete [] s;
      delete [] h;
    }



    fprintf(stderr, "getSequence(full) vs getSequence(part) & random access\n");
    {
      char  *p = new char [128 * 1024 * 1024];

      for (u32bit tl=0; tl<2 * SF->getNumberOfSequences(); tl++) {
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
    }


    delete SF;
  }


  ////////////////////////////////////////


  fprintf(stderr, "seqCache.\n");
  {
    seqCache  *SC = new seqCache(argv[1]);
    delete SC;
  }


  fprintf(stderr, "seqInCore stream.\n");
  {
    seqCache  *SC = new seqCache(argv[1], 0, true);
    seqInCore *IC = SC->getSequenceInCore();

    while (IC) {
      delete IC;
      IC = SC->getSequenceInCore();
    }
  }


  fprintf(stderr, "seqInCore iteration.\n");
  {
    seqCache  *SC = new seqCache(argv[1], 0, true);

    for (u32bit sid=0; sid<SC->getNumberOfSequences(); sid++) {
      seqInCore *IC = SC->getSequenceInCore(sid);
      delete IC;
    }
    delete SC;
  }


  fprintf(stderr, "seqInCore cache random access.\n");
  {
    seqCache  *SC = new seqCache(argv[1], nums/10, true);

    for (u32bit sid=0; sid<SC->getNumberOfSequences(); sid++) {
      u32bit    sid = mtRandom32(mtctx) % SC->getNumberOfSequences();

      //  Cache owns these.
      seqInCore *IC = SC->getSequenceInCore(sid);
    }
    delete SC;
  }


  ////////////////////////////////////////

  fprintf(stderr, "seqStream\n");
  {
    seqStream *S = new seqStream(argv[1]);
    delete S;
  }

  ////////////////////////////////////////

  fprintf(stderr, "merStream\n");
  {
    merStream *M = new merStream(new kMerBuilder(22),
                                 new seqStream(argv[1]),
                                 true, true);
    delete M;
  }

  ////////////////////////////////////////


  for (u32bit sid=0; sid<nums; sid++) {
    delete [] allh[sid];
    delete [] alls[sid];
  }

  delete [] allh;
  delete [] alls;

  return(0);
}

