
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

class alignStats {
public:
  alignStats() {
    startTime       = getTime();
    reportThreshold = 0;

    clear();
  };
  ~alignStats() {
  };

  void        clear(void) {
    nSkipped = 0;
    nPassed  = 0;
    nFailed  = 0;

    nFailExtA = 0;
    nFailExtB = 0;
    nFailExt  = 0;

    nExtendedA = 0;
    nExtendedB = 0;

    nPartial  = 0;
    nDovetail = 0;

    nExt5a = 0;
    nExt3a = 0;
    nExt5b = 0;
    nExt3b = 0;
  };


  alignStats &operator+=(alignStats &that) {
    nSkipped += that.nSkipped;
    nPassed  += that.nPassed;
    nFailed  += that.nFailed;

    nFailExtA += that.nFailExtA;
    nFailExtB += that.nFailExtB;
    nFailExt  += that.nFailExt;

    nExtendedA += that.nExtendedA;
    nExtendedB += that.nExtendedB;

    nPartial  += that.nPartial;
    nDovetail += that.nDovetail;

    nExt5a += that.nExt5a;
    nExt3a += that.nExt3a;
    nExt5b += that.nExt5b;
    nExt3b += that.nExt3b;

    return(*this);
  };

  void   reportStatus(void) {

    if (nPassed + nFailed < reportThreshold)
      return;

    reportThreshold += 10000;

    fprintf(stderr, "Tested %9" F_U64P " olaps -- Skipped %8.4f%% -- Passed %8.4f%% -- %8.2f olaps/sec\n",
            nPassed + nFailed,
            100.0 * nSkipped / (nPassed + nFailed),
            100.0 * nPassed  / (nPassed + nFailed),
            (nPassed + nFailed) / (getTime() - startTime));
  };

  void    reportFinal(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, " -- %" F_U64P " overlaps processed.\n", nPassed + nFailed);
    fprintf(stderr, " -- %" F_U64P " skipped.\n", nSkipped);
    fprintf(stderr, " -- %" F_U64P " failed %" F_U64P " passed (%.4f%%).\n", nFailed, nPassed, 100.0 * nPassed / (nPassed + nFailed));
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P " failed initial alignment, allowing A to extend\n", nFailExtA);
    fprintf(stderr, " -- %" F_U64P " failed initial alignment, allowing B to extend\n", nFailExtB);
    fprintf(stderr, " -- %" F_U64P " failed initial alignment\n", nFailExt);
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P " partial overlaps (of any quality)\n", nPartial);
    fprintf(stderr, " -- %" F_U64P " dovetail overlaps (before extensions, of any quality)\n", nDovetail);
    fprintf(stderr, " --\n");
    fprintf(stderr, " -- %" F_U64P "/%" F_U64P " A read dovetail extensions\n", nExt5a, nExt3a);
    fprintf(stderr, " -- %" F_U64P "/%" F_U64P " B read dovetail extensions\n", nExt5b, nExt3b);
  };

  double        startTime;
  uint64        reportThreshold;

  uint64        nSkipped;
  uint64        nPassed;
  uint64        nFailed;

  uint64        nFailExtA;
  uint64        nFailExtB;
  uint64        nFailExt;

  uint64        nExtendedA;
  uint64        nExtendedB;

  uint64        nPartial;
  uint64        nDovetail;

  uint64        nExt5a;
  uint64        nExt3a;
  uint64        nExt5b;
  uint64        nExt3b;
};

