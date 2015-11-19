
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    kmer/meryl/mapMers.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-OCT-18 to 2014-APR-11
 *      are Copyright 2006-2008,2013-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-JUL-22
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bio++.H"
#include "seqCache.H"
#include "merStream.H"
#include "libmeryl.H"
#include "existDB.H"

#define  OP_NONE    0
#define  OP_STATS   1
#define  OP_REGIONS 2
#define  OP_DETAILS 3

int
main(int argc, char **argv) {
  uint32    merSize    = 16;
  char     *merylFile  = 0L;
  char     *fastaFile  = 0L;
  bool      beVerbose  = false;
  uint32    loCount    = 0;
  uint32    hiCount    = ~uint32ZERO;
  uint32    operation  = OP_NONE;

  //  For OP_STATS

  uint32         Clen = 0;
  uint32         Cmax = 4 * 1024 * 1024;
  uint32        *C    = new uint32 [Cmax];

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merSize = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-mers") == 0) {
      merylFile = argv[++arg];

    } else if (strcmp(argv[arg], "-seq") == 0) {
      fastaFile = argv[++arg];

    } else if (strcmp(argv[arg], "-v") == 0) {
      beVerbose = true;

    } else if (strcmp(argv[arg], "-lo") == 0) {
      loCount = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-hi") == 0) {
      hiCount = strtouint32(argv[++arg]);

    } else if (strcmp(argv[arg], "-stats") == 0) {
      operation = OP_STATS;

    } else if (strcmp(argv[arg], "-regions") == 0) {
      operation = OP_REGIONS;

    } else if (strcmp(argv[arg], "-details") == 0) {
      operation = OP_DETAILS;

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
    }
    arg++;
  }

  if ((operation == OP_NONE) || (merylFile == 0L) || (fastaFile == 0L)) {
    fprintf(stderr, "usage: %s [-stats | -regions | -details] -m mersize -mers mers -seq fasta > output\n", argv[0]);
    exit(1);
  }

#if 0
  existDB       *E = NULL;

  if (fileExists("junk.existDB")) {
    fprintf(stderr, "loading from junk.existDB\n");
    E = new existDB("junk.existDB");
    fprintf(stderr, "loaded\n");
  } else {
    exit(1);
    E = new existDB(merylFile, merSize, existDBcounts, loCount, hiCount);
    E->saveState("junk.existDB");
  }
#endif

  existDB  *E = new existDB(merylFile, merSize, existDBcounts, loCount, hiCount);
  seqCache *F = new seqCache(fastaFile);

  fprintf(stderr, "Begin.\n");


  for (uint32 Sid=0; Sid < F->getNumberOfSequences(); Sid++) {
    seqInCore  *S  = F->getSequenceInCore(Sid);
    merStream  *MS = new merStream(new kMerBuilder(merSize),
                                   new seqStream(S->sequence(), S->sequenceLength()),
                                   true, true);

    //  with counts, report mean, mode, median, min, max for each frag.
    if (operation == OP_STATS) {

      Clen = 0;

      while (MS->nextMer())
        C[Clen++] = E->count(MS->theFMer()) + E->count(MS->theRMer());

      uint64         mean     = uint64ZERO;
      uint64         min      = ~uint64ZERO;
      uint64         max      = uint64ZERO;
      uint64         hist[16]  = { 0 };

      //  Histogram values are powers of two, e.g., <=1, <=2, <=4, <=8, <=16, <=32, <=64, <=128, <=256, <=512, <=1024, <=4096, <=8192, <=328768

      for (uint32 i=0; i<Clen; i++) {
        mean += C[i];

        if (min > C[i])
          min = C[i];
        if (max < C[i])
          max = C[i];

        hist[ logBaseTwo64(C[i]) ]++;
      }

      if (Clen > 0) {
        mean /= Clen;

      } else {
        mean = uint64ZERO;
        min  = uint64ZERO;
        max  = uint64ZERO;
      }

      fprintf(stdout,
              "%s\t"
              uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"
              uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"
              uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\n",
              S->header(),
              mean, min, max,
              hist[ 0], hist[ 1], hist[ 2], hist[ 3], hist[ 4], hist[ 5], hist[ 6], hist[ 7],
              hist[ 8], hist[ 9], hist[10], hist[11], hist[12], hist[13], hist[14], hist[15]);
    }


    //  without counts, reports regions with mer coverage.
    //  Orientation tells us nothing, since the mers are probably canonical
    if (operation == OP_REGIONS) {
      uint64    beg = ~uint64ZERO;
      uint64    end = ~uint64ZERO;
      uint64    pos = ~uint64ZERO;

      uint64    numCovReg = 0;
      uint64    lenCovReg = 0;

      while (MS->nextMer()) {
        if (E->exists(MS->theFMer()) || E->exists(MS->theRMer())) {
          pos = MS->thePositionInSequence();

          if (beg == ~uint64ZERO)
            beg = end = pos;

          if (pos <= end + merSize) {
            end = pos;
          } else {
            fprintf(stdout, "%s\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\n", S->header(), beg, end+merSize, end+merSize - beg);
            numCovReg++;
            lenCovReg += end+merSize - beg;
            beg = end = pos;
          }
        } else {
          fprintf(stdout, "%s\t"uint64FMT"\tuncovered\n", S->header(), MS->thePositionInSequence());
        }
      }

      if (beg != ~uint64ZERO)
        fprintf(stdout, "%s\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\n", S->header(), beg, end+merSize, end+merSize - beg);

      fprintf(stderr, "numCovReg: "uint64FMT"\n", numCovReg);
      fprintf(stderr, "lenCovReg: "uint64FMT"\n", lenCovReg);
    }



    if (operation == OP_DETAILS) {
      char  merString[256];

      while (MS->nextMer()) {
        uint64 beg = MS->thePositionInSequence();
        uint64 end = beg + merSize;
        uint64 fnt = E->count(MS->theFMer());
        uint64 rnt = E->count(MS->theRMer());

        fprintf(stdout, "%s\t%s\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\t"uint64FMT"\n",
                S->header(),
                MS->theFMer().merToString(merString),
                beg,
                end,
                fnt,
                rnt,
                fnt + rnt);
      }
    }


    delete MS;
    delete S;
  }


  delete F;
  delete E;
}
