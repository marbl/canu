
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
 *  Modifications by:
 *
 *    Brian P. Walenz from 2006-MAY-09 to 2014-APR-11
 *      are Copyright 2006-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-OCT-07
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "util++.H"
#include "bio++.H"
#include "existDB.H"

#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

int
main(int argc, char **argv) {
  char  *merFile   = 0L;
  char  *queryFile = 0L;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merFile = argv[++arg];
    } else if (strcmp(argv[arg], "-q") == 0) {
      queryFile = argv[++arg];
    } else {
      fprintf(stderr, "Unknown arg '%s'\n", argv[arg]);
    }
    arg++;
  }

  existDB      *E = new existDB(merFile, 22, existDBnoFlags, 0, ~uint32ZERO);
  seqCache     *Q = new seqCache(queryFile);
  seqInCore    *S = Q->getSequenceInCore();

  intervalList<uint64>  IL;
  speedCounter          SC(" %8f frags (%8.5f frags/sec)\r", 1, 1000, true);

  while (S) {
    merStream     *MS = new merStream(new kMerBuilder(22),
                                      new seqStream(S->sequence(), S->sequenceLength()),
                                      true, true);

    IL.clear();

    while (MS->nextMer()) {
      if (E->exists(MS->theFMer())) {
        IL.add(MS->thePositionInSequence(), 22);
      }
    }

    IL.merge();

    if (IL.sumOfLengths() > 0) {
      fprintf(stdout, "%5.2f\n",
              100.0 * IL.sumOfLengths() / (double)S->sequenceLength());
    }

    delete MS;
    delete S;

    SC.tick();

    S = Q->getSequenceInCore();
  }

  delete Q;
  delete E;

  return(0);
}

