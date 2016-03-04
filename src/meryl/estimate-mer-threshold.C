
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
 *    src/AS_MER/estimate-mer-threshold.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2008-DEC-11 to 2014-APR-11
 *      are Copyright 2008-2010,2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "libmeryl.H"


int
main(int argc, char **argv) {
  char              *gkpPath = 0L;
  char              *merCountsFile = 0L;

  merylStreamReader *MF  = 0L;

  uint32             maxCount = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merCountsFile = argv[++arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((merCountsFile == 0L) || (err)) {
    fprintf(stderr, "usage: %s -m mercounts\n", argv[0]);
    fprintf(stderr, "  -m mercounts    file of mercounts\n");
    exit(1);
  }

  MF = new merylStreamReader(merCountsFile);

  //  Examine the counts, pick a reasonable upper limit.

  uint64  totalUsefulDistinct = MF->numberOfDistinctMers() - MF->numberOfUniqueMers();
  uint64  totalUsefulAll      = MF->numberOfTotalMers()    - MF->numberOfUniqueMers();
  uint64  distinct            = 0;
  uint64  total               = 0;
  uint32  Xcoverage           = 8;

  fprintf(stderr, "distinct: "F_U64"\n", MF->numberOfDistinctMers());
  fprintf(stderr, "unique:   "F_U64"\n", MF->numberOfUniqueMers());
  fprintf(stderr, "total:    "F_U64"\n", MF->numberOfTotalMers());

  //  Pass 0: try to deduce the X coverage we have.  The
  //  pattern we should see in mer counts is an initial spike
  //  for unique mers (these contain errors), then a drop into
  //  a valley, and a bump at the X coverage.
  //
  //  .
  //  .      ...
  //  ..  ..........
  //  .................
  //
  //  If this pattern is not found, we fallback to the default
  //  guess of 8x coverage.
  //

  uint32  i  = 0;
  uint32  iX = 0;

  fprintf(stderr, "distinct: "F_U64"\n", MF->numberOfDistinctMers());
  fprintf(stderr, "unique:   "F_U64"\n", MF->numberOfUniqueMers());
  fprintf(stderr, "total:    "F_U64"\n", MF->numberOfTotalMers());

  fprintf(stderr, "Xcoverage zero 1 0 "F_U64"\n", MF->histogram(1));

  for (i=2; (i < MF->histogramLength()) && (MF->histogram(i-1) > MF->histogram(i)); i++)
    fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(i-1), MF->histogram(i));

  iX = i - 1;

  for (; i < MF->histogramLength(); i++) {
    if (MF->histogram(iX) < MF->histogram(i)) {
      fprintf(stderr, "Xcoverage incr "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
      iX = i;
    } else {
      //fprintf(stderr, "Xcoverage drop "F_U32" "F_U64" "F_U64"\n", i, MF->histogram(iX), MF->histogram(i));
    }
  }

  fprintf(stderr, "Guessed X coverage is "F_U32"\n", iX);

  Xcoverage = iX;

  //  Pass 1: look for a reasonable limit, using %distinct and %total.
  //
  for (i=2; (i < MF->histogramLength()) && (maxCount == 0); i++) {
    distinct += MF->histogram(i);
    total    += MF->histogram(i) * i;

    //  If we cover 99% of all the distinct mers, that's reasonable.
    //
    if ((distinct / (double)totalUsefulDistinct) > 0.99)
      maxCount = i;

    //  If we're a somewhat high count, and we're covering 2/3
    //  of the total mers, assume that there are lots of
    //  errors (or polymorphism) that are preventing us from
    //  covering many distinct mers.
    //
    if ((i > 25 * Xcoverage) && ((total / (double)totalUsefulAll) > (2.0 / 3.0)))
      maxCount = i;
  }

  fprintf(stderr, "Set maxCount to "F_U32", which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
          i, 100.0 * distinct / totalUsefulDistinct, 100.0 * total / totalUsefulAll);


  //  Pass 2: if the limit is relatively small compared to our
  //  guessed Xcoverage, and %total is high, keep going to
  //  close 75% of the gap in total coverage.  So if the TC is
  //  90%, we'd keep going until TC is 97.5%.
  //
  //  If we're WAY low compared to X coverage, close the gap
  //  too, but not as much.  This only happens if we're
  //  covering 99% of the distinct, so we're already in good
  //  shape.  The genome doesn't appear to be very repetitive.
  //
  if (((maxCount <  5 * Xcoverage)) ||
      ((maxCount < 50 * Xcoverage) && (total / (double)totalUsefulAll > 0.90))) {
    double  closeAmount = 0.75;

    if (total / (double)totalUsefulAll <= 0.90)
      closeAmount = 0.5;

    //  No, really.  This is just 0.75 * (1-TC) + TC
    double  desiredTC = closeAmount + (1 - closeAmount) * total / (double)totalUsefulAll;

    for (; (i < MF->histogramLength()) && (total / (double)totalUsefulAll < desiredTC); i++) {
      distinct += MF->histogram(i);
      total    += MF->histogram(i) * i;
    }

    maxCount = i;

    fprintf(stderr, "Reset maxCount to "F_U32", which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
            maxCount, 100.0 * distinct / totalUsefulDistinct, 100.0 * total / totalUsefulAll);
  }

  fprintf(stdout, F_U32"\n", maxCount);

  return(0);
}
