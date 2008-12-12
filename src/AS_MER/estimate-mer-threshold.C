
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2006-2008, J. Craig Venter Institute
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

const char *mainid = "$Id: estimate-mer-threshold.C,v 1.1 2008-12-12 03:44:03 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_OVS_overlapStore.h"
}

#include "AS_MER_gkpStore_to_FastABase.H"

#include "bio++.H"
#include "sweatShop.H"
#include "positionDB.H"
#include "libmeryl.H"


int
main(int argc, char **argv) {
  char              *gkpPath = 0L;
  char              *merCountsFile = 0L;

  GateKeeperStore   *gkp = 0L;
  merylStreamReader *MF  = 0L;

  u32bit             maxCount = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpPath = argv[++arg];

    } else if (strcmp(argv[arg], "-m") == 0) {
      merCountsFile = argv[++arg];

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((gkpPath == 0L) || (merCountsFile == 0L) || (err)) {
    fprintf(stderr, "usage: %s -g gkpStore -m mercounts\n", argv[0]);
    fprintf(stderr, "  -g gkpStore     path to the gkpStore\n");
    fprintf(stderr, "  -m mercounts    file of mercounts\n");
    exit(1);
  }

  gkpStoreFile::registerFile();

  MF = new merylStreamReader(merCountsFile);

  //  Examine the counts, pick a reasonable upper limit.

  uint64  totalUsefulDistinct = MF->numberOfDistinctMers() - MF->numberOfUniqueMers();
  uint64  totalUsefulAll      = MF->numberOfTotalMers()    - MF->numberOfUniqueMers();
  uint64  distinct            = 0;
  uint64  total               = 0;
  uint32  Xcoverage           = 8;
  uint32  i                   = 0;

  fprintf(stderr, "distinct: "u64bitFMT"\n", MF->numberOfDistinctMers());
  fprintf(stderr, "unique:   "u64bitFMT"\n", MF->numberOfUniqueMers());
  fprintf(stderr, "total:    "u64bitFMT"\n", MF->numberOfTotalMers());

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
  fprintf(stderr, "Xcoverage zero %d %d %d\n", 1, 0, MF->histogram(1));

  for (i=2; (i < MF->histogramLength()) && (MF->histogram(i-1) > MF->histogram(i)); i++) {
    fprintf(stderr, "Xcoverage drop %d %d %d\n", i, MF->histogram(i-1), MF->histogram(i));
  }

  for (; (i < MF->histogramLength()) && (MF->histogram(i-1) < MF->histogram(i)); i++) {
    fprintf(stderr, "Xcoverage incr %d %d %d\n", i, MF->histogram(i-1), MF->histogram(i));
    Xcoverage = i;
  }

  fprintf(stderr, "Xcoverage done %d %d %d\n", i, MF->histogram(i-1), MF->histogram(i));

  fprintf(stderr, "Guessed X coverage is %d\n", Xcoverage);

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

  fprintf(stderr, "Set maxCount to "u32bitFMT", which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
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

    fprintf(stderr, "Reset maxCount to "u32bitFMT", which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
            maxCount, 100.0 * distinct / totalUsefulDistinct, 100.0 * total / totalUsefulAll);
  }

  fprintf(stdout, u32bitFMT"\n", maxCount);

  return(0);
}
