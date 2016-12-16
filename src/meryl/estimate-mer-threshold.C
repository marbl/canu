
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

#include "splitToWords.H"

//  Try to deduce the X coverage we have.  The pattern we should see in mer counts is an initial
//  spike for unique mers (these contain errors), then a drop into a valley, and a bump at the X
//  coverage.
//
//  .
//  .      ...
//  ..  ..........
//  .................
//
double
guessCoverage(uint32 *hist, uint32 histLen) {
  uint32  i = 2;

  while ((i < histLen) && (hist[i-1] > hist[i]))
    i++;

  uint32  iX = i - 1;

  while (i < histLen) {
    if (hist[iX] < hist[i])
      iX = i;
    i++;
  }

  fprintf(stderr, "Guessed X coverage is " F_U32 "\n", iX);

  return(iX);
}



void
loadHistogram(merylStreamReader *MF,
              uint64 &nDistinct,
              uint64 &nUnique,
              uint64 &nTotal,
              uint32 &histLen, uint32* &hist) {

  nDistinct = MF->numberOfDistinctMers();
  nUnique   = MF->numberOfUniqueMers();
  nTotal    = MF->numberOfTotalMers();

  histLen   = MF->histogramLength();
  hist      = new uint32 [histLen];

  for (uint32 hh=0; hh<histLen; hh++)
    hist[hh] = MF->histogram(hh);
}



void
loadHistogram(FILE *HF,
              uint64 &nDistinct,
              uint64 &nUnique,
              uint64 &nTotal,
              uint32 &histLen, uint32* &hist) {
  char    L[1024];
  uint32  histMax;

  nDistinct = 0;
  nUnique   = 0;
  nTotal    = 0;

  histLen   = 0;
  histMax   = 1048576;
  hist      = new uint32 [histMax];

  memset(hist, 0, sizeof(uint32) * histMax);

  fgets(L, 1024, HF);

  while (!feof(HF)) {
    splitToWords  W(L);

    uint32  h = W(0);
    uint32  c = W(1);

    while (h >= histMax)
      resizeArray(hist, histLen, histMax, histMax * 2, resizeArray_copyData | resizeArray_clearNew);

    hist[h] = c;

    histLen = (histLen < h) ? h : histLen;

    fgets(L, 1024, HF);
  }

  histLen++;

  nUnique = hist[1];

  for (uint32 hh=0; hh<histLen; hh++) {
    nDistinct += hist[hh];
    nTotal    += hist[hh] * hh;
  }
}




int
main(int argc, char **argv) {
  char              *gkpPath = 0L;
  char              *merCountsFile = 0L;
  char              *histogramFile = 0L;

  double             expectedCoverage = 0;
  double             guessedCoverage = 0;

  argc = AS_configure(argc, argv);

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-m") == 0) {
      merCountsFile = argv[++arg];

    } else if (strcmp(argv[arg], "-h") == 0) {
      histogramFile = argv[++arg];

    } else if (strcmp(argv[arg], "-c") == 0) {
      expectedCoverage = atof(argv[++arg]);

    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if (((merCountsFile == NULL) && (histogramFile == NULL)) || (err)) {
    fprintf(stderr, "usage: %s [-c coverage] [-m mercounts] [-h histogram]\n", argv[0]);
    fprintf(stderr, "INPUTS: (exactly one)\n");
    fprintf(stderr, "  -m mercounts    file of mercounts from meryl\n");
    fprintf(stderr, "  -h histogram    histogram from meryl\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "  -c coverage     expected coverage of reads\n");
    exit(1);
  }

  uint64  nDistinct = 0;
  uint64  nUnique   = 0;
  uint64  nTotal    = 0;

  uint32   histLen  = 0;
  uint32  *hist     = NULL;

  if (merCountsFile) {
    merylStreamReader *MF = new merylStreamReader(merCountsFile);
    loadHistogram(MF, nDistinct, nUnique, nTotal, histLen, hist);
    delete MF;
  }

  if (histogramFile) {
    FILE *HF = fopen(histogramFile, "r");
    loadHistogram(HF, nDistinct, nUnique, nTotal, histLen, hist);
    fclose(HF);
  }

  //  Examine the counts, pick a reasonable upper limit.

  fprintf(stderr, "RAW MER COUNTS:\n");
  fprintf(stderr, "  distinct: %12" F_U64P " (different kmer sequences)\n", nDistinct);
  fprintf(stderr, "  unique:   %12" F_U64P " (single-copy kmers)\n",        nUnique);
  fprintf(stderr, "  total:    %12" F_U64P " (kmers in sequences)\n",       nTotal);
  fprintf(stderr, "\n");

  guessedCoverage = (expectedCoverage > 0) ? expectedCoverage : guessCoverage(hist, histLen);



  //  Pass 1: look for a reasonable limit, using %distinct and %total.
  //
  uint64  totalUsefulDistinct = nDistinct - nUnique;
  uint64  totalUsefulAll      = nTotal    - nUnique;
  uint64  distinct            = 0;
  uint64  total               = 0;

  uint32  maxCount = 0;
  uint32  extCount = 0;

  uint32  kk = 2;

  //  If we cover 99% of all the distinct mers, that's reasonable.
  //
  //  If we're a somewhat high count, and we're covering 2/3 of the total mers, assume that there
  //  are lots of errors (or polymorphism) that are preventing us from covering many distinct mers.
  //

  for (; kk < histLen; kk++) {
    distinct += hist[kk];
    total    += hist[kk] * kk;

    if (((distinct / (double)totalUsefulDistinct) > 0.9975) ||
        (((total   / (double)totalUsefulAll)      > 0.6667) && (kk > 50 * guessedCoverage))) {
      maxCount = kk;
      break;
    }
  }

  fprintf(stderr, "Set maxCount to " F_U32 " (" F_U32 " kmers), which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
          maxCount,
          hist[maxCount],
          100.0 * distinct / totalUsefulDistinct,
          100.0 * total    / totalUsefulAll);

  //  Compute an average number of kmers around this count.

  int32  min = (maxCount - 25 < 2)       ? 2       : maxCount - 25;
  int32  max = (maxCount + 25 > histLen) ? histLen : maxCount + 26;
  int64  avg = 0;
  int64  tot = 0;

  for (int32 ii=min; ii<max; ii++) {
    avg += ii * hist[ii];
    tot += ii * hist[ii];
  }

  avg /= 51;

  fprintf(stderr, "Average number of kmers between count " F_S32 " and " F_S32 " is " F_S64 "\n", min, max, avg);

  //  Scan forward until we find a big gap in kmers.  This lets us ignore wildly over represented
  //  kmers in small assemblies.

  while (max < histLen) {
    if (tot == 0)
      break;

    tot -= min * hist[min];  //  Subtract out the last min, move ahead one.
    tot += max * hist[max];  //  Add in the next max, move ahead one.

    min++;
    max++;
  }

  uint32  limit;

  if (tot > 0) {
    limit = histLen;
    fprintf(stderr, "No break in kmer coverage found.\n");
  } else {
    limit = min;
    fprintf(stderr, "Found break in kmer coverage at %d\n", limit);
  }

  //  Scan forward until we see a 10x increase in the number of kmers, OR, until we cover 90% of the sequence and are at sufficient coverage.

  for (; kk < limit; kk++) {
    if (avg * 10 < kk * hist[kk])
      break;

    if (((double)total / totalUsefulAll >= 0.9) && (kk > 50 * guessedCoverage))
      break;

    distinct += hist[kk];
    total    += hist[kk] * kk;

    if (hist[kk] > 0)
      extCount = kk;
  }

  if (extCount > 0)
    maxCount = extCount;

  fprintf(stderr, "Set maxCount to " F_U32 " (" F_U32 " kmers), which will cover %.2f%% of distinct mers and %.2f%% of all mers.\n",
          maxCount,
          hist[maxCount],
          100.0 * distinct / totalUsefulDistinct,
          100.0 * total    / totalUsefulAll);

  fprintf(stdout, F_U32"\n", maxCount);

  return(0);
}
