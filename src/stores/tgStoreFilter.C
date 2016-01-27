
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
 *    src/AS_BAT/markRepeatUnique.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-MAR-31 to 2014-APR-15
 *      are Copyright 2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2014-APR-14
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-OCT-09 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-30
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "gkStore.H"
#include "tgStore.H"

#include "intervalList.H"

#include <algorithm>

using namespace std;



//  Stats on repeat labeling of input unitigs.
//
class ruLabelStat {
public:
  ruLabelStat() {
    num = 0;
    len = 0;
  };

  void     operator+=(uint64 len_) {
    num++;
    len += len_;
  };

  uint32   num;
  uint64   len;
};

ruLabelStat  repeat_LowReads;
ruLabelStat  repeat_LowCovStat;
ruLabelStat  repeat_Short;
ruLabelStat  repeat_SingleSpan;
ruLabelStat  repeat_LowCov;

ruLabelStat  repeat_IsSingleton;
ruLabelStat  repeat_IsUnique;
ruLabelStat  repeat_IsRepeat;





intervalList<int32> *
computeCoverage(tgTig *tig) {
  intervalList<int32>    IL;

  for (uint32 ii=0; ii<tig->numberOfChildren(); ii++) {
    tgPosition  *pos = tig->getChild(ii);

    IL.add(pos->min(), pos->max() - pos->min());
  }

  return(new intervalList<int32>(IL));
}







int
main(int argc, char **argv) {
  char             *gkpName      = NULL;
  char             *tigName      = NULL;
  int32             tigVers      = -1;

  int64             genomeSize   = 0;

  char              outName[FILENAME_MAX];
  char             *outPrefix = NULL;

  FILE             *outLOG = NULL;
  FILE             *outSTA = NULL;

  //  From the original description of these values (CGB_UNIQUE_CUTOFF):
  //    A threshold value for Gene's coverage statistic.  Values above this value have never been
  //    known to associated with unitigs with fragments that are not contiguous in the genome.
  double            cgbUniqueCutoff           = 10.0;
  double            cgbDefinitelyUniqueCutoff = 10.0;

  double            singleReadMaxCoverage     = 1.0;    //  Reads covering more than this will demote the unitig
  uint32            lowCovDepth               = 2;
  double            lowCovFractionAllowed     = 1.0;
  uint32            tooLong                   = UINT32_MAX;
  uint32            tooShort                  = 1000;
  uint32            minReads                  = 2;
  double            minPopulous               = 0;

  uint32            bgnID = 0;
  uint32            endID = 0;
  uint32            maxID = 0;

  gkStore          *gkpStore = NULL;
  tgStore          *tigStore = NULL;
  tgStoreType       tigMode  = tgStoreModify;


  argc = AS_configure(argc, argv);

  fprintf(stderr, "this is obsolete.  do not use.\n");
  exit(1);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      tigMode  = tgStoreReadOnly;


    } else if (strcmp(argv[arg], "-j") == 0) {
      cgbUniqueCutoff = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-k") == 0) {
      cgbDefinitelyUniqueCutoff = atof(argv[++arg]);


    } else if (strcmp(argv[arg], "-span") == 0) {
      singleReadMaxCoverage = atof(argv[++arg]);  //  If a single read spans more than this fraction of unitig, it is demoted

    } else if (strcmp(argv[arg], "-lowcov") == 0) {
      lowCovDepth            = atoi(argv[++arg]);  //  Coverage below this is too low
      lowCovFractionAllowed  = atof(argv[++arg]);  //  If unitig has more than this fraction low coverage, it is demoted

    } else if (strcmp(argv[arg], "-reads") == 0) {
      arg++;
      minReads    = atoi(argv[arg]);  //  If unitig has fewer than this number of reads it is demoted
      minPopulous = atof(argv[arg]);

      if (minPopulous > 1.0)
        minPopulous = 0;

    } else if (strcmp(argv[arg], "-long") == 0) {
      tooLong = atoi(argv[++arg]);  //  Unitigs longer than this cannot be demoted

    } else if (strcmp(argv[arg], "-short") == 0) {
      tooShort = atoi(argv[++arg]);  //  Unitigs shorter than this are demoted


    } else {
      err++;
    }

    arg++;
  }

  if (gkpName == NULL)
    err++;
  if (tigName == NULL)
    err++;

  if (err) {
    fprintf(stderr, "usage: %s -g gkpStore -t tigStore version\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G <G>       Mandatory, path G to a gkpStore directory.\n");
    fprintf(stderr, "  -T <T> <v>   Mandatory, path T to a tigStore, and version V.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -j J         Unitig is not unique if astat is below J (cgbUniqueCutoff)\n");
    fprintf(stderr, "  -k K         (unused) (cgbDefinitelyUniqueCutoff)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -span F      Unitig is not unique if a single read spans more than fraction F (default 1.0) of unitig\n");
    fprintf(stderr, "  -lowcov D F  Unitig is not unique if fraction F (default 1.0) of unitig is below read depth D (default 2)\n");
    fprintf(stderr, "  -reads R     Unitig is not unique if unitig has fewer than R (default 2) reads\n");
    fprintf(stderr, "               If R is fractional, the least populous unitigs containing fraction R of reads are marked as repeat\n");
    fprintf(stderr, "               Example: unitigs with 9 or fewer reads contain 10%% of the reads.  -reads 0.10 would mark these are repeat.\n");
    fprintf(stderr, "  -long L      Unitig is unique if unitig is at least L (default unlimited) bases long\n");
    fprintf(stderr, "  -short S     Unitig is not unique if unitig is shorter than S (default 1000) bases long\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -o <name>    Prefix for output files.\n");
    fprintf(stderr, "  -n           Do not update the tigStore.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Algorithm:  The first rule to trigger will mark the unitig.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  1)  A unitig with a single read is NOT unique.\n");
    fprintf(stderr, "  2)  A unitig with fewer than R (-reads) reads is NOT unique.\n");
    fprintf(stderr, "  3)  A unitig with a single read spanning fraction F (-span) of the unitig is NOT unique.\n");
    fprintf(stderr, "  4)  A unitig longer than L (-length) bases IS unique.\n");
    fprintf(stderr, "  5)  A unitig with astat less than J (-j) is NOT unique.\n");
    fprintf(stderr, "  6)  A unitig with fraction F below coverage D (-lowcov) is NOT unique.\n");
    fprintf(stderr, "  7)  A unitig shorter than S (-short) bases long is NOT unique.\n");
    fprintf(stderr, "  8)  Otherwise, the unitig IS unique.\n");

    if (gkpName == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (tigName == NULL)
      fprintf(stderr, "No input tigStore (-T option) supplied.\n");

    if (outPrefix == NULL)
      fprintf(stderr, "No output prefix (-o option) supplied.\n");

    exit(1);
  }

  gkpStore     = gkStore::gkStore_open(gkpName, gkStore_readOnly);
  tigStore     = new tgStore(tigName, tigVers, tgStoreReadOnly);

  if (endID == 0)
    endID = tigStore->numTigs();

  maxID = tigStore->numTigs();

  errno = 0;

  sprintf(outName, "%s.log", outPrefix);

  outLOG = fopen(outName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

  fprintf(outLOG, "tigID\trho\tcovStat\tarrDist\n");

  sprintf(outName, "%s.stats", outPrefix);

  outSTA = fopen(outName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outName, strerror(errno)), exit(1);

  fprintf(stderr, "Command Line options:\n");
  fprintf(stderr, "  singleReadMaxCoverage    %f\n", singleReadMaxCoverage);
  fprintf(stderr, "  lowCoverage              %u coverage %f fraction\n", lowCovDepth, lowCovFractionAllowed);
  if (minPopulous > 0)
    fprintf(stderr, "  minReads                 recomputed based on least populous fraction %.4f\n", minPopulous);
  else
    fprintf(stderr, "  minReads                 %u\n", minReads);
  fprintf(stderr, "  tooLong                  %u\n", tooLong);
  fprintf(stderr, "  tooShort                 %u\n", tooShort);

  //
  //  Load fragment data
  //

  fprintf(stderr, "Loading fragment data.\n");

  double      globalRate     = 0;

  bool       *isNonRandom = new bool   [gkpStore->gkStore_getNumReads() + 1];
  uint32     *fragLength  = new uint32 [gkpStore->gkStore_getNumReads() + 1];

  for (uint32 fi=1; fi <= gkpStore->gkStore_getNumReads(); fi++) {
    gkRead     *read = gkpStore->gkStore_getRead(fi);
    gkLibrary  *libr = gkpStore->gkStore_getLibrary(read->gkRead_libraryID());

    isNonRandom[fi] = libr->gkLibrary_isNonRandom();
    fragLength[fi]  = read->gkRead_sequenceLength();
  }

  //
  //  Compute:
  //    global coverage histogram
  //    global single-read fraction covered - number of bases in unitigs with specific fraction covered
  //    global number of reads per unitig
  //
  //  This is using GAPPED lengths, because they're faster, and we don't need the actual ungapped positions.

  fprintf(stderr, "Generating statistics.\n");

  uint32    covHistogramMax      = 1048576;
  uint32   *covHistogram         = new uint32   [covHistogramMax];
  uint32  **utgCovHistogram      = new uint32 * [maxID];
  uint32   *utgCovData           = new uint32   [maxID * lowCovDepth];

  memset(covHistogram,    0, sizeof(uint32) * covHistogramMax);
  memset(utgCovHistogram, 0, sizeof(uint32) * maxID);
  memset(utgCovData,      0, sizeof(uint32) * maxID * lowCovDepth);

  uint32   *singleReadCoverageHistogram   = new uint32 [1001];
  double   *singleReadCoverage            = new double [maxID];

  memset(singleReadCoverageHistogram, 0, sizeof(uint32) * 1001);
  memset(singleReadCoverage,          0, sizeof(double) * maxID);

  uint32    numReadsPerUnitigMax = maxID + 1;
  uint32   *numReadsPerUnitig    = new uint32 [numReadsPerUnitigMax];

  memset(numReadsPerUnitig, 0, sizeof(uint32) * numReadsPerUnitigMax);

  for (uint32 uu=0; uu<maxID; uu++) {
    tgTig  *tig    = tigStore->loadTig(uu);
    uint32  tigLen = tig->length(true);

    if (tig == NULL)
      continue;

    utgCovHistogram[tig->tigID()] = utgCovData + tig->tigID() * lowCovDepth;

    if (tig->numberOfChildren() == 1)
      continue;

    //  Global coverage histogram.

    intervalList<int32>  *ID = computeCoverage(tig);

    for (uint32 ii=0; ii<ID->numberOfIntervals(); ii++) {
      if (ID->depth(ii) < lowCovDepth)
        utgCovHistogram[tig->tigID()][ID->depth(ii)] += ID->hi(ii) - ID->lo(ii) + 1;

      if (ID->depth(ii) < covHistogramMax)
        covHistogram[ID->depth(ii)] += ID->hi(ii) - ID->lo(ii) + 1;
    }

    delete ID;  ID = NULL;

    //  Single read max fraction covered.

    uint32  covMax = 0;
    uint32  cov;

    for (uint32 ff=0; ff<tig->numberOfChildren(); ff++) {
      tgPosition *pos = tig->getChild(ff);

      cov = 1000 * (pos->max() - pos->min()) / tigLen;

      if (covMax < cov)
        covMax = cov;
    }

    singleReadCoverageHistogram[covMax]++;

    singleReadCoverage[tig->tigID()] = covMax / 1000.0;

    //  Number of reads per unitig

    numReadsPerUnitig[uu] = tig->numberOfChildren();

    //fprintf(stderr, "unitig %u covMax %f\n", tig->tigID(), covMax / 1000.0);

    tigStore->unloadTig(uu);
  }

  //
  //  Analyze our collected data, decide on some thresholds.
  //

  fprintf(stderr, "Analyzing statistics.\n");

  if ((minPopulous > 0.0) &&
      (minPopulous < 1.0)) {
    uint32    maxReadsPerUnitig = 0;

    for (uint32 uu=0; uu<maxID; uu++)
      maxReadsPerUnitig = MAX(maxReadsPerUnitig, numReadsPerUnitig[uu] + 1);

    uint32   *totReadsPerNumReads  = new uint32 [maxReadsPerUnitig];
    uint32    totReads             = 0;

    memset(totReadsPerNumReads, 0, sizeof(uint32) * maxReadsPerUnitig + 1);

    for (uint32 uu=0; uu<maxID; uu++) {
      totReadsPerNumReads[numReadsPerUnitig[uu]] += numReadsPerUnitig[uu];
      totReads                                   += numReadsPerUnitig[uu];
    }

    double   xx = 0.0;

    for (minReads=0; xx < minPopulous; minReads++)
      xx += (double)totReadsPerNumReads[minReads] / totReads;

    minReads--;
    xx -= (double)totReadsPerNumReads[minReads] / totReads;

    fprintf(stderr, "  minReads                 %u based on least populous unitigs with fraction %.4f of reads (wanted fraction %.4f)\n",
            minReads, xx, minPopulous);

    uint32  hist[101] = {0};
    double  vals[101] = {0};

    xx = 0.0;

    for (uint32 ff=0; ff<maxReadsPerUnitig; ff++) {
      xx += (double)totReadsPerNumReads[ff] / totReads;
      assert(xx <= 1.0 + 1e9 * ff);  //  For rounding issues
      hist[(int32)(100 * xx)] = ff;
      vals[(int32)(100 * xx)] = xx;
    }
    for (uint32 xx=0; xx<101; xx++)
      if (hist[xx] > 0)
        fprintf(outSTA, "minReads %u with fraction %.04f of reads\n", hist[xx], vals[xx]);

    delete [] totReadsPerNumReads;
  }


  //
  //  Apply the thresholds to unitigs.  The first half of these are the historical CGW rules.
  //

  fprintf(stderr, "Processing unitigs %u to %u.\n", bgnID, endID);

  for (uint32 uu=bgnID; uu<endID; uu++) {
    tgTig  *tig = tigStore->loadTig(uu);

    if (tig == NULL) {
      fprintf(outLOG, "unitig %d not present\n", uu);
      continue;
    }

    //  This uses UNGAPPED lengths, because they make more sense to humans.

    uint32        tigLen = tig->length(false);

    uint32  lowCovBases = 0;
    for (uint32 ll=0; ll<lowCovDepth; ll++)
      lowCovBases += utgCovHistogram[tig->tigID()][ll];

    bool          isUnique    = true;
    bool          isSingleton = false;


    if (tig->numberOfChildren() == 1) {
      fprintf(outLOG, "unitig %d not unique -- singleton\n",
              tig->tigID());
      isUnique    = false;
      isSingleton = true;
    }

    else if (tig->numberOfChildren() < minReads) {
      fprintf(outLOG, "unitig %d not unique -- %u reads, need at least %d\n",
              tig->tigID(), tig->numberOfChildren(), minReads);
      repeat_LowReads += tigLen;
      isUnique = false;
    }

    else if (singleReadCoverage[tig->tigID()] > singleReadMaxCoverage) {
      fprintf(outLOG, "unitig %d not unique -- single read spans fraction %f of unitig (>= %f)\n",
              tig->tigID(),
              singleReadCoverage[tig->tigID()],
              singleReadMaxCoverage);
      repeat_SingleSpan += tigLen;
      isUnique = false;
    }

    else if (tigLen >= tooLong) {
      fprintf(outLOG, "unitig %d IS unique -- too long to be repeat, %u > allowed %u\n",
              tig->tigID(),
              tigLen, tooLong);
      isUnique = true;
    }

    else if (tigStore->getCoverageStat(tig->tigID()) < cgbUniqueCutoff) {
      fprintf(outLOG, "unitig %d not unique -- coverage stat %f, needs to be at least %f\n",
              tig->tigID(), tigStore->getCoverageStat(tig->tigID()), cgbUniqueCutoff);
      repeat_LowCovStat += tigLen;
      isUnique = false;
    }

    else if ((double)lowCovBases / tigLen > lowCovFractionAllowed) {
      fprintf(outLOG, "unitig %d not unique -- too many low coverage bases, %u out of %u bases, fraction %f > allowed %f\n",
              tig->tigID(),
              lowCovBases, tigLen,
              (double)lowCovBases / tigLen, lowCovFractionAllowed);
      repeat_LowCov += tigLen;
      isUnique = false;
    }

    //  This was an attempt to not blindly call all short unitigs as non-unique.  It didn't work so
    //  well in initial limited testing.  The threshold is arbitrary; older versions used
    //  cgbDefinitelyUniqueCutoff.  If used, be sure to disable the real check after this!
#if 0
    else if ((tigStore->getCoverageStat(tig->tigID()) < cgbUniqueCutoff * 10) &&
             (tigLen < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH)) {
      fprintf(outLOG, "unitig %d not unique -- length %d too short, need to be at least %d AND coverage stat %d must be larger than %d\n",
              tig->tigID(), tigLen, CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH,
              tigStore->getCoverageStat(tig->tigID()), cgbUniqueCutoff * 10);
      repeat_Short += tigLen;
      isUnique = false;
    }
#endif

    else if (tigLen < tooShort) {
      fprintf(outLOG, "unitig %d not unique -- length %d too short, need to be at least %d\n",
              tig->tigID(), tigLen, tooShort);
      repeat_Short += tigLen;
      isUnique = false;
    }

    else {
      fprintf(outLOG, "unitig %d not repeat -- no test failed\n", tig->tigID());
    }

    //
    //  Allow flag to override the rules and force it to be unique or repeat.  AKA, toggling.
    //

    if (isUnique) {
      repeat_IsUnique += tigLen;
#warning NOT setSuggestUnique
      //tigStore->setSuggestUnique(tig->tigID());
      tigStore->setSuggestRepeat(tig->tigID(), false);

    } else if (isSingleton) {
      repeat_IsSingleton += tigLen;
#warning NOT setSuggestUnique
      //tigStore->setSuggestUnique(tig->tigID(), false);
      tigStore->setSuggestRepeat(tig->tigID());

    } else {
      repeat_IsRepeat += tigLen;
#warning NOT setSuggestUnique
      //tigStore->setSuggestUnique(tig->tigID(), false);
      tigStore->setSuggestRepeat(tig->tigID());
    }
  }

  fprintf(outSTA, "classification     number of unitigs    total length\n");
  fprintf(outSTA, "  unique:          %17"F_U32P"  %14"F_U64P"\n", repeat_IsUnique.num,     repeat_IsUnique.len);
  fprintf(outSTA, "  singleton:       %17"F_U32P"  %14"F_U64P"\n", repeat_IsSingleton.num,  repeat_IsSingleton.len);
  fprintf(outSTA, "  repeat:          %17"F_U32P"  %14"F_U64P"\n", repeat_IsRepeat.num,     repeat_IsRepeat.len);
  fprintf(outSTA, "    too few reads: %17"F_U32P"  %14"F_U64P"\n", repeat_LowReads.num,     repeat_LowReads.len);
  fprintf(outSTA, "    low cov stat:  %17"F_U32P"  %14"F_U64P"\n", repeat_LowCovStat.num,   repeat_LowCovStat.len);
  fprintf(outSTA, "    too short:     %17"F_U32P"  %14"F_U64P"\n", repeat_Short.num,        repeat_Short.len);
  fprintf(outSTA, "    spanning read: %17"F_U32P"  %14"F_U64P"\n", repeat_SingleSpan.num,   repeat_SingleSpan.len);
  fprintf(outSTA, "    low coverage:  %17"F_U32P"  %14"F_U64P"\n", repeat_LowCov.num,       repeat_LowCov.len);

  fclose(outLOG);
  fclose(outSTA);

  delete [] isNonRandom;
  delete [] fragLength;

  gkpStore->gkStore_close();

  delete tigStore;

  exit(0);
}
