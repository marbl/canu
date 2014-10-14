
/**************************************************************************
 * Copyright (C) 2014, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id:  $";

#include "AS_global.H"
#include "AS_PER_gkpStore.H"
#include "MultiAlign.H"
#include "MultiAlignStore.H"

#include "AS_CGB_histo.H"

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
ruLabelStat  repeat_MicroHet;
ruLabelStat  repeat_SingleSpan;
ruLabelStat  repeat_LowCov;

ruLabelStat  repeat_IsSingleton;
ruLabelStat  repeat_IsUnique;
ruLabelStat  repeat_IsRepeat;





intervalList<int32> *
computeCoverage(MultiAlignT *ma) {
  uint32                 maNum = GetNumIntMultiPoss(ma->f_list);
  intervalList<int32>    IL;

  for (uint32 ii=0; ii<maNum; ii++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, ii);

    int32   bgn = MIN(imp->position.bgn, imp->position.end);
    int32   end = MAX(imp->position.bgn, imp->position.end);
    int32   len = end - bgn;

    IL.add(bgn, len);
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

  //  MicroHet probability is actually the probability of the sequence being UNIQUE, based on
  //  microhet considerations.  Falling below threshhold makes something a repeat.
  double            cgbApplyMicrohetCutoff    = -1;     //  Basically turns it off, unless enabled
  double            cgbMicrohetProb           = 1.e-5;  //  Scores less than this are considered repeats

  double            cgbUniqueCutoff           = CGB_UNIQUE_CUTOFF;
  double            cgbDefinitelyUniqueCutoff = CGB_UNIQUE_CUTOFF;
  double            singleReadMaxCoverage     = 1.0;    //  Reads covering more than this will demote the unitig
  uint32            lowCovDepth               = 2;
  double            lowCovFractionAllowed     = 1.0;
  uint32            tooLong                   = UINT32_MAX;
  uint32            tooShort                  = 1000;
  uint32            minReads                  = 2;
  double            minPopulous               = 0;

  AS_IID            bgnID = 0;
  AS_IID            endID = 0;
  AS_IID            maxID = 0;

  MultiAlignStore  *tigStore = NULL;
  gkStore          *gkpStore = NULL;

  bool              doUpdate = true;
  const bool        isUnitig = true; // Required for tigStore accessor

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outPrefix = argv[++arg];

    } else if (strcmp(argv[arg], "-n") == 0) {
      doUpdate = false;


    } else if (strcmp(argv[arg], "-e") == 0) {
      cgbMicrohetProb = atof(argv[++arg]);

    } else if (strcmp(argv[arg], "-i") == 0) {
      cgbApplyMicrohetCutoff = atof(argv[++arg]);

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
    fprintf(stderr, "  -g <G>       Mandatory, path G to a gkpStore directory.\n");
    fprintf(stderr, "  -t <T> <v>   Mandatory, path T to a tigStore, and version V.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -e P         Microhet probability (default 1e-5)\n");
    fprintf(stderr, "  -i C         Microhet cutoff (default -1)\n");
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
    fprintf(stderr, "  6)  A unitig with microhet probability P (-e) and astat below C (-i) is NOT unique.\n");
    fprintf(stderr, "  7)  A unitig with fraction F below coverage D (-lowcov) is NOT unique.\n");
    fprintf(stderr, "  8)  A unitig shorter than S (-short) bases long is NOT unique.\n");
    fprintf(stderr, "  9)  Otherwise, the unitig IS unique.\n");

    if (gkpName == NULL)
      fprintf(stderr, "No gatekeeper store (-g option) supplied.\n");

    if (tigName == NULL)
      fprintf(stderr, "No input tigStore (-t option) supplied.\n");

    if (outPrefix == NULL)
      fprintf(stderr, "No output prefix (-o option) supplied.\n");

    exit(1);
  }

  gkpStore     = new gkStore(gkpName, false, false);
  tigStore     = new MultiAlignStore(tigName, tigVers, 0, 0, doUpdate, doUpdate, false);

  if (endID == 0)
    endID = tigStore->numUnitigs();

  maxID = tigStore->numUnitigs();

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

  bool       *isNonRandom = new bool   [gkpStore->gkStore_getNumFragments() + 1];
  uint32     *fragLength  = new uint32 [gkpStore->gkStore_getNumFragments() + 1];

  gkStream   *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment  fr;

  while(fs->next(&fr)) {
    uint32 iid = fr.gkFragment_getReadIID();

    isNonRandom[iid] = fr.gkFragment_getIsNonRandom();
    fragLength[iid]  = fr.gkFragment_getClearRegionLength(AS_READ_CLEAR_OBTCHIMERA);

    if ((iid % 10000000) == 0)
      fprintf(stderr, "Loading fragment information %9d out of %9d\n", iid, gkpStore->gkStore_getNumFragments());
  }
  
  delete fs;

  //
  //  Compute:
  //    global coverage histogram
  //    global single-read fraction covered - number of bases in unitigs with specific fraction covered
  //    global number of reads per unitig
  //

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
    MultiAlignT  *ma    = tigStore->loadMultiAlign(uu, isUnitig);

    if (ma == NULL)
      continue;

    utgCovHistogram[ma->maID] = utgCovData + ma->maID * lowCovDepth;

    if (GetNumIntMultiPoss(ma->f_list) == 1)
      continue;

    //  This MUST use gapped lengths, otherwise we'd need to translate all the read coords from
    //  gapped to ungapped.

    uint32        maLen = GetMultiAlignLength(ma);
    uint32        maNum = GetNumIntMultiPoss(ma->f_list);

    assert(ma->data.num_frags == GetNumIntMultiPoss(ma->f_list));

    //  Global coverage histogram.

    intervalList<int32>  *ID = computeCoverage(ma);

    for (uint32 ii=0; ii<ID->numberOfIntervals(); ii++) {
      if (ID->depth(ii) < lowCovDepth)
        utgCovHistogram[ma->maID][ID->depth(ii)] += ID->hi(ii) - ID->lo(ii) + 1;

      if (ID->depth(ii) < covHistogramMax)
        covHistogram[ID->depth(ii)] += ID->hi(ii) - ID->lo(ii) + 1;
    }

    //  Single read max fraction covered.

    uint32  covMax = 0;
    uint32  cov;

    for (uint32 ff=0; ff<maNum; ff++) {
      IntMultiPos  *frg = GetIntMultiPos(ma->f_list, ff);

      if (frg->position.bgn < frg->position.end)
        cov = 1000 * (frg->position.end - frg->position.bgn) / maLen;
      else
        cov = 1000 * (frg->position.bgn - frg->position.end) / maLen;

      if (covMax < cov)
        covMax = cov;
    }

    singleReadCoverageHistogram[covMax]++;

    singleReadCoverage[ma->maID] = covMax / 1000.0;

    //  Number of reads per unitig

    numReadsPerUnitig[uu] = maNum;

    //fprintf(stderr, "unitig %u covMax %f\n", ma->maID, covMax / 1000.0);
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

  fprintf(stderr, "Processing unitigs.\n");

  for (uint32 uu=bgnID; uu<endID; uu++) {
    MultiAlignT  *ma = tigStore->loadMultiAlign(uu, isUnitig);

    if (ma == NULL) {
      fprintf(outLOG, "unitig %d not present\n", uu);
      continue;
    }

    //  This uses UNGAPPED lengths, because they make more sense to humans.

    uint32        maLen = GetMultiAlignUngappedLength(ma);
    uint32        maNum = GetNumIntMultiPoss(ma->f_list);

    uint32  lowCovBases = 0;
    for (uint32 ll=0; ll<lowCovDepth; ll++)
      lowCovBases += utgCovHistogram[ma->maID][ll];

    bool          isUnique    = true;
    bool          isSingleton = false;


    if (maNum == 1) {
      fprintf(outLOG, "unitig %d not unique -- singleton\n",
              ma->maID);
      isUnique    = false;
      isSingleton = true;
    }

    else if (maNum < minReads) {
      fprintf(outLOG, "unitig %d not unique -- %u reads, need at least %d\n",
              ma->maID, maNum, minReads);
      repeat_LowReads += maLen;
      isUnique = false;
    }

    else if (singleReadCoverage[ma->maID] > singleReadMaxCoverage) {
      fprintf(outLOG, "unitig %d not unique -- single read spans fraction %f of unitig (>= %f)\n",
              ma->maID,
              singleReadCoverage[ma->maID],
              singleReadMaxCoverage);
      repeat_SingleSpan += maLen;
      isUnique = false;
    }

    else if (maLen >= tooLong) {
      fprintf(outLOG, "unitig %d IS unique -- too long to be repeat, %u > allowed %u\n",
              ma->maID,
              maLen, tooLong);
      isUnique = true;
    }

    else if (tigStore->getUnitigCoverageStat(ma->maID) < cgbUniqueCutoff) {
      fprintf(outLOG, "unitig %d not unique -- coverage stat %d, needs to be at least %f\n",
              ma->maID, tigStore->getUnitigCoverageStat(ma->maID), cgbUniqueCutoff);
      repeat_LowCovStat += maLen;
      isUnique = false;
    }

    else if ((tigStore->getUnitigMicroHetProb(ma->maID) < cgbMicrohetProb) &&
             (tigStore->getUnitigCoverageStat(ma->maID) < cgbApplyMicrohetCutoff)) {
      fprintf(outLOG, "unitig %d not unique -- low microhetprob %f (< %f) and low coverage stat %d (< %f)\n",
              ma->maID,
              tigStore->getUnitigMicroHetProb(ma->maID), cgbMicrohetProb,
              tigStore->getUnitigCoverageStat(ma->maID), cgbApplyMicrohetCutoff);
      repeat_MicroHet += maLen;
      isUnique = false;
    }

    else if ((double)lowCovBases / maLen > lowCovFractionAllowed) {
      fprintf(outLOG, "unitig %d not unique -- too many low coverage bases, %u out of %u bases, fraction %f > allowed %f\n",
              ma->maID,
              lowCovBases, maLen,
              (double)lowCovBases / maLen, lowCovFractionAllowed);
      repeat_LowCov += maLen;
      isUnique = false;
    }

    //  This was an attempt to not blindly call all short unitigs as non-unique.  It didn't work so
    //  well in initial limited testing.  The threshold is arbitrary; older versions used
    //  cgbDefinitelyUniqueCutoff.  If used, be sure to disable the real check after this!
#if 0
    else if ((tigStore->getUnitigCoverageStat(ma->maID) < cgbUniqueCutoff * 10) &&
             (maLen < CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH)) {
      fprintf(outLOG, "unitig %d not unique -- length %d too short, need to be at least %d AND coverage stat %d must be larger than %d\n",
              ma->maID, maLen, CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH,
              tigStore->getUnitigCoverageStat(ma->maID), cgbUniqueCutoff * 10);
      repeat_Short += maLen;
      isUnique = false;
    }
#endif

    else if (maLen < tooShort) {
      fprintf(outLOG, "unitig %d not unique -- length %d too short, need to be at least %d\n",
              ma->maID, maLen, tooShort);
      repeat_Short += maLen;
      isUnique = false;
    }

    else {
      fprintf(outLOG, "unitig %d not repeat -- no test failed\n", ma->maID);
    }

    //
    //  Allow flag to override the rules and force it to be unique or repeat.  AKA, toggling.
    //

    if (isUnique) {
      repeat_IsUnique += maLen;
      tigStore->setUnitigSuggestUnique(ma->maID);
      tigStore->setUnitigSuggestRepeat(ma->maID, false);

    } else if (isSingleton) {
      repeat_IsSingleton += maLen;
      tigStore->setUnitigSuggestUnique(ma->maID, false);
      tigStore->setUnitigSuggestRepeat(ma->maID);

    } else {
      repeat_IsRepeat += maLen;
      tigStore->setUnitigSuggestUnique(ma->maID, false);
      tigStore->setUnitigSuggestRepeat(ma->maID);
    }
  }


  fprintf(stderr, "\n");
  //fprintf(stderr, "Processed %d unitigs with %d fragments.\n", 0, 0);
  fprintf(stderr, "\n");
  fprintf(stderr, "classification     number of unitigs    total length\n");
  fprintf(stderr, "  unique:          %17"F_U32P"  %14"F_U64P"\n", repeat_IsUnique.num,     repeat_IsUnique.len);
  fprintf(stderr, "  singleton:       %17"F_U32P"  %14"F_U64P"\n", repeat_IsSingleton.num,  repeat_IsSingleton.len);
  fprintf(stderr, "  repeat:          %17"F_U32P"  %14"F_U64P"\n", repeat_IsRepeat.num,     repeat_IsRepeat.len);
  fprintf(stderr, "    too few reads: %17"F_U32P"  %14"F_U64P"\n", repeat_LowReads.num,     repeat_LowReads.len);
  fprintf(stderr, "    low cov stat:  %17"F_U32P"  %14"F_U64P"\n", repeat_LowCovStat.num,   repeat_LowCovStat.len);
  fprintf(stderr, "    too short:     %17"F_U32P"  %14"F_U64P"\n", repeat_Short.num,        repeat_Short.len);
  fprintf(stderr, "    microhet:      %17"F_U32P"  %14"F_U64P"\n", repeat_MicroHet.num,     repeat_MicroHet.len);
  fprintf(stderr, "    spanning read: %17"F_U32P"  %14"F_U64P"\n", repeat_SingleSpan.num,   repeat_SingleSpan.len);
  fprintf(stderr, "    low coverage:  %17"F_U32P"  %14"F_U64P"\n", repeat_LowCov.num,       repeat_LowCov.len);



  fclose(outLOG);
  fclose(outSTA);

  delete [] isNonRandom;
  delete [] fragLength;

  delete gkpStore;
  delete tigStore;
}
