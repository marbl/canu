
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: splitUnitigs.C,v 1.2 2011-12-19 02:20:06 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "MultiAlign.h"
#include "MultiAlignStore.h"
#include "MultiAlignment_CNS.h"
#include "MultiAlignment_CNS_private.h"

#include <map>
#include <set>
#include <vector>

using namespace std;

#define READ_TRIM_BASES          (AS_OVERLAP_MIN_LEN / 2 - 1)
#define MAX_SEQUENCE_COVERAGE     1
#define MIN_BAD_CLONE_COVERAGE    3
#define MAX_GOOD_CLONE_COVERAGE   1
#define LN_2                      0.693147

#define CGW_CUTOFF 5

//gkStore          *gkpStore     = NULL;
//MultiAlignStore  *tigStore     = NULL;
AS_IID           *matePair     = NULL;
AS_IID           *library      = NULL;

CNS_Options options = { CNS_OPTIONS_SPLIT_ALLELES_DEFAULT,
                        CNS_OPTIONS_MIN_ANCHOR_DEFAULT,
                        CNS_OPTIONS_DO_PHASING_DEFAULT };



struct splitInterval {
  int32  bgn;
  int32  end;
  bool   isGood;
};

inline
void
incrementInterval(uint32 *coverage,
                  int32   minpos,
                  int32   maxpos) {
  for (uint32 i=minpos; i<maxpos; i++)
    coverage[i]++;
}


static
void
createReadCoverageMap(uint32       *rc,
                      MultiAlignT  *ma,
                      uint32        maLen,
                      uint32        minSplit) {
  uint32 n = GetNumIntMultiPoss(ma->f_list);

  for (uint32 i=0; i<n; i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    int32 minPos = MIN(imp->position.bgn, imp->position.end) + READ_TRIM_BASES;
    int32 maxPos = MAX(imp->position.bgn, imp->position.end) - READ_TRIM_BASES;

    //  Undo the offsets if we're near the end of a unitig where we expect coverage to be low.

#if 1
    if ((i < 5) || (maxPos < minSplit)) {
      minPos -= READ_TRIM_BASES;
      maxPos += READ_TRIM_BASES;
    }

    if ((i >= n-5) || (minPos + minSplit > maLen)) {
      minPos -= READ_TRIM_BASES;
      maxPos += READ_TRIM_BASES;
    }

    minPos = MAX(0, minPos);
    maxPos = MIN(maLen, maxPos);
#endif

    incrementInterval(rc, minPos, maxPos);
  }
}


static
void
createCloneCoverageMapExternal(uint32       *gcc,
                               uint32       *bcc,
                               IntMultiPos  *imp,
                               MultiAlignT  *ma,
                               uint32        maLen) {
  bool       frgforward = (imp->position.bgn < imp->position.end);
  int32      frgmin     = MIN(imp->position.bgn, imp->position.end);
  int32      frgmax     = MAX(imp->position.bgn, imp->position.end);

  gkLibrary *lib       = gkpStore->gkStore_getLibrary(library[imp->ident]);
  uint32     liborient = lib->orientation;
  int32      distMin   = lib->mean - CGW_CUTOFF * lib->stddev;
  int32      distMax   = lib->mean + CGW_CUTOFF * lib->stddev;

  switch (liborient) {
    case AS_READ_ORIENT_UNKNOWN:
      break;

    case AS_READ_ORIENT_INNIE:
      if ((frgforward == true)  && (frgmin < maLen - distMax))   incrementInterval(bcc, MAX(0, frgmin),           MIN(maLen, frgmin + distMax));
      if ((frgforward == false) && (frgmin >         distMax))   incrementInterval(bcc, MAX(0, frgmin - distMax), MIN(maLen, frgmax));
      break;

    case AS_READ_ORIENT_OUTTIE:
      if ((frgforward == false) && (frgmin < maLen - distMax))   incrementInterval(bcc, MAX(0, frgmin),           MIN(maLen, frgmin + distMax));
      if ((frgforward == true)  && (frgmin >         distMax))   incrementInterval(bcc, MAX(0, frgmin - distMax), MIN(maLen, frgmax));
      break;

    case AS_READ_ORIENT_NORMAL:
      assert(0);
      break;

    case AS_READ_ORIENT_ANTINORMAL:
      assert(0);
      break;

    default:
      break;
  }
}


static
void
createCloneCoverageMapInternal(uint32       *gcc,
                               uint32       *bcc,
                               IntMultiPos  *imp,
                               IntMultiPos  *mmp,
                               MultiAlignT  *ma,
                               uint32        maLen) {
  bool       frgforward = (imp->position.bgn < imp->position.end);
  int32      frgmin     = MIN(imp->position.bgn, imp->position.end);
  int32      frgmax     = MAX(imp->position.bgn, imp->position.end);

  bool       matforward = (mmp->position.bgn < mmp->position.end);
  int32      matmin     = MIN(mmp->position.bgn, mmp->position.end);
  int32      matmax     = MAX(mmp->position.bgn, mmp->position.end);

  assert(frgmin < matmin);

  gkLibrary *lib       = gkpStore->gkStore_getLibrary(library[imp->ident]);
  uint32     liborient = lib->orientation;
  int32      distMin   = lib->mean - CGW_CUTOFF * lib->stddev;
  int32      distMax   = lib->mean + CGW_CUTOFF * lib->stddev;

  int32  distance = matmax - frgmin;

  if (frgforward == matforward)
    distance = 0;

  //  Well, this turned out simpler than expected.  We assume the fragments come in ordered (imp
  //  before mmp - that's the assert above).  If the fragments are in the same orientation we set
  //  the distance to zero.
  //
  //  Then, for innie oriented fragments, the pair is good if the first fragment is forward and the
  //  distance is correct.  For the bad pairs, the bad region depends only on the orientation of the
  //  fragment.

  switch (liborient) {
    case AS_READ_ORIENT_UNKNOWN:
      break;

      case AS_READ_ORIENT_INNIE:
        if ((frgforward == true) && (distMin < distance) && (distance < distMax)) {
          incrementInterval(gcc, frgmin, matmax);
        } else {
          if (frgforward == true)   incrementInterval(bcc, MAX(0, frgmin),           MIN(maLen, frgmin + distMax));
          if (matforward == true)   incrementInterval(bcc, MAX(0, matmin),           MIN(maLen, matmin + distMax));

          if (frgforward == false)  incrementInterval(bcc, MAX(0, frgmin - distMax), MIN(maLen, frgmax));
          if (matforward == false)  incrementInterval(bcc, MAX(0, matmin - distMax), MIN(maLen, matmax));
        }
        break;

      case AS_READ_ORIENT_OUTTIE:
        if ((frgforward == false) && (distMin < distance) && (distance < distMax)) {
          incrementInterval(gcc, frgmin, matmax);
        } else {
          if (frgforward == false)  incrementInterval(bcc, MAX(0, frgmin),           MIN(maLen, frgmin + distMax));
          if (matforward == false)  incrementInterval(bcc, MAX(0, matmin),           MIN(maLen, matmin + distMax));

          if (frgforward == true)   incrementInterval(bcc, MAX(0, frgmin - distMax), MIN(maLen, frgmax));
          if (matforward == true)   incrementInterval(bcc, MAX(0, matmin - distMax), MIN(maLen, matmax));
        }
        break;

      case AS_READ_ORIENT_NORMAL:
        assert(0);
        break;

      case AS_READ_ORIENT_ANTINORMAL:
        assert(0);
        break;

      default:
        break;
  }
}


static
void
createCloneCoverageMap(uint32       *gcc,
                       uint32       *bcc,
                       MultiAlignT  *ma,
                       uint32        maLen) {
  
  map<AS_IID, uint32>  fragIdx;
  uint32               numFrag = GetNumIntMultiPoss(ma->f_list);

  for (uint32 i=0; i<numFrag; i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);

    fragIdx[imp->ident] = i;
  }


  for (uint32 i=0; i<numFrag; i++) {
    IntMultiPos *imp = GetIntMultiPos(ma->f_list, i);
    AS_IID       mmi = matePair[imp->ident];

    if (mmi == 0)
      //  Not mated
      continue;

    if (fragIdx.find(mmi) == fragIdx.end()) {
      createCloneCoverageMapExternal(gcc, bcc, imp, ma, maLen);

    } else {
      IntMultiPos *mmp = GetIntMultiPos(ma->f_list, fragIdx[mmi]);

      if (MIN(imp->position.bgn, imp->position.end) < MIN(mmp->position.bgn, mmp->position.end))
        //  Compute clone coverage for the pair of fragments once - we'll see both fragments
        //  as 'imp', but we only compute for one ordering
        createCloneCoverageMapInternal(gcc, bcc, imp, mmp, ma, maLen);
    }
  }
}





int
main(int argc, char **argv) {
  char      *gkpName   = NULL;
  char      *tigName   = NULL;
  int32      tigVers   = -1;

  AS_IID     bgnID     = 0;
  AS_IID     onlID     = UINT32_MAX;
  AS_IID     endID     = UINT32_MAX;

  int32      minLength = INT32_MAX;
  int32      minSplit  = INT32_MAX;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-t") == 0) {
      tigName = argv[++arg];
      tigVers = atoi(argv[++arg]);

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
    fprintf(stderr, "  -g         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -t         Mandatory path to a tigStore (can exist or not).\n");

    if (gkpName == NULL)
      fprintf(stderr, "No gatekeeper store (-G option) supplied.\n");

    if (tigName == NULL)
      fprintf(stderr, "No output tigStore (-T option) supplied.\n");

    exit(1);
  }


  //  The gatekeeper store MUST have an updated insert size.  We DO NOT recompute it here.

  gkpStore     = new gkStore(gkpName, FALSE, FALSE);

  tigStore     = new MultiAlignStore(tigName, tigVers, 0, 0, TRUE, FALSE, FALSE);

  matePair     = new AS_IID [gkpStore->gkStore_getNumFragments() + 1];
  library      = new AS_IID [gkpStore->gkStore_getNumFragments() + 1];


  gkStream   *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment  fr;

  while(fs->next(&fr)) {
    uint32 iid = fr.gkFragment_getReadIID();

    matePair[iid] = fr.gkFragment_getMateIID();
    library[iid]  = fr.gkFragment_getLibraryIID();

    if ((iid % 10000000) == 0)
      fprintf(stderr, "Loading fragment information %9d out of %9d\n", iid, gkpStore->gkStore_getNumFragments());
  }
  
  delete fs;

  //  Determine minimums

  for (uint32 i=1; i<gkpStore->gkStore_getNumLibraries() + 1; i++) {
    gkLibrary *lib  = gkpStore->gkStore_getLibrary(i);

    fprintf(stderr, "LIB %d %f +- %f\n", i, lib->mean, lib->stddev);

    if (lib->mean > 0) {
      minLength = MIN(minLength, lib->mean + CGW_CUTOFF * lib->stddev);
      minSplit  = MIN(minSplit,  lib->mean - CGW_CUTOFF * lib->stddev);
    }
  }

  fprintf(stderr, "minLength = %d\n", minLength);
  fprintf(stderr, "minSplit = %d\n", minSplit);

  //  Over every unitig, analyze and maybe split

  bgnID = 0;
  endID = tigStore->numUnitigs();

  if (onlID < endID) {
    bgnID = onlID;
    endID = onlID + 1;
  }
  
  fprintf(stderr, "Analyzing unitig for b="F_U32" to e="F_U32"\n", bgnID, endID);

  uint32   rcMax = 1024 * 1024;
  uint32   rcBgn = 0;
  uint32   rcEnd = 0;
  uint32  *rc    = new uint32 [rcMax];
  uint32  *gcc   = new uint32 [rcMax];
  uint32  *bcc   = new uint32 [rcMax];

  for (uint32 i=bgnID; i<endID; i++) {
    MultiAlignT  *maOrig = tigStore->loadMultiAlign(i, TRUE);

    if (maOrig == NULL)
      continue;

    uint32        maLen  = GetMultiAlignLength(maOrig);

    if (maLen < minLength)
      continue;

    while (rcMax < maLen) {
      rcMax *= 2;
      delete [] rc;     rc  = new uint32 [rcMax];
      delete [] gcc;    gcc = new uint32 [rcMax];
      delete [] bcc;    bcc = new uint32 [rcMax];
    }

    memset(rc, 0, sizeof(uint32) * maLen);

    createReadCoverageMap(rc, maOrig, maLen, minSplit);

    //  Pick out a thick region in the middle to analyze - ignore low coverage at the end of the
    //  unitig.
    //
    //  We tried ignoring the ends with no good clone coverage, but found plenty of examples where
    //  we'd want to trim off those ends.  Example: no good clone coverage, a hump of read coverage
    //  that drops back to one, and bad clone coverage.
    //
    //  The original would skip single coverage regions only.  This wasn't working as a single
    //  contain, or overlap would make it stop.  We'd then split on the next 1-coverage area,
    //  trimming off a single read (or two or three) that likely isn't a problem.  There was some
    //  magic in the original that somehow skipped these bad splits; I couldn't find it.

    uint32 minBase = READ_TRIM_BASES;
    uint32 maxBase = maLen - READ_TRIM_BASES;
    uint32 curBase = 0;

    //while ((minBase < maLen - READ_TRIM_BASES) && (gcc[minBase] == 0))
    //  minBase++;
    while ((minBase < maLen - READ_TRIM_BASES) && (rc[minBase] <= 2))
      minBase++;

    //while ((maxBase > READ_TRIM_BASES) && (gcc[maxBase] == 1))
    //  maxBase--;
    while ((maxBase > READ_TRIM_BASES) && (rc[maxBase] <= 2))
      maxBase--;

    //  Find a first candidate interval

    for (curBase=minBase; curBase<maxBase; curBase++)
      if (rc[curBase] <= MAX_SEQUENCE_COVERAGE)
        break;

    if (curBase >= maxBase)
      //  No interval, all good!
      continue;

    //  Read coverage was bad, what is mate coverage doing?

    memset(gcc, 0, sizeof(uint32) * maLen);
    memset(bcc, 0, sizeof(uint32) * maLen);

    createCloneCoverageMap(gcc, bcc, maOrig, maLen);

    //  Identify potential chimeric intervals

    bool                   inInterval = false;
    splitInterval          interval;
    vector<splitInterval>  intervals;

    interval.bgn    = 0;
    interval.end    = curBase;
    interval.isGood = true;


    for (uint32 chkBase=curBase; chkBase<maxBase; chkBase++) {
      if ((rc[chkBase]  <= MAX_SEQUENCE_COVERAGE) &&
          (bcc[chkBase] >= MIN_BAD_CLONE_COVERAGE) &&
          (gcc[chkBase] <= MAX_GOOD_CLONE_COVERAGE)) {
        //  A bad part of town.
        if (interval.isGood == true) {
          //  Switching from good to bad.
          intervals.push_back(interval);
          interval.bgn    = chkBase;
          interval.end    = chkBase;
          interval.isGood = false;
        } else {
          //  Extending a bad
          interval.end    = chkBase;
        }
      } else {
        //  Upper middle class, nice cars, etc.
        if (interval.isGood == false) {
          //  Moving on up!
          intervals.push_back(interval);
          interval.bgn    = chkBase;
          interval.end    = chkBase;
          interval.isGood = true;
        } else {
          //  Extending a good
          interval.end    = chkBase;
        }
      }
    }

    //  Add the last interval, and maybe a final good one to span the unitig.
    if (interval.isGood == true) {
      interval.end = maLen;
      intervals.push_back(interval);
    } else {
      intervals.push_back(interval);
      interval.bgn = maxBase;
      interval.end = maLen;
      interval.isGood = true;
      intervals.push_back(interval);
    }

    if (intervals.size() == 0)
      fprintf(stderr, "Found no intervals for unitig %d\n", maOrig->maID);

    if (intervals.size() <= 1)
      continue;

#if 1
    {
      char  N[FILENAME_MAX];
      FILE *F;

      sprintf(N, "splitUnitigs-%08d.dat", maOrig->maID);
      F = fopen(N, "w");
      for (uint32 i=0; i<maLen; i++)
        fprintf(F, "%d\t%u\t%u\t%u\n", i, rc[i], bcc[i], gcc[i]);
      fclose(F);

      sprintf(N, "splitUnitigs-%08d.gp", maOrig->maID);
      F = fopen(N, "w");
      fprintf(F, "set terminal png\n");
      fprintf(F, "set output \"splitUnitigs-%08d.png\"\n", maOrig->maID);
      fprintf(F, "plot \"splitUnitigs-%08d.dat\" using 1:2 with lines title \"RC\", \"splitUnitigs-%08d.dat\" using 1:3 with lines title \"BCC\", \"splitUnitigs-%08d.dat\" using 1:4 with lines title \"GCC\"\n",
            maOrig->maID, maOrig->maID, maOrig->maID);
      fclose(F);

      sprintf(N, "gnuplot < splitUnitigs-%08d.gp", maOrig->maID);
      system(N);
    }
#endif

    //  Be nice and report the intervals

    for (uint32 i=0; i<intervals.size(); i++)
      fprintf(stderr, "unitig %d interval %2d %d,%d %s\n",
              maOrig->maID, i, intervals[i].bgn, intervals[i].end, intervals[i].isGood ? "good" : "bad");

    //  Merge bad intervals that are separated by a tiny good interval?  Nah, just do it on the fly
    //  when splitting.

    //  The splitting works in multiple passes.  The first pass, any fragment that touches a bad
    //  interval is moved to the corresponding bad unitig.  We don't know exactly where the chimeric
    //  break is, so this is the best we can do to isolate it (other than shattering).  The second
    //  pass moves fragments from the original layout to new unitigs as long as they're contiguous.

    //  One stupid case needs to be cleaned up.  The first (few) fragments in a good region could be
    //  contained by the fragment in the bad region, and not actially connected to the good region.
    //  This is handled using maBgn/maEnd when a fragment is added to a good region.

    MultiAlignT  **maNew = new MultiAlignT * [intervals.size()];
    int32         *maBgn = new int32 [intervals.size()];
    int32         *maEnd = new int32 [intervals.size()];

    for (uint32 i=0; i<intervals.size(); i++) {
      maNew[i] = CreateEmptyMultiAlignT();
      maBgn[i] = maLen;
      maEnd[i] = 0;
    }

    uint32 n = GetNumIntMultiPoss(maOrig->f_list);

    for (uint32 i=0; i<n; i++) {
      IntMultiPos *imp    = GetIntMultiPos(maOrig->f_list, i);

      if (imp->ident == 0)
        continue;

      int32       minpos = MIN(imp->position.bgn, imp->position.end);
      int32       maxpos = MAX(imp->position.bgn, imp->position.end);
      int32       dest   = INT32_MAX;

      //  Search for a destination bad unitig for this fragment.
      for (uint32 ii=0; ii<intervals.size(); ii++)
        if ((intervals[ii].isGood == false) && (minpos <= intervals[ii].end) && (intervals[ii].bgn <= maxpos)) {
          dest = ii;
          break;
        }

      //  If the fragment touched a bad interval, dest is set, and we move the fragment to that new unitig.
      if (dest < intervals.size()) {
        AppendVA_IntMultiPos(maNew[dest]->f_list, imp);
        maBgn[dest] = MIN(maBgn[dest], minpos);
        maEnd[dest] = MAX(maEnd[dest], maxpos);
        imp->ident = 0;
        continue;
      }

      //  Otherwise, search for a destination good unitig.
      for (uint32 ii=0; ii<intervals.size(); ii++)
        if ((intervals[ii].isGood == true) && (minpos <= intervals[ii].end) && (intervals[ii].bgn <= maxpos)) {
          dest = ii;
          break;
        }


      //  If the new fragment doesn't overlap with what is already in this interval, and the stuff in this
      //  interval does overlap with the previous interval, move it all there.
      if ((dest > 0) &&
          (maEnd[dest] - AS_OVERLAP_MIN_LEN <= minpos) &&
          (maBgn[dest]                      <= maEnd[dest-1] - AS_OVERLAP_MIN_LEN)) {

        fprintf(stderr, "Fixing contains.\n");

        for (uint32 iii=0; iii<GetNumIntMultiPoss(maNew[dest]->f_list); iii++) {
          IntMultiPos *ttt = GetIntMultiPos(maNew[dest]->f_list, iii);

          AppendVA_IntMultiPos(maNew[dest-1]->f_list, ttt);

          fprintf(stderr, " prev %d,%d -- %d %d,%d (no overlap to new %d,%d)\n",
                  maBgn[dest], maEnd[dest],
                  ttt->ident, ttt->position.bgn, ttt->position.end,
                  minpos, maxpos);

          maBgn[dest-1] = MIN(ttt->position.bgn, maBgn[dest-1]);
          maBgn[dest-1] = MIN(ttt->position.end, maBgn[dest-1]);

          maEnd[dest-1] = MAX(ttt->position.bgn, maEnd[dest-1]);
          maEnd[dest-1] = MAX(ttt->position.end, maEnd[dest-1]);
        }

        ResetToRange_VA(maNew[dest]->f_list, 0);
        maBgn[dest] = maLen;
        maEnd[dest] = 0;
      }


      //  If it touched a good interval, move it there.
      if (dest < intervals.size()) {
        AppendVA_IntMultiPos(maNew[dest]->f_list, imp);
        maBgn[dest] = MIN(maBgn[dest], minpos);
        maEnd[dest] = MAX(maEnd[dest], maxpos);
        imp->ident = 0;
        continue;
      }

      //  We should never get here.  The original unitig should be covered completely by intervals.
      assert(0);
    }


    //  Run these through consensus (if the original had a consensus sequence) and add to the store.

    for (uint32 i=0; i<intervals.size(); i++) {
      if (GetNumIntMultiPoss(maNew[i]->f_list) == 0)
        //  Possibly a tiny good interval between two bad intervals.  Not sure how this would occur
        //  though.
        continue;

      maNew[i]->maID = tigStore->numUnitigs();

      fprintf(stderr, "Creating new unitig %d with "F_SIZE_T" fragments\n",
              maNew[i]->maID, GetNumIntMultiPoss(maNew[i]->f_list));

      if (GetNumchars(maOrig->consensus) > 1)
        if (MultiAlignUnitig(maNew[i], gkpStore, &options, NULL) == false)
          fprintf(stderr, "  Unitig %d FAILED.\n", maNew[i]->maID);

      tigStore->insertMultiAlign(maNew[i], TRUE, FALSE);

      DeleteMultiAlignT(maNew[i]);
    }

    delete [] maNew;

    //  Now mark the original unitig as deleted.

    tigStore->deleteMultiAlign(maOrig->maID, true);
  }

  delete [] rc;
  delete [] gcc;
  delete [] bcc;

  delete [] matePair;
  delete [] library;

  delete gkpStore;
  delete tigStore;
}
