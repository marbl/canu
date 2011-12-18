
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

const char *mainid = "$Id: splitUnitigs.C,v 1.1 2011-12-18 08:13:22 brianwalenz Exp $";

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
#define MAX_GOOD_CLONE_COVERAGE   0
#define LN_2                      0.693147


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
  int32      distMin   = lib->mean - 3 * lib->stddev;
  int32      distMax   = lib->mean + 3 * lib->stddev;

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
  int32      distMin   = lib->mean - 3 * lib->stddev;
  int32      distMax   = lib->mean + 3 * lib->stddev;

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

  uint32     minLength = UINT32_MAX;
  uint32     minSplit  = UINT32_MAX;

  argc = AS_configure(argc, argv);

  int err = 0;
  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-G") == 0) {
      gkpName = argv[++arg];

    } else if (strcmp(argv[arg], "-T") == 0) {
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
    fprintf(stderr, "usage: %s -G gkpStore -T tigStore version\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -G         Mandatory path to a gkpStore.\n");
    fprintf(stderr, "  -T         Mandatory path to a tigStore (can exist or not).\n");

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

    if (lib->mean > 0) {
      minLength = MIN(minLength, lib->mean + 3 * lib->stddev);
      minSplit  = MIN(minSplit,  lib->mean - 3 * lib->stddev);
    }
  }

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
    uint32        maLen  = GetMultiAlignLength(maOrig);

    if (maOrig == NULL)
      continue;

    if (maLen < minLength)
      continue;

    //

    while (rcMax < maLen) {
      rcMax *= 2;
      delete [] rc;     rc  = new uint32 [rcMax];
      delete [] gcc;    gcc = new uint32 [rcMax];
      delete [] bcc;    bcc = new uint32 [rcMax];
    }

    memset(rc, 0, sizeof(uint32) * maLen);

    createReadCoverageMap(rc, maOrig, maLen, minSplit);

    //  Pick out a thick region in the middle to analyze - ignore low coverage
    //  at the end of the unitig

    uint32 minBase = READ_TRIM_BASES;
    uint32 maxBase = maLen - READ_TRIM_BASES;
    uint32 curBase = 0;

    while ((minBase < maLen - READ_TRIM_BASES) && (rc[minBase] <= 1))
      minBase++;

    while ((maxBase > READ_TRIM_BASES) && (rc[maxBase] <= 1))
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

    MultiAlignT  **maNew = new MultiAlignT * [intervals.size()];

    for (uint32 i=0; i<intervals.size(); i++)
      maNew[i] = CreateEmptyMultiAlignT();

    uint32 n = GetNumIntMultiPoss(maOrig->f_list);

    for (uint32 i=0; i<n; i++) {
      IntMultiPos *imp    = GetIntMultiPos(maOrig->f_list, i);

      if (imp->ident == 0)
        continue;

      uint32       minpos = MIN(imp->position.bgn, imp->position.end);
      uint32       maxpos = MAX(imp->position.bgn, imp->position.end);
      uint32       dest   = UINT32_MAX;

      //  Search for a destination bad unitig for this fragment.
      for (uint32 ii=0; ii<intervals.size(); ii++)
        if ((intervals[ii].isGood == false) && (minpos <= intervals[ii].end) && (intervals[ii].bgn <= maxpos)) {
          dest = ii;
          break;
        }

      //  If the fragment touched a bad interval, dest is set, and we move the fragment to that new unitig.
      if (dest < intervals.size()) {
        AppendVA_IntMultiPos(maNew[dest]->f_list, imp);
        imp->ident = 0;
        continue;
      }

      //  Otherwise, search for a destination good unitig.
      for (uint32 ii=0; ii<intervals.size(); ii++)
        if ((intervals[ii].isGood == true) && (minpos <= intervals[ii].end) && (intervals[ii].bgn <= maxpos)) {
          dest = ii;
          break;
        }

      //  If it touched a good interval, move it there.
      if (dest < intervals.size()) {
        AppendVA_IntMultiPos(maNew[dest]->f_list, imp);
        imp->ident = 0;
        continue;
      }

      //  We should never get here.  The original unitig should be covered completely by intervals.
      assert(0);
    }

    //  Run these through consensus, and add to the store.

    for (uint32 i=0; i<intervals.size(); i++) {
      if (GetNumIntMultiPoss(maNew[i]->f_list) == 0)
        continue;

      maNew[i]->maID = tigStore->numUnitigs();

      fprintf(stderr, "Creating new unitig %d with "F_SIZE_T" fragments\n",
              maNew[i]->maID, GetNumIntMultiPoss(maNew[i]->f_list));

      if (MultiAlignUnitig(maNew[i], gkpStore, &options, NULL) == false)
        fprintf(stderr, "  Unitig %d FAILED.\n", maNew[i]->maID);

      //  Add the new unitig to the store (the asserts are duplicated in cgw too).

      tigStore->insertMultiAlign(maNew[i], TRUE, FALSE);

      DeleteMultiAlignT(maNew[i]);
    }

    delete [] maNew;
  }

  delete [] rc;
  delete [] gcc;
  delete [] bcc;

  delete [] matePair;
  delete [] library;

  delete gkpStore;
  delete tigStore;
}
