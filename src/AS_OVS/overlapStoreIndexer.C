
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

const char *mainid = "$Id: overlapStoreIndexer.C,v 1.7 2012-06-02 08:36:24 brianwalenz Exp $";

#include "AS_global.h"

#include "AS_PER_gkpStore.h"

#include "overlapStore.h"

#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"
#include "AS_OBT_acceptableOverlap.h"

#include <ctype.h>
#include <unistd.h>  //  sysconf()

#include <vector>

using namespace std;


#define AS_OVS_CURRENT_VERSION  2


bool
testIndex(char *ovlName, bool doFixes) {
  char name[FILENAME_MAX];
  FILE *I = NULL;
  FILE *F = NULL;

  sprintf(name, "%s/idx", ovlName);

  errno = 0;
  I = fopen(name, "r");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s' for reading: %s\n", name, strerror(errno)), exit(1);

  fprintf(stderr, "TESTING '%s'\n", name);

  if (doFixes) {
    sprintf(name, "%s/idx.fixed", ovlName);

    errno = 0;
    F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s' for writing: %s\n", name, strerror(errno)), exit(1);

    fprintf(stderr, "WITH FIXES TO '%s'\n", name);
  }

  OverlapStoreOffsetRecord  O;

  uint32  curIID = 0;
  uint32  minIID = UINT32_MAX;
  uint32  maxIID = 0;

  uint32  nErrs = 0;

  while (1 == AS_UTL_safeRead(I, &O, "offset", sizeof(OverlapStoreOffsetRecord), 1)) {
    bool  maxIncreases   = (maxIID < O.a_iid);
    bool  errorDecreased = ((O.a_iid < curIID));
    bool  errorGap       = ((O.a_iid > 0) && (curIID + 1 != O.a_iid));

    if (O.a_iid < minIID)
      minIID = O.a_iid;

    if (maxIncreases)
      maxIID = O.a_iid;

    if (errorDecreased)
      fprintf(stderr, "ERROR: index decreased from "F_U32" to "F_U32"\n", curIID, O.a_iid), nErrs++;
    else if (errorGap)
      fprintf(stderr, "ERROR: gap between "F_U32" and "F_U32"\n", curIID, O.a_iid), nErrs++;

    if ((maxIncreases == true) && (errorGap == false)) {
      if (doFixes)
        AS_UTL_safeWrite(F, &O, "offset", sizeof(OverlapStoreOffsetRecord), 1);

    } else if (O.numOlaps > 0) {
      fprintf(stderr, "ERROR: lost overlaps a_iid "F_U32" fileno "F_U32" offset "F_U32" numOlaps "F_U32"\n",
              O.a_iid, O.fileno, O.offset, O.numOlaps);
    }

    curIID = O.a_iid;
  }

  fprintf(stderr, "minIID "F_U32"\n", minIID);
  fprintf(stderr, "maxIID "F_U32"\n", maxIID);

  fclose(I);

  if (F)
    fclose(F);

  return(nErrs == 0);
}


int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  uint32          maxJob       = 0;
  uint32          cntJob       = 0;

  bool            deleteIntermediates = false;

  bool            doExplicitTest = false;
  bool            doFixes        = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      maxJob = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-f") == 0) {
      doFixes = true;

    } else if (strcmp(argv[arg], "-t") == 0) {
      doExplicitTest = true;
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-delete") == 0) {
      deleteIntermediates = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if ((maxJob == 0) && (doExplicitTest == false))
    err++;
  if (err) {
    fprintf(stderr, "usage: %s ...\n", argv[0]);
    fprintf(stderr, "  -o x.ovlStore    path to overlap store to build the final index for\n");
    fprintf(stderr, "  -F s             number of slices used in bucketizing/sorting\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -t x.ovlStore    explicitly test a previously constructed index\n");
    fprintf(stderr, "  -f               when testing, also create a new 'idx.fixed' which might resolve\n");
    fprintf(stderr, "                   rare problems\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -delete          remove intermediate files when the index is successfully created\n");
    exit(1);
  }

  if (doExplicitTest == true) {
    bool passed = testIndex(ovlName, doFixes);

    exit((passed == true) ? 0 : 1);
  }

  //  Check that all segments are present.  Every segment should have an ovs file.

  for (uint32 i=1; i<=maxJob; i++) {
    uint32  complete = 0;

    sprintf(name, "%s/%04d", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" data not present  (%s)\n", i, name);

    sprintf(name, "%s/%04d.ovs", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" store not present (%s)\n", i, name);

    sprintf(name, "%s/%04d.idx", ovlName, i);
    if (AS_UTL_fileExists(name, FALSE, FALSE) == true)
      complete++;
    else
      fprintf(stderr, "ERROR: Segment "F_U32" index not present (%s)\n", i, name);

    if (complete == 3)
      cntJob++;
  }

  if (cntJob != maxJob) {
    fprintf(stderr, "ERROR: Expected "F_U32" segments, only found "F_U32".\n", maxJob, cntJob);
    exit(1);
  }


  OverlapStoreInfo    ovspiece;
  OverlapStoreInfo    ovs;

	ovs.ovsMagic              = 1;
	ovs.ovsVersion            = AS_OVS_CURRENT_VERSION;
  ovs.numOverlapsPerFile    = 1024 * 1024 * 1024 / sizeof(OVSoverlapINT);
  ovs.smallestIID           = UINT_MAX;
  ovs.largestIID            = 0;
  ovs.numOverlapsTotal      = 0;
  ovs.highestFileIndex      = 0;
	ovs.maxReadLenInBits      = AS_READ_MAX_NORMAL_LEN_BITS;

  OverlapStoreOffsetRecord missing;

  missing.a_iid     = 0;
  missing.fileno    = 1;
  missing.offset    = 0;
  missing.numOlaps  = 0;

  //  Open the new master index output file

  sprintf(name, "%s/idx", ovlName);

  errno = 0;
  FILE  *idx = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  ovs.highestFileIndex = maxJob;

  //  Special case, we need an empty index for the zeroth fragment.

  AS_UTL_safeWrite(idx, &missing, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);

  //  Process each

  for (uint32 i=1; i<=maxJob; i++) {
    sprintf(name, "%s/%04d.ovs", ovlName, i);

    fprintf(stderr, "Processing '%s'\n", name);

    if (AS_UTL_fileExists(name, FALSE, FALSE) == false) {
      fprintf(stderr, "ERROR: file '%s' not found.\n", name);
      exit(1);
    }

    {
      errno = 0;
      FILE *F = fopen(name, "r");
      if (errno)
        fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
      AS_UTL_safeRead(F, &ovspiece, "ovspiece", sizeof(OverlapStoreInfo), 1);
      fclose(F);
    }

    //  Add empty index elements for missing overlaps

    if (ovspiece.numOverlapsTotal == 0) {
      fprintf(stderr, "  No overlaps found.\n");
      continue;
    }

    assert(ovspiece.smallestIID <= ovspiece.largestIID);

    if (ovs.largestIID + 1 < ovspiece.smallestIID)
      fprintf(stderr, "  Adding empty records for fragments "F_U64" to "F_U64"\n",
              ovs.largestIID + 1, ovspiece.smallestIID - 1);

    while (ovs.largestIID + 1 < ovspiece.smallestIID) {
      missing.a_iid     = ovs.largestIID + 1;
      //missing.fileno    = set elsewhere
      //missing.offset    = set elsewhere
      //missing.numOlaps  = 0;

      AS_UTL_safeWrite(idx, &missing, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);

      ovs.largestIID++;
    }

    //  Copy index elements for existing overlaps.  While copying, update the supposed position
    //  of any fragments with no overlaps.  Without doing this, accessing the store beginning
    //  or ending at such a fragment will fail.

    {
      sprintf(name, "%s/%04d.idx", ovlName, i);

      errno = 0;
      FILE  *F = fopen(name, "r");
      if (errno)
        fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

      uint32                     recsLen = 0;
      uint32                     recsMax = 1024 * 1024;
      OverlapStoreOffsetRecord  *recs    = new OverlapStoreOffsetRecord [recsMax];

      recsLen = AS_UTL_safeRead(F, recs, "recs", sizeof(OverlapStoreOffsetRecord), recsMax);

      if (recsLen > 0) {
        if (ovs.largestIID + 1 != recs[0].a_iid)
          fprintf(stderr, "ERROR: '%s' starts with iid "F_U32", but store only up to "F_U64"\n",
                  name, recs[0].a_iid, ovs.largestIID);
        assert(ovs.largestIID + 1 == recs[0].a_iid);
      }

      while (recsLen > 0) {
        missing.fileno = recs[recsLen-1].fileno;  //  Update location of missing stuff.
        missing.offset = recs[recsLen-1].offset;

				AS_UTL_safeWrite(idx, recs, "recs", sizeof(OverlapStoreOffsetRecord), recsLen);
        recsLen = AS_UTL_safeRead(F, recs, "recs", sizeof(OverlapStoreOffsetRecord), recsMax);
      }

      delete [] recs;

      fclose(F);
    }

    //  Update

    ovs.smallestIID = MIN(ovs.smallestIID, ovspiece.smallestIID);
    ovs.largestIID  = MAX(ovs.largestIID,  ovspiece.largestIID);

    ovs.numOverlapsTotal += ovspiece.numOverlapsTotal;

    fprintf(stderr, "  Now finished with fragments "F_U64" to "F_U64" -- "F_U64" overlaps.\n",
            ovs.smallestIID, ovs.largestIID, ovs.numOverlapsTotal);
  }

  fclose(idx);


  //  Dump the new store info file

  {
    sprintf(name, "%s/ovs", ovlName);

    errno = 0;
    FILE  *F = fopen(name, "w");
    if (errno)
      fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &ovs, "ovs", sizeof(OverlapStoreInfo), 1);

    fclose(F);
  }

  //  Diagnostics.

  fprintf(stderr, "new store from "F_U64" to "F_U64" with "F_U64" overlaps.\n",
          ovs.smallestIID,
          ovs.largestIID,
          ovs.numOverlapsTotal);

  if (testIndex(ovlName, false) == false) {
    fprintf(stderr, "ERROR: index failed tests.\n");
    exit(1);
  }

  //  Remove intermediates.  For the buckets, we keep going until there are 10 in a row not present.
  //  During testing, on a microbe using 2850 buckets, some buckets were empty.

  if (deleteIntermediates) {
  for (uint32 i=1; i<=maxJob; i++) {
    char name[FILENAME_MAX];

    sprintf(name, "%s/%04u.idx", ovlName, i);
    fprintf(stderr, "unlink %s\n", name);
    AS_UTL_unlink(name);

    sprintf(name, "%s/%04u.ovs", ovlName, i);
    fprintf(stderr, "unlink %s\n", name);
    AS_UTL_unlink(name);
  }

  for (uint32 missing=0, i=1; missing<10; i++) {
    char name[FILENAME_MAX];

    sprintf(name, "%s/bucket%04d", ovlName, i);
    if (AS_UTL_fileExists(name, TRUE, FALSE) == FALSE) {
      missing++;
      continue;
    }

    missing = 0;

    sprintf(name, "%s/bucket%04d/sliceSizes", ovlName, i);
    fprintf(stderr, "unlink %s\n", name);
    AS_UTL_unlink(name);

    sprintf(name, "%s/bucket%04d", ovlName, i);
    fprintf(stderr, "rmdir %s\n", name);
    rmdir(name);
  }
  }

  exit(0);
}

