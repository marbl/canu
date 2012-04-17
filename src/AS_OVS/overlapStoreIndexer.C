
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

const char *mainid = "$Id: overlapStoreIndexer.C,v 1.3 2012-04-17 04:04:51 brianwalenz Exp $";

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

#define DELETE_INTERMEDIATE

int
main(int argc, char **argv) {
  char           *ovlName      = NULL;
  uint32          maxJob       = 0;
  uint32          cntJob       = 0;

  bool            deleteIntermediates = false;

  char            name[FILENAME_MAX];

  argc = AS_configure(argc, argv);

  int err=0;
  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-o") == 0) {
      ovlName = argv[++arg];

    } else if (strcmp(argv[arg], "-F") == 0) {
      maxJob = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-delete") == 0) {
      deleteIntermediates = true;

    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
    }

    arg++;
  }
  if (ovlName == NULL)
    err++;
  if (maxJob == 0)
    err++;
  if (err) {
    exit(1);
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

  missing.offset    = 0;
  missing.numOlaps  = 0;

  //  Open the new master index output file

  sprintf(name, "%s/idx", ovlName);

  errno = 0;
  FILE  *idx = fopen(name, "w");
  if (errno)
    fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  ovs.highestFileIndex = maxJob;

  //  Process each

  for (uint32 i=1; i<=maxJob; i++) {
    sprintf(name, "%s/%04d.ovs", ovlName, i);

    fprintf(stderr, "Processing '%s'\n", name);

    //if (AS_UTL_fileExists(name, FALSE, FALSE) == false)
    //  error;

    {
      errno = 0;
      FILE *F = fopen(name, "r");
      if (errno)
        fprintf(stderr, "ERROR: Failed to open '%s': %s\n", name, strerror(errno)), exit(1);
      AS_UTL_safeRead(F, &ovspiece, "ovspiece", sizeof(OverlapStoreInfo), 1);
      fclose(F);
    }

    //  Add empty index elements for missing overlaps

    if (ovs.largestIID > 0) {
      fprintf(stderr, "  Adding empty records for fragments "F_U64" to "F_U64"\n",
              ovs.largestIID + 1, ovspiece.smallestIID);

      while (ovs.largestIID + 1 < ovspiece.smallestIID) {
        missing.a_iid     = ovs.largestIID + 1;
				missing.fileno    = 0;  //  POSSIBLY WRONG
				missing.offset    = 0;  //  POSSIBLY WRONG
				missing.numOlaps  = 0;

				AS_UTL_safeWrite(idx, &missing, "AS_OVS_writeOverlapToStore offset", sizeof(OverlapStoreOffsetRecord), 1);

				ovs.largestIID++;
      }
    }

    //  Copy index elements for existing overlaps

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

      while (recsLen > 0) {
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

