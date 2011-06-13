
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

const char *mainid = "$Id: overlap_partition.C,v 1.1 2011-06-13 03:18:37 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_fileIO.h"
#include "AS_MSG_pmesg.h"
#include "AS_GKP_include.h"


//  Reads gkpStore, outputs four files:
//    ovlbat - batch names
//    ovljob - job names
//    ovlopt - overlapper options
//    ovlinf - (unused)

uint32  batchMax = 1000;


void
outputJob(FILE   *BAT,
          FILE   *JOB,
          FILE   *OPT,
          uint32  hashBeg,
          uint32  hashEnd,
          uint32  refBeg,
          uint32  refEnd,
          uint32  maxNumFrags,
          uint32  maxLength,
          uint32 &batchSize,
          uint32 &batchName,
          uint32 &jobName) {

  fprintf(BAT, "%03"F_U32P"\n", batchName);
  fprintf(JOB, "%06"F_U32P"\n", jobName);

  if (maxNumFrags == 0)
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd);
  else
    fprintf(OPT, "-h "F_U32"-"F_U32" -r "F_U32"-"F_U32" --hashstrings "F_U32" --hashdatalen "F_U32"\n",
            hashBeg, hashEnd, refBeg, refEnd, maxNumFrags, maxLength);

  refBeg = refEnd + 1;

  batchSize++;

  if (batchSize >= batchMax) {
    batchSize = 0;
    batchName++;
  }

  jobName++;
}



void
partitionFrags(gkStore    *gkp,
               FILE       *BAT,
               FILE       *JOB,
               FILE       *OPT,
               uint32      ovlHashBlockSize,
               uint32      ovlRefBlockSize) {
  uint32  hashBeg = 1;
  uint32  hashEnd = 0;

  uint32  refBeg = 1;
  uint32  refEnd = 0;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32     numFrags = gkp->gkStore_getNumFragments();

  while (hashBeg < numFrags) {
    hashEnd = hashBeg + ovlHashBlockSize - 1;

    if (hashEnd > numFrags)
      hashEnd = numFrags;

    refBeg = 1;
    refEnd = 0;

    while (refBeg < hashEnd) {
      refEnd = refBeg + ovlRefBlockSize - 1;

      if (refEnd > numFrags)
        refEnd = numFrags;

      outputJob(BAT, JOB, OPT, hashBeg, hashEnd, refBeg, refEnd, 0, 0, batchSize, batchName, jobName);

      refBeg = refEnd + 1;
    }

    hashBeg = hashEnd + 1;
  }
}





void
partitionLength(gkStore    *gkp,
                FILE       *BAT,
                FILE       *JOB,
                FILE       *OPT,
                uint32      ovlHashBlockLength,
                uint32      ovlRefBlockSize) {
  uint32  hashBeg = 1;
  uint32  hashEnd = 0;

  uint32  refBeg = 1;
  uint32  refEnd = 0;

  uint32  batchSize = 0;
  uint32  batchName = 1;
  uint32  jobName   = 1;

  uint32     numFrags = gkp->gkStore_getNumFragments();

  gkFragment  fr;
  gkStream   *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  while (hashBeg < numFrags) {
    uint32  len    = 0;

    assert(hashEnd == hashBeg - 1);

    do {
      if (fs->next(&fr) == false)
        break;

      hashEnd = fr.gkFragment_getReadIID();

      //  Non deleted fragments contribute one byte per untrimmed base
      if (fr.gkFragment_getIsDeleted() == false)
        len += fr.gkFragment_getClearRegionLength();

      //  Even deleted fragments contribute one byte (the terminating zero)
      len += 1;

    } while (len < ovlHashBlockLength);

    assert(hashEnd <= numFrags);

    refBeg = 1;
    refEnd = 0;

    while (refBeg < hashEnd) {
      refEnd = refBeg + ovlRefBlockSize - 1;

      if (refEnd > numFrags)
        refEnd = numFrags;

      outputJob(BAT, JOB, OPT, hashBeg, hashEnd, refBeg, refEnd, hashEnd - hashBeg, len, batchSize, batchName, jobName);

      refBeg = refEnd + 1;
    }

    hashBeg = hashEnd + 1;
  }

  delete fs;
}





int
main(int argc, char **argv) {
  char            *gkpStoreName        = NULL;
  gkStore         *gkpStore            = NULL;

  char            *outputPrefix        = NULL;
  char             outputName[FILENAME_MAX];

  uint32           ovlHashBlockLength  = 0;
  uint32           ovlHashBlockSize    = 0;
  uint32           ovlRefBlockSize     = 0;

  int arg = 1;
  int err = 0;

  while (arg < argc) {
    if        (strcmp(argv[arg], "-g") == 0) {
      gkpStoreName = argv[++arg];

    } else if (strcmp(argv[arg], "-bl") == 0) {
      ovlHashBlockLength = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-bs") == 0) {
      ovlHashBlockSize   = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-rs") == 0) {
      ovlRefBlockSize    = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-o") == 0) {
      outputPrefix = argv[++arg];

    } else {
      err++;
    }

    arg++;
  }
  if (err) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    exit(1);
  }

  if ((ovlHashBlockLength > 0) && (ovlHashBlockSize > 0))
    fprintf(stderr, "ERROR:  At most one of -bl and -bs can be non-zero.\n"), exit(1);

  gkStore   *gkp      = new gkStore(gkpStoreName, FALSE, FALSE, true);

  errno = 0;

  sprintf(outputName, "%s/ovlbat", outputPrefix);
  FILE *BAT = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s/ovljob", outputPrefix);
  FILE *JOB = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  sprintf(outputName, "%s/ovlopt", outputPrefix);
  FILE *OPT = fopen(outputName, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s': %s\n", outputName, strerror(errno)), exit(1);

  if (ovlHashBlockLength == 0)
    partitionFrags(gkp, BAT, JOB, OPT, ovlHashBlockSize, ovlRefBlockSize);
  else
    partitionLength(gkp, BAT, JOB, OPT, ovlHashBlockLength, ovlRefBlockSize);

  fclose(BAT);
  fclose(JOB);
  fclose(OPT);

  delete gkp;

  exit(0);
}
