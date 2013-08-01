
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
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

const char *mainid = "$Id: initialTrim.C,v 1.31 2012-03-05 09:09:17 brianwalenz Exp $";

//  Read a fragStore, does quality trimming based on quality scores,
//  intersects the quality trim with a vector trim, and updates the
//  original clear range in the store.

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.H"
#include "AS_PER_gkpStore.H"

#include "trim.H"

int
main(int argc, char **argv) {
  char   *gkpName             = 0L;
  FILE   *logFile             = 0L;
  bool    beVerbose           = false;
  bool    doUpdate            = true;

  argc = AS_configure(argc, argv);

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-log", 2) == 0) {
      arg++;
      logFile = stderr;
      if (strcmp(argv[arg], "-") != 0) {
        errno=0;
        logFile = fopen(argv[arg], "w");
        if (errno)
          fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);
      }

    } else if (strncmp(argv[arg], "-frg", 2) == 0) {
      gkpName = argv[++arg];

    } else if (strncmp(argv[arg], "-v", 2) == 0) {
      beVerbose = true;

    } else if (strncmp(argv[arg], "-n", 2) == 0) {
      doUpdate = false;

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (!gkpName)) {
    fprintf(stderr, "usage: %s [-q quality] [-update] [-replace] [-log logfile] -frg some.gkpStore\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  -q quality    Find quality trim points using 'quality' as the base.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -log X        Report the iid, original trim and new quality trim\n");
    fprintf(stderr, "  -frg F        Operate on this gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -v            Be uselessly verbose (for debugging)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
    fprintf(stderr, "    iid originalBegin originalEnd newBegin newEnd\n");
    fprintf(stderr, "    iid origBegin origEnd qualBegin qualEnd vecBeg vecEnd newBegin newEnd\n");
    exit(1);
  }

  gkStore  *gkpStore = new gkStore(gkpName, FALSE, TRUE);

  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTINITIAL);
  gkpStore->gkStore_metadataCaching(true);

  gkFragment    fr;
  gkLibrary    *lr = NULL;

  uint32        qltL = 0;
  uint32        qltR = 0;
  uint32        vecL = 0;
  uint32        vecR = 0;
  uint32        finL = 0;
  uint32        finR = 0;

  //  One per fragment
  uint32        stat_alreadyDeleted = 0;
  uint32        stat_noTrimming     = 0;
  uint32        stat_merBased       = 0;
  uint32        stat_flowBased      = 0;
  uint32        stat_qualityBased   = 0;
  uint32        stat_NOLIBRARY      = 0;  //  DOUBLE COUNTS

  //  One per fragment
  uint32        stat_noVecClr       = 0;
  uint32        stat_noHQnonVec     = 0;
  uint32        stat_HQtrim5        = 0;
  uint32        stat_HQtrim3        = 0;
  uint32        stat_LQtrim5        = 0;
  uint32        stat_LQtrim3        = 0;
  uint32        stat_tooShort       = 0;

  //  Used to be a library parameter.
  double        minQuality = qual.lookupNumber(12);

  if (logFile)
    fprintf(logFile, "iid\torigL\torigR\tqltL\tqltR\tfinalL\tfinalR\tvecL\tvecR\tdeleted?\n");

  for (AS_IID iid=1; iid<=gkpStore->gkStore_getNumFragments(); iid++) {
    if (beVerbose)
      fprintf(stderr, "Loading fragment "F_IID".\n", iid);
    gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    lr = gkpStore->gkStore_getLibrary(fr.gkFragment_getLibraryIID());

    if (fr.gkFragment_getIsDeleted()) {
      stat_alreadyDeleted++;
      continue;
    }

    if ((lr) && (lr->doTrim_initialNone == true)) {
      stat_noTrimming++;
      continue;
    }

    if ((lr) && (lr->doTrim_initialMerBased == true)) {
      stat_merBased++;
      continue;
    }

    if ((lr) && (lr->doTrim_initialFlowBased == true)) {
      stat_flowBased++;
      qltL = 0;
      qltR = fr.gkFragment_getSequenceLength();

    } else if ((lr) && (lr->doTrim_initialQualityBased)) {
      stat_qualityBased++;
      doTrim(&fr, minQuality, qltL, qltR);

    } else {
      //  We should be aborting here, read not in a library!
#warning GET RID OF THESE READS NOT IN A LIBRARY
      stat_qualityBased++;
      stat_NOLIBRARY++;
      doTrim(&fr, minQuality, qltL, qltR);
    }

    finL = qltL;
    finR = qltR;

    //  Intersect with the vector clear range.
    fr.gkFragment_getClearRegion(vecL, vecR, AS_READ_CLEAR_VEC);

    if (vecL > vecR) {
      //  No vector clear defined.
      stat_noVecClr++;
      finL = qltL;
      finR = qltR;

    } else if ((vecL > finR) || (vecR < finL)) {
      //  No intersection; trust nobody.
      stat_noHQnonVec++;
      finL = 0;
      finR = 0;

    } else {
      //  They intersect.  Pick the largest begin and the smallest end

      if (finL < vecL) {
        stat_HQtrim5++;
        finL = vecL;
      } else {
        stat_LQtrim5++;
      }
      if (vecR < finR) {
        stat_HQtrim3++;
        finR = vecR;
      } else {
        stat_LQtrim3++;
      }
    }

    //  Update the clear ranges

    if ((finL + AS_READ_MIN_LEN) > finR)
      stat_tooShort++;

    if (beVerbose)
      fprintf(stderr, "Updating fragment "F_IID".\n", iid);
    if (doUpdate) {
      fr.gkFragment_setClearRegion(finL, finR, AS_READ_CLEAR_OBTINITIAL);
      gkpStore->gkStore_setFragment(&fr);
    }

    if ((finL + AS_READ_MIN_LEN) > finR) {
      if (beVerbose)
        fprintf(stderr, "Deleting fragment "F_IID" (mate "F_IID").\n",
                iid, fr.gkFragment_getMateIID());
      if (doUpdate)
        gkpStore->gkStore_delFragment(iid);
    }

    if (logFile)
      fprintf(logFile, F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"%s\n",
              iid,
              fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_CLR),
              fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_CLR),
              qltL,
              qltR,
              fr.gkFragment_getClearRegionBegin(AS_READ_CLEAR_VEC),
              fr.gkFragment_getClearRegionEnd  (AS_READ_CLEAR_VEC),
              finL,
              finR,
              ((finL + AS_READ_MIN_LEN) > finR) ? " (deleted)" : "");
  }

  delete gkpStore;

  fprintf(stdout, "Fragments trimmed using:\n");
  fprintf(stdout, "  alreadyDeleted   "F_U32"\n", stat_alreadyDeleted);
  fprintf(stdout, "  noTrimming       "F_U32"\n", stat_noTrimming);
  fprintf(stdout, "  merBased         "F_U32"\n", stat_merBased);
  fprintf(stdout, "  flowBased        "F_U32"\n", stat_flowBased);
  fprintf(stdout, "  qualityBased     "F_U32"\n", stat_qualityBased);
  
  fprintf(stdout, "Fragments trimmed using:\n");
  fprintf(stdout, "  NOLIBRARY        "F_U32"\n", stat_NOLIBRARY);

  fprintf(stdout, "Trimming result:\n");
  fprintf(stdout, " no vector clear range known:  "F_U32" (trimed to quality clear)\n", stat_noVecClr);
  fprintf(stdout, " no HQ non-vector sequence:    "F_U32" (deleted)\n", stat_noHQnonVec);
  fprintf(stdout, " HQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_HQtrim5 + stat_HQtrim3);
  fprintf(stdout, "     5' end:                   "F_U32"\n", stat_HQtrim5);
  fprintf(stdout, "     3' end:                   "F_U32"\n", stat_HQtrim3);
  fprintf(stdout, " LQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_LQtrim5 + stat_LQtrim3);
  fprintf(stdout, "     5' end:                   "F_U32"\n", stat_LQtrim5);
  fprintf(stdout, "     3' end:                   "F_U32"\n", stat_LQtrim3);
  fprintf(stdout, " final clear range too short:  "F_U32" (deleted)\n", stat_tooShort);

  return(0);
}
