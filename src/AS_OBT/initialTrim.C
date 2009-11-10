
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

const char *mainid = "$Id: initialTrim.C,v 1.25 2009-11-10 05:45:55 brianwalenz Exp $";

//  Read a fragStore, does quality trimming based on quality scores,
//  intersects the quality trim with a vector trim, and updates the
//  original clear range in the store.

#include <stdio.h>
#include <stdlib.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"

#include "util++.H"
#include "trim.H"
#include "constants.H"

int
main(int argc, char **argv) {
  char   *gkpName             = 0L;
  FILE   *logFile             = 0L;

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
    fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
    fprintf(stderr, "    iid originalBegin originalEnd newBegin newEnd\n");
    fprintf(stderr, "    uid,iid origBegin origEnd qualBegin qualEnd vecBeg vecEnd newBegin newEnd\n");
    exit(1);
  }

  gkStore  *gkpStore = new gkStore(gkpName, FALSE, TRUE);

  gkpStore->gkStore_metadataCaching(true);
  gkpStore->gkStore_enableClearRange(AS_READ_CLEAR_OBTINITIAL);

  gkFragment    fr;
  gkLibrary    *lr = NULL;

  uint32        qltL = 0;
  uint32        qltR = 0;
  uint32        vecL = 0;
  uint32        vecR = 0;
  uint32        finL = 0;
  uint32        finR = 0;

  uint32        stat_immutable      = 0;
  uint32        stat_donttrim       = 0;
  uint32        stat_alreadyDeleted = 0;
  uint32        stat_noVecClr       = 0;
  uint32        stat_noHQnonVec     = 0;
  uint32        stat_HQtrim5        = 0;
  uint32        stat_HQtrim3        = 0;
  uint32        stat_LQtrim5        = 0;
  uint32        stat_LQtrim3        = 0;
  uint32        stat_tooShort       = 0;

  if (logFile)
    fprintf(logFile, "uid,iid\torigL\torigR\tqltL\tqltR\tvecL\tvecR\tfinalL\tfinalR\tdeleted?\n");

  for (uint32 iid=1; iid<=gkpStore->gkStore_getNumFragments(); iid++) {
    gkpStore->gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    if (fr.gkFragment_getLibraryIID() != 0) {
      lr = gkpStore->gkStore_getLibrary(fr.gkFragment_getLibraryIID());
    }

    if (fr.gkFragment_getIsDeleted()) {
      stat_alreadyDeleted++;
      continue;
    }

    if ((lr) && (lr->doNotOverlapTrim)) {
      stat_immutable++;
      continue;
    }

    if ((lr) && (lr->doNotQVTrim)) {
      stat_donttrim++;
      qltL = 0;
      qltR = fr.gkFragment_getSequenceLength();
    } else {
      double  minQuality = qual.lookupNumber(12);
      if ((lr) && (lr->goodBadQVThreshold > 0))
        minQuality = qual.lookupNumber(lr->goodBadQVThreshold);

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

    fr.gkFragment_setClearRegion(finL, finR, AS_READ_CLEAR_OBTINITIAL);

    gkpStore->gkStore_setFragment(&fr);

    if ((finL + AS_READ_MIN_LEN) > finR)
      gkpStore->gkStore_delFragment(iid);

    if (logFile)
      fprintf(logFile, "%s,"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"%s\n",
              AS_UID_toString(fr.gkFragment_getReadUID()),
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

  fprintf(stdout, "Fragments with:\n");
  fprintf(stdout, " no changes allowed:           "F_U32"\n", stat_immutable);
  fprintf(stdout, " no QV trim allowed:           "F_U32"\n", stat_donttrim);
  fprintf(stdout, " already deleted               "F_U32"\n", stat_alreadyDeleted);
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
