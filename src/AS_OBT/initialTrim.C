
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

//  Read a fragStore, does quality trimming based on quality scores,
//  intersects the quality trim with a vector trim, and updates the
//  original clear range in the store.

#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <ctype.h>
//#include <math.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
}

#include "util++.H"
#include "trim.H"
#include "constants.H"

int
main(int argc, char **argv) {
  bool    doUpdate            = false;
  char   *gkpName             = 0L;
  FILE   *logFile             = 0L;

  int arg = 1;
  int err = 0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-update", 2) == 0) {
      doUpdate = true;

    } else if (strncmp(argv[arg], "-log", 2) == 0) {
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
    fprintf(stderr, "  -update       Update the clear range in the fragStore.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -log X        Report the iid, original trim and new quality trim\n");
    fprintf(stderr, "  -frg F        Operate on this gkpStore\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
    fprintf(stderr, "    iid originalBegin originalEnd newBegin newEnd\n");
    fprintf(stderr, "    uid,iid origBegin origEnd qualBegin qualEnd vecBeg vecEnd newBegin newEnd\n");
    exit(1);
  }

  GateKeeperStore  *gkpStore = openGateKeeperStore(gkpName, doUpdate);
  if (gkpStore == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpName);
    exit(1);
  }

  gkpStore->frg = convertStoreToMemoryStore(gkpStore->frg);

  uint32        firstElem = getFirstElemFragStore(gkpStore);
  uint32        lastElem  = getLastElemFragStore(gkpStore) + 1;

  fragRecord    fr = {0};

  uint32        qltL = 0;
  uint32        qltR = 0;
  uint32        finL = 0;
  uint32        finR = 0;

  uint32        stat_immutable      = 0;
  uint32        stat_donttrim       = 0;
  uint32        stat_noVecClr       = 0;
  uint32        stat_noHQnonVec     = 0;
  uint32        stat_HQtrim5        = 0;
  uint32        stat_HQtrim3        = 0;
  uint32        stat_LQtrim5        = 0;
  uint32        stat_LQtrim3        = 0;
  uint32        stat_tooShort       = 0;

  if (logFile)
    fprintf(logFile, "uid,iid\torigL\torigR\tqltL\tqltR\tvecL\tvecR\tfinalL\tfinalR\tdeleted?\n");

  for (uint32 iid=firstElem; iid<lastElem; iid++) {

    getFrag(gkpStore, iid, &fr, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);

    GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkpStore, getFragRecordLibraryIID(&fr));

    //  Bail now if we've been told to not modify this read.  We do
    //  not print a message in the log.
    //
    if        ((gklr) && (gklr->doNotOverlapTrim)) {
      stat_immutable++;
    } else {

      if ((gklr) && (gklr->doNotQVTrim)) {
        stat_donttrim++;
        qltL = 0;
        qltR = getFragRecordSequenceLength(&fr);
      } else {
        double  minQuality = qual.lookupNumber(12);
        if ((gklr) && (gklr->goodBadQVThreshold > 0))
          minQuality = qual.lookupNumber(gklr->goodBadQVThreshold);

        doTrim(&fr, minQuality, qltL, qltR);
      }

      finL = qltL;
      finR = qltR;

      //  Intersect with the vector clear range, if it exists

      if (fr.gkfr.hasVectorClear == false) {
        //  no vector clear known, use the quality clear
        stat_noVecClr++;

      } else if ((fr.gkfr.clearBeg[AS_READ_CLEAR_VEC] > finR) || (fr.gkfr.clearEnd[AS_READ_CLEAR_VEC] < finL)) {
        //  don't intersect; trust nobody
        stat_noHQnonVec++;
        finL = 0;
        finR = 0;

      } else {
        //  They intersect.  Pick the largest begin and the smallest end

        if (finL < fr.gkfr.clearBeg[AS_READ_CLEAR_VEC]) {
          stat_HQtrim5++;
          finL = fr.gkfr.clearBeg[AS_READ_CLEAR_VEC];
        } else {
          stat_LQtrim5++;
        }
        if (fr.gkfr.clearEnd[AS_READ_CLEAR_VEC] < finR) {
          stat_HQtrim3++;
          finR = fr.gkfr.clearEnd[AS_READ_CLEAR_VEC];
        } else {
          stat_LQtrim3++;
        }
      }

      //  Update the clear ranges

      if ((finL + AS_FRAG_MIN_LEN) > finR)
        stat_tooShort++;

      if (doUpdate) {
        setFragRecordClearRegion(&fr, finL, finR, AS_READ_CLEAR_OBTINI);
        setFrag(gkpStore, iid, &fr);

        if ((finL + AS_FRAG_MIN_LEN) > finR)
          delFrag(gkpStore, iid);
      }
    }

    if (logFile)
      fprintf(logFile, "%s,"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"%s\n",
              AS_UID_toString(getFragRecordUID(&fr)),
              iid,
              getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_ORIG),
              getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_ORIG),
              qltL,
              qltR,
              fr.gkfr.clearBeg[AS_READ_CLEAR_VEC],
              fr.gkfr.clearEnd[AS_READ_CLEAR_VEC],
              finL,
              finR,
              ((finL + AS_FRAG_MIN_LEN) > finR) ? " (deleted)" : "");
  }

  closeGateKeeperStore(gkpStore);

  fprintf(stderr, "Fragments with:\n");
  fprintf(stderr, " no changes allowed:           "F_U32"\n", stat_immutable);
  fprintf(stderr, " no QV trim allowed:           "F_U32"\n", stat_donttrim);
  fprintf(stderr, " no vector clear range known:  "F_U32" (trimed to quality clear)\n", stat_noVecClr);
  fprintf(stderr, " no HQ non-vector sequence:    "F_U32" (deleted)\n", stat_noHQnonVec);
  fprintf(stderr, " HQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_HQtrim5 + stat_HQtrim3);
  fprintf(stderr, "     5' end:                   "F_U32"\n", stat_HQtrim5);
  fprintf(stderr, "     3' end:                   "F_U32"\n", stat_HQtrim3);
  fprintf(stderr, " LQ vector trimmed:            "F_U32" (trimmed to intersection)\n", stat_LQtrim5 + stat_LQtrim3);
  fprintf(stderr, "     5' end:                   "F_U32"\n", stat_LQtrim5);
  fprintf(stderr, "     3' end:                   "F_U32"\n", stat_LQtrim3);
  fprintf(stderr, " final clear range too short:  "F_U32" (deleted)\n", stat_tooShort);
}
