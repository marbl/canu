
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
}

#include "util++.H"
#include "trim.H"
#include "constants.H"

//  Read a fragStore, does quality trimming based on quality scores,
//  modifies the original clear range in the store.
//
//  Optionally intersects the quality trim with a vector trim.
//
//  Optionally does NOT modify a list of fragment UIDs

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-q quality] [-update] [-replace] [-log logfile] -frg some.gkpStore\n", name);
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
}



int
main(int argc, char **argv) {
  double  minQuality          = qual.lookupNumber(20);
  bool    doUpdate            = false;
  char   *gkpStore            = 0L;
  FILE   *logFile             = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-quality", 2) == 0) {
      minQuality = qual.lookupNumber(atoi(argv[++arg]));

    } else if (strncmp(argv[arg], "-update", 2) == 0) {
      doUpdate = true;

    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);
    } else if (strncmp(argv[arg], "-frg", 2) == 0) {
      gkpStore = argv[++arg];

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (!gkpStore) {
    usage(argv[0]);
    exit(1);
  }

  GateKeeperStore  *gkp = openGateKeeperStore(gkpStore, doUpdate);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStore);
    exit(1);
  }

  uint32        firstElem = getFirstElemFragStore(gkp);
  uint32        lastElem  = getLastElemFragStore(gkp) + 1;

  fragRecord   *fr = new_fragRecord();

  uint32        qltL = 0;
  uint32        qltR = 0;
  uint32        vecL = 0;
  uint32        vecR = 0;

  uint32        stat_notPresent  = 0;
  uint32        stat_noIntersect = 0;
  uint32        stat_change      = 0;
  uint32        stat_noChange    = 0;
  uint32        stat_immutable   = 0;

  if (logFile)
    fprintf(logFile, "uid,iid\torigL\torigR\tqltL\tqltR\tvecL\tvecR\tfinalL\tfinalR\tdeleted?\n");

  for (uint32 iid=firstElem; iid<lastElem; iid++) {

    getFrag(gkp, iid, fr, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);

    GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkp, getFragRecordLibraryIID(fr));

    //  Bail now if we've been told to not modify this read.  We do
    //  not print a message in the log.
    //
    if ((gklr) && (gklr->doNotOverlapTrim)) {
      stat_immutable++;
      continue;
    }

    doTrim(fr, minQuality, qltL, qltR);

    vecL = qltL;
    vecR = qltR;

    //  Intersect with the vector clear range, if it exists
    //
    if (fr->gkfr.hasVectorClear == false) {
      //  iid not present in our input list, do nothing.
      stat_notPresent++;
    } else if ((fr->gkfr.clearBeg[AS_READ_CLEAR_VEC] > vecR) || (fr->gkfr.clearEnd[AS_READ_CLEAR_VEC] < vecL)) {
      //  don't intersect; trust the quality clear
      stat_noIntersect++;
    } else {
      //  They intersect.  Pick the largest begin and the smallest end

      bool changed = false;

      if (vecL < fr->gkfr.clearBeg[AS_READ_CLEAR_VEC]) {
        changed = true;
        vecL = fr->gkfr.clearBeg[AS_READ_CLEAR_VEC];
      }
      if (fr->gkfr.clearEnd[AS_READ_CLEAR_VEC] < vecR) {
        changed = true;
        vecR = fr->gkfr.clearEnd[AS_READ_CLEAR_VEC];
      }

      if (changed)
        stat_change++;
      else
        stat_noChange++;
    }


    if (logFile) {
      fprintf(logFile, F_U64","F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"%s\n",
              getFragRecordUID(fr),
              iid,
              getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_ORIG),
              getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_ORIG),
              qltL,
              qltR,
              fr->gkfr.clearBeg[AS_READ_CLEAR_VEC],
              fr->gkfr.clearEnd[AS_READ_CLEAR_VEC],
              vecL,
              vecR,
              ((vecL + AS_FRAG_MIN_LEN) > vecR) ? " (deleted)" : "");
    }

    if (doUpdate) {
      setFragRecordClearRegion(fr, vecL, vecR, AS_READ_CLEAR_OBTINI);
      setFrag(gkp, iid, fr);

      if ((vecL + AS_FRAG_MIN_LEN) > vecR)
        delFrag(gkp, iid);
    }
  }

  closeGateKeeperStore(gkp);

  fprintf(stderr, "Fragments with no vector clear:  "F_U32"\n", stat_notPresent);
  fprintf(stderr, "Fragments with no intersection:  "F_U32"\n", stat_noIntersect);
  fprintf(stderr, "Fragments with vector trimmed:   "F_U32"\n", stat_change);
  fprintf(stderr, "Fragments low quality vector:    "F_U32"\n", stat_noChange);
  fprintf(stderr, "Fragments marked immutable:      "F_U32"\n", stat_immutable);
}
