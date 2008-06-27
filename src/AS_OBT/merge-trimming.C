
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

#include "trim.H"
#include "constants.H"

//  Reads the output of sort-overlap-trim, does the actual trim-point
//  decision, updates the frgStore.

uint32         lineMax = 128 * 1024;
char          *line    = 0L;

void
readLine(FILE *F) {
  fgets(line, lineMax, F);
  line[lineMax-1] = 0;
  assert(strlen(line) < (lineMax-1));
}

class mode5 {
public:
  mode5() {
    memset(_histo, 0, sizeof(uint32) * AS_READ_MAX_LEN);
    _mode5 = 999999999;
  };
  ~mode5() {
  };

  void  add(uint32 x) {
    if (x < 2048)
      _histo[x]++;
  };

  void   compute(void) {
    _mode5 = 1;
    for (uint32 i=1; i<2048; i++)
      if (_histo[i] > _histo[_mode5])
        _mode5 = i;
  };

  uint32 get(void) {
    return(_mode5);
  };

private:
  uint32   _histo[AS_READ_MAX_LEN];
  uint32   _mode5;
};

mode5 *
findModeOfFivePrimeMode(GateKeeperStore *gkp, char *ovlFile) {
  mode5         *modes = new mode5 [getNumGateKeeperLibraries(gkp) + 1];
  fragRecord     fr;

  errno = 0;
  FILE *O = fopen(ovlFile, "r");
  if (errno)
    fprintf(stderr, "Can't open overlap-trim file %s: %s\n", ovlFile, strerror(errno)), exit(1);

  readLine(O);
  while (!feof(O)) {
    splitToWords W(line);

    AS_IID  id = atoi(W[0]);
    getFrag(gkp, id, &fr, FRAG_S_INF);

    AS_IID  lb = getFragRecordLibraryIID(&fr);

    modes[lb].add(atoi(W[4]) + getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBTINI));

    readLine(O);
  }

  fclose(O);

  for (uint32 i=0; i<=getNumGateKeeperLibraries(gkp); i++)
    modes[i].compute();

  return(modes);
}



int
main(int argc, char **argv) {
  uint32   stats[32]         = {0};
  FILE    *O                 = 0L;
  FILE    *logFile           = 0L;
  FILE    *staFile           = 0L;
  char    *frgStore          = 0L;
  char    *ovlFile           = 0L;
  bool     doModify          = true;  //  Make this false for testing

  line = new char [lineMax];

  if (argc < 5) {
    fprintf(stderr, "usage: %s [-log log] -frg frgStore -ovl overlap-consolidated\n", argv[0]);
    fprintf(stderr, "  -ovl o                Read consolidated overlaps from here.\n");
    fprintf(stderr, "  -log x                Write a record of changes to 'x', summary statistics to 'x.stats'\n");
    fprintf(stderr, "  -frg f                'f' is our frag store\n");
    exit(1);
  }

  int arg=1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
    } else if (strncmp(argv[arg], "-ovl", 2) == 0) {
      ovlFile = argv[++arg];
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);

      char *staName = new char [strlen(argv[arg]) + 16];
      sprintf(staName, "%s.stats", argv[arg]);

      staFile = fopen(staName, "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the statistics: %s\n", staName, strerror(errno)), exit(1);

      delete [] staName;
    }
    arg++;
  }


  //  Open the frgStore, prepare for reading fragments
  //
  GateKeeperStore *gkp = openGateKeeperStore(frgStore, doModify);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStore);
    exit(1);
  }

  gkp->frg = convertStoreToMemoryStore(gkp->frg);

  uint32      firstElem = getFirstElemFragStore(gkp);
  uint32      lastElem  = getLastElemFragStore(gkp) + 1;
  fragRecord  fr;


  //  Open the overlap file
  //
  errno = 0;
  O = fopen(ovlFile, "r");
  if (errno)
    fprintf(stderr, "Can't open overlap-trim file %s: %s\n", ovlFile, strerror(errno)), exit(1);


  double  minQuality = qual.lookupNumber(20);


  //  Find the mode of the 5'mode.  Read the whole overlap-trim file,
  //  counting the 5'mode
  //
  mode5 *modes = findModeOfFivePrimeMode(gkp, ovlFile);


  //  Stream through the overlap-trim file, applying some rules to
  //  decide on the correct trim points

  uint64  lid = 0;
  uint64  iid = 0;

  readLine(O);
  while (!feof(O)) {
    splitToWords  W(line);
    iid    = atoi(W[0]);


    //  Report the frags that had no overlap -- we update the clear region
    //  to the Q20 clear region.  If the two clear regions do not intersect,
    //  we delete the fragment.
    //
    lid++;
    while (lid < iid) {
      getFrag(gkp, lid, &fr, FRAG_S_INF | FRAG_S_QLT);

      uint32 qltL0 = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBTINI);
      uint32 qltR0 = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBTINI);
      uint32 qltL1 = 0;
      uint32 qltR1 = 0;
      AS_UID uid   = getFragRecordUID(&fr);

      GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkp, getFragRecordLibraryIID(&fr));

      //  If not already deleted, update the clear.  Updating the
      //  clear on deleted fragments usually results in the log
      //  showing merge-trimming deleted the fragment, which is
      //  incorrect.

      if (getFragRecordIsDeleted(&fr) == 0) {
        if ((gklr) && (gklr->doNotQVTrim)) {
          //  Just leave it as is.  It should show up as a singleton.
          //qltL1 = 0;
          //qltR1 = 0;
        } else {
          doTrim(&fr, minQuality, qltL1, qltR1);
        }

        //  Pick the bigger of the L's and the lesser of the R's.  If
        //  L<R still, then the Q0 and Q1 trimming intersect, and we
        //  should use that intersection for the clear range.
        //  Otherwise, delete the fragment.

        uint32 l = (qltL0 < qltL1) ? qltL1 : qltL0;
        uint32 r = (qltR0 < qltR1) ? qltR0 : qltR1;

        if (l + AS_FRAG_MIN_LEN < r) {
          if (doModify) {
            setFragRecordClearRegion(&fr, l, r, AS_READ_CLEAR_OBT);
            setFrag(gkp, lid, &fr);
          }

          if (logFile)
            fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" (no overlaps)\n",
                    AS_UID_toString(uid), lid, qltL0, qltR0, l, r);
        } else {
          //  What?  No intersect...too small?  Delete it!
          //
          if (doModify)
            delFrag(gkp, lid);

          if (logFile)
            if (l < r)
              fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" (no overlaps, intersection too short, deleted)\n",
                      AS_UID_toString(uid), lid, qltL0, qltR0, qltL1, qltR1);
            else
              fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" (no overlaps, no intersection, deleted)\n",
                      AS_UID_toString(uid), lid, qltL0, qltR0, qltL1, qltR1);
        }
      }

      lid++;
    }



    //  Read the fragment from the store, compute the quality trim
    //  points.  All the values from overlap are off by the original
    //  clear range, we add it back in as we decode the string.
    //
    uint32 qltL = 0;
    uint32 qltR = 0;

    getFrag(gkp, iid, &fr, FRAG_S_INF | FRAG_S_QLT);

    uint32 qltLQ1 = getFragRecordClearRegionBegin(&fr, AS_READ_CLEAR_OBTINI);
    uint32 qltRQ1 = getFragRecordClearRegionEnd  (&fr, AS_READ_CLEAR_OBTINI);
    AS_UID uid    = getFragRecordUID(&fr);
    AS_IID lib    = getFragRecordLibraryIID(&fr);

    GateKeeperLibraryRecord  *gklr = getGateKeeperLibrary(gkp, lib);

    //  Only proceed if we're mutable.
    //
    if ((gklr) && (gklr->doNotOverlapTrim)) {
      if (logFile)
        fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" (immutable)\n",
                AS_UID_toString(uid), lid, qltLQ1, qltRQ1, qltLQ1, qltRQ1);
    } else {
      uint32 min5   = atoi(W[1]) + qltLQ1;
      uint32 minm5  = atoi(W[2]) + qltLQ1;
      uint32 minm5c = atoi(W[3]);
      uint32 mode5  = atoi(W[4]) + qltLQ1;
      uint32 mode5c = atoi(W[5]);
      uint32 max3   = atoi(W[6]) + qltLQ1;
      uint32 maxm3  = atoi(W[7]) + qltLQ1;
      uint32 maxm3c = atoi(W[8]);
      uint32 mode3  = atoi(W[9]) + qltLQ1;
      uint32 mode3c = atoi(W[10]);

      //  Adjust the mode and max/min(m) counts and values if the
      //  min/max or min/max(m) values are within OBT_MODE_WIGGLE
      //  (currently, 5 bp) of the mode

      if ((min5 != mode5) && ((min5 + OBT_MODE_WIGGLE) >= mode5)){
        if (minm5 < mode5) {
          if (min5 < minm5) {
            mode5c++;
          }
          mode5c += minm5c;
          minm5c = mode5c;
        } else if (minm5 == mode5) {
          mode5c++;
          minm5c++;
          assert(minm5c == mode5c);
        } else {
          assert(0);
        }
        min5 = mode5;
        minm5 = mode5;
      } else if ((minm5 != mode5) && ((minm5 + OBT_MODE_WIGGLE) >= mode5)){
        mode5c += minm5c;
        minm5c = mode5c;
        minm5 = mode5;
      }

      if ((max3 != mode3) && ((max3 - OBT_MODE_WIGGLE) <= mode3)){
        if (maxm3 > mode3) {
          if (max3 > maxm3) {
            mode3c++;
          }
          mode3c += maxm3c;
          maxm3c = mode3c;
        } else if (maxm3 == mode3) {
          mode3c++;
          maxm3c++;
          assert(maxm3c == mode3c);
        } else {
          assert(0);
        }
        max3 = mode3;
        maxm3 = mode3;
      } else if ((maxm3 != mode3) && ((maxm3 - OBT_MODE_WIGGLE) <= mode3)){
        mode3c += maxm3c;
        maxm3c = mode3c;
        maxm3 = mode3;
      }

#define BADIDEA
#ifdef BADIDEA
      if ((gklr) && (gklr->doNotQVTrim)) {
        qltLQ1 = 0;
        qltRQ1 = getFragRecordSequenceLength(&fr);

        qltL   = min5;
        if (mode5c > 4) qltL = mode5;
        if (minm5c > 8) qltL = minm5;

        qltR   = max3;
        if (mode3c > 4) qltR = mode3;
        if (maxm3c > 8) qltR = maxm3;
      } else {
        doTrim(&fr, minQuality, qltL, qltR);
      }
#else
      doTrim(&fr, minQuality, qltL, qltR);
#endif

      uint32 left  = 0;
      uint32 right = 0;

      //
      //  The sort makes sure that all overlaps are "good quality".
      //

      //  0) if the left quality trim is less than the mode-of-the-mode,
      //     reset it.
      //
      if (qltL < qltLQ1) {
        stats[15]++;
        qltL = qltLQ1;
      }
      if (qltR > qltRQ1) {
        stats[16]++;
        qltR = qltRQ1;
      }
      if (qltL < modes[lib].get()) {
        stats[0]++;
        qltL = modes[lib].get();
      }
      if (qltR < qltL) {
        stats[17]++;
        qltR = qltL;
      }

      //  1) if the quality-range is < 100, or
      //  2) if the quality-range is < 200 and the intersection with the
      //     overlap-range is < 100
      //  be more conservative
      //
      if ((qltL + OBT_CQ_LENGTH > qltR) ||
          ((qltL + OBT_CQO_LENGTH > qltR) &&
           ((min5 + OBT_CQO_OVERLAP > qltR) ||
            (qltL + OBT_CQO_OVERLAP > max3)))) {
        stats[1]++;  //  short quality

#if 0
        if ((qltL  + OBT_CQ_SHORT > qltR) ||
            (mode5 + OBT_CQ_SHORT > qltR) ||
            (qltL  + OBT_CQ_SHORT > mode3)) {
#else
          // test granger email re ixodes
        if ((qltL  + OBT_CQ_SHORT > qltR)) {
#endif
          stats[2]++;
          left  = 0;
          right = 0;
        } else {
          stats[3]++;
          left = mode5;
          right = mode3;
        }
      } else {
        stats[4]++;

        if ((minm5 < qltL + OBT_QLT_CLOSE_5) && (minm5c > 1)) {
          stats[6]++;
          left = minm5;
        } else if ((mode5 < qltL + OBT_QLT_CLOSE_5) && (mode5c > 0)) {
          stats[5]++;
          left = mode5;
        } else if ((mode5c > 0) && ((min5 <= qltL) || (min5 < qltL + OBT_QLT_FAR_5))) {
          stats[11]++;
          left = min5;
        } else {
          stats[7]++;
          left = qltL;
        }

        if ((maxm3 >= qltR) && (maxm3c > 1)) {
          stats[8]++;
          right = maxm3;
        } else if ((mode3 == max3) && (mode3 == maxm3) && (mode3c > 1) &&
                   ((mode3 >= qltR) || (qltR < mode3 + OBT_QLT_MODE3))) {
          stats[9]++;
          right = mode3;
        } else if ((maxm3c > 1) && (maxm3 < qltR) && (max3 > qltR) && (max3 < maxm3 + OBT_QLT_CLOSE_MAXM3)) {
          stats[14]++;
          right = maxm3;
        } else if ((mode3c > 0) && ((max3 >= qltR) || (qltR < max3 + OBT_QLT_CLOSE_MAX3))) {
          stats[12]++;
          right = max3;
        } else {
          stats[10]++;
          right = qltR;
        }

        if ((left == 0) && (right == 0)) {
          stats[18]++;
          fprintf(stderr, "INVALID CLEAR from OVL:\t"F_U64"\t"F_U32"\t"F_U32"\t->\t"F_U32"\t"F_U32"\t--\t%s\n",
                  iid, qltL, qltR, left, right,
                  line);
        }
      }

      //  If after all that, we have an invalid range, throw it out.
      //
      //  What has happened (once) is for a fragment:
      //
      //                OVL                QUAL
      //      ------|---------|---------|---------|---------
      //              --------
      //               -------
      //                ------
      //
      //  (that is, lots of 5' alignment in the overlap, none in 3', no
      //  intersection between OVL and QUAL)
      //
      //  We picked the 5'mode (the right end of OVL) and the left
      //  quality (left end of QUAL), which isn't a valid region.
      //
      //  XXX: We should check if the quality of the OVL region is
      //  decent, and then use that if it is also large.
      //
      if ((left + right > 0) && ((left + AS_FRAG_MIN_LEN) > right)) {
        stats[13]++;
#if 0
        fprintf(stderr, "INVALID CLEAR:\t"F_U64"\t"F_U32"\t"F_U32"\t->\t"F_U32"\t"F_U32"\t--\t%s\n",
                iid, qltL, qltR, left, right,
                line);
#endif
        left  = 0;
        right = 0;
      }

      //  Report the iid, quality trim points, and final trimming.  You
      //  can use the input to this executable to get the overlap trim
      //  points.  The original trim points are elsewhere.
      //

      if ((left == 0) && (right == 0)) {
        stats[19]++;

        if (logFile)
          fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32" (deleted, too short)\n",
                  AS_UID_toString(uid), iid, qltL, qltR, left, right);

        if (doModify)
          delFrag(gkp, iid);
      } else {
        stats[20]++;

        if (logFile)
          fprintf(logFile, "%s,"F_U64"\t"F_U32"\t"F_U32"\t"F_U32"\t"F_U32"\n",
                  AS_UID_toString(uid), iid, qltL, qltR, left, right);

        if (doModify) {
          setFragRecordClearRegion(&fr, left, right, AS_READ_CLEAR_OBT);
          setFrag(gkp, iid, &fr);
        }
      }
    }

    lid = iid;

    readLine(O);
  }

  closeGateKeeperStore(gkp);
  fclose(O);

  //  Report statistics
  //

  fprintf(staFile, F_U32":\treset qltL to mode-of-5'mode\n", stats[0]);
  fprintf(staFile, F_U32":\treset qltL to vector left\n", stats[15]);
  fprintf(staFile, F_U32":\treset qltR to vector right\n", stats[16]);
  fprintf(staFile, F_U32":\treset qltR to qltL due to inconsistency\n", stats[17]);
  fprintf(staFile, F_U32":\tshort quality:\n", stats[1]);
  fprintf(staFile, F_U32":\t  very short quality < %d or very short in common < %d, discard frag\n", stats[2], OBT_CQ_SHORT, OBT_CQ_SHORT);
  fprintf(staFile, F_U32":\t  short quality use overlap modes\n", stats[3]);
  fprintf(staFile, F_U32":\tuse the min/max/mode:\n", stats[4]);
  fprintf(staFile, F_U32":\t  use mode (5')\n", stats[5]);
  fprintf(staFile, F_U32":\t  use min>1 (5')\n", stats[6]);
  fprintf(staFile, F_U32":\t  use quality (5')\n", stats[7]);
  fprintf(staFile, F_U32":\t  use max>1 (3')\n", stats[8]);
  fprintf(staFile, F_U32":\t  use mode (3')\n", stats[9]);
  fprintf(staFile, F_U32":\t  use quality (3')\n", stats[10]);
  fprintf(staFile, F_U32":\t  use min (5')\n", stats[11]);
  fprintf(staFile, F_U32":\t  use max (3')\n", stats[12]);
  fprintf(staFile, F_U32":\t  use max>1 close to max (3')\n", stats[14]);
  fprintf(staFile, F_U32":\tinvalid clear after merging overlaps (should be 0)\n", stats[18]);
  fprintf(staFile, F_U32":\tshort or inconsistent\n", stats[13]);
  fprintf(staFile, F_U32":\tdeleted fragment due to zero clear\n", stats[19]);
  fprintf(staFile, F_U32":\tupdated fragment clear range\n", stats[20]);

  fclose(logFile);
  fclose(staFile);

  return(0);
}
