#include "trim.H"
#include "constants.H"
#include "maps.H"

using namespace std;
#include <set>

//  Reads the output of sort-overlap-trim, does the actual trim-point
//  decision, updates the frgStore.
//
//  There is a great deal of overlap with qualityTrim -- we need to do
//  exactly the same stuff as done there, and then some.



//  XXX:  This is an ugly hack, someone should make it read variable
//  length lines.  Possilby already in libutil....
u32bit         lineMax = 128 * 1024;
char          *line    = 0L;


void
readLine(FILE *F) {
  fgets(line, lineMax, F);
  line[lineMax-1] = 0;
  assert(strlen(line) < (lineMax-1));
}


u32bit
findModeOfFivePrimeMode(FragStoreHandle fs, char *name) {
  u32bit         mode5;
  u32bit         iid;
  u32bit         histo[2048] = {0};


  ReadStructp    rd = new_ReadStruct();
  unsigned int   cl = 0;
  unsigned int   cr = 0;

  errno = 0;
  FILE *O = fopen(name, "r");
  if (errno)
    fprintf(stderr, "Can't open overlap-trim file %s: %s\n", name, strerror(errno)), exit(1);

  readLine(O);
  while (!feof(O)) {
    splitToWords W(line);
    iid   = strtou32bit(W[0], 0L);
    mode5 = strtou32bit(W[4], 0L);

    //  Grab the fragment to get the clear range, so we can find the real mode5
    //
    getFragStore(fs, iid, FRAG_S_ALL, rd);
    getClearRegion_ReadStruct(rd, &cl, &cr, READSTRUCT_ORIGINAL);

    mode5 += cl;

    if (mode5 < 2047)
      histo[mode5]++;

    readLine(O);
  }

  fclose(O);

  //  find the largest value -- but ignore zero!
  //
  mode5 = 1;
  for (u32bit i=1; i<2048; i++)
    if (histo[i] > histo[mode5])
      mode5 = i;

  return(mode5);
}




int
main(int argc, char **argv) {
  u32bit   stats[32]         = {0};
  FILE    *O                 = 0L;
  FILE    *logFile           = 0L;
  FILE    *staFile           = 0L;
  char    *frgStore          = 0L;
  char    *ovlFile           = 0L;
  bool     doModify          = true;  //  Make this false for testing
  char    *immutableFileName = 0L;

  line = new char [lineMax];

  if (argc < 5) {
    fprintf(stderr, "usage: %s [-immutable uidlist] [-log log] -frg frgStore -ovl overlap-trim\n", argv[0]);
    fprintf(stderr, "  -immutable uidlist    Never, ever modify these fragments.\n");
    fprintf(stderr, "  -ovl o                Read consolidated overlaps from here.n");
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
    } else if (strncmp(argv[arg], "-immutable", 2) == 0) {
      immutableFileName = argv[++arg];
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);

      char staName[strlen(argv[arg]) + 16];
      sprintf(staName, "%s.stats", argv[arg]);

      staFile = fopen(staName, "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the statistics: %s\n", staName, strerror(errno)), exit(1);
    }
    arg++;
  }




  ////////////////////////////////////////
  //
  //  Open the frgStore, prepare for reading fragments
  //
  FragStoreHandle fs = openFragStore(frgStore, doModify ? "r+" : "r");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open fragStore %s!\n", frgStore);
    exit(1);
  }
  u32bit   firstElem = getFirstElemFragStore(fs);
  u32bit   lastElem  = getLastElemFragStore(fs) + 1;
  ReadStructp       rd = new_ReadStruct();


  ////////////////////////////////////////
  //
  //  Open the overlap file
  //
  errno = 0;
  O = fopen(ovlFile, "r");
  if (errno)
    fprintf(stderr, "Can't open overlap-trim file %s: %s\n", ovlFile, strerror(errno)), exit(1);


  double  minQuality = qual.lookupNumber(20);


  ////////////////////////////////////////
  //
  //  Build a list of the immutable fragments
  //
  vectorMap  immutable;
  immutable.readImmutableMap(immutableFileName);

  

  ////////////////////////////////////////
  //
  //  Find the mode of the 5'mode.  Read the whole overlap-trim file, counting the 5'mode
  //
  u32bit modemode5 = findModeOfFivePrimeMode(fs, ovlFile);

  if (staFile)
    fprintf(staFile, "Mode of the 5'-mode is "u32bitFMT"\n", modemode5);


  //  Stream through the overlap-trim file, applying some rules to
  //  decide on the correct trim points

  u64bit  lid = 0;
  u64bit  iid = 0;

  readLine(O);
  while (!feof(O)) {
    splitToWords  W(line);
    iid    = strtou64bit(W[0], 0L);


    //  Report the frags that had no overlap -- we update the clear region
    //  to the Q20 clear region.  If the two clear regions do not intersect,
    //  we delete the fragment.
    //
    lid++;
    while (lid < iid) {
      getFragStore(fs, lid, FRAG_S_ALL, rd);

      u32bit  qltL0, qltR0, qltL1, qltR1;

      unsigned int l=0;
      unsigned int r=0;
      getClearRegion_ReadStruct(rd, &l, &r, READSTRUCT_ORIGINAL);
      qltL0 = l;
      qltR0 = r;

      doTrim(rd, minQuality, qltL1, qltR1);

      //  Pick the bigger of the L's and the lesser of the R's.  If
      //  L<R still, then the Q0 and Q1 trimming intersect, and we
      //  should use that intersection for the clear range.
      //  Otherwise, delete the fragment.

      if (l < qltL1)  l = qltL1;
      if (r > qltR1)  r = qltR1;

      if (l < r) {
        if (doModify) {
          setClearRegion_ReadStruct(rd, l, r, READSTRUCT_OVL);
          if (setFragStore(fs, lid, rd)) {
            fprintf(stderr, "setFragStore() failed.\n");
            exit(1);
          }
        }

        if (logFile)
          fprintf(logFile, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT" (no overlaps)\n",
                  lid, qltL0, qltR0, l, r);
      } else {
        //  What?  No intersect?  Delete it!
        if (doModify)
          deleteFragStore(fs, lid);

        if (logFile)
          fprintf(logFile, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT" (no overlaps, no intersection, deleted)\n",
                  lid, qltL0, qltR0, qltL1, qltR1);
      }

      lid++;
    }



    //  Read the fragment from the store, compute the quality trim
    //  points.  All the values from overlap are off by the original
    //  clear range, we add it back in as we decode the string.
    //
    u32bit qltL = 0;
    u32bit qltR = 0;

    getFragStore(fs, iid, FRAG_S_ALL, rd);
    unsigned int l=0;
    unsigned int r=0;
    getClearRegion_ReadStruct(rd, &l, &r, READSTRUCT_ORIGINAL);
    u32bit qltLQ1 = l;
    u32bit qltRQ1 = r;

    u64bit uid = 0;
    getAccID_ReadStruct(rd, &uid);

    //  Only proceed if we're mutable.
    //
    if (immutable.exists(uid)) {
      if (logFile)
        fprintf(logFile, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT" (immutable)\n",
                lid, qltLQ1, qltRQ1, qltLQ1, qltRQ1);
    } else {
      u32bit min5   = strtou32bit(W[1], 0L) + qltLQ1;
      u32bit minm5  = strtou32bit(W[2], 0L) + qltLQ1;
      u32bit minm5c = strtou32bit(W[3], 0L);
      u32bit mode5  = strtou32bit(W[4], 0L) + qltLQ1;
      u32bit mode5c = strtou32bit(W[5], 0L);
      u32bit max3   = strtou32bit(W[6], 0L) + qltLQ1;
      u32bit maxm3  = strtou32bit(W[7], 0L) + qltLQ1;
      u32bit maxm3c = strtou32bit(W[8], 0L);
      u32bit mode3  = strtou32bit(W[9], 0L) + qltLQ1;
      u32bit mode3c = strtou32bit(W[10], 0L);

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


      doTrim(rd, minQuality, qltL, qltR);

      u32bit left  = 0;
      u32bit right = 0;

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
      if (qltL < modemode5) {
        stats[0]++;
        qltL = modemode5;
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
        if ((qltL  + OBT_CQ_SHORT > qltR) ||
            (mode5 + OBT_CQ_SHORT > qltR) ||
            (qltL  + OBT_CQ_SHORT > mode3)) {
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
          fprintf(stderr, "INVALID CLEAR from OVL: "u64bitFMTW(8)" "u32bitFMTW(3)"-"u32bitFMTW(4)" -> "u32bitFMTW(3)"-"u32bitFMTW(4)" -- %s\n",
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
      if ((left + right > 0) && ((left + OBT_MIN_LENGTH) > right)) {
        stats[13]++;
#if 0
        fprintf(stderr, "INVALID CLEAR: "u64bitFMTW(8)" "u32bitFMTW(3)"-"u32bitFMTW(4)" -> "u32bitFMTW(3)"-"u32bitFMTW(4)" -- %s\n",
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
          fprintf(logFile, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT" (deleted, too short)\n",
                  iid, qltL, qltR, left, right);

        if (doModify)
          deleteFragStore(fs, iid);
      } else {
        stats[20]++;

        if (logFile)
          fprintf(logFile, u64bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
                  iid, qltL, qltR, left, right);

        if (doModify) {
          setClearRegion_ReadStruct(rd, left, right, READSTRUCT_OVL);
          if (setFragStore(fs, iid, rd)) {
            fprintf(stderr, "setFragStore() failed.\n");
            exit(1);
          }
        }
      }
    }

    lid = iid;

    readLine(O);
  }

  closeFragStore(fs);
  fclose(O);

  //  Report statistics
  //

  fprintf(staFile, u32bitFMTW(8)": reset qltL to mode-of-5'mode\n", stats[0]);
  fprintf(staFile, u32bitFMTW(8)": reset qltL to vector left\n", stats[15]);
  fprintf(staFile, u32bitFMTW(8)": reset qltR to vector right\n", stats[16]);
  fprintf(staFile, u32bitFMTW(8)": reset qltR to qltL due to inconsistency\n", stats[17]);
  fprintf(staFile, u32bitFMTW(8)": short quality:\n", stats[1]);
  fprintf(staFile, u32bitFMTW(8)":   very short quality < %d or very short in common < %d, discard frag\n", stats[2], OBT_CQ_SHORT, OBT_CQ_SHORT);
  fprintf(staFile, u32bitFMTW(8)":   short quality use overlap modes\n", stats[3]);
  fprintf(staFile, u32bitFMTW(8)": use the min/max/mode:\n", stats[4]);
  fprintf(staFile, u32bitFMTW(8)":   use mode (5')\n", stats[5]);
  fprintf(staFile, u32bitFMTW(8)":   use min>1 (5')\n", stats[6]);
  fprintf(staFile, u32bitFMTW(8)":   use quality (5')\n", stats[7]);
  fprintf(staFile, u32bitFMTW(8)":   use max>1 (3')\n", stats[8]);
  fprintf(staFile, u32bitFMTW(8)":   use mode (3')\n", stats[9]);
  fprintf(staFile, u32bitFMTW(8)":   use quality (3')\n", stats[10]);
  fprintf(staFile, u32bitFMTW(8)":   use min (5')\n", stats[11]);
  fprintf(staFile, u32bitFMTW(8)":   use max (3')\n", stats[12]);
  fprintf(staFile, u32bitFMTW(8)":   use max>1 close to max (3')\n", stats[14]);
  fprintf(staFile, u32bitFMTW(8)": invalid clear after merging overlaps (should be 0)\n", stats[18]);
  fprintf(staFile, u32bitFMTW(8)": short or inconsistent\n", stats[13]);
  fprintf(staFile, u32bitFMTW(8)": deleted fragment due to zero clear\n", stats[19]);
  fprintf(staFile, u32bitFMTW(8)": updated fragment clear range\n", stats[20]);

  fclose(logFile);
  fclose(staFile);

  return(0);
}
