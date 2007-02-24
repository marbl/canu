#include "trim.H"
#include "maps.H"
#include "constants.H"

//  Read a fragStore, does quality trimming based on quality scores,
//  modifies the original clear range in the store.
//
//  Optionally intersects the quality trim with a vector trim.
//
//  Optionally does NOT modify a list of fragment UIDs

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-q quality] [-v vectorlist] [-i uidlist] [-update] [-replace] [-log logfile] -frg some.gkpStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -q quality    Find quality trim points using 'quality' as the base.\n");
  fprintf(stderr, "  -v vector     Intersect the quality trim with a vector trim.\n");
  fprintf(stderr, "  -i uidlist    Never, ever modify these fragments.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -update       Update the clear range in the fragStore.\n");
  fprintf(stderr, "  -replace      Replace non ACGT with random ACGT with low quality (unimplemented).\n");
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
  char   *vectorFileName      = 0L;
  char   *immutableFileName   = 0L;
  bool    doUpdate            = false;
  bool    doReplace           = false;
  char   *gkpStore            = 0L;
  FILE   *logFile             = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-quality", 2) == 0) {
      minQuality = qual.lookupNumber(strtou32bit(argv[++arg], 0L));
    } else if (strncmp(argv[arg], "-vector", 2) == 0) {
      vectorFileName = argv[++arg];
    } else if (strncmp(argv[arg], "-immutable", 2) == 0) {
      immutableFileName = argv[++arg];

    } else if (strncmp(argv[arg], "-update", 2) == 0) {
      doUpdate = true;
    } else if (strncmp(argv[arg], "-replace", 2) == 0) {
      doReplace = true;

    } else if (strncmp(argv[arg], "-memory", 2) == 0) {
      ++arg;
      //batchMemory = strtou32bit(argv[arg], 0L) * 1024 * 1024;
      //batchSize   = batchMemory / sizeof();

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

  //srand48(time(NULL));

  vectorMap  m;
  m.readVectorMap(vectorFileName);
  m.readImmutableMap(immutableFileName);

  //  Open the store
  //
  GateKeeperStore  *gkp = openGateKeeperStore(gkpStore, doUpdate);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStore);
    exit(1);
  }

  u32bit        firstElem = getFirstElemFragStore(gkp);
  u32bit        lastElem  = getLastElemFragStore(gkp) + 1;

  fragRecord   *fr = new_fragRecord();

  u32bit        qltL = 0;
  u32bit        qltR = 0;
  u32bit        vecL = 0;
  u32bit        vecR = 0;

  u32bit        stat_notPresent  = 0;
  u32bit        stat_noIntersect = 0;
  u32bit        stat_change      = 0;
  u32bit        stat_noChange    = 0;
  u32bit        stat_immutable   = 0;

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    u64bit uid = 0;

    getFrag(gkp, elem, fr, FRAG_S_INF | FRAG_S_SEQ | FRAG_S_QLT);
    uid = getFragRecordUID(fr);

    //  Bail now if we've been told to not modify this read.  We do
    //  not print a message in the log.
    //
    if (m.exists(uid) && (m[uid].immutable == 1)) {
      stat_immutable++;
      continue;
    }

    doTrim(fr, minQuality, qltL, qltR);
    vecL = qltL;
    vecR = qltR;

    //  Intersect with the vector clear range, if it exists
    //
    if (!m.exists(uid) || (m[uid].hasVec == 0)) {
      //  uid not present in our input list, do nothing.
      stat_notPresent++;
    } else if ((m[uid].vecL > vecR) || (m[uid].vecR < vecL)) {
      //  don't intersect, we've already set vecL and vecR appropriately.
      stat_noIntersect++;
    } else {
      //  They intersect.  Pick the largest begin and the smallest end

      bool changed = false;

      if (vecL < m[uid].vecL) {
        changed = true;
        vecL = m[uid].vecL;
      }
      if (m[uid].vecR < vecR) {
        changed = true;
        vecR = m[uid].vecR;
      }

      if (changed)
        stat_change++;
      else
        stat_noChange++;
    }


    if (logFile) {
      unsigned int      clrBeg = getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_ORIG);
      unsigned int      clrEnd = getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_ORIG);

      fprintf(logFile, u64bitFMT","u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"%s\n",
              uid, elem,
              (u32bit)clrBeg, (u32bit)clrEnd,
              qltL, qltR,
              m[uid].vecL, m[uid].vecR,
              vecL, vecR,
              ((vecL + OBT_MIN_LENGTH) > vecR) ? " (deleted)" : "");

    }

    if (doReplace) {
    }

    if (doUpdate) {
      setFragRecordClearRegion(fr, vecL, vecR, AS_READ_CLEAR_OBTINI);
      setFrag(gkp, elem, fr);

      if ((vecL + OBT_MIN_LENGTH) > vecR)
        delFrag(gkp, elem);
    }
  }

  closeGateKeeperStore(gkp);

  fprintf(stderr, "Fragments with no vector clear:  "u32bitFMT"\n", stat_notPresent);
  fprintf(stderr, "Fragments with no intersection:  "u32bitFMT"\n", stat_noIntersect);
  fprintf(stderr, "Fragments with vector trimmed:   "u32bitFMT"\n", stat_change);
  fprintf(stderr, "Fragments low quality vector:    "u32bitFMT"\n", stat_noChange);
  fprintf(stderr, "Fragments marked immutable:      "u32bitFMT"\n", stat_immutable);
}
