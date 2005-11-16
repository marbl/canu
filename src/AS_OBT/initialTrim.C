#include "trim.H"

using namespace std;
#include <map>

//  Read a fragStore, does quality trimming based on quality scores,
//  modifies the original clear range in the store.
//
//  Optionally intersects the quality trim with a vector trim.
//
//  Optionally does NOT modify a list of fragment UIDs
//


void
usage(char *name) {
  fprintf(stderr, "usage: %s [-q quality] [-v vectorlist] [-i uidlist] [-update] [-replace] [-log logfile] -frg some.frgStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -q quality    Find quality trim points using 'quality' as the base.\n");
  fprintf(stderr, "  -v vector     Intersect the quality trim with a vector trim.\n");
  fprintf(stderr, "  -i uidlist    Never, ever modify these fragments.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -update       Update the clear range in the fragStore.\n");
  fprintf(stderr, "  -replace      Replace non ACGT with random ACGT with low quality (unimplemented).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -log X        Report the iid, original trim and new quality trim\n");
  fprintf(stderr, "  -frg F        Operate on this frgStore\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
  fprintf(stderr, "    iid originalBegin originalEnd newBegin newEnd\n");
  fprintf(stderr, "    uid,iid origBegin origEnd qualBegin qualEnd vecBeg vecEnd newBegin newEnd\n");
}


class vectorInfo {
public:
  vectorInfo() {
    vecL      = 0;
    vecR      = 0;
    hasVec    = 0;
    immutable = 0;
  };

  u32bit  vecL      : 12;
  u32bit  vecR      : 12;
  u32bit  hasVec    : 1;
  u32bit  immutable : 1;
};


bool
readVectorIntersection(char *vectorFileName, map<u64bit, vectorInfo> &m) {
  bool fatal = false;

  if (vectorFileName == 0L)
    return(false);

  errno = 0;
  FILE *intFile = fopen(vectorFileName, "r");
  if (errno)
    fprintf(stderr, "Can't open '%s': %s\n", vectorFileName, strerror(errno)), exit(1);

  char intLine[1024] = {0};
  fgets(intLine, 1024, intFile);
  while (!feof(intFile)) {
    chomp(intLine);
    splitToWords  W(intLine);
    if ((W[0] == 0L) || (W[1] == 0L) || (W[2] == 0L) || (W[3] != 0L)) {
      fprintf(stderr, "readVectorIntersection()-- Invalid line '%s'\n", intLine);
      fatal = true;
    } else {

      //  These look base-based!

      u64bit  uid  = strtou64bit(W[0], 0L);
      u32bit  intl = strtou32bit(W[1], 0L) - 1;
      u32bit  intr = strtou32bit(W[2], 0L) - 1;

      //  Silently swap, if needed

      if (intl > intr) {
        u32bit  s = intl;
        intl = intr;
        intr = s;
      }

      vectorInfo v;
      v.vecL      = intl;
      v.vecR      = intr;
      v.hasVec    = 1;
      v.immutable = 0;

      m[uid] = v;
    }

    fgets(intLine, 1024, intFile);
  }

  fclose(intFile);

  return(fatal);
}




bool
readImmutable(char *immutableFileName, map<u64bit, vectorInfo> &m) {
  bool fatal = false;

  if (immutableFileName == 0L)
    return(false);

  errno = 0;
  FILE *intFile = fopen(immutableFileName, "r");
  if (errno)
    fprintf(stderr, "Can't open '%s': %s\n", immutableFileName, strerror(errno)), exit(1);

  char intLine[1024] = {0};
  fgets(intLine, 1024, intFile);
  while (!feof(intFile)) {
    chomp(intLine);
    splitToWords  W(intLine);
    if ((W[0] == 0L) || (W[1] != 0L)) {
      fprintf(stderr, "readImmutable()-- Invalid line '%s'\n", intLine);
      fatal = true;
    } else {
      u64bit  uid  = strtou64bit(W[0], 0L);
      m[uid].immutable = 1;
    }

    fgets(intLine, 1024, intFile);
  }

  fclose(intFile);

  return(fatal);
}



int
main(int argc, char **argv) {
  double  minQuality          = qual.lookupNumber(20);
  char   *vectorFileName      = 0L;
  char   *immutableFileName   = 0L;
  bool    doUpdate            = false;
  bool    doReplace           = false;
  char   *frgStore            = 0L;
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
      frgStore = argv[++arg];

    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if (!frgStore) {
    usage(argv[0]);
    exit(1);
  }

  //srand48(time(NULL));

  map<u64bit, vectorInfo> m;

  bool fatal = false;
  fatal |= readVectorIntersection(vectorFileName, m);
  fatal |= readImmutable(immutableFileName, m);
  if (fatal) {
    fprintf(stderr, "%s: I encountered errors reading the input.  Fatal error.\n", argv[0]);
    exit(1);
  }

  //  Open the store
  //
  FragStoreHandle   fs = openFragStore(frgStore, doUpdate ? "r+" : "r");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  u32bit        firstElem = getFirstElemFragStore(fs);
  u32bit        lastElem  = getLastElemFragStore(fs) + 1;

  ReadStructp   rd = new_ReadStruct();

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

    getFragStore(fs, elem, FRAG_S_ALL, rd);
    getAccID_ReadStruct(rd, &uid);

    if (m[uid].immutable == 1) {
      stat_immutable++;
      continue;
    }

    doTrim(rd, minQuality, qltL, qltR);
    vecL = qltL;
    vecR = qltR;

    if (m[uid].hasVec == 0) {
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
      unsigned int      clrBeg = 0;
      unsigned int      clrEnd = 0;
      getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_ORIGINAL);

      fprintf(logFile, u64bitFMT","u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\n",
              uid, elem,
              (u32bit)clrBeg, (u32bit)clrEnd,
              qltL, qltR,
              m[uid].vecL, m[uid].vecR,
              vecL, vecR);
    }

    if (doReplace) {
    }

    if (doUpdate) {
      setClearRegion_ReadStruct(rd, vecL, vecR, READSTRUCT_ORIGINAL);
    }

    if (doReplace || doUpdate) {
      if (setFragStore(fs, elem, rd)) {
        fprintf(stderr, "setFragStore() failed.\n");
        exit(1);
      }
    }
  }

  closeFragStore(fs);

  fprintf(stderr, "Fragments with no vector clear:  "u32bitFMT"\n", stat_notPresent);
  fprintf(stderr, "Fragments with no intersection:  "u32bitFMT"\n", stat_noIntersect);
  fprintf(stderr, "Fragments with vector trimmed:   "u32bitFMT"\n", stat_change);
  fprintf(stderr, "Fragments low quality vector:    "u32bitFMT"\n", stat_noChange);
  fprintf(stderr, "Fragments marked immutable:      "u32bitFMT"\n", stat_immutable);
}







