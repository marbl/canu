#include "trim.H"

using namespace std;
#include <map>

//  Intersects the Q1 trim in a frag store with the one supplied.  Updates the trimming, writes a report.

void
usage(char *name) {
  fprintf(stderr, "usage: %s [options] -intersect file -frg some.frgStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -quiet     Don't update the clear range in the fragStore.\n");
  fprintf(stderr, "  -update    Update the clear range in the fragStore (the default).\n");
  fprintf(stderr, "  -memory    Use X MB of memory per batch -- greatly improves performance!\n");
  fprintf(stderr, "  -intersect Intersect our quality trim range with the range in 'file'\n");
  fprintf(stderr, "             Format: 'UID low high'\n");
  fprintf(stderr, "             NOTE:  Assume base-based positions!\n");
  fprintf(stderr, "  -frg       Operate on this frgStore\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
  fprintf(stderr, "    iid  originalBegin originalEnd  intersectBegin intersectEnd  newBegin newEnd\n");

}

int
main(int argc, char **argv) {
  bool    doUpdate     = true;
  char   *intFileName  = 0L;
  char   *frgStore     = 0L;
  FILE   *logFile      = 0L;
  u32bit  batchMemory  = 0;  //  XXX: unimplemented
  u32bit  batchSize    = 0;

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-intersect", 2) == 0) {
      intFileName = argv[++arg];
    } else if (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
    } else if (strncmp(argv[arg], "-update", 2) == 0) {
      doUpdate = true;
    } else if (strncmp(argv[arg], "-quiet", 2) == 0) {
      doUpdate = false;
    } else if (strncmp(argv[arg], "-memory", 2) == 0) {
      ++arg;
      //batchMemory = strtou32bit(argv[arg], 0L) * 1024 * 1024;
      //batchSize   = batchMemory / sizeof();
    } else if (strncmp(argv[arg], "-log", 2) == 0) {
      errno=0;
      logFile = fopen(argv[++arg], "w");
      if (errno)
        fprintf(stderr, "Failed to open %s for writing the log: %s\n", argv[arg], strerror(errno)), exit(1);
    } else {
      fprintf(stderr, "Invalid option: '%s'\n", argv[arg]);
      usage(argv[0]);
      exit(1);
    }
    arg++;
  }

  if ((!frgStore) || (!intFileName)) {
    usage(argv[0]);
    exit(1);
  }

  //  Open the store
  //
  FragStoreHandle   fs = openFragStore(frgStore, doUpdate ? "r+" : "r");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(fs);
  u32bit   lastElem  = getLastElemFragStore(fs) + 1;

  ReadStructp       rd = new_ReadStruct();

  u32bit            qltL = 0;
  u32bit            qltR = 0;


  //  Read the intersection into a map.  Being lazy, we're actually
  //  using a map, and twice.
  //
  map<u64bit, u32bit>   intL, intR;
  errno = 0;
  FILE *intFile = fopen(intFileName, "r");
  if (errno)
    fprintf(stderr, "Can't open '%s': %s\n", intFileName, strerror(errno)), exit(1);
  //speedCounter  *C = new speedCounter("%7.2f uids -- %5.2f Muids/second\r",
  //                                    1000000.0, 0x7fff, true);
  //C->enableLiner();
  char intLine[1024] = {0};
  fgets(intLine, 1024, intFile);
  while (!feof(intFile)) {
    chomp(intLine);
    splitToWords  W(intLine);
    if ((W[0] == 0L) || (W[1] == 0L) || (W[2] == 0L)) {
      fprintf(stderr, "WARNING: Invalid line '%s'\n", intLine);
    } else {
      u64bit  uid  = strtou64bit(W[0], 0L);
      u32bit  intl = strtou32bit(W[1], 0L);
      u32bit  intr = strtou32bit(W[2], 0L);

      //  Silently swap, if needed
      if (intl > intr) {
        u32bit  s = intl;
        intl = intr;
        intr = s;
      }

      //  These look base-based!
      intL[uid] = intl - 1;
      intR[uid] = intr - 1;
    }

    //C->tick();
    fgets(intLine, 1024, intFile);
  }
  fclose(intFile);
  //delete C;

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFragStore(fs, elem, FRAG_S_ALL, rd);

    //  Read the existing trim from the fragstore.
    unsigned int      clrBeg = 0;
    unsigned int      clrEnd = 0;
    getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_ORIGINAL);
    qltL = clrBeg;
    qltR = clrEnd;

    //  Read the new trim to intersect with
    u64bit uid = 0;
    getAccID_ReadStruct(rd, &uid);

    //  Do the intersection, remembering if anything changed.
    //
    //  There are several error cases:
    //  1) the two trims don't intersect.
    //  2) the new trim is invalid (length zero)
    //  3) the new trim is backwards -- we've already fixed this
    //
    bool  changed = false;

    if ((intL.count(uid) == 0) || (intR.count(uid) == 0)) {
      //  uid not present in our input list, do nothing.
    } else if ((intL[uid] > clrEnd) || (intR[uid] < clrBeg)) {
      //  don't intersect, we've already set qltL and qltR appropriately.
    } else {
      //  They intersect.  Pick the largest begin and the smallest end
      if (qltL < intL[uid]) {
        changed = true;
        qltL = intL[uid];
      }
      if (qltR > intR[uid]) {
        changed = true;
        qltR = intR[uid];
      }
    }

    if (changed) {
      sprintf(intLine, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t^\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\n",
              elem,
              (u32bit)clrBeg,
              (u32bit)clrEnd,
              intL[uid],
              intR[uid],
              qltL,
              qltR);
      if (logFile)
        fprintf(logFile, intLine);
      else
        fprintf(stdout, intLine);

      if (doUpdate) {
        setClearRegion_ReadStruct(rd, qltL, qltR, READSTRUCT_ORIGINAL);
        if (setFragStore(fs, elem, rd)) {
          fprintf(stderr, "setFragStore() failed.\n");
          exit(1);
        }
      }
    }
  }

  closeFragStore(fs);
}
