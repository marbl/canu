#include "trim.H"

//  Read a fragStore, does quality trimming based on quality scores,
//  modifies the original clear range in the store.

void
usage(char *name) {
  fprintf(stderr, "usage: %s [-q quality] [-trim] [-replace] -frg some.frgStore\n", name);
  fprintf(stderr, "\n");
  fprintf(stderr, "  -q         Find quality trim points using 'quality' as the base.\n");
  fprintf(stderr, "  -update    Update the clear range in the fragStore.\n");
  fprintf(stderr, "  -replace   Replace non ACGT with random ACGT with low quality.\n");
  fprintf(stderr, "  -frg       Operate on this frgStore\n");
  fprintf(stderr, "  -log       Report the iid, original trim and new quality trim\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  -replace is NOT IMPLEMENTED.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  A report of the trimming is printed to stdout:\n");
  fprintf(stderr, "    iid originalBegin originalEnd newBegin newEnd\n");

}

int
main(int argc, char **argv) {
  double  minQuality = qual.lookupNumber(20);
  bool    doUpdate   = false;
  bool    doReplace  = false;
  char   *frgStore   = 0L;
  FILE   *logFile    = 0L;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-q") == 0) {
      ++arg;
      minQuality = qual.lookupNumber(strtou32bit(argv[arg], 0L));
    } else if (strncmp(argv[arg], "-update", 2) == 0) {
      doUpdate = true;
    } else if (strncmp(argv[arg], "-replace", 2) == 0) {
      doReplace = true;
    } else if (strncmp(argv[arg], "-frg", 2) == 0) {
      frgStore = argv[++arg];
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

  if (!frgStore) {
    usage(argv[0]);
    exit(1);
  }

  srand48(time(NULL));

  //  Open the store
  //
  FragStoreHandle   fs = openFragStore(frgStore, "rw+");
  if (fs == NULLSTOREHANDLE) {
    fprintf(stderr, "Failed to open %s\n", frgStore);
    exit(1);
  }

  u32bit   firstElem = getFirstElemFragStore(fs);
  u32bit   lastElem  = getLastElemFragStore(fs) + 1;

  ReadStructp       rd = new_ReadStruct();

  u32bit            qltL = 0;
  u32bit            qltR = 0;

  for (u32bit elem=firstElem; elem<lastElem; elem++) {
    getFragStore(fs, elem, FRAG_S_ALL, rd);

    doTrim(rd, minQuality, qltL, qltR);

    if (logFile) {
      unsigned int      clrBeg = 0;
      unsigned int      clrEnd = 0;
      getClearRegion_ReadStruct(rd, &clrBeg, &clrEnd, READSTRUCT_ORIGINAL);

      fprintf(logFile, u32bitFMT"\t"u32bitFMT"\t"u32bitFMT"\t->\t"u32bitFMT"\t"u32bitFMT"\n",
              elem,
              (u32bit)clrBeg,
              (u32bit)clrEnd,
              qltL,
              qltR);
    }

    if (doReplace) {
    }

    if (doUpdate) {
      setClearRegion_ReadStruct(rd, qltL, qltR, READSTRUCT_ORIGINAL);
    }

    if (doReplace || doUpdate) {
      if (setFragStore(fs, elem, rd)) {
        fprintf(stderr, "setFragStore() failed.\n");
        exit(1);
      }
    }
  }

  closeFragStore(fs);
}
