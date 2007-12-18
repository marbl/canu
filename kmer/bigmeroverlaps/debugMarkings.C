#include "util++.H"
#include "bio++.H"

#include "markFile.H"

//  Writes a bunch of not really useful stuff to stdout

int
main(int argc, char **argv) {
  char           *fName = 0L;
  char           *mName = 0L;
  u32bit          merSize = 100;
  bool            show = false;

  int arg = 1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-fasta") == 0) {
      fName = argv[++arg];
    } else if (strcmp(argv[arg], "-markfile") == 0) {
      mName = argv[++arg];
    } else if (strcmp(argv[arg], "-mersize") == 0) {
      merSize = strtou32bit(argv[++arg], 0L);
    } else if (strcmp(argv[arg], "-show") == 0) {
      show = true;
    } else {
      fprintf(stderr, "ERROR: unknown option '%s'\n", argv[arg]);
      exit(1);
    }

    arg++;
  }

  if ((fName == 0L) || (mName == 0L)) {
    fprintf(stderr, "usage: %s -fasta f.fasta -markfile marks [-mersize m] [-show]\n", argv[0]);
    exit(1);
  }

  FastAstream  *FS = new FastAstream(fName);
  merStream    *MS = new merStream(merSize, FS);

  markFile   *distinct = new markFile(mName);

  u64bit      nMers = distinct->numberOfMers();
  u64bit      nMark = distinct->numberMarked();

  u64bit      nMersActual = 0;
  u64bit      nMarkActual = 0;

  u64bit      merNumber = 0;

  char        merstring[1024];

  while (MS->nextMer()) {
    u64bit  dist = distinct->get(merNumber);

    if (dist)
      nMarkActual++;

    if (show) {
      if (dist)
        fprintf(stdout, "DIST\t%s\n", MS->theFMer().merToString(merstring));
      else
        fprintf(stdout, "dupl\t%s\n", MS->theFMer().merToString(merstring));
    }

    nMersActual++;
    merNumber++;
  }

  fprintf(stderr, "merNumber: "u64bitFMT"\n", merNumber);
  fprintf(stderr, "nMers: Found "u64bitFMT", file reported "u64bitFMT"\n", nMersActual, nMers);
  fprintf(stderr, "nMark: Found "u64bitFMT", file reported "u64bitFMT"\n", nMarkActual, nMark);

  return(0);
}
