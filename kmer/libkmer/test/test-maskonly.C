#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

//  Tests a positionDB when using an existDB for masking.
//
//  existDB can be either include or exclude
//  positionDB can use include, exclude or threshold
//

#define MERSIZE 14

int
main(int argc, char **argv) {
  existDB      *include;
  existDB      *exclude;
  positionDB   *full;
  positionDB   *incl;
  positionDB   *excl;
  positionDB   *thrs;

  if (argc != 4) {
    fprintf(stderr, "usage: %s seq.fasta mask.fasta incl.fasta\n", argv[0]);
    exit(1);
  }

  char *seqName = argv[1];
  char *mskName = argv[2];
  char *incName = argv[3];

  fprintf(stderr, "BUILDING EXCLUDE\n");
  exclude = new existDB(mskName, MERSIZE, existDBnoFlags, u32bitZERO, ~u32bitZERO);

  fprintf(stderr, "BUILDING INCLUDE\n");
  include = new existDB(incName, MERSIZE, existDBnoFlags, u32bitZERO, ~u32bitZERO);

  seqStream *F = new seqStream(seqName, true);
  merStream *T = new merStream(new kMerBuilder(MERSIZE), F);

  fprintf(stderr, "BUILDING FULL\n");
  full = new positionDB(T, MERSIZE, 0,      0L,      0L, 0L, 0, 0, true);
  full->saveState("junk-full");
  delete full;

  fprintf(stderr, "BUILDING INCL\n");
  incl = new positionDB(T, MERSIZE, 0,      0L, include, 0L, 0, 0, true);
  incl->saveState("junk-incl");
  delete incl;

  fprintf(stderr, "BUILDING EXCL\n");
  excl = new positionDB(T, MERSIZE, 0, exclude,      0L, 0L, 0, 0, true);
  excl->saveState("junk-excl");
  delete excl;

  fprintf(stderr, "BUILDING THRS\n");
  thrs = new positionDB(T, MERSIZE, 0,      0L,      0L, 0L, 1, 0, true);
  thrs->saveState("junk-thrs");
  delete thrs;

  full = new positionDB("junk-full");
  incl = new positionDB("junk-incl");
  excl = new positionDB("junk-excl");
  thrs = new positionDB("junk-thrs");

  char    themer[1000];
  u32bit  mernum = 0;

  u32bit  err = 0;

  //  Check everything looks ok
  T->rewind();
  while (T->nextMer()) {

    if (!full->existsExact(T->theFMer())) {
      fprintf(stderr, "Didn't find mer "u32bitFMT" %s in full.\n", mernum, T->theFMer().merToString(themer));
      err++;
    }

    if (include->exists(T->theFMer())) {
      if (!incl->existsExact(T->theFMer())) {
        fprintf(stderr, "Didn't find mer "u32bitFMT" %s in incl.\n", mernum, T->theFMer().merToString(themer));
        err++;
      }
    } else {
      if (incl->existsExact(T->theFMer())) {
        fprintf(stderr, "Found extra mer "u32bitFMT" %s in incl.\n", mernum, T->theFMer().merToString(themer));
        err++;
      }
    }

    if (exclude->exists(T->theFMer())) {
      if (excl->existsExact(T->theFMer())) {
        fprintf(stderr, "Found extra mer "u32bitFMT" %s in excl.\n", mernum, T->theFMer().merToString(themer));
        err++;
      }
    } else {
      if (!excl->existsExact(T->theFMer())) {
        fprintf(stderr, "Didn't find mer "u32bitFMT" %s in excl.\n", mernum, T->theFMer().merToString(themer));
        err++;
      }
    }

    mernum++;
  }
  
  delete T;
  delete F;

  exit(err > 0);
}
