#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "positionDB.H"
#include "fasta.H"

//  Driver for the positionDB creation.  Reads a sequence.fasta, builds
//  a positionDB for the mers in the file, and then writes the internal
//  structures to disk.
//
//  The positionDB constructor is smart enough to read either a pre-built
//  image or a regular multi-fasta file.


char const *usage =
"usage: %s [-mersize s] [-merskip k] [-tablesize t]\n"
"         default: s=20, k=0, t=auto\n"
"\n"
"       The default operation is to create an image:\n"
"         [build opts] <sequence.fasta> <datafile>\n"
"\n"
"       To dump information about an image:\n"
"         -dump datafile\n"
"\n"
"       To run sanity tests:\n"
"         -buildonly [build opts] sequence.fasta\n"
"           --  just builds a table and exits\n"
"\n"
"         -existence [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               in sequence.fasta cannot be found\n"
"\n"
"         -extra [build opts] sequence.fasta\n"
"           --  builds (or reads) a table reports if any mers\n"
"               NOT in sequence.fasta are be found\n"
"\n";

int
main(int argc, char **argv) {
  u32bit    mersize = 20;
  u32bit    merskip = 0;
  u32bit    tblsize = 0;

  if (argc < 3) {
fprintf(stderr, usage, argv[0]);
    exit(1);
  }

  int arg = 1;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-mersize", 6) == 0) {
      arg++;
      mersize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-merskip", 6) == 0) {
      arg++;
      merskip = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-tablesize", 2) == 0) {
      arg++;
      tblsize = atoi(argv[arg]);
    } else if (strncmp(argv[arg], "-dump", 2) == 0) {
      positionDB *e = new positionDB(argv[argc-1], false);
      e->printState(stdout);
      delete e;
      exit(0);
    }
    arg++;
  }

  fprintf(stderr, "Building table with merSize %u, merSkip %u and table size %u\n", mersize, merskip, tblsize);

  u64bit    seqlen  = 0;
  u32bit    ts      = 25;

  fprintf(stderr, "Opening '%s', finding length of sequences\n", argv[argc-2]);
  FastAWrapper  *F = new FastAWrapper(argv[argc-2]);
  F->openIndex();
  for (u32bit i=0; i<F->getNumberOfSequences(); i++)
    seqlen += F->sequenceLength(i);

  if (seqlen < 64 * 1024 * 1024) ts = 24;
  if (seqlen < 16 * 1024 * 1024) ts = 23;
  if (seqlen <  4 * 1024 * 1024) ts = 22;
  if (seqlen <  2 * 1024 * 1024) ts = 21;
  if (seqlen <  1 * 1024 * 1024) ts = 20;

  if (tblsize == 0)
    tblsize = ts;

  if (ts != tblsize)
    fprintf(stderr, "For sequence length %lu, suggested tableSize is %u, you want %u\n", seqlen, ts, tblsize);

  positionDB  *e = new positionDB(0L, argv[argc-2], mersize, merskip, tblsize, true);

  fprintf(stderr, "Writing positionDB state to '%s'\n", argv[argc-1]);
  e->saveState(argv[argc-1]);

  delete e;
  exit(0);
}
