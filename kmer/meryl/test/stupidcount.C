#include "bio++.H"

//  Reads a sequence file, outputs a list of the mers in it.  You can
//  then pipe this to unix sort and uniq to do a mercount.  You
//  probably don't want to count large things this way...

int
main(int argc, char **argv) {
  char    *seqName = 0L;
  u32bit   merSize = 20;

  int arg=1;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-s") == 0) {
      seqName = argv[++arg];
    } else if (strcmp(argv[arg], "-m") == 0) {
      merSize = strtou32bit(argv[++arg], 0L);
    }
    arg++;
  }

  if (seqName == 0L) {
    fprintf(stderr, "usage: %s [-m mersize] -s seqfile.fasta\n", argv[0]);
    exit(1);
  }

  seqStream       *CS = new seqStream(seqName, true);
  merStream       *MS = new merStream(new kMerBuilder(merSize), CS);
  char             str[1024];

  while (MS->nextMer())
    fprintf(stdout, "%s\n", MS->theFMer().merToString(str));

  delete MS;
  delete CS;

  exit(0);
}
