#include "bio++.H"

int
main(int argc, char **argv) {
  char     *filename = 0L;
  u32bit    windowsize = 1000;
  u32bit    skipsize   = 1000;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-f") == 0) {
      filename = argv[++arg];
    } else if (strcmp(argv[arg], "-w") == 0) {
      windowsize = atoi(argv[++arg]);
    } else if (strcmp(argv[arg], "-s") == 0) {
      skipsize = atoi(argv[++arg]);
    } else {
      fprintf(stderr, "unknown option '%s'\n", argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err) || (filename == 0L)) {
    fprintf(stderr, "usage: %s -f some.fasta\n", argv[0]);
    exit(1);
  }

  seqFile    *F = openSeqFile(filename);
  seqInCore  *S;

  while ((S = F->getSequenceInCore()) != 0L) {
    char       *seq = S->sequence();
    u32bit      len = S->sequenceLength();
    u32bit      GCt = 0;
    u32bit      GC  = 0;
    u32bit      i   = 0;

    for (i=0; i<windowsize; i++) {
      if ((seq[i] == 'G') || (seq[i] == 'g') ||
          (seq[i] == 'C') || (seq[i] == 'c')) {
        GC++;
        GCt++;
      }
    }

    while (i < len) {
      u32bit  avepos = (i - 1) - (windowsize - 1) / 2;
      if ((avepos % skipsize) == 0)
        fprintf(stdout, u32bitFMT"\t%.4f\n",
                avepos,
                100.0 * GC / (double)windowsize);

      if ((seq[i-windowsize] == 'G') || (seq[i-windowsize] == 'g') ||
          (seq[i-windowsize] == 'C') || (seq[i-windowsize] == 'c'))
        GC--;

      if ((seq[i] == 'G') || (seq[i] == 'g') ||
          (seq[i] == 'C') || (seq[i] == 'c')) {
        GC++;
        GCt++;
      }

      i++;
    }

    fprintf(stderr, "overall GC %.4f\n", 100.0 * GCt / (double)len);

    delete S;
  }

  delete F;
}
