
#include "bio++.H"

//  Reads a fastafile, samples randomly and outputs colorspace tags.

#define MS 25

int
main(int argc, char **argv) {

  if (argv[1] == 0L) {
    fprintf(stderr, "usage: %s N some.fasta\n", argv[0]);
    exit(1);
  }

  int         N = atoi(argv[1]);

  seqFile    *F = openSeqFile(argv[2]);
  seqInCore  *s = F->getSequenceInCore();

  u32bit  pos = 0;
  u32bit  len = s->sequenceLength();

  char    seq[1024] = {0};

  mt_s   *mtctx = mtInit(time(0));

  for (u32bit i=0; i<N; i++) {
    pos = mtRandom32(mtctx) % (len - MS);

    char  n = 'N';
    char  l = 'N';

#if 1
    for (u32bit x=0; x<MS; x++) {
      n = s->sequence()[pos++];
      seq[x] = baseToColor[l][n];
      l = n;
    }
#else
    pos += MS;
    for (u32bit x=0; x<MS; x++) {
      n = s->sequence()[pos--];
      seq[x] = baseToColor[l][n];
      l = n;
    }
#endif

    //  Insert an error.

    if (0) {
      u32bit e = mtRandom32(mtctx) % MS;
      seq[e] = seq[e] + 1;
      if (seq[e] > '4')
        seq[e] = '0';
    }

    //seq[10] = 3 - seq[10];

    fprintf(stdout, ">idx"u32bitFMT"pos"u32bitFMT"\nN%s\n", i, pos - MS, seq);
  }
  
}
