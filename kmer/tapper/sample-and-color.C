
#include "bio++.H"

//  Reads a fastafile, samples randomly and outputs colorspace tags.

//  Max of 33 for now.
#define MS 23

int
main(int argc, char **argv) {

  if (argv[1] == 0L) {
    fprintf(stderr, "usage: %s N some.fasta\n", argv[0]);
    exit(1);
  }

  u32bit      N = atoi(argv[1]);

  seqFile    *F = openSeqFile(argv[2]);
  seqInCore  *s = F->getSequenceInCore();

  u32bit  pos = 0;
  u32bit  len = s->sequenceLength();

  char    cor[1024] = {0};
  char    seq[1024] = {0};

  char    acgt[4] = {'A', 'C', 'G', 'T'};

  mt_s   *mtctx = mtInit(time(0));

  for (u32bit i=0; i<N; i++) {
    pos = mtRandom32(mtctx) % (len - MS);

    char  n = acgt[mtRandom32(mtctx) % 4];
    char  l = n;

    cor[0] = n;
    seq[0] = n;

    bool   doForward = (mtRandom32(mtctx) & 0x1000) == 0x1000;

    doForward = false;

    if (doForward) {
      //  Forward
      u32bit sp = pos;
      for (u32bit x=1; x<=MS; x++) {
        n = s->sequence()[sp++];
        cor[x] = n;
        seq[x] = baseToColor[l][n];
        l = n;
      }
    } else {
      //  Reverse
      u32bit sp = pos + MS - 1;
      for (u32bit x=1; x<=MS; x++) {
        n = complementSymbol[s->sequence()[sp--]];
        cor[x] = n;
        seq[x] = baseToColor[l][n];
        l = n;
      }
    }

    //  Insert 0 to 3 errors.

    char     errors[256] = {0};
    char     errort[256] = {0};
    u32bit   nerr = 0;

    for (u32bit xx=0; xx<nerr; xx++) {
      u32bit e = mtRandom32(mtctx) % MS + 1;
      char   o = seq[e];
      seq[e] = seq[e] + 1;
      if (seq[e] > '3')
        seq[e] = '0';
      sprintf(errort, "_%c-%c@%d", o, seq[e], e);
      strcat(errors, errort);
    }

    //seq[ 1]++;   if (seq[ 1] > '3')   seq[ 1] = '0';
    //seq[33]++;   if (seq[33] > '3')   seq[33] = '0';
    seq[10]++;   if (seq[10] > '3')   seq[10] = '0';

    fprintf(stdout, ">i"u32bitFMT"_p"u32bitFMT"_%s%s\n%s\n", i, pos, cor+1, errors, seq);
  }
  
}
