#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

#define TAG_LEN_MAX   32

#include "tapperTag.H"

//  Convert reads from ASCI to tapper binary.

int
main(int argc, char **argv) {
  char  *tagseq  = 0L;
  char  *tagqlt  = 0L;
  char  *tagout  = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-tags") == 0) {
      tagout  = argv[++arg];
    } else if (strcmp(argv[arg], "-tagseq") == 0) {
      tagseq  = argv[++arg];
    } else if (strcmp(argv[arg], "-tagqlt") == 0) {
      tagqlt  = argv[++arg];
    } else {
      err++;
    }
    arg++;
  }
  if ((err) || (tagout == 0L) || (tagseq == 0L) || (tagqlt == 0L)) {
    fprintf(stderr, "usage: %s -tagseq xx.csfasta -tagqlt xx.qual -tagout xx.tapperTags\n", argv[0]);
    exit(1);
  }

  FILE *seq = fopen(tagseq, "r");
  FILE *qlt = fopen(tagseq, "r");

  tapperTag      *TT = new tapperTag;
  tapperTagFile  *TF = new tapperTagFile(tagout);
  speedCounter   *CT = new speedCounter(" %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

  char    seqhdr[1024];
  char    seqseq[1024];
  char    qlthdr[1024];
  char    qltseq[1024];
  u64bit  qltnum[1024];

  u64bit  UID = 1;

  splitToWords  S;

  while (!feof(seq) && !feof(qlt)) {
    fgets(seqhdr, 1024, seq);  chomp(seqhdr);
    fgets(seqseq, 1024, seq);  chomp(seqseq);

    fgets(qlthdr, 1024, qlt);  chomp(qlthdr);
    fgets(qltseq, 1024, qlt);  chomp(qltseq);

    S.split(qltseq);

    if (strcmp(seqhdr, qlthdr) != 0)
      fprintf(stderr, "WARNING:  Got seq '%s' and qlt '%s'\n", seqhdr, qlthdr);

    for (u32bit i=0; i<S.numWords(); i++)
      qltnum[i] = atoi(S[i]);

    TT->encode(UID++, seqseq, qltnum);
    TF->put(TT);

    CT->tick();
  }

  delete CT;
  delete TF;
  delete TT;

  fclose(seq);
  fclose(qlt);
}
