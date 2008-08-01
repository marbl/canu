#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

#define TAG_LEN_MAX   32

#include "tapperTag.H"

//  Convert reads from ASCI to tapper binary.
//
//  ASSUMPTIONS
//
//  1) User is smart enough to give the correct set of mated files.
//  Code doesn't check that an F tag goes with an R tag, just that the
//  tag coordinates agree.  It is possible to mate an F to an F if the
//  wrong inputs are given.
//
//  2) Tag coords are 16-bit integers.  File UIDs are 16-bit integers.
//


int
tapperTagCompare(const void *a, const void *b) {
  tapperTag const *A = (tapperTag const *)a;
  tapperTag const *B = (tapperTag const *)b;
  if (A->tagID()  < B->tagID()) return(-1);
  return(A->tagID() != B->tagID());
}


bool
readTag(u64bit fileUID, FILE *seq, FILE *qlt, tapperTag *T) {
  static char    seqhdr[1024];
  static char    seqseq[1024];
  static char    qlthdr[1024];
  static char    qltseq[1024];
  static u64bit  qltnum[1024];
  static splitToWords  S;

  seqhdr[0] = 0;
  seqseq[0] = 0;
  qlthdr[0] = 0;
  qltseq[0] = 0;

  if (feof(seq) || feof(qlt))
    return(false);

  fgets(seqhdr, 1024, seq);
  while (seqhdr[0] == '#')
    fgets(seqhdr, 1024, seq);
  fgets(seqseq, 1024, seq);

  fgets(qlthdr, 1024, qlt);
  while (qlthdr[0] == '#')
    fgets(qlthdr, 1024, qlt);
  fgets(qltseq, 1024, qlt);

  if ((seqhdr[0] == 0) || (qlthdr[0] == 0))
    return(false);

  chomp(seqhdr);
  chomp(seqseq);
  chomp(qlthdr);
  chomp(qltseq);

  if (strcmp(seqhdr, qlthdr) != 0)
    fprintf(stderr, "WARNING:  Got seq '%s' and qlt '%s'\n", seqhdr, qlthdr);

  //  Assumes the header is >461_28_1918_F3
  //  -- copies it to the left by one to remove the >
  //  -- the loop below doesn't move the zero-terminator
  //  -- resulting string is "461 28 1918 F33"
  //
  for (u32bit i=1; seqhdr[i]; i++) {
    if (seqhdr[i] == '_')
      seqhdr[i] = ' ';
    seqhdr[i-1] = seqhdr[i];
  }

  S.split(seqhdr);

  u64bit UID;

  UID   = fileUID;
  UID <<= 16;
  UID  |= strtou64bit(S[0], 0L) & u64bitMASK(16);
  UID <<= 16;
  UID  |= strtou64bit(S[1], 0L) & u64bitMASK(16);
  UID <<= 16;
  UID  |= strtou64bit(S[2], 0L) & u64bitMASK(16);

  S.split(qltseq);

#warning quality values are getting fudged here
  for (u32bit i=0; i<S.numWords(); i++) {
    qltnum[i] = (S[i][0] == '-') ? 0 : strtou64bit(S[i], 0L);
    if (qltnum[i] > 31)
      qltnum[i] = 31;
  }

  T->encode(UID, seqseq, qltnum);

#define TEST_ENCODING
#ifdef TEST_ENCODING
  {
    char    seqtst[1024];
    u64bit  qlttst[1024];
    u64bit  tst = T->decode(seqtst, qlttst);
    u32bit  len = strlen(seqtst);
    u32bit  fail = 0;
    u64bit  qltsum=0, tstsum=0;

#if 0
    //  We don't encode QV precisely
    for (u32bit l=0; l<len; l++) {
      qltsum += qltnum[l];
      tstsum += qlttst[l];
      if ((seqseq[l] != seqtst[l]) || (qltnum[l] != qlttst[l]))
        fail++;
    }
#endif

    if ((tst != UID) || (fail)) {
      fprintf(stderr, "FAIL:  ("u64bitHEX",%s,"u64bitFMT") != ("u64bitHEX",%s,"u64bitFMT")\n",
              UID, seqseq, qltsum,
              tst, seqtst, tstsum);
      for (u32bit l=0; l<len; l++)
        fprintf(stderr, "  %2d -- "u64bitFMT" "u64bitFMT"\n", l, qltnum[l], qlttst[l]);
    }
  }

#endif

  return(true);
}



void
dumpTagFile(char *tagfile) {
  tapperTagFile  *TF = new tapperTagFile(tagfile);
  tapperTag       a, b;
  u64bit          ida, idb;
  char            seqa[265], seqb[256];
  char            quaa[256], quab[256];
  u64bit          qvsa[256], qvsb[256];
  u32bit          i;

  if (TF->metaData()->isPairedTagFile()) {
    while (TF->get(&a, &b)) {
      ida = a.decode(seqa, qvsa);
      idb = b.decode(seqb, qvsb);
      for (i=0; seqa[i+1]; i++)
        quaa[i] = qvsa[i] + '0';
      for (i=0; seqb[i+1]; i++)
        quab[i] = qvsb[i] + '0';
      fprintf(stdout, u64bitHEX"\t%s/%s\t"u64bitHEX"\t%s/%s\n",
              ida, seqa, quaa, idb, seqb, quab);
    }
  } else {
    while (TF->get(&a)) {
      ida = a.decode(seqa, qvsa);
      for (i=0; seqa[i+1]; i++)
        quaa[i] = qvsa[i] + '0';
      fprintf(stdout, u64bitHEX"\t%s/%s\n",
              ida, seqa, quaa);
    }
  }

  delete TF;
}



int
main(int argc, char **argv) {
  char  *prefix  = 0L;

  u64bit  tagfuid = 0,   tagruid = 0;
  char   *tagfseq = 0L, *tagrseq  = 0L;
  char   *tagfqlt = 0L, *tagrqlt  = 0L;

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-tagout") == 0) {
      prefix   = argv[++arg];

    } else if (strcmp(argv[arg], "-tags") == 0) {
      tagfuid  = strtou64bit(argv[++arg], 0L);
      tagfseq  = argv[++arg];
      tagfqlt  = argv[++arg];

    } else if (strcmp(argv[arg], "-ftags") == 0) {
      tagfuid  = strtou64bit(argv[++arg], 0L);
      tagfseq  = argv[++arg];
      tagfqlt  = argv[++arg];
    } else if (strcmp(argv[arg], "-rtags") == 0) {
      tagruid  = strtou64bit(argv[++arg], 0L);
      tagrseq  = argv[++arg];
      tagrqlt  = argv[++arg];

    } else if (strcmp(argv[arg], "-dump") == 0) {
      dumpTagFile(argv[++arg]);
      exit(0);

    } else {
      err++;
    }
    arg++;
  }
  if ((tagfseq == 0L) || (tagfqlt == 0L))  err++;
  if ((tagfseq != 0L) && (tagfqlt == 0L))  err++;
  if ((tagfseq == 0L) && (tagfqlt != 0L))  err++;
  if ((err) || (prefix == 0L)) {
    fprintf(stderr, "usage: %s -tagout prefix  -tags fileUID xx.csfasta xx.qual\n", argv[0]);
    fprintf(stderr, "usage: %s -tagout prefix -ftags fileUID ff.csfasta ff.qual -rtags fileUID rr.csfasta rr.qual\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "unmated tags will be placed in 'prefix.frag.tapperTags'\n");
    fprintf(stderr, "  mated tags will be placed in 'prefix.mate.tapperTags'\n");
    exit(1);
  }

  u64bit  numTagsF = 0, maxTagsF = 0;
  u64bit  numTagsR = 0, maxTagsR = 0;

  tapperTag      *TF = 0L;
  tapperTag      *TR = 0L;

  //
  //  Suck in all the F tags.
  //
  if (tagfseq) {
    FILE *fseq = fopen(tagfseq, "r");
    FILE *fqlt = fopen(tagfqlt, "r");

    speedCounter *CT = new speedCounter(" reading F tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

    maxTagsF = sizeOfFile(tagfseq) / 44 + 1000000;
    TF       = new tapperTag [maxTagsF];

    while (readTag(tagfuid, fseq, fqlt, TF + numTagsF)) {
      numTagsF++;
      if (numTagsF >= maxTagsF)
        fprintf(stderr, "Too many F tags.  Boom.\n"), exit(1);
      CT->tick();
    }
    delete CT;

    fclose(fseq);
    fclose(fqlt);

    maxTagsF = numTagsF;
    numTagsF = 0;
  }

  //
  //  Suck in all the R tags.
  //
  if (tagrseq) {
    FILE *rseq = fopen(tagrseq, "r");
    FILE *rqlt = fopen(tagrqlt, "r");

    speedCounter *CT = new speedCounter(" reading R tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

    maxTagsR = sizeOfFile(tagrseq) / 44 + 1000000;
    TR       = new tapperTag [maxTagsR];;

    while (readTag(tagruid, rseq, rqlt, TR + numTagsR)) {
      numTagsR++;
      if (numTagsR >= maxTagsR)
        fprintf(stderr, "Too many R tags.  Boom.\n"), exit(1);
      CT->tick();
    }
    delete CT;

    fclose(rseq);
    fclose(rqlt);

    maxTagsR = numTagsR;
    numTagsR = 0;
  }

  //
  //  Sort them.
  //
  qsort_mt(TF, maxTagsF, sizeof(tapperTag), tapperTagCompare, 4, 4 * 1024 * 1024);
  qsort_mt(TR, maxTagsR, sizeof(tapperTag), tapperTagCompare, 4, 4 * 1024 * 1024);

  //
  //  Merge to find pairs, output.
  //
  char            fragout[FILENAME_MAX];
  char            mateout[FILENAME_MAX];

  sprintf(fragout, "%s.frag.tapperTags", prefix);
  sprintf(mateout, "%s.mate.tapperTags", prefix);

  tapperTagFile  *TOfrag = 0L;
  tapperTagFile  *TOmate = 0L;

  speedCounter *CF = new speedCounter(" writing frag tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);
  speedCounter *CM = new speedCounter(" writing mate tags %7.0f sequences -- %5.0f sequences/second\r", 1.0, 0x1ffff, true);

  while ((numTagsF < maxTagsF) && (numTagsR < maxTagsR)) {
    u64bit   fID = TF[numTagsF].tagID() & u64bitMASK(48);
    u64bit   rID = TR[numTagsR].tagID() & u64bitMASK(48);

    if (fID == rID) {
      if (TOmate == 0L)
        TOmate = new tapperTagFile(mateout);
      TOmate->put(TF + numTagsF, TR + numTagsR);
      numTagsF++;
      numTagsR++;
      CM->tick();
    } else if (fID < rID) {
      if (TOfrag == 0L)
        TOfrag = new tapperTagFile(fragout);
      TOfrag->put(TF + numTagsF);
      numTagsF++;
      CF->tick();
    } else {
      if (TOfrag == 0L)
        TOfrag = new tapperTagFile(fragout);
      TOfrag->put(TR + numTagsR);
      numTagsR++;
      CF->tick();
    }
  }
  while (numTagsF < maxTagsF) {
    if (TOfrag == 0L)
      TOfrag = new tapperTagFile(fragout);
    TOfrag->put(TF + numTagsF);
    numTagsF++;
    CF->tick();
  }
  while (numTagsR < maxTagsR) {
    if (TOfrag == 0L)
      TOfrag = new tapperTagFile(fragout);
    TOfrag->put(TR + numTagsR);
    numTagsR++;
    CF->tick();
  }

  delete CF;
  delete CM;

  delete TOmate;
  delete TOfrag;

  delete [] TR;
  delete [] TF;
}
