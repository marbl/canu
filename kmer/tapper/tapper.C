#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"


#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperHit.H"
#include "tapperComputation.H"


void*
tapperReader(void *G) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = 0L;

 again:
  seqInCore *a = g->QF->getSequenceInCore();
  seqInCore *b = g->QF->getSequenceInCore();

  if ((a == 0L) || (b == 0L))
    return(0L);

  s = new tapperComputation(a, b);

  delete a;
  delete b;

  if ((s->tag1size != s->tag2size) ||
      (s->tag1size != g->tagSize) ||
      (s->tag2size != g->tagSize)) {
    fprintf(stderr, "tag size mismatch for '%s' ("u32bitFMT") and '%s' ("u32bitFMT"); expected "u32bitFMT".\n",
            a->header(), s->tag1size, b->header(), s->tag2size, g->tagSize);
    delete s;
    goto again;
  }

  return(s);
}


void
tapperWriter(void *G, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = (tapperComputation *)S;

  s->writeHits(stdout);

  delete s;
}



inline
u32bit
tapperWorker_addHits(u64bit *hits, u32bit hitsLen, u64bit *posn, u64bit posnLen, tapperGlobalData *g, bool ta2, bool rev) {
  u64bit  ta2rev = u64bitZERO;

  if (ta2)
    ta2rev |= 0x02;
  if (rev)
    ta2rev |= 0x01;

  for (u32bit i=0; i<posnLen; i++) {
    u64bit  pos = posn[i];
    u64bit  seq = g->SS->sequenceNumberOfPosition(pos);

    pos -= g->SS->startOf(seq);
    seq  = g->SS->IIDOf(seq);

    //  seqIdx | positionInSeq | isTag2 | isReversed

    hits[hitsLen++] = (seq << 33) | (pos << 2) | ta2rev;
  }

  return(hitsLen);
}


int
u64bitcompare(void const *a, void const *b) {
  u64bit A = *((u64bit const *)a);
  u64bit B = *((u64bit const *)b);
  if (A < B)  return(-1);
  if (A > B)  return(1);
  return(0);
}


void
tapperWorker(void *G, void *T, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperThreadData  *t = (tapperThreadData  *)T;
  tapperComputation *s = (tapperComputation *)S;

  t->posn1fLen = 0;
  t->posn1rLen = 0;
  t->posn2fLen = 0;
  t->posn2rLen = 0;

  g->PS->getUpToNMismatches(s->tag1f, 3, t->posn1f, t->posn1fMax, t->posn1fLen);
  g->PS->getUpToNMismatches(s->tag1r, 3, t->posn1r, t->posn1rMax, t->posn1rLen);

  g->PS->getUpToNMismatches(s->tag2f, 3, t->posn2f, t->posn2fMax, t->posn2fLen);
  g->PS->getUpToNMismatches(s->tag2r, 3, t->posn2r, t->posn2rMax, t->posn2rLen);

  //  Now need to make sense of this stuff.
  //
  //  1.  convert from positions to seqId/pos
  //  2.  build list of: refId[31] refPos[31]  tag2?  rev?
  //  3.  sort numerically (and so by refId and refPos).
  //  4.  scan the list
  //      - report tags that map to any scaffold without their mate
  //      - rewrite the list to have scaffolds with mated tags

  u32bit  hitsLen = t->posn1fLen + t->posn1rLen + t->posn2fLen + t->posn2rLen;

  if (hitsLen == 0)
    return;

  fprintf(stderr, "Found "u64bitFMT","u64bitFMT" hits for %s and "u64bitFMT","u64bitFMT" for %s\n",
          t->posn1fLen, t->posn1rLen, s->tag1name,
          t->posn2fLen, t->posn2rLen, s->tag2name);

  u64bit *hits    = new u64bit [hitsLen];

  hitsLen = 0;
  hitsLen = tapperWorker_addHits(hits, hitsLen, t->posn1f, t->posn1fLen, g, false, false);
  hitsLen = tapperWorker_addHits(hits, hitsLen, t->posn1r, t->posn1rLen, g, false, true);
  hitsLen = tapperWorker_addHits(hits, hitsLen, t->posn2f, t->posn2fLen, g, true,  false);
  hitsLen = tapperWorker_addHits(hits, hitsLen, t->posn2r, t->posn2rLen, g, true,  true);

#warning use a real sort
  qsort(hits, hitsLen, sizeof(u64bit), u64bitcompare);

  for (u32bit i=0; i<hitsLen; i++) {
    u64bit  so  = (hits[i] >> 33) & u64bitMASK(31);
    u64bit  po  = (hits[i] >> 2)  & u64bitMASK(31);
    u64bit  ta2 = (hits[i] >> 1)  & 0x01;
    u64bit  rev = (hits[i] >> 0)  & 0x01;

    fprintf(stderr, "hit: "u64bitFMT" @ "u64bitFMT" tag2:"u64bitFMT" rev:"u64bitFMT"\n", so, po, ta2, rev);
  }

  //
  //  Do something useful with the hits, like verify them against the ACGT reference.
  //
  

  //  Stuff all the good hits into tappeComputation's hits entry.

  for (u32bit i=0; i<hitsLen; i++) {
    tapperHit  h = {0};

    u64bit  so  = (hits[i] >> 33) & u64bitMASK(31);
    u64bit  po  = (hits[i] >> 2)  & u64bitMASK(31);
    u64bit  ta2 = (hits[i] >> 1)  & 0x01;
    u64bit  rev = (hits[i] >> 0)  & 0x01;

    h.seqIdx     = so;

    if (ta2) {
      h.setOrientation(0, 0, 1, rev);
      h.tag2pos = po;
    } else {
      h.setOrientation(1, rev, 0, 0);
      h.tag1pos = po;
    }

    s->addHit(h);
  }

  delete [] hits;
}



int
main(int argc, char **argv) {
  tapperGlobalData  *g = new tapperGlobalData();

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-genomic") == 0) {
      g->genName = argv[++arg];
    } else if (strcmp(argv[arg], "-color") == 0) {
      g->colName = argv[++arg];
    } else if (strcmp(argv[arg], "-queries") == 0) {
      g->qryName = argv[++arg];

    } else if (strcmp(argv[arg], "-maxerror") == 0) {
      g->maxError = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-numthreads") == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strcmp(argv[arg], "-verbose") == 0) {
      g->beVerbose = true;
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err > 0) || (g->genName == 0L) || (g->colName == 0L) || (g->qryName == 0L)) {
    fprintf(stderr, "usage: %s -genomic g.fasta -color tmp.colors -queries q.fasta\n", argv[0]);
    exit(1);
  }

  g->initialize();

  sweatShop *ss = new sweatShop(tapperReader, tapperWorker, tapperWriter);

  ss->loaderQueueSize(10240);
  ss->writerQueueSize(10240);

  ss->numberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new tapperThreadData(g));

  ss->run(g, true);  //  true == verbose

  delete g;
  delete ss;

  fprintf(stderr, "\nSuccess!  Bye.\n");
  return(0);
}
