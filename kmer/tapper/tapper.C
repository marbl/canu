#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"


#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperHit.H"
#include "tapperComputation.H"

//  To make this work for single reads, add a flag to global that says
//  we're doing single reads.  Modify tapperReader to read a single
//  read.
//
//  Use a whole new tapperWorker.  It's simple compared to the mated
//  one (basically done now).
//
//  Tapper output hopefully will work for both single and paired
//  reads.


void*
tapperReaderMated(void *G) {
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


void*
tapperReaderSingle(void *G) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = 0L;

 again:
  seqInCore *a = g->QF->getSequenceInCore();

  if (a == 0L)
    return(0L);

#warning need to get rid of duplicate a
  s = new tapperComputation(a, a);

  delete a;

  if (s->tag1size != g->tagSize) {
    fprintf(stderr, "tag size mismatch for '%s' ("u32bitFMT"); expected "u32bitFMT".\n",
            a->header(), s->tag1size, g->tagSize);
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



//  The big value of this function is to convert from a chained
//  position to a (seqID,pos), and save the hit onto a bitpacked list.
//  This list can then be numerically sorted to order all hits.  Of
//  course, we could have just sorted the original chained positions.
//
//  It saves (seqID,pos,isTag2,isReverse)
//
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



//  Compose the colors from beg to end.
char
composeColors(char *colors, u32bit beg, u32bit end) {
  char c = colors[beg];

  for (u32bit x=beg; x<end; x++)
    c = baseToColor[c][colors[x]];

  return(c);
}



//  Reencode the genomic sequence, using the same reference start as the tag.
//
//  tag          T3 1 0 2 1 1   encoded seeding with T
//  tag          A C C T G T    decoded
//
//  genomic A 3 2 2 1 0 2 3 1   encoded seeding with A, and at the wrong spot, notice the tag matches
//  genomic  T C T G G A T G C  real bases, we just have this around
//  genomic      T1 0 2 3 1 3   reencode, seeding with a T at the same start as the tag
//
//  After reencoding, there is no agreement with the colorspace tag.
//
//  genomic  T C A C C T G T C  the correct sequence
//  genomic   2 1 1 0 2 1 1 2   in color, note the first 3 in the tag doesn't match
//  genomic     TA C C T G T C  reencoding with the first tag letter
//  genomic     T31 0 2 1 1 2   now a match
//
//  Thus, if we locally reencode the genomic, we should end up
//  with a color space sequence with the allowed number of
//  errors.  At this point we can correct any color encoding
//  issues (isolated changes) or detect any real base changes.
//
//  Count the number of differences between the tag and the reencoded genomic
//
u32bit
alignToReference(tapperGlobalData *g,
                 u64bit  so,
                 u64bit  po,
                 char   *tag, char *tagACGT, u32bit len,
                 char   *ref, char *refACGT) {

  char       *seq       = g->GS->getSequenceInCore(so)->sequence();
  u32bit      ti = 0;
  u32bit      si = po;
  u32bit      errs = 0;  //  number of errors
  u32bit      errp[64];  //  location of the errors

  //  Analyze color[] and t[], correct differences, call base changes.
  //  Return an ACGT sequence for the tag.

  //  Compose the colors together.  At points where the compositions
  //  disagree, the base at that point is different.  The composition
  //  tells us how to transform the reference letter to the base at
  //  this position, in one step.
  //
  //  If our final composed value is different, then either we end on
  //  a SNP, or we have an error somewhere.  The choice here is
  //  arbitrary, and made depending on where that error is.

  ref[0] = tag[0];
  ref[1] = baseToColor[ref[0]][seq[po]];

  if (ref[1] != tag[1]) {
    errp[errs] = 1;
    errs++;
  }

  for (ti=2, si=po+1; ti<len; ti++, si++) {
    ref[ti] = baseToColor[seq[si-1]][seq[si]];

    if (tag[ti] != ref[ti]) {
      errp[errs] = ti;
      errs++;
    }
  }

  //  Note that errp[] is actaully 1-based; the first position is never
  //  an error; it's the reference base.

  if ((errs == 0) || (g->maxError < errs))
    return(errs);

  //  Compose the colors.  If these are the same, we call the two
  //  colors 'consistent' - the ACGT reads begin and end with the same
  //  base.
  //
  u32bit  tagc = composeColors(tag, 1, len);
  u32bit  refc = composeColors(ref, 1, len);

  //  Copy the tag colors so we can correct single color flaws to
  //  generate the final ACGT.
  //
  char        cor[64];

  strcpy(cor, tag);



  if (errs == 1) {

    //  len-1 because we lose one letter when converting from color to acgt.
    //  errp[] is 1-based, so == len-1 is the last base in the tag.
    //
    if ((errp[0]) == (len - 1)) {
      //  Exactly at the end, believe it's the start of a SNP.  For a
      //  SNP to be at the start, we'll need two color changes.  One
      //  color change at the start is a color error.  The end is
      //  different, because we just don't know what the next color
      //  would have been.  The beginning has the reference base to
      //  anchor it.
      //
    } else {
      //  Correct it.
      cor[errp[0]] = ref[errp[0]];
    }
  }


  //  If there are two errors, and the whole string is consistent,
  //  we'll call it a single mismatch block, if it is small.
  //  Otherwise, we'll fix the two errors.
  //
  if (errs == 2) {
    if ((tagc == refc) && (errp[1] - errp[0] < 4)) {
      //  MNP of size 4.
    } else if (tagc == refc) {
      //  Big MNP?  What to do?  Correct, I guess.
      cor[errp[0]] = ref[errp[0]];
      cor[errp[1]] = ref[errp[1]];
    } else {
      //  Correct 'em.
      cor[errp[0]] = ref[errp[0]];
      cor[errp[1]] = ref[errp[1]];
    }
  }


  //  If there are three errors...
  //
  //  String consistent, and MNP block is small, let it through.
  //
  //  String not consistent, we need to see if the first two or last
  //  two errors are consistent.  If so, declare the third error is a
  //  color error and correct it.
  //
  if (errs == 3) {
    if        ((tagc == refc) && (errp[2] - errp[0] < 5)) {
      //  MNP of size 5
    } else if (tagc == refc) {
      //  Big MNP?  What to do?  Correct, I Guess.
      cor[errp[0]] = ref[errp[0]];
      cor[errp[1]] = ref[errp[1]];
      cor[errp[2]] = ref[errp[2]];
    } else if (composeColors(tag, 0, errp[2]) == composeColors(ref, 0, errp[2])) {
      //  First two ok, fix the third.
      cor[errp[2]] = ref[errp[2]];
    } else if (composeColors(tag, errp[0]+1, len) == composeColors(ref, errp[0]+1, len)) {
      //  Last two ok, fix the first.
      cor[errp[0]] = ref[errp[0]];
    } else {
      //  Nothing consistent, fix all of 'em.
      cor[errp[0]] = ref[errp[0]];
      cor[errp[1]] = ref[errp[1]];
      cor[errp[2]] = ref[errp[2]];
    }
  }


  if (errs == 4) {
    fprintf(stderr, "Four errors detected.  Code doesn't know what to do.\n");
    exit(1);
  }

  tagACGT[0] = baseToColor[cor[0]][cor[1]];
  refACGT[0] = baseToColor[ref[0]][ref[1]];
  for (ti=1, si=2; si<len; ti++, si++) {
    tagACGT[ti] = baseToColor[tagACGT[ti-1]][cor[si]];
    refACGT[ti] = baseToColor[refACGT[ti-1]][ref[si]];
    if (tagACGT[ti] != refACGT[ti]) {
      tagACGT[ti] = toUpper[tagACGT[ti]];
      refACGT[ti] = toUpper[refACGT[ti]];
    }
  }
  tagACGT[len-1] = 0;
  refACGT[len-1] = 0;

  fprintf(stderr, "tag:%s/%s\nref:%s/%s\n", tag, tagACGT, ref, refACGT);

  
  return(errs);
}



void
tapperWorkerMated(void *G, void *T, void *S) {
#if 0
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

  if (t->posn1fLen + t->posn1rLen + t->posn2fLen + t->posn2rLen == 0)
    return;

  if (t->hits1Max < t->posn1fLen + t->posn1rLen) {
    t->hits1Max = t->posn1fLen + t->posn1rLen;
    delete [] t->hits1;
    t->hits1 = new u64bit [t->posn1fLen + t->posn1rLen];
  }
  t->hits1Len = 0;
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1f, t->posn1fLen, g, false, false);
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1r, t->posn1rLen, g, false, true);

  if (t->hits2Max < t->posn2fLen + t->posn2rLen) {
    t->hits2Max = t->posn2fLen + t->posn2rLen;
    delete [] t->hits2;
    t->hits2 = new u64bit [t->posn2fLen + t->posn2rLen];
  }
  t->hits2Len = 0;
  t->hits2Len = tapperWorker_addHits(t->hits2, t->hits2Len, t->posn2f, t->posn2fLen, g, true,  false);
  t->hits2Len = tapperWorker_addHits(t->hits2, t->hits2Len, t->posn2r, t->posn2rLen, g, true,  true);

#warning score hits
  //  Then score hits individually, throwing out any that are junk.
  //  Probably should do that before we addHits?

  //  Sort so we can do mate rescue.
#warning use a real sort
  qsort(t->hits1, t->hits1Len, sizeof(u64bit), u64bitcompare);
  qsort(t->hits2, t->hits2Len, sizeof(u64bit), u64bitcompare);
#endif
}


void
tapperWorkerSingle(void *G, void *T, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperThreadData  *t = (tapperThreadData  *)T;
  tapperComputation *s = (tapperComputation *)S;

  t->posn1fLen = 0;
  t->posn1rLen = 0;

  g->PS->getUpToNMismatches(s->tag1f, g->maxError, t->posn1f, t->posn1fMax, t->posn1fLen);
  g->PS->getUpToNMismatches(s->tag1r, g->maxError, t->posn1r, t->posn1rMax, t->posn1rLen);

  if (t->posn1fLen + t->posn1rLen == 0)
    return;

  if (t->hits1Max < t->posn1fLen + t->posn1rLen) {
    t->hits1Max = t->posn1fLen + t->posn1rLen;
    delete [] t->hits1;
    t->hits1 = new u64bit [t->posn1fLen + t->posn1rLen];
  }
  t->hits1Len = 0;
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1f, t->posn1fLen, g, false, false);
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1r, t->posn1rLen, g, false, true);

  for (u32bit i=0; i<t->hits1Len; i++) {
    tapperHit  h   = {0};
    u64bit     so  = (t->hits1[i] >> 33) & u64bitMASK(31);  //  Sequence it hits
    u64bit     po  = (t->hits1[i] >> 2)  & u64bitMASK(31);  //  Position in sequence
    u64bit     ta2 = (t->hits1[i] >> 1)  & 0x01;            //  isTag2
    u64bit     rev = (t->hits1[i] >> 0)  & 0x01;            //  isReversed

    //  We're ignoring the first base, so adjust the position to
    //  include it.
    po--;

    //
    //  Do something useful with the hits, like verify them against the ACGT reference.
    //

    u32bit mm = 0;
    char   tagACGT[64];
    char   refACGT[64];

    if (ta2) {
      mm = alignToReference(g, so, po, (rev) ? s->tag2rseq : s->tag2fseq, tagACGT, s->tag2size+2, s->refseq, refACGT);
    } else {
      mm = alignToReference(g, so, po, (rev) ? s->tag1rseq : s->tag1fseq, tagACGT, s->tag1size+2, s->refseq, refACGT);
    }

    //fprintf(stderr, "hit: "u64bitFMT" @ "u64bitFMT" tag2:"u64bitFMT" rev:"u64bitFMT"\n", so, po, ta2, rev);

    if (mm <= g->maxError) {
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
  }
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

    } else if (strcmp(argv[arg], "-mated") == 0) {
      g->isMated = true;

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

  sweatShop *ss = new sweatShop((g->isMated) ? tapperReaderMated : tapperReaderSingle,
                                (g->isMated) ? tapperWorkerMated : tapperWorkerSingle,
                                tapperWriter);

  ss->loaderQueueSize(10240);
  ss->writerQueueSize(10240);

  ss->numberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new tapperThreadData(g));

  ss->run(g, g->beVerbose);

  delete g;
  delete ss;

  fprintf(stderr, "\nSuccess!  Bye.\n");
  return(0);
}
