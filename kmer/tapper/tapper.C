#include "bio++.H"
#include "existDB.H"
#include "positionDB.H"

#define TAG_LEN_MAX   32

#include "tapperTag.H"
#include "tapperGlobalData.H"
#include "tapperHit.H"
#include "tapperThreadData.H"
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
  return(0L);
}


void*
tapperReaderSingle(void *G) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = 0L;
  tapperTag          a;

  if (g->TF->get(&a))
    s = new tapperComputation(&a, 0L);

  return(s);
}





void
tapperWriter(void *G, void *S) {
  //tapperGlobalData  *g = (tapperGlobalData  *)G;
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

    //  Search ignores first letter, align needs it.  This makes for a
    //  very special case, 0, which isn't a full match.
    if (pos > 0) {
      pos--;
      hits[hitsLen++] = (seq << 33) | (pos << 2) | ta2rev;
    }
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
//
inline
char
composeColors(char *colors, u32bit beg, u32bit end) {
  char c = colors[beg];

  for (u32bit x=beg; x<end; x++)
    c = baseToColor[c][colors[x]];

  return(c);
}


//  Returns true if the the i and j errors result in a consistent
//  encoding, and they're not too far away.  Consistent in that the
//  sequence before agrees and the sequence after agrees.
//
inline
bool
isConsistent(char     *ref, char     *tag,
             u32bit    i,   u32bit    j) {
  return(composeColors(ref, i, j) == composeColors(tag, i, j));
}



//  Analyze tag[] and ref[], correct differences, call base changes.
//  Return an ACGT sequence for the tag.
//
//  Compose the colors together.  At points where the compositions
//  disagree, the base at that point is different.  The composition
//  tells us how to transform the reference letter to the base at
//  this position, in one step.
//
//  If our final composed value is different, then either we end on
//  a SNP, or we have an error somewhere.  The choice here is
//  arbitrary, and made depending on where that error is.
//
void
tapperHit::alignToReference(tapperGlobalData *g,
                            u32bit  so_in,
                            u32bit  po_in,
                            char   *tag_in, u32bit len_in) {

  //  This function is NOT a bottleneck.  Don't bother optimizing.

  u32bit      errs = 0;               //  number of errors
  u32bit      errp[TAG_LEN_MAX];      //  location of the errors
  char       _tagCOREC[TAG_LEN_MAX];  //  For holding corrected color calls, only to generate ACGT align

  _seqIdx = so_in;
  _seqPos = po_in;
  _tagIdx = 0;

  //  _rev: Yeah, we assume ASCII and UNIX newlines all over the
  //  place.  A forward read starts with a reference base; reverse
  //  reads have a number here.

  _pad                = 0;
  _len                = len_in;
  _rev                = (tag_in[0] < 'A') ? true : false;
  _rank               = 0x000000000000ffffllu;

  _basesMismatch      = len_in;  //  Set at end

  _colorMismatch      = 0;       //  Set when parsing errors
  _colorInconsistent  = 0;       //  Set when parsing errors

  _tagCOLOR[0] = 0;
  _tagCOREC[0] = 0;
  _refCOLOR[0] = 0;

  _tagACGT[0] = 0;
  _refACGT[0] = 0;


  //  Copy the tag.
  //
  //  A bit of devilish trickery to make a reverse read look like a
  //  forward read - we locally reverse the reference and read,
  //  process as if the reverse read is a forward read, then clean up
  //  at the end.
  //
  {
    strncpy(_tagCOLOR, tag_in, _len);

    if (_rev) {
      reverseString(_tagCOLOR, _len);
      _tagCOLOR[0] = complementSymbol[_tagCOLOR[0]];
    }

    strncpy(_tagCOREC, _tagCOLOR, _len);

    _tagCOLOR[_len] = 0;
    _tagCOREC[_len] = 0;
  }


  //  Copy the reference and convert the genomic sequence to
  //  color space using the reference base of the read.
  //
  {
    char     *seq = g->GS->getSequenceInCore(so_in)->sequence();

    strncpy(_refACGT, seq + po_in, _len-1);
    _refACGT[_len-1] = 0;

    if (_rev)
      reverseComplementSequence(_refACGT, _len - 1);

    _refCOLOR[0] = _tagCOLOR[0];  //  ALWAYS the reference encoding base, as long as we copy the tag first.
    _refCOLOR[1] = baseToColor[_refCOLOR[0]][_refACGT[0]];

    for (u32bit ti=2; ti<_len; ti++)
      _refCOLOR[ti] = baseToColor[_refACGT[ti-2]][_refACGT[ti-1]];

    _refCOLOR[_len] = 0;
  }

  //fprintf(stderr, "tag: %s ref: %s %s\n", _tagCOLOR, _refCOLOR, _refACGT);

  //  Count the number of color space errors
  //
  //  Note that errp[] is actaully 1-based; the first position is
  //  never an error; it's the reference base.

  for (u32bit ti=1; ti<_len; ti++) {
    if (_tagCOLOR[ti] != _refCOLOR[ti]) {
      errp[errs] = ti;
      errs++;
    }
  }

  //
  //  The following if blocks correct single color errors using very
  //  complicated rules.
  //


  if        (errs == 0) {
    _colorMismatch     = 0;
    _colorInconsistent = 0;

  } else if (errs == 1) {
    //  Always corrected, just to get an ACGT alignment.  We can't
    //  tell if the color mismatch is an error, or if the error is
    //  adjacent to the mismatch, which would have resulted in a valid
    //  SNP.
    _colorMismatch     = 0;
    _colorInconsistent = 1;
    _tagCOREC[errp[0]] = _refCOLOR[errp[0]];

  } else if (errs == 2) {
    bool  ok21  = isConsistent(_refCOLOR, _tagCOLOR, 1, _len) && (errp[1] - errp[0] < 4);

    if (ok21) {
      //  MNP of size 4.
      _colorMismatch     = 2;
      _colorInconsistent = 0;
    } else {
      //  Correct 'em.
      _colorMismatch     = 0;
      _colorInconsistent = 2;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
      _tagCOREC[errp[1]] = _refCOLOR[errp[1]];
    }

  } else if (errs == 3) {
    bool   ok21 = isConsistent(_refCOLOR, _tagCOLOR,         1, errp[2]) && (errp[1] - errp[0] < 4);
    bool   ok22 = isConsistent(_refCOLOR, _tagCOLOR, errp[0]+1, _len)    && (errp[2] - errp[1] < 4);

    bool   ok31 = isConsistent(_refCOLOR, _tagCOLOR,         1, _len)    && (errp[2] - errp[0] < 5);

    if        (ok31) {
      //  MNP of size 5
      _colorMismatch     = 3;
      _colorInconsistent = 0;
    } else if (ok21) {
      //  First two ok, fix the third.
      _colorMismatch     = 2;
      _colorInconsistent = 1;
      _tagCOREC[errp[2]] = _refCOLOR[errp[2]];
    } else if (ok22) {
      //  Last two ok, fix the first.
      _colorMismatch     = 2;
      _colorInconsistent = 1;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
    } else {
      //  Nothing consistent, fix all of 'em.
      _colorMismatch     = 0;
      _colorInconsistent = 3;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
      _tagCOREC[errp[1]] = _refCOLOR[errp[1]];
      _tagCOREC[errp[2]] = _refCOLOR[errp[2]];
    }

  } else if (errs == 4) {
    bool   ok21 = isConsistent(_refCOLOR, _tagCOLOR,         1, errp[2]) && (errp[1] - errp[0] < 4);
    bool   ok22 = isConsistent(_refCOLOR, _tagCOLOR, errp[0]+1, errp[2]) && (errp[2] - errp[1] < 4);
    bool   ok23 = isConsistent(_refCOLOR, _tagCOLOR, errp[1]+1, _len)    && (errp[3] - errp[2] < 4);

    bool   ok31 = isConsistent(_refCOLOR, _tagCOLOR,         1, errp[3]) && (errp[2] - errp[0] < 5);
    bool   ok32 = isConsistent(_refCOLOR, _tagCOLOR, errp[0]+1, _len)    && (errp[3] - errp[1] < 5);

    bool   ok41 = isConsistent(_refCOLOR, _tagCOLOR,         1, _len)    && (errp[3] - errp[0] < 6);

    //  With two exceptions, exactly one of the ok's will be true.
    //  The exceptions are:
    //
    //  a) ok21 and ok23 will imply ok41.  However there is nothing to
    //  correct here.  We just need to make sure that we stop
    //  processing rules on ok41.
    //
    //  b) ok41 and ok22.  Not sure if this can ever happen, but like
    //  case a, we're ok if we stop after ok41.
    //

    if        (ok41) {
      //  MNP of size 6
      _colorMismatch     = 4;
      _colorInconsistent = 0;
    } else if (ok31) {
      //  First three ok, fix the last one.
      _colorMismatch     = 3;
      _colorInconsistent = 1;
      _tagCOREC[errp[3]] = _refCOLOR[errp[3]];
    } else if (ok32) {
      //  Last three ok, fix the first one.
      _colorMismatch     = 3;
      _colorInconsistent = 1;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
    } else if (ok21) {
      //  First two ok, fix the last two.
      _colorMismatch     = 2;
      _colorInconsistent = 2;
      _tagCOREC[errp[2]] = _refCOLOR[errp[2]];
      _tagCOREC[errp[3]] = _refCOLOR[errp[3]];
    } else if (ok22) {
      //  Middle two ok, fix the outties.
      _colorMismatch     = 2;
      _colorInconsistent = 2;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
      _tagCOREC[errp[3]] = _refCOLOR[errp[3]];
    } else if (ok23) {
      //  Last two ok, fix the first two.
      _colorMismatch     = 2;
      _colorInconsistent = 2;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
      _tagCOREC[errp[3]] = _refCOLOR[errp[3]];
    } else {
      //  Nothing consistent, fix all of 'em.
      _colorMismatch     = 0;
      _colorInconsistent = 4;
      _tagCOREC[errp[0]] = _refCOLOR[errp[0]];
      _tagCOREC[errp[1]] = _refCOLOR[errp[1]];
      _tagCOREC[errp[2]] = _refCOLOR[errp[2]];
      _tagCOREC[errp[3]] = _refCOLOR[errp[3]];
    }
  } else if (errs == 5) {
    //fprintf(stderr, "Five errors detected.  Code doesn't know what to do.\n");
    _colorMismatch     = 0;
    _colorInconsistent = 5;
  } else if (errs == 6) {
    //fprintf(stderr, "Six errors detected.  Code doesn't know what to do.\n");
    _colorMismatch     = 0;
    _colorInconsistent = 6;
  } else {
    //fprintf(stderr, "Wow, you got a lot of errors.  Code doesn't know what to do.\n");
    _colorMismatch     = 0;
    _colorInconsistent = errs;
  }

  //  Compute alignments of corrected color strings.

  _basesMismatch = 0;

  _tagACGT[0] = baseToColor[_tagCOREC[0]][_tagCOREC[1]];
  _refACGT[0] = baseToColor[_refCOLOR[0]][_refCOLOR[1]];
  for (u32bit ti=1; ti<_len; ti++) {
    _tagACGT[ti] = baseToColor[_tagACGT[ti-1]][_tagCOREC[ti+1]];
    _refACGT[ti] = baseToColor[_refACGT[ti-1]][_refCOLOR[ti+1]];
  }
  _tagACGT[_len-1] = 0;
  _refACGT[_len-1] = 0;

  for (u32bit ti=0; ti<_len-1; ti++) {
    if (_tagACGT[ti] != _refACGT[ti]) {
      _basesMismatch++;

      _tagACGT[ti] = toUpper[_tagACGT[ti]];
      _refACGT[ti] = toUpper[_refACGT[ti]];
    }
  }

  if (_rev) {
    //  Undo the tag and ref reversals.
    _tagCOLOR[0] = complementSymbol[_tagCOLOR[0]];
    reverseString(_tagCOLOR, _len);

    _tagCOREC[0] = complementSymbol[_tagCOREC[0]];
    reverseString(_tagCOREC, _len);

    _refCOLOR[0] = complementSymbol[_refCOLOR[0]];
    reverseString(_refCOLOR, _len);

    //  Reverse complement the alignments

    reverseComplementSequence(_tagACGT, _len-1);
    reverseComplementSequence(_refACGT, _len-1);

    //  Adjust the error positions...once we start caring about positions.
  }

  return;
}



void
tapperWorkerMated(void *G, void *T, void *S) {
}


void
tapperWorkerSingle(void *G, void *T, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperThreadData  *t = (tapperThreadData  *)T;
  tapperComputation *s = (tapperComputation *)S;
  tapperHit          h;

  //
  //  Get the hits.
  //

  t->posn1fLen = 0;
  t->posn1rLen = 0;

  g->PS->getUpToNMismatches(s->tag1f, g->maxColorError, t->posn1f, t->posn1fMax, t->posn1fLen);
  g->PS->getUpToNMismatches(s->tag1r, g->maxColorError, t->posn1r, t->posn1rMax, t->posn1rLen);

  if (t->posn1fLen + t->posn1rLen == 0)
    return;

  //  True, only the mated worker really needs to put hits into the
  //  hits array, but we need to call tapperWorker_addHits to decode
  //  the positions.

  if (t->hits1Max < t->posn1fLen + t->posn1rLen) {
    t->hits1Max = t->posn1fLen + t->posn1rLen;
    delete [] t->hits1;
    t->hits1 = new u64bit [t->posn1fLen + t->posn1rLen];
  }

  t->hits1Len = 0;
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1f, t->posn1fLen, g, false, false);
  t->hits1Len = tapperWorker_addHits(t->hits1, t->hits1Len, t->posn1r, t->posn1rLen, g, false, true);

  //
  //  Do something useful with the hits, like verify them against the ACGT reference.
  //

  for (u32bit i=0; i<t->hits1Len; i++) {
    u64bit so  = (t->hits1[i] >> 33) & u64bitMASK(31);  //  Sequence it hits
    u64bit po  = (t->hits1[i] >> 2)  & u64bitMASK(31);  //  Position in sequence
    u64bit rev = (t->hits1[i] >> 0)  & 0x01;            //  isReversed

    h.alignToReference(g, so, po,
                       (rev) ? s->tag1rseq : s->tag1fseq,
                       s->tag1size+2);

    if ((h.numberOfBaseMismatches()       <= g->maxBaseError)  &&
        (h.numberOfColorMismatches()      <= g->maxColorError) &&
        (h.numberOfColorInconsistencies() <= g->maxColorError))
      s->addHit(h);
  }

  //
  //  Sort the hits by score, apply some heuristics to output only
  //  those we care about.
  //

  s->sortHits();
}



int
main(int argc, char **argv) {
  tapperGlobalData  *g = new tapperGlobalData();

  fprintf(stderr, "tapperHit: %d bytes\n", (int)sizeof(tapperHit));

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strcmp(argv[arg], "-genomic") == 0) {
      g->genName = argv[++arg];
    } else if (strcmp(argv[arg], "-color") == 0) {
      g->colName = argv[++arg];
    } else if (strcmp(argv[arg], "-queries") == 0) {
      g->qryName = argv[++arg];

    } else if (strcmp(argv[arg], "-maxcolorerror") == 0) {
      g->maxColorError = strtou32bit(argv[++arg], 0L);

    } else if (strcmp(argv[arg], "-maxbaseerror") == 0) {
      g->maxBaseError = strtou32bit(argv[++arg], 0L);

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

  ss->loaderQueueSize(65536);
  ss->writerQueueSize(16384);

  ss->numberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new tapperThreadData(g));

  ss->run(g, g->beVerbose);

  delete g;
  delete ss;

  fprintf(stderr, "\nSuccess!  Bye.\n");
  return(0);
}
