#include "tapperTag.H"
#include "tapperHit.H"
#include "tapperGlobalData.H"
#include "tapperThreadData.H"
#include "tapperComputation.H"

#undef VERBOSEWORKER

//  Very expensive.  Compare the obvious O(n^2) happy mate finding
//  algorithm against the O(n) algorithm.
//
#undef DEBUG_MATES


void*
tapperReader(void *G) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = 0L;
  tapperTag          a, b;

  if (g->TF->metaData()->isPairedTagFile()) {
    if (g->TF->get(&a, &b))
      s = new tapperComputation(&a, &b);
  } else {
    if (g->TF->get(&a))
      s = new tapperComputation(&a, 0L);
  }

  return(s);
}



void
tapperWriter(void *G, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperComputation *s = (tapperComputation *)S;
  tapperResultIndex  result;

  //  Build the result index.

  result._tag1id = s->tag1id;
  result._tag2id = s->tag2id;

  result._maxColrMismatchMapped = g->maxColorError;
  result._maxBaseMismatchMapped = g->maxBaseError;

  result._mean   = g->TF->metaData()->mean();
  result._stddev = g->TF->metaData()->stddev();

  if (s->resultFragmentLen > g->repeatThreshold) {
    result._numFragmentDiscarded = s->resultFragmentLen;
    result._numFragment          = 0;
  } else {
    result._numFragmentDiscarded = 0;
    result._numFragment          = s->resultFragmentLen;
  }

  result._numSingleton = s->resultSingletonLen;
  result._numMated     = s->resultMatedLen;
  result._numTangled   = s->resultTangledLen;

  result._pad1 = 0;
  result._pad2 = 0;

  //  Now write.

  g->outIndex->    putRecord(&result);

  g->outFragment-> putRecord(s->resultFragment,     result._numFragment);
  g->outMated->    putRecord(s->resultMated,        result._numMated);
  g->outSingleton->putRecord(s->resultSingleton,    result._numSingleton);
  g->outTangled->  putRecord(s->resultTangled,      result._numTangled);

  g->outAlignQual->putRecord(s->alignQualHistogram, s->alignQualHistogramLen);

  delete s;
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

  //  so_in and po_in are the sequence iid and position in that
  //  sequence where the tag maps.
  //
  //  tag_in is the full tag with reference base at start or end,
  //  either T010203010331 or 01031031033G.  len_in is the length of
  //  the COLOR CALLS in tag_in, NOT the strlen of it.
  //

  u32bit      errs = 0;               //  number of errors
  u32bit      errp[TAG_LEN_MAX];      //  location of the errors
  char       _tagCOREC[TAG_LEN_MAX];  //  For holding corrected color calls, only to generate ACGT align

  _seqIdx = so_in;
  _seqPos = po_in;
  _tagIdx = 0;

  //  _rev: Yeah, we assume ASCII and UNIX newlines all over the
  //  place.  A forward read starts with a reference base; reverse
  //  reads have a number here.
  //
  //  _len -- the length of the tag + reference base.
  //       -- number of color calls / ACGT + 1.
  //
  _pad                = 0;
  _len                = len_in + 1;
  _rev                = (tag_in[0] < 'A') ? true : false;

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
  //  at the end.  See tapperComputation.H for what is in the reverse
  //  tag.
  //
  {
    if (_rev) {
      for (u32bit i=0, j=_len-1; i<_len; i++, j--)
        _tagCOLOR[i] = _tagCOREC[i] = tag_in[j];
      _tagCOLOR[0] = _tagCOREC[0] = complementSymbol[_tagCOLOR[0]];
    } else {
      for (u32bit i=0; i<_len; i++)
        _tagCOLOR[i] = _tagCOREC[i] = tag_in[i];
    }
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
      reverseComplementSequence(_refACGT, _len-1);

    _refCOLOR[0] = _tagCOLOR[0];  //  ALWAYS the reference encoding base, as long as we copy the tag first.
    _refCOLOR[1] = baseToColor[_refCOLOR[0]][_refACGT[0]];

    for (u32bit ti=2; ti<_len; ti++)
      _refCOLOR[ti] = baseToColor[_refACGT[ti-2]][_refACGT[ti-1]];

    _refCOLOR[_len] = 0;
  }

  //fprintf(stderr, "tag: %s %s ref: %s %s\n", tag_in, _tagCOLOR, _refCOLOR, _refACGT);

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

  //fprintf(stderr, "tag: %s %s ref: %s %s "u32bitFMT" "u32bitFMT" "u32bitFMT"\n",
  //        tag_in, _tagCOLOR, _refCOLOR, _refACGT, _basesMismatch, _colorMismatch, _colorInconsistent);

  return;
}




//  The big value of this function is to convert from a chained
//  position to a (seqID,pos), and save the hit onto a bitpacked list.
//  This list can then be numerically sorted to order all hits.  Of
//  course, we could have just sorted the original chained positions.
//
//  It saves (seqID,pos,isTag2,isReverse)
//
inline
void
tapperWorker_addHits(u64bit    *posn, u64bit posnLen,
                     tapperGlobalData   *g,
                     tapperComputation  *s,
                     bool                rev,
                     bool                tag1) {
  tapperHit  h;
  char      *tagseq;
  u32bit     taglen;

  if (tag1) {
    tagseq = (rev) ? s->tag1rseq : s->tag1fseq;
    taglen = s->tag1size;
  } else {
    tagseq = (rev) ? s->tag2rseq : s->tag2fseq;
    taglen = s->tag2size;
  }

  for (u32bit i=0; i<posnLen; i++) {
    u64bit  pos = posn[i];
    u64bit  seq = g->SS->sequenceNumberOfPosition(pos);

    pos -= g->SS->startOf(seq);
    seq  = g->SS->IIDOf(seq);

    //  Search ignores first letter, align needs it.  This makes for a
    //  very special case, 0, which isn't a full match.

    if (pos > 0) {
      pos--;

      h.alignToReference(g, seq, pos, tagseq, taglen);

      if ((h.numberOfBaseMismatches()                                     <= g->maxBaseError)  &&
          (h.numberOfColorMismatches() + h.numberOfColorInconsistencies() <= g->maxColorError)) {
        s->addHit(g, h, tag1);
      }
    }
  }
}


void
tapperWorker(void *G, void *T, void *S) {
  tapperGlobalData  *g = (tapperGlobalData  *)G;
  tapperThreadData  *t = (tapperThreadData  *)T;
  tapperComputation *s = (tapperComputation *)S;

  //
  //  Get the hits.
  //

#ifdef VERBOSEWORKER
  fprintf(stderr, "GET HITS %s %s.\n", s->tag1fseq, s->tag2fseq);
#endif

  t->posn1fLen = t->posn1rLen = t->posn2fLen = t->posn2rLen = 0;

  if (s->tag1size > 0) {
    g->PS->getUpToNMismatches(s->tag1f, g->maxColorError, t->posn1f, t->posn1fMax, t->posn1fLen);
    g->PS->getUpToNMismatches(s->tag1r, g->maxColorError, t->posn1r, t->posn1rMax, t->posn1rLen);
  }

  if (s->tag2size > 0) {
    g->PS->getUpToNMismatches(s->tag2f, g->maxColorError, t->posn2f, t->posn2fMax, t->posn2fLen);
    g->PS->getUpToNMismatches(s->tag2r, g->maxColorError, t->posn2r, t->posn2rMax, t->posn2rLen);
  }

  //  Quit if nothing there.

  if (t->posn1fLen + t->posn1rLen + t->posn2fLen + t->posn2rLen == 0)
    return;

#ifdef VERBOSEWORKER
  fprintf(stderr, " raw hits: "u64bitFMT" "u64bitFMT" "u64bitFMT" "u64bitFMT"\n",
          t->posn1fLen, t->posn1rLen, t->posn2fLen, t->posn2rLen);
#endif

  //
  //  Align to reference to get rid of the 3/4 false hits.
  //

#ifdef VERBOSEWORKER
  fprintf(stderr, "ALIGN TO REFERENCE.\n");
#endif

  tapperWorker_addHits(t->posn1f, t->posn1fLen, g, s, false, true);
  tapperWorker_addHits(t->posn1r, t->posn1rLen, g, s, true,  true);

  tapperWorker_addHits(t->posn2f, t->posn2fLen, g, s, false, false);
  tapperWorker_addHits(t->posn2r, t->posn2rLen, g, s, true,  false);

  //  Quit if nothing there.

  if (s->tag1hitsLen + s->tag2hitsLen == 0)
    return;

  //
  //  If mated, tease out any valid mate relationships and build the
  //  results.  If fragment, just build.
  //

#ifdef VERBOSEWORKER
  fprintf(stderr, "REPORT.\n");
#endif

  //  OUTPUT CASE 1 - nothing.
  if ((s->tag1size == 0) && (s->tag2size == 0)) {
    assert(0);

    //  OUTPUT CASE 2 - unmated fragments
  } else if ((s->tag1size > 0) && (s->tag2size == 0)) {
    s->resultFragment    = new tapperResultFragment [s->tag1hitsLen];
    s->resultFragmentLen = s->tag1hitsLen;

    for (u32bit i=0; i<s->tag1hitsLen; i++) {
      s->resultFragment[i]._seq = s->tag1hits[i]._seqIdx;
      s->resultFragment[i]._pos = s->tag1hits[i]._seqPos;

      s->resultFragment[i]._bits = 0;

      s->resultFragment[i]._qual._tag1valid             = 1;
      s->resultFragment[i]._qual._tag1basesMismatch     = s->tag1hits[i]._basesMismatch;
      s->resultFragment[i]._qual._tag1colorMismatch     = s->tag1hits[i]._colorMismatch;
      s->resultFragment[i]._qual._tag1colorInconsistent = s->tag1hits[i]._colorInconsistent;
      s->resultFragment[i]._qual._tag1rev               = s->tag1hits[i]._rev;
    }

    //  OUTPUT CASE 3 - unmated fragments (but wrong set, should always be in tag1)
  } else if ((s->tag1size == 0) && (s->tag2size > 0)) {
    assert(0);

    //  OUTPUT CASE 4 - mated fragments
  } else if ((s->tag1size > 0) && (s->tag2size > 0)) {
    if (t->tangle == 0L)
      t->tangle = new intervalList [g->GS->fasta()->getNumberOfSequences()];

    if ((t->numHappiesMax < s->tag1hitsLen) || (t->numHappiesMax < s->tag2hitsLen)) {
      delete [] t->tag1happies;
      delete [] t->tag1mate;
      delete [] t->tag2happies;
      delete [] t->tag2mate;

      t->numHappiesMax = MAX(s->tag1hitsLen, s->tag2hitsLen) + 16 * 1024;

      t->tag1happies = new u32bit [t->numHappiesMax];
      t->tag1mate    = new u32bit [t->numHappiesMax];
      t->tag2happies = new u32bit [t->numHappiesMax];
      t->tag2mate    = new u32bit [t->numHappiesMax];
    }

#ifdef VERBOSEWORKER
    fprintf(stderr, "  Found "u32bitFMT" and "u32bitFMT" hits.\n", s->tag1hitsLen, s->tag2hitsLen);
#endif

    //  Sort by position.
    s->sortHitsByPosition();

    u32bit  mean   = g->TF->metaData()->mean();
    u32bit  stddev = g->TF->metaData()->stddev();

    tapperHit *t1h = s->tag1hits;
    tapperHit *t2h = s->tag2hits;

    //  Pass zero, clear.  Tangles are cleared below.
    //
    for (u32bit a=0; a<s->tag1hitsLen; a++)
      t->tag1happies[a] = 0;
    for (u32bit b=0; b<s->tag2hitsLen; b++)
      t->tag2happies[b] = 0;

    //  Pass one.  Count the number of times each fragment is in a
    //  happy relationship.
    //
    {
#ifdef DEBUG_MATES
      u32bit  debug_numHappies = 0;
      u64bit  debug_happyCheck = 0;

      for (u32bit a=0; a<s->tag1hitsLen; a++) {
        for (u32bit b=0; b<s->tag2hitsLen; b++) {
          if (t1h[a].happy(t2h[b], mean, stddev)) {
            debug_numHappies += 1;
            debug_happyCheck += t1h[a]._seqPos ^ t2h[b]._seqPos;
          }
        }
      }
#endif

      u32bit  bbaserev = 0;
      u32bit  bbasefor = 0;

      for (u32bit a=0; a<s->tag1hitsLen; a++) {

        //  Both lists of hits are sorted by position.  For each tag1 (a)
        //  hit, we first advance the bbase to the first hit that is
        //  within the proper distance before the a tag.  Then scan forward
        //  until the b tag is too far away to be mated.

        u32bit b = 0;

        if (t1h[a]._rev == true) {
          while ((bbaserev < s->tag2hitsLen) && (t1h[a].mateTooFarBefore(t2h[bbaserev], mean, stddev)))
            bbaserev++;
          b = bbaserev;
        } else {
          while ((bbasefor < s->tag2hitsLen) && (t1h[a].mateTooFarBefore(t2h[bbasefor], mean, stddev)))
            bbasefor++;
          b = bbasefor;
        }

        //  Now, until the b read is too far away to be mated, check
        //  for happiness and do stuff.

        for (; (b<s->tag2hitsLen) && (t1h[a].mateTooFarAfter(t2h[b], mean, stddev) == false); b++) {
          if (t1h[a].happy(t2h[b], mean, stddev)) {

#ifdef DEBUG_MATES
            debug_numHappies -= 1;
            debug_happyCheck -= t1h[a]._seqPos ^ t2h[b]._seqPos;
#endif

            //  Count.
            t->tag1happies[a]++;
            t->tag2happies[b]++;

            //  Add the previous mate pair if we just became tangled.
            //  It is possible for both to be == 2, but in that case,
            //  we've already added the previous mate pair.
            if ((t->tag1happies[a] == 2) && (t->tag2happies[b] == 1)) {
              u32bit c  = t->tag1mate[a];
              u32bit mn = MIN(t1h[a]._seqPos,               t2h[c]._seqPos);
              u32bit mx = MAX(t1h[a]._seqPos + s->tag1size, t2h[c]._seqPos + s->tag2size);

              t->tangle[t1h[a]._seqIdx].add(mn, mx-mn);
            }

            if ((t->tag1happies[a] == 1) && (t->tag2happies[b] == 2)) {
              u32bit c  = t->tag2mate[b];
              u32bit mn = MIN(t1h[c]._seqPos,               t2h[b]._seqPos);
              u32bit mx = MAX(t1h[c]._seqPos + s->tag1size, t2h[b]._seqPos + s->tag2size);

              t->tangle[t1h[c]._seqIdx].add(mn, mx-mn);
            }

            //  Finally, add the current mate pair to the tangle.
            if ((t->tag1happies[a] >= 2) || (t->tag2happies[b] >= 2)) {
              u32bit mn = MIN(t1h[a]._seqPos,               t2h[b]._seqPos);
              u32bit mx = MAX(t1h[a]._seqPos + s->tag1size, t2h[b]._seqPos + s->tag2size);

              t->tangle[t1h[a]._seqIdx].add(mn, mx-mn);
            }

            //  Remember the mate; only valid if tag1happies[a] and
            //  tag2happies[b] both == 1.
            t->tag1mate[a] = b;
            t->tag2mate[b] = a;
          }
        }
      }

#ifdef DEBUG_MATES
      if ((debug_numHappies != 0) || (debug_happyCheck != 0)) {
        FILE *df = fopen("tapper.DEBUG_MATES.err", "w");

        fprintf(df, "numHappies: "u64bitFMT"\n", debug_numHappies);
        fprintf(df, "happyCheck: "u64bitFMT"\n", debug_happyCheck);

        for (u32bit a=0; a<s->tag1hitsLen; a++)
          fprintf(df, "a="u32bitFMT" ori=%c pos="u32bitFMT","u32bitFMT"\n",
                  a, t1h[a]._rev ? 'r' : 'f', t1h[a]._seqIdx, t1h[a]._seqPos);

        for (u32bit b=0; b<s->tag2hitsLen; b++)
          fprintf(df, "b="u32bitFMT" ori=%c pos="u32bitFMT","u32bitFMT"\n",
                  b, t2h[b]._rev ? 'r' : 'f', t2h[b]._seqIdx, t2h[b]._seqPos);

        u32bit  bbaserev = 0;
        u32bit  bbasefor = 0;

        for (u32bit a=0; a<s->tag1hitsLen; a++) {
          u32bit b = 0;

          if (t1h[a]._rev == true) {
            while ((bbaserev < s->tag2hitsLen) && (t1h[a].mateTooFarBefore(t2h[bbaserev], mean, stddev))) {
              fprintf(df, "rev bbaserev <- "u32bitFMT" + 1\n", bbaserev);
              bbaserev++;
            }
            b = bbaserev;
          } else {
            while ((bbasefor < s->tag2hitsLen) && (t1h[a].mateTooFarBefore(t2h[bbasefor], mean, stddev))) {
              fprintf(df, "rev bbasefor <- "u32bitFMT" + 1\n", bbasefor);
              bbasefor++;
            }
            b = bbasefor;
          }

          for (; (b<s->tag2hitsLen) && (t1h[a].mateTooFarAfter(t2h[b], mean, stddev) == false); b++) {
            fprintf(df, "test a="u32bitFMT" b="u32bitFMT"\n", a, b);
            if (t1h[a].happy(t2h[b], mean, stddev)) {
              fprintf(df, "HAPPY CLEVER     a="u32bitFMT" b="u32bitFMT"\n", a, b);
            }
          }
        }

        for (u32bit a=0; a<s->tag1hitsLen; a++) {
          for (u32bit b=0; b<s->tag2hitsLen; b++) {
            if (t1h[a].happy(t2h[b], mean, stddev)) {
              fprintf(df, "HAPPY EXHAUSTIVE a="u32bitFMT" b="u32bitFMT"\n", a, b);
            }
          }
        }

        fclose(df);
      }
      assert(debug_numHappies == 0);
      assert(debug_happyCheck == 0);
#endif


#ifdef VERBOSEWORKER
      fprintf(stderr, "  Paired.\n");
#endif
    }


    //  Allocate space for the outputs.  We can kind of guess how much
    //  to grab.  Not perfect.  Can do a lot better.
    //
    s->resultFragment   = new tapperResultFragment   [s->tag1hitsLen + s->tag2hitsLen];
    s->resultMated      = new tapperResultMated      [MIN(s->tag1hitsLen, s->tag2hitsLen)];
    s->resultSingleton  = new tapperResultSingleton  [s->tag1hitsLen + s->tag2hitsLen];
    s->resultTangled    = new tapperResultTangled    [MIN(s->tag1hitsLen, s->tag2hitsLen)];


    //  Pass two.  For anything with zero happies, emit to the
    //  singleton file.  For anything with a pair of single happies,
    //  emit to the happy mate file.
    //
    for (u32bit a=0; a<s->tag1hitsLen; a++) {
      if (t->tag1happies[a] == 0) {
        if (s->tag1hits[a].happyNearEnd(true, mean, stddev, g->GS->fasta()->sequenceLength(s->tag1hits[a]._seqIdx))) {
          tapperResultSingleton *f = s->resultSingleton + s->resultSingletonLen++;

          f->_seq = s->tag1hits[a]._seqIdx;
          f->_pos = s->tag1hits[a]._seqPos;

          f->_bits = 0;

          f->_qual._tag1valid             = 1;
          f->_qual._tag1basesMismatch     = s->tag1hits[a]._basesMismatch;
          f->_qual._tag1colorMismatch     = s->tag1hits[a]._colorMismatch;
          f->_qual._tag1colorInconsistent = s->tag1hits[a]._colorInconsistent;
          f->_qual._tag1rev               = s->tag1hits[a]._rev;
        } else {
          tapperResultFragment *f = s->resultFragment + s->resultFragmentLen++;

          f->_seq = s->tag1hits[a]._seqIdx;
          f->_pos = s->tag1hits[a]._seqPos;

          f->_bits = 0;

          f->_qual._tag1valid             = 1;
          f->_qual._tag1basesMismatch     = s->tag1hits[a]._basesMismatch;
          f->_qual._tag1colorMismatch     = s->tag1hits[a]._colorMismatch;
          f->_qual._tag1colorInconsistent = s->tag1hits[a]._colorInconsistent;
          f->_qual._tag1rev               = s->tag1hits[a]._rev;
        }
      }
    }

    for (u32bit b=0; b<s->tag2hitsLen; b++) {
      if (t->tag2happies[b] == 0) {
        if (s->tag2hits[b].happyNearEnd(false, mean, stddev, g->GS->fasta()->sequenceLength(s->tag2hits[b]._seqIdx))) {
          tapperResultSingleton *f = s->resultSingleton + s->resultSingletonLen++;

          f->_seq = s->tag2hits[b]._seqIdx;
          f->_pos = s->tag2hits[b]._seqPos;

          f->_bits = 0;

          f->_qual._tag2valid             = 1;
          f->_qual._tag2basesMismatch     = s->tag2hits[b]._basesMismatch;
          f->_qual._tag2colorMismatch     = s->tag2hits[b]._colorMismatch;
          f->_qual._tag2colorInconsistent = s->tag2hits[b]._colorInconsistent;
          f->_qual._tag2rev               = s->tag2hits[b]._rev;
        } else {
          tapperResultFragment *f = s->resultFragment + s->resultFragmentLen++;

          f->_seq = s->tag2hits[b]._seqIdx;
          f->_pos = s->tag2hits[b]._seqPos;

          f->_bits = 0;

          f->_qual._tag2valid             = 1;
          f->_qual._tag2basesMismatch     = s->tag2hits[b]._basesMismatch;
          f->_qual._tag2colorMismatch     = s->tag2hits[b]._colorMismatch;
          f->_qual._tag2colorInconsistent = s->tag2hits[b]._colorInconsistent;
          f->_qual._tag2rev               = s->tag2hits[b]._rev;
        }
      }
    }

    for (u32bit a=0; a<s->tag1hitsLen; a++) {
      u32bit b = t->tag1mate[a];

      if ((t->tag1happies[a] == 1) && (t->tag2happies[b] == 1)) {
        tapperResultMated     *m = s->resultMated + s->resultMatedLen++;

        assert(t->tag1mate[a] == b);
        assert(t->tag2mate[b] == a);

        m->_seq  = s->tag1hits[a]._seqIdx;
        m->_pos1 = s->tag1hits[a]._seqPos;
        m->_pos2 = s->tag2hits[b]._seqPos;

        m->_bits = 0;

        m->_qual._tag1valid             = 1;
        m->_qual._tag1basesMismatch     = s->tag1hits[a]._basesMismatch;
        m->_qual._tag1colorMismatch     = s->tag1hits[a]._colorMismatch;
        m->_qual._tag1colorInconsistent = s->tag1hits[a]._colorInconsistent;
        m->_qual._tag1rev               = s->tag1hits[a]._rev;

        m->_qual._tag2valid             = 1;
        m->_qual._tag2basesMismatch     = s->tag2hits[b]._basesMismatch;
        m->_qual._tag2colorMismatch     = s->tag2hits[b]._colorMismatch;
        m->_qual._tag2colorInconsistent = s->tag2hits[b]._colorInconsistent;
        m->_qual._tag2rev               = s->tag2hits[b]._rev;
      }
    }

    //  Pass three.  Emit and then clear the tangles.
    //
    {
      u32bit simax = g->GS->fasta()->getNumberOfSequences();

      for (u32bit si=0; si<simax; si++) {

        if (t->tangle[si].numberOfIntervals() > 0) {
          t->tangle[si].merge();

          for (u32bit ti=0; ti<t->tangle[si].numberOfIntervals(); ti++) {
            tapperResultTangled   *x = s->resultTangled + s->resultTangledLen++;

            x->_tag1count = 0;
            x->_tag2count = 0;

            x->_seq = si;

            x->_bgn = t->tangle[si].lo(ti);
            x->_end = t->tangle[si].hi(ti);
          }

          //  This is persistent; clear it for the next mate pair.
          t->tangle[si].clear();
        }
      }
    }
  }
}



int
main(int argc, char **argv) {
  tapperGlobalData  *g = new tapperGlobalData();

  fprintf(stderr, "sizeof(tapperResultIndex) --       %d\n", sizeof(tapperResultIndex));
  fprintf(stderr, "sizeof(tapperResultQV) --          %d\n", sizeof(tapperResultQV));
  fprintf(stderr, "sizeof(tapperResultHistogram) --   %d\n", sizeof(tapperResultHistogram));
  fprintf(stderr, "sizeof(tapperResultFragment) --    %d\n", sizeof(tapperResultFragment));
  fprintf(stderr, "sizeof(tapperResultMated) --       %d\n", sizeof(tapperResultMated));
  fprintf(stderr, "sizeof(tapperResultSingleton) --   %d\n", sizeof(tapperResultSingleton));
  fprintf(stderr, "sizeof(tapperResultTangled) --     %d\n", sizeof(tapperResultTangled));
  fprintf(stderr, "sizeof(tapperHit) --               %d\n", sizeof(tapperHit));
  fprintf(stderr, "sizeof(tapperTag) --               %d\n", sizeof(tapperTag));

  int arg=1;
  int err=0;
  while (arg < argc) {
    if        (strncmp(argv[arg], "-genomic", 2) == 0) {
      g->genName = argv[++arg];
    } else if (strncmp(argv[arg], "-queries", 2) == 0) {
      g->qryName = argv[++arg];
    } else if (strncmp(argv[arg], "-prefix", 2) == 0) {
      g->outName = argv[++arg];

    } else if (strncmp(argv[arg], "-begin", 2) == 0) {
      g->bgnRead = strtou32bit(argv[++arg], 0L);
    } else if (strncmp(argv[arg], "-end", 2) == 0) {
      g->endRead = strtou32bit(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-repeatthreshold", 2) == 0) {
      g->repeatThreshold = strtou32bit(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-maxcolorerror", 5) == 0) {
      g->maxColorError = strtou32bit(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-maxbaseerror", 5) == 0) {
      g->maxBaseError = strtou32bit(argv[++arg], 0L);

    } else if (strncmp(argv[arg], "-maxmemory", 5) == 0) {
      g->maxMemory = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-numthreads", 2) == 0) {
      g->numThreads = atoi(argv[++arg]);

    } else if (strncmp(argv[arg], "-verbose", 2) == 0) {
      g->beVerbose = true;
    } else {
      fprintf(stderr, "%s: unknown option '%s'\n", argv[0], argv[arg]);
      err++;
    }
    arg++;
  }
  if ((err > 0) || (g->genName == 0L) || (g->qryName == 0L) || (g->outName == 0L)) {
    fprintf(stderr, "usage: %s [opts]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "  MANDATORY\n");
    fprintf(stderr, "          -genomic genomic.fasta\n");
    fprintf(stderr, "          -queries tags.tapperTags\n");
    fprintf(stderr, "          -prefix  output-prefix\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  OPTIONAL\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "          -begin b               Start aligning at read b\n");
    fprintf(stderr, "          -end   e               Stop aligning at read e\n");
    fprintf(stderr, "                                 NOTE!  When mapping mated reads, this will\n");
    fprintf(stderr, "                                 start/stop at matepair b/2 and e/2.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "          -repeatthreshold x     Do not report fragment alignments for tags\n");
    fprintf(stderr, "                                 with more than x alignments.  Singletons, mated\n");
    fprintf(stderr, "                                 tags and are still reported and computed using\n");
    fprintf(stderr, "                                 all alignments.  The default is "u32bitFMT".\n", g->repeatThreshold);
    fprintf(stderr, "\n");
    fprintf(stderr, "          -maxcolorerror  n\n");
    fprintf(stderr, "          -maxbaseerror   n\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "          -maxmemory      m (MB)\n");
    fprintf(stderr, "          -numthreads     n\n");
    fprintf(stderr, "          -verbose\n");
            
    exit(1);
  }

  g->initialize();

  sweatShop *ss = new sweatShop(tapperReader, tapperWorker, tapperWriter);

  ss->loaderQueueSize(2000);
  ss->writerQueueSize(30000);  //  TESTING!  Should be 1000-ish!

  ss->numberOfWorkers(g->numThreads);

  for (u32bit w=0; w<g->numThreads; w++)
    ss->setThreadData(w, new tapperThreadData(g));

  ss->run(g, g->beVerbose);

  delete g;
  delete ss;

  fprintf(stderr, "\nSuccess!  Bye.\n");
  return(0);
}
