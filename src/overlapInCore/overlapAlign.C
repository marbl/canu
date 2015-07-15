#include "overlapAlign.H"

#include "kMer.H"
#include "merStream.H"

#include "Display_Alignment.H"

#undef  DEBUG_ALGORTHM         //  Some details.
#undef  DEBUG_HITS             //  Lots of details.

#undef  SEED_NON_OVERLAPPING   //  Allow mismatches in seeds
#define SEED_OVERLAPPING


overlapAlign::overlapAlign(pedAlignType   alignType,
                           double         maxErate,
                           int32          merSize) {
  _alignType        = alignType;
  _maxErate         = maxErate;
  _merSizeInitial   = merSize;

  _aID      = UINT32_MAX;
  _aStr     = NULL;
  _aLen     = 0;
  _aLoOrig  = 0;
  _aHiOrig  = 0;

  _bID      = UINT32_MAX;
  _bStr     = NULL;
  _bLen     = 0;
  _bLoOrig  = 0;
  _bHiOrig  = 0;

  _bFlipped = false;
  _bRevMax  = 0;
  _bRev     = NULL;

  _editDist = new prefixEditDistance(_alignType, _maxErate);

  _minDiag  = 0;
  _maxDiag  = 0;

  _merSize  = 0;

  //  Initialize Constants

  for (uint32 ii=0; ii<256; ii++) {
    acgtToBit[ii] = 0x00;
    acgtToVal[ii] = 0x03;
  }

  acgtToBit['a'] = acgtToBit['A'] = 0x00;  //  Bit encoding of ACGT
  acgtToBit['c'] = acgtToBit['C'] = 0x01;
  acgtToBit['g'] = acgtToBit['G'] = 0x02;
  acgtToBit['t'] = acgtToBit['T'] = 0x03;

  acgtToVal['a'] = acgtToVal['A'] = 0x00;  //  Word is valid if zero
  acgtToVal['c'] = acgtToVal['C'] = 0x00;
  acgtToVal['g'] = acgtToVal['G'] = 0x00;
  acgtToVal['t'] = acgtToVal['T'] = 0x00;

  merMask[0] = 0x0000000000000000llu;
  merMask[1] = 0x0000000000000003llu;

  for (uint32 ii=2; ii<33; ii++)
    merMask[ii] = (merMask[ii-1] << 2) | 0x03;

  assert(merMask[ 6] == 0x0000000000000fffllu);
  assert(merMask[17] == 0x00000003ffffffffllu);
  assert(merMask[26] == 0x000fffffffffffffllu);
  assert(merMask[32] == 0xffffffffffffffffllu);
}



overlapAlign::~overlapAlign() {
  delete    _editDist;
  delete [] _bRev;
}



void
overlapAlign::initialize(uint32 aID, char *aStr, int32 aLen, int32 aLo, int32 aHi,
                         uint32 bID, char *bStr, int32 bLen, int32 bLo, int32 bHi, bool bFlipped) {

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "overlapAlign::initialize()--  A %5u %5u - %5u  B %5u %5u - %5u\n", aID, aLo, aHi, bID, bLo, bHi);
#endif

  _aID      = aID;
  _aStr     = aStr;
  _aLen     = aLen;
  _aLoOrig  = aLo;
  _aHiOrig  = aHi;

  _bID      = bID;
  _bStr     = bStr;
  _bLen     = bLen;
  _bLoOrig  = bLo;
  _bHiOrig  = bHi;

  _bFlipped  = bFlipped;

  if (_bFlipped == true) {
    if (_bRevMax < _bLen) {
      delete [] _bRev;
      _bRevMax = _bLen + 1000;
      _bRev    = new char [_bRevMax];
    }

    memcpy(_bRev, bStr, sizeof(char) * (_bLen + 1));

    reverseComplementSequence(_bRev, _bLen);

    _bStr = _bRev;

    _bLoOrig = _bLen - bLo;  //  Now correct for the reverse complemented sequence
    _bHiOrig = _bLen - bHi;
  }

#ifdef DEBUG_ALGORITHM
  //fprintf(stderr, "overlapAlign:initialize()--  A %5u %5u - %5u  B %5u %5u - %5u\n", aID, _aLoOrig, _aHiOrig, bID, _bLoOrig, _bHiOrig);
#endif

  assert(_aLoOrig < _aHiOrig);
  assert(_bLoOrig < _bHiOrig);

  //_editDist doesn't need to be cleared.

  _minDiag = 0;
  _maxDiag = 0;

  _merSize = 0;

  _aMap.clear();
  _bMap.clear();

  _rawhits.clear();
  _hits.clear();

  _bestResult.clear();
}




//  Set the min/max diagonal we will accept seeds for.  It's just the min/max diagonal for the
//  two endpoints extended by half the erate times the align length.
//
//  Returns false if there is no chance of an overlap - the computed alignment length is too small.
//
bool
overlapAlign::findMinMaxDiagonal(int32  minLength) {

  _minDiag = 0;
  _maxDiag = 0;

  int32  aALen = _aHiOrig - _aLoOrig;
  int32  bALen = _bHiOrig - _bLoOrig;

  int32  alignLen = (aALen < bALen) ? bALen : aALen;

  if (alignLen < minLength)
    return(false);

  int32  bgnDiag = _aLoOrig - _bLoOrig;
  int32  endDiag = _aHiOrig - _bHiOrig;

  if (bgnDiag < endDiag) {
    _minDiag = bgnDiag - _maxErate * alignLen / 2;
    _maxDiag = endDiag + _maxErate * alignLen / 2;

  } else {
    _minDiag = endDiag - _maxErate * alignLen / 2;
    _maxDiag = bgnDiag + _maxErate * alignLen / 2;
  }

  //  For very very short overlaps (mhap kindly reports 4 bp overlaps) reset the min/max
  //  to something a little more permissive.
  //if (_minDiag > -5)  _minDiag = -5;
  //if (_maxDiag <  5)  _maxDiag =  5;

  if (_minDiag > _maxDiag)
    fprintf(stderr, "ERROR: _minDiag=%d >= _maxDiag=%d\n", _minDiag, _maxDiag);
  assert(_minDiag <= _maxDiag);

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "overlapAlign::findMinMaxDiagonal()--  min %d max %d -- span %d -- alignLen %d\n", _minDiag, _maxDiag, _maxDiag-_minDiag, alignLen);
#endif

  return(true);
}




void
overlapAlign::fastFindMersA(bool dupIgnore) {

  uint64   mer = 0x0000000000000000llu;
  uint64   val = 0xffffffffffffffffllu;

  uint64   mask = merMask[_merSize];

  //  Restict the mers we seed with the those in the overlap (plus a little wiggle room).  If this
  //  isn't done, overlaps in repeats are sometimes lost.  In the A read, mers will drop out (e.g.,
  //  if a repeat at both the 5' and 3' end, and we have an overlap on one end only).  In the B
  //  read, the mers can now be on the wrong diagonal.

  int32  bgn = _aLoOrig - _merSize - _merSize;
  int32  end = _aHiOrig + _merSize;

  if (bgn < 0)
    bgn = 0;

  //  Create mers.  Since 'val' was initialized as invalid until the first _merSize things
  //  are pushed on, no special case is needed to load the mer.  It costs us two extra &'s
  //  and the test for saving the valid mer while we initialize.

  for (int32 seqpos=bgn; (seqpos < end) && (_aStr[seqpos] != 0); seqpos++) {
    mer <<= 2;
    val <<= 2;

    mer |= acgtToBit[_aStr[seqpos]];
    val |= acgtToVal[_aStr[seqpos]];

    mer &= merMask[_merSize];
    val &= merMask[_merSize];

    if (val != 0x0000000000000000)
      //  Not a valid mer.
      continue;

    //  +1 - consider a 1-mer.  The first time through we have a valid mer, but seqpos == 0.
    //  To get an aMap position of zero (the true position) we need to add one.

    if (_aMap.find(mer) != _aMap.end()) {
      if (dupIgnore == true)
        _aMap[mer] = INT32_MAX;  //  Duplicate mer, now ignored!
    } else {
      _aMap[mer] = seqpos + 1 - _merSize;
    }
  }

  //fprintf(stderr, "Found %u hits in A at mersize %u dupIgnore %u t %u %u\n", _aMap.size(), _merSize, dupIgnore, t[0], t[1]);
}


void
overlapAlign::fastFindMersB(bool dupIgnore) {

  uint64   mer = 0x0000000000000000llu;
  uint64   val = 0xffffffffffffffffllu;

  uint64   mask = merMask[_merSize];

  //  Like the A read, we limit to mers in the overlap region.

  int32  bgn = _bLoOrig - _merSize - _merSize;
  int32  end = _bHiOrig + _merSize;

  if (bgn < 0)
    bgn = 0;

  //  Create mers.  Since 'val' was initialized as invalid until the first _merSize things
  //  are pushed on, no special case is needed to load the mer.  It costs us two extra &'s
  //  and the test for saving the valid mer while we initialize.

  for (int32 seqpos=bgn; (seqpos < end) && (_bStr[seqpos] != 0); seqpos++) {
    mer <<= 2;
    val <<= 2;

    mer |= acgtToBit[_bStr[seqpos]];
    val |= acgtToVal[_bStr[seqpos]];

    mer &= merMask[_merSize];
    val &= merMask[_merSize];

    if (val != 0x0000000000000000)
      //  Not a valid mer.
      continue;

    if (_aMap.find(mer) == _aMap.end())
      //  Not in the A sequence, don't care.
      continue;

    int32  apos = _aMap[mer];
    int32  bpos = seqpos + 1 - _merSize;

    if (apos == INT32_MAX)
      //  Exists too many times in aSeq, don't care.
      continue;

    if ((apos - bpos < _minDiag) ||
        (apos - bpos > _maxDiag))
      //  Too different.
      continue;

    if (_bMap.find(mer) != _bMap.end()) {
      if (dupIgnore == true)
        _bMap[mer] = INT32_MAX;  //  Duplicate mer, now ignored!
    } else {
      _bMap[mer] = bpos;
    }
  }

  //fprintf(stderr, "Found %u hits in B at mersize %u dupIgnore %u t %u %u %u %u %u\n", _bMap.size(), _merSize, dupIgnore, t[0], t[1], t[2], t[3], t[4]);
}





//  Find seeds - hash the kmer and position from the first read, then lookup
//  each kmer in the second read.  For unique hits, save the diagonal.  Then what?
//  If the diagonal is too far from the expected diagonal (based on the overlap),
//  ignore the seed.
//
//  If dupIgnore == true, ignore duplicates.  Otherwise, use the first occurrence
//
//  Returns false if no mers found.
//
bool
overlapAlign::findSeeds(bool dupIgnore) {

  _merSize = _merSizeInitial;

  //  Find mers in A

 findMersAgain:

  fastFindMersA(dupIgnore);

  if (_aMap.size() == 0) {
    _aMap.clear();
    _bMap.clear();

    _merSize--;

    if ((_merSize < 8) && (dupIgnore == true)) {
      _merSize   = _merSizeInitial + 2;
      dupIgnore =  false;
    }

    if (_merSize >= 8)
      goto findMersAgain;
  }

  //  Find mers in B

  fastFindMersB(dupIgnore);

  if (_bMap.size() == 0) {
    _aMap.clear();
    _bMap.clear();

    _merSize--;

    if ((_merSize < 8) && (dupIgnore == true)) {
      _merSize   = _merSizeInitial + 2;
      dupIgnore =  false;
    }

    if (_merSize >= 8)
      goto findMersAgain;
  }

  //  Still zero?  Didn't find any unique seeds anywhere.

  if (_bMap.size() == 0) {
#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "overlapAlign::findSeeds()--  No seeds found.\n");
#endif
    return(false);
  }

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "overlapAlign::findSeeds()--  Found %u seeds.\n", _bMap.size());
#endif
  return(true);
}


//  For unique mers in B, if the mer is also unique in A, add a hit.

bool
overlapAlign::findHits(void) {

  for (map<uint64,int32>::iterator bit=_bMap.begin(); bit != _bMap.end(); bit++) {
    uint64  kmer = bit->first;
    int32   bpos = bit->second;

    if (bpos == INT32_MAX)
      //  Exists too many times in bSeq, don't care about it.
      continue;

    int32  apos = _aMap[kmer];

    assert(apos != INT32_MAX);        //  Should never get a bMap if the aMap isn't set

    if ((apos - bpos < _minDiag) ||
        (apos - bpos > _maxDiag))
      fprintf(stderr, "kmer "F_X64" apos - bpos = %d  _minDiag = %d  _maxDiag = %d\n",
              kmer, apos-bpos, _minDiag, _maxDiag);
    assert(apos - bpos >= _minDiag);  //  ...these too.
    assert(apos - bpos <= _maxDiag);
    
    _rawhits.push_back(exactMatch(apos, bpos, _merSize));
  }

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "overlapAlign::findHits()--  Found %u hits.\n", _rawhits.size());
#endif

  return(true);
}


bool
overlapAlign::chainHits(void) {

  //  Sort by aPos (actually by length, then by aPos, but length is constant here).

  sort(_rawhits.begin(), _rawhits.end());

#ifdef DEBUG_HITS
  for (uint32 rr=0; rr<_rawhits.size(); rr++)
    if (_rawhits[rr].aBgn - _rawhits[rr].bBgn == 2289)
      fprintf(stderr, "HIT: %d - %d diag %d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn);
#endif

  //  Chain the _hits.  This chains hits that are on the same diagonal and contiguous.

  _hits.push_back(_rawhits[0]);

  for (uint32 rr=1; rr<_rawhits.size(); rr++) {
    uint32  hh    = _hits.size() - 1;
    bool    merge = false;

    assert(_rawhits[rr-1].aBgn < _rawhits[rr].aBgn);

#ifdef SEED_NON_OVERLAPPING
    //  Allows non-overlapping (mismatch only) hits - but breaks seeding of alignment
    //           --------- end of next hit -----------   ----- end of existing hit -----
    int32   da = _rawhits[rr].aBgn + _rawhits[rr].tLen - _hits[hh].aBgn - _hits[hh].tLen;
    int32   db = _rawhits[rr].bBgn + _rawhits[rr].tLen - _hits[hh].bBgn - _hits[hh].tLen;

    assert(da > 0);

    if (da == db)
      merge = true;
#endif

#ifdef SEED_OVERLAPPING
    //  Requires overlapping hits - full block of identity
    //           - start of next -   ------- end of existing -------
    int32   da = _rawhits[rr].aBgn - _hits[hh].aBgn - _hits[hh].tLen;
    int32   db = _rawhits[rr].bBgn - _hits[hh].bBgn - _hits[hh].tLen;

    if ((da < 0) && (da == db))
      merge = true;
#endif

    if (merge) {
      //fprintf(stderr, "MERGE HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
      _hits[hh].tLen = _rawhits[rr].aBgn + _rawhits[rr].tLen - _hits[hh].aBgn;
    } else {
      //fprintf(stderr, "NEW   HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
      _hits.push_back(_rawhits[rr]);
    }
  }

  //  Sort by longest

  sort(_hits.begin(), _hits.end());

#ifdef DEBUG_HITS
  for (uint32 hh=0; hh<_hits.size(); hh++) {
    fprintf(stderr, "hit %02u %5d-%5d diag %d len %3u\n",
            hh,
            _hits[hh].aBgn, _hits[hh].bBgn,
            _hits[hh].aBgn - _hits[hh].bBgn,
            _hits[hh].tLen);
  }
#endif

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "overlapAlign::chainHits()--  Found %u chains of hits.\n", _hits.size());
#endif

  return(true);
}



bool
overlapAlign::processHits(void) {

  //  The expected worst score, in numer of matches.
  double   expectedScore = (1 - _maxErate) * min(_aHiOrig - _aLoOrig, _bHiOrig - _bLoOrig);

  //  Scratch space for finding alignments

  int32  aLo=0, aHi=0;
  int32  bLo=0, bHi=0;

  //  Go!

  for (uint32 hh=0; hh<_hits.size(); hh++) {
    Match_Node_t  match;

    match.Start  = _hits[hh].aBgn;    //  Begin position in a
    match.Offset = _hits[hh].bBgn;    //  Begin position in b
    match.Len    = _hits[hh].tLen;    //  tLen can include mismatches if alternate scoring is used!
    match.Next   = 0;                 //  Not used here

#ifdef SEED_NON_OVERLAPPING
    match.Offset = _merSize;  //  Really should track this in the hits, oh well.
#endif

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "\n");
    fprintf(stderr, "overlapAlign::processHits()-- Extend_Alignment Astart %d Bstart %d length %d\n", match.Start, match.Offset, match.Len);
#endif

    int32           errors  = 0;
    pedOverlapType  olapType = _editDist->Extend_Alignment(&match,         //  Initial exact match, relative to start of string
                                                           _aStr, _aLen,
                                                           _bStr, _bLen,
                                                           aLo,   aHi,    //  Output: Regions which the match extends
                                                           bLo,   bHi,
                                                           errors);

    aHi++;  //  Add one to the end point because Extend_Alignment returns the base-based coordinate.
    bHi++;

    int32  olapLen   = min(aHi - aLo, bHi - bLo);
    double olapQual  = (double)errors / olapLen;
    double olapScore = olapLen * (1 - olapQual);

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "overlapAlign::processHits()-- hit %2u at a=%5d b=%5d -- ORIG A %6u %5d-%5d (%5d) %s B %6u %5d-%5d (%5d)  %.4f -- REALIGN %s A %5d-%5d  B %5d-%5d  errors %4d  erate %6.4f = %6u / %6u deltas %d %d %d %d %d\n",
            hh, _hits[hh].aBgn, _hits[hh].bBgn,
            _aID, _aLoOrig, _aHiOrig, _aLen,
            _bFlipped ? "<-" : "->",
            _bID, _bLoOrig, _bHiOrig, _bLen,
            erate(),
            toString(olapType),
            aLo,
            aHi,
            (_bFlipped == false) ? bLo : bHi,
            (_bFlipped == false) ? bHi : bLo,
            errors,
            olapQual, errors, olapLen,
            _editDist->Left_Delta[0],
            _editDist->Left_Delta[1],
            _editDist->Left_Delta[2],
            _editDist->Left_Delta[3],
            _editDist->Left_Delta[4]);
#endif

    //  Is this a better overlap than what we have?

    if (_bestResult.score() <= olapScore) {
#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "overlapAlign::processHits()--  - save best! - score %f previous %f expected %f\n", olapScore, _bestResult.score(), expectedScore);
#endif

      _bestResult.save(aLo, aHi, bLo, bHi, olapLen, olapQual, olapType, _editDist->Left_Delta_Len, _editDist->Left_Delta);
    }

    //  If pretty crappy, keep looking.

    if (_bestResult.score() < 0.5 * expectedScore) {
#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "overlapAlign::processHits()--  - too short: score %f < 0.5 * expected = %f\n",
              _bestResult.score(), 0.5 * expectedScore);
#endif
      continue;
    }

    //  If this IS a dovetail, we're done in all cases.

    if (olapType == pedDovetail) {
#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "overlapAlign::processHits()-- DOVETAIL return - score %f expected %f\n", _bestResult.score(), expectedScore);
#endif
      return(true);
    }

    //  Is this still a decent overlap?  Continue on to the next seed.  Decent if olapScore
    //  is at least 1/2 of the bestScore.

    if (0.5 * _bestResult.score() < olapScore) {
#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "overlapAlign::processHits()--  - decent score, keep looking: score %f > 0.5 * expected = %f\n",
              _bestResult.score(), 0.5 * expectedScore);
#endif
      continue;
    }

    //  Nope, this overlap is crap.  Assume that the rest of the seeds are crap too and give up.

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "overlapAlign::processHits()-- REST_CRAP return - score %f expected %f\n", _bestResult.score(), expectedScore);
#endif
    return(_bestResult.score() >= 0.5 * expectedScore);
  }

  //  We ran out of seeds to align.  No overlap found.

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "overlapAlign::processHits()-- NO_SEEDS return - score %f expected %f\n", _bestResult.score(), expectedScore);
#endif
  return(_bestResult.score() >= 0.5 * expectedScore);
}



void
overlapAlign::display(bool withAlign) {

  fprintf(stderr, "A %5u - %5u %s B %5u - %5u  olap length %d erate %6.4f type %s\n",
          _bestResult._aLo, _bestResult._aHi,
          _bFlipped ? "<--" : "-->",
          _bestResult._bLo, _bestResult._bHi,
          length(), erate(), toString(type()));

  if (withAlign)
    Display_Alignment(astr() + _bestResult._aLo, _bestResult._aHi - _bestResult._aLo,
                      bstr() + _bestResult._bLo, _bestResult._bHi - _bestResult._bLo,
                      delta(),
                      deltaLen());
}
