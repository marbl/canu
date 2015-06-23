#include "overlapAlign.H"

#include "kMer.H"
#include "merStream.H"


//  Set the min/max diagonal we will accept seeds for.  It's just the min/max diagonal for the
//  two endpoints extended by half the erate times the align length.
//
//  Returns false if there is no chance of an overlap - the computed alignment length is too small.
//
bool
overlapAlign::findMinMaxDiagonal(int32  minLength) {

  _minDiag = 0;
  _maxDiag = 0;

  int32  aALen = _aHi - _aLo;
  int32  bALen = _bHi - _bLo;

  int32  alignLen = (aALen < bALen) ? bALen : aALen;

  if (alignLen < minLength)
    return(false);

  int32  bgnDiag = _aLo - _bLo;
  int32  endDiag = _aHi - _bHi;

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

  //fprintf(stderr, "findMinMaxDiagonal()--  min %d max %d -- span %d -- alignLen %d\n", _minDiag, _maxDiag, _maxDiag-_minDiag, alignLen);

  return(true);
}




void
overlapAlign::initializeConstants(void) {

  if (merMask[32] == 0xffffffffffffffffllu)
    return;

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


void
overlapAlign::fastFindMersA(bool dupIgnore) {

  uint64   mer = 0x0000000000000000llu;
  uint64   val = 0xffffffffffffffffllu;

  uint64   mask = merMask[_merSize];

  //  Restict the mers we seed with the those in the overlap (plus a little wiggle room).  If this
  //  isn't done, overlaps in repeats are sometimes lost.  In the A read, mers will drop out (e.g.,
  //  if a repeat at both the 5' and 3' end, and we have an overlap on one end only).  In the B
  //  read, the mers can now be on the wrong diagonal.

  int32  bgn = _aLo - _merSize - _merSize;
  int32  end = _aHi + _merSize;

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

  assert(_bLo < _bHi);

  int32  bgn = _bLo - _merSize - _merSize;
  int32  end = _bHi + _merSize;

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

  initializeConstants();

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

  if (_bMap.size() == 0)
    return(false);

  return(true);
}


//  For unique mers in B, if the mer is also unique in A, add a hit.

void
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
}


void
overlapAlign::chainHits(void) {

  //  Sort by aPos (actually by length, then by aPos, but length is constant here).

  sort(_rawhits.begin(), _rawhits.end());

#if 0
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

#if 0
    //  Allows non-overlapping (mismatch only) hits - but breaks seeding of alignment
    //           --------- end of next hit -----------   ----- end of existing hit -----
    int32   da = _rawhits[rr].aBgn + _rawhits[rr].tLen - _hits[hh].aBgn - _hits[hh].tLen;
    int32   db = _rawhits[rr].bBgn + _rawhits[rr].tLen - _hits[hh].bBgn - _hits[hh].tLen;

    assert(da > 0);

    if (da == db)
      merge = true;
#endif

#if 1
    //  Requires overlapping hits - full block of identity
    //           - start of next -   ------- end of existing -------
    int32   da = _rawhits[rr].aBgn - _hits[hh].aBgn - _hits[hh].tLen;
    int32   db = _rawhits[rr].bBgn - _hits[hh].bBgn - _hits[hh].tLen;

    if ((da < 0) && (da == db))
      merge = true;
#endif

    if (merge) {
#if 0
      fprintf(stderr, "MERGE HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
#endif
      _hits[hh].tLen = _rawhits[rr].aBgn + _rawhits[rr].tLen - _hits[hh].aBgn;
    } else {
#if 0
      fprintf(stderr, "NEW   HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
#endif
      _hits.push_back(_rawhits[rr]);
    }
  }

  //  Sort by longest

  sort(_hits.begin(), _hits.end());

#if 0
  for (uint32 hh=0; hh<_hits.size(); hh++) {
    fprintf(stderr, "hit %02u %5d-%5d diag %d len %3u\n",
            hh,
            _hits[hh].aBgn, _hits[hh].bBgn,
            _hits[hh].aBgn - _hits[hh].bBgn,
            _hits[hh].tLen);
  }
#endif
}



bool
overlapAlign::processHits(void) {

  //  Recompute.

  for (uint32 hh=0; hh<_hits.size(); hh++) {
    Match_Node_t  match;

    match.Start  = _hits[hh].aBgn;    //  Begin position in a
    match.Offset = _hits[hh].bBgn;    //  Begin position in b
    match.Len    = _hits[hh].tLen;    //  tLen can include mismatches if alternate scoring is used!
    match.Next   = 0;                 //  Not used here

    int32      errors  = 0;
    Overlap_t  ovltype = _editDist->Extend_Alignment(&match,         //  Initial exact match, relative to start of string
                                                     _aStr, _aLen,
                                                     _bStr, _bLen,
                                                     _aLo,  _aHi,    //  Output: Regions which the match extends
                                                     _bLo,  _bHi,
                                                     errors,
                                                     _partialOverlaps);

    _aHi++;  //  Add one to the end point because Extend_Alignment returns the base-based coordinate.
    _bHi++;

    _olapLen  = min(_aHi - _aLo, _bHi - _bLo);
    _olapQual = (double)errors / _olapLen;

    if (_olapLen < 40)
      ovltype = NONE;

    if (_bFlipped == true) {
      _bLo = _bLen - _bLo;  //  Now correct for the original forward sequence
      _bHi = _bLen - _bHi;  //  Done early just for the print below
    }

#if 0
    fprintf(stdout, "hit %2u a %5d b %5d -- ORIG A %6u %5d-%5d %s B %6u %5d-%5d  %.4f -- REALIGN type %d A %5d-%5d (%5d)  B %5d-%5d (%5d)  errors %4d  erate %6.4f = %6u / %6u deltas %d %d %d %d %d%s\n",
            hh, _hits[hh].aBgn, _hits[hh].bBgn,
            _aID, _aLoOrig, _aHiOrig,
            _bFlipped ? "<-" : "->",
            _bID, _bLoOrig, _bHiOrig,
            erate(),
            ovltype,
            _aLo, _aHi, _aLen,
            _bLo, _bHi, _bLen,
            errors,
            _olapQual, errors, _olapLen,
            _editDist->Left_Delta[0],
            _editDist->Left_Delta[1],
            _editDist->Left_Delta[2],
            _editDist->Left_Delta[3],
            _editDist->Left_Delta[4],
            ((ovltype != DOVETAIL) && (_partialOverlaps == false)) ? "  FAILED" : "");
#endif

    //  If not a good overlap, keep searching for one.  If partial overlaps, ovltype isn't set, and
    //  we need to do something else to reject.

    if ((ovltype != DOVETAIL) && (_partialOverlaps == false))
      continue;

    //  A good overlap.  Stop searching and return what we found.

    assert(_editDist->Right_Delta_Len == 0);

    _deltaLen = _editDist->Left_Delta_Len;

    if (_deltaMax < _deltaLen) {
      delete [] _delta;
      _deltaMax = _deltaLen + 4096;
      _delta    = new int32 [_deltaMax];
    }

    memcpy(_delta, _editDist->Left_Delta, sizeof(int32) * _deltaLen);

    return(true);
  }

  //  We ran out of seeds to align.  No overlap found.
  return(false);
}

