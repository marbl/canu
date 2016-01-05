
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/overlapInCore/overlapAlign.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2015-JUN-17 to 2015-AUG-25
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-13
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "NDalign.H"

#include "kMer.H"
#include "merStream.H"

#include "stddev.H"

#include "Display_Alignment.H"

#undef  DEBUG_ALGORITHM         //  Some details.
#undef  DEBUG_HITS              //  Lots of details (chainHits())

#undef  SEED_NON_OVERLAPPING    //  Allow mismatches in seeds
#define SEED_OVERLAPPING

#define DISPLAY_WIDTH         250


NDalign::NDalign(pedAlignType   alignType,
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

  _editDist = new NDalgorithm(_alignType, _maxErate);

  _minDiag  = 0;
  _maxDiag  = 0;

  _merSize  = 0;

  _topDisplay = NULL;
  _botDisplay = NULL;
  _resDisplay = NULL;

  //  Initialize Constants

  for (uint32 ii=0; ii<256; ii++) {
    acgtToBit[ii] = 0x00;
    acgtToVal[ii] = 0x03;
  }

  //  Do NOT include lowercase letters here.  These are 'don't care' matches
  //  in consensus, and are usually matches we want to be skipping.
  //
  //  ATGtTTgaTGACC vs ATGtTTgaTGACC
  //  ATG-TT--TGACC    ATGTTT---GACC
  //
  //  While, yes, the second form is simpler (and would be correct had the last T been
  //  optional) the first is what consensus is expecting to find.

#if 0
  acgtToBit['a'] = acgtToBit['A'] = 0x00;  //  Bit encoding of ACGT
  acgtToBit['c'] = acgtToBit['C'] = 0x01;
  acgtToBit['g'] = acgtToBit['G'] = 0x02;
  acgtToBit['t'] = acgtToBit['T'] = 0x03;

  acgtToVal['a'] = acgtToVal['A'] = 0x00;  //  Word is valid if zero
  acgtToVal['c'] = acgtToVal['C'] = 0x00;
  acgtToVal['g'] = acgtToVal['G'] = 0x00;
  acgtToVal['t'] = acgtToVal['T'] = 0x00;
#else
  acgtToBit['A'] = 0x00;  //  Bit encoding of ACGT
  acgtToBit['C'] = 0x01;
  acgtToBit['G'] = 0x02;
  acgtToBit['T'] = 0x03;

  acgtToVal['A'] = 0x00;  //  Word is valid if zero
  acgtToVal['C'] = 0x00;
  acgtToVal['G'] = 0x00;
  acgtToVal['T'] = 0x00;
#endif

  merMask[0] = 0x0000000000000000llu;
  merMask[1] = 0x0000000000000003llu;

  for (uint32 ii=2; ii<33; ii++)
    merMask[ii] = (merMask[ii-1] << 2) | 0x03;

  assert(merMask[ 6] == 0x0000000000000fffllu);
  assert(merMask[17] == 0x00000003ffffffffllu);
  assert(merMask[26] == 0x000fffffffffffffllu);
  assert(merMask[32] == 0xffffffffffffffffllu);
}



NDalign::~NDalign() {
  delete    _editDist;
  delete [] _bRev;

  delete [] _topDisplay;
  delete [] _botDisplay;
  delete [] _resDisplay;
}



void
NDalign::initialize(uint32 aID, char *aStr, int32 aLen, int32 aLo, int32 aHi,
                    uint32 bID, char *bStr, int32 bLen, int32 bLo, int32 bHi, bool bFlipped) {

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "NDalign::initialize()--  A %u %d-%d  B %u %d-%d\n", aID, aLo, aHi, bID, bLo, bHi);
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

  _hitr = UINT32_MAX;

  _bestResult.clear();
}




//  Set the min/max diagonal we will accept seeds for.  It's just the min/max diagonal for the
//  two endpoints extended by half the erate times the align length.
//
//  Returns false if there is no chance of an overlap - the computed alignment length is too small.
//
//  The adjusts change the sequence start and end position (bgn + bgnAdj, end - endAdj),
//  and are intended to compensate for supplying sequences larger than absolutely necessary
//  (for when the exact boundaries aren't known).

bool
NDalign::findMinMaxDiagonal(int32  minLength,
                                 uint32 AbgnAdj, uint32 AendAdj,
                                 uint32 BbgnAdj, uint32 BendAdj) {

  _minDiag = 0;
  _maxDiag = 0;

  int32  aALen = (_aHiOrig - AendAdj) - (_aLoOrig + AbgnAdj);
  int32  bALen = (_bHiOrig - BendAdj) - (_bLoOrig + BbgnAdj);

  int32  alignLen = (aALen < bALen) ? bALen : aALen;

  if (alignLen < minLength)
    return(false);

  int32  bgnDiag = (_aLoOrig + AbgnAdj) - (_bLoOrig + BbgnAdj);
  int32  endDiag = (_aHiOrig - AendAdj) - (_bHiOrig - BendAdj);

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
    fprintf(stderr, "NDalign::findMinMaxDiagonal()-- ERROR: _minDiag=%d >= _maxDiag=%d\n", _minDiag, _maxDiag);
  assert(_minDiag <= _maxDiag);

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "NDalign::findMinMaxDiagonal()--  A: _aLoOrig %d + AbgnAdj %d -- _aHiOrig %d - AendAdj %d\n", _aLoOrig, AbgnAdj, _aHiOrig, AendAdj);
  fprintf(stderr, "NDalign::findMinMaxDiagonal()--  B: _bLoOrig %d + BbgnAdj %d -- _bHiOrig %d - BendAdj %d\n", _bLoOrig, BbgnAdj, _bHiOrig, BendAdj);
  fprintf(stderr, "NDalign::findMinMaxDiagonal()--  min %d max %d -- span %d -- alignLen %d\n", _minDiag, _maxDiag, _maxDiag-_minDiag, alignLen);
#endif

  return(true);
}




void
NDalign::fastFindMersA(bool dupIgnore) {

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
NDalign::fastFindMersB(bool dupIgnore) {

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
NDalign::findSeeds(bool dupIgnore) {

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
    fprintf(stderr, "NDalign::findSeeds()--  No seeds found.\n");
#endif
    return(false);
  }

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "NDalign::findSeeds()--  Found %u seeds.\n", _bMap.size());
#endif
  return(true);
}


//  For unique mers in B, if the mer is also unique in A, add a hit.

bool
NDalign::findHits(void) {

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
      fprintf(stderr, "NDalign::findHits()-- kmer "F_X64" apos - bpos = %d  _minDiag = %d  _maxDiag = %d\n",
              kmer, apos-bpos, _minDiag, _maxDiag);
    assert(apos - bpos >= _minDiag);  //  ...these too.
    assert(apos - bpos <= _maxDiag);

    _rawhits.push_back(exactMatch(apos, bpos, _merSize));
  }

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "NDalign::findHits()--  Found %u hits.\n", _rawhits.size());
#endif

  return(true);
}


bool
NDalign::chainHits(void) {

  //  Sort by aPos (actually by length, then by aPos, but length is constant here).

  sort(_rawhits.begin(), _rawhits.end());

#ifdef DEBUG_HITS
  for (uint32 rr=0; rr<_rawhits.size(); rr++)
    if (_rawhits[rr].aBgn - _rawhits[rr].bBgn == 2289)
      fprintf(stderr, "NDalign::chainHits()-- HIT: %d - %d diag %d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn);
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
      //fprintf(stderr, "NDalign::chainHits()-- MERGE HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
      _hits[hh].tLen = _rawhits[rr].aBgn + _rawhits[rr].tLen - _hits[hh].aBgn;
    } else {
      //fprintf(stderr, "NDalign::chainHits()-- NEW   HIT: %d - %d diag %d  da=%d db=%d\n", _rawhits[rr].aBgn, _rawhits[rr].bBgn, _rawhits[rr].aBgn - _rawhits[rr].bBgn, da, db);
      _hits.push_back(_rawhits[rr]);
    }
  }

  //  Sort by longest

  sort(_hits.begin(), _hits.end());

#ifdef DEBUG_HITS
  for (uint32 hh=0; hh<_hits.size(); hh++) {
    fprintf(stderr, "NDalign::chainHits()-- hit %02u %5d-%5d diag %d len %3u\n",
            hh,
            _hits[hh].aBgn, _hits[hh].bBgn,
            _hits[hh].aBgn - _hits[hh].bBgn,
            _hits[hh].tLen);
  }
#endif

#ifdef DEBUG_ALGORITHM
  fprintf(stderr, "NDalign::chainHits()--  Found %u chains of hits.\n", _hits.size());
#endif

  return(true);
}



bool
NDalign::makeNullHit(void) {

  _hits.push_back(exactMatch(0, 0, 0));

  return(true);
}



bool
NDalign::processHits(void) {

  //  If the first time here, set the hit iterator to zero, otherwise move to the next one.
  //  And then return if there are no more hits to iterate over.

  if (_hitr == UINT32_MAX)
    _hitr = 0;
  else
    _hitr++;

  if (_hitr >= _hits.size())
    return(false);

  //  While hits, process them.
  //
  //  If a good hit is found, return, leaving hitr as is.  The next time we enter this function,
  //  we'll increment hitr and process the next hit.  If no good hit is found, we iterate the loop
  //  until a good one is found, or we run out of hits.

  for (; _hitr < _hits.size(); _hitr++) {
    Match_Node_t  match;

    match.Start  = _hits[_hitr].aBgn;    //  Begin position in a
    match.Offset = _hits[_hitr].bBgn;    //  Begin position in b
    match.Len    = _hits[_hitr].tLen;    //  tLen can include mismatches if alternate scoring is used!
    match.Next   = 0;                 //  Not used here


#ifdef SEED_NON_OVERLAPPING
    match.Offset = _merSize;  //  Really should track this in the hits, oh well.
#endif

#ifdef DEBUG_ALGORITHM
    fprintf(stderr, "NDalign::processHits()--\n");
    fprintf(stderr, "NDalign::processHits()-- Extend_Alignment Astart %d Bstart %d length %d\n", match.Start, match.Offset, match.Len);
    fprintf(stderr, "NDalign::processHits()--\n");

    fprintf(stderr, ">A len=%u\n", _aLen);
    fwrite(_aStr, sizeof(char), _aLen, stderr);
    fprintf(stderr, "\n");

    fprintf(stderr, ">B len=%u\n", _bLen);
    fwrite(_bStr, sizeof(char), _bLen, stderr);
    fprintf(stderr, "\n");
#endif

    int32  aLo=0, aHi=0;
    int32  bLo=0, bHi=0;

    char   origScore[1024];

    pedOverlapType  olapType = _editDist->Extend_Alignment(&match,         //  Initial exact match, relative to start of string
                                                           _aStr, _aLen,
                                                           _bStr, _bLen,
                                                           aLo,   aHi,    //  Output: Regions which the match extends
                                                           bLo,   bHi);

    aHi++;  //  Add one to the end point because Extend_Alignment returns the base-based coordinate.
    bHi++;

    //  Is this a better overlap than what we have?  Save it and update statistics.

    if (((score() <  _editDist->score())) ||
        ((score() <= _editDist->score()) && (length() > ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2))) {

#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "NDalign::processHits()--\n");
      fprintf(stderr, "NDalign::processHits()-- %s\n", (length() > 0) ? "Save better alignment" : "First alignment");
      fprintf(stderr, "NDalign::processHits()--\n");

      sprintf(origScore, "NDalign::processHits()--  OLD length %u erate %f score %u (%d-%d %d-%d)\n",
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
#endif

      _bestResult.save(aLo, aHi, bLo, bHi, _editDist->score(), olapType, _editDist->Left_Delta_Len, _editDist->Left_Delta);

      display("NDalign::processHits()-- ", false);

      _bestResult.setErate(1.0 - (double)(_matches + _gapmatches) / (length() - _freegaps));

#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "%sNDalign::processHits()--  NEW length %u erate %f score %u (%d-%d %d-%d)\n",
              origScore,
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
#endif

    } else {
      olapType = pedBothBranch;

#ifdef DEBUG_ALGORITHM
      fprintf(stderr, "NDalign::processHits()--  DON'T save alignment - OLD length %u erate %f score %u (%d-%d %d-%d) ",
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
      fprintf(stderr, "NDalign::processHits()--  NEW length %u score %u coords %u-%u %u-%u\n",
              ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2,
              _editDist->score(),
              aLo, aHi, bLo, bHi);
#endif
    }

    //  If a dovetail, we're done.  Let the client figure out if the quality is good.

    if (olapType == pedDovetail)
      return(true);

  }  //  Over all seeds.

  //  No more seeds to align.  Did we save an alignment?

  return(score() > 0);
}


#define ABS(x)  (((x) < 0) ? -(x) : (x))


//  A simple scan for blocks of gaps.  The delta values here are relative to the last gap,
//  so a block of too much gap will have N delta values with small values.
//
bool
NDalign::scanDeltaForBadness(bool verbose, bool showAlign) {
  double  ema       = 0.0;
  double  alpha     = 0.001;  //  Smaller will increase the averaging length

  int32   badBlocks = 0;

  for (uint32 ii=0; _resDisplay[ii]; ii++) {
    if (_resDisplay[ii] == '^')
      ema = computeExponentialMovingAverage(alpha, ema, 1.0);
    else
      ema = computeExponentialMovingAverage(alpha, ema, 0.0);

    if (ema > 0.25)
      badBlocks++;
  }

  if ((verbose == true) && (badBlocks > 0)) {
    fprintf(stderr, "NDalign::scanForDeltaBadness()--  Potential bad alignment: found %d bad blocks (alpha %f)\n",
            badBlocks, alpha);

    if (showAlign == true)
      display("NDalign::scanForDeltaBadness()-- ", true);
  }

  return(badBlocks > 0);
}



void
NDalign::realignForward(bool displayAlgorithm, bool displayAlign) {
  Match_Node_t  match;

  match.Start  = abgn();    //  Begin position in a
  match.Offset = bbgn();    //  Begin position in b
  match.Len    = 0;
  match.Next   = 0;         //  Not used here

  int32  aLo=0, aHi=0;
  int32  bLo=0, bHi=0;

  char   origScore[1024];

  if (displayAlgorithm)
    fprintf(stderr, "NDalign::realignForward()--\n");

  pedOverlapType  olapType = _editDist->Extend_Alignment(&match,        //  Initial exact match, relative to start of string
                                                         _aStr, _aLen,
                                                         _bStr, _bLen,
                                                         aLo,   aHi,    //  Output: Regions which the match extends
                                                         bLo,   bHi);

  aHi++;  //  Add one to the end point because Extend_Alignment returns the base-based coordinate.
  bHi++;

  //  Is this a better overlap than what we have?

  if (((score() <  _editDist->score())) ||
      ((score() <= _editDist->score()) && (length() > ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2))) {

    if (displayAlgorithm) {
      fprintf(stderr, "NDalign::realignForward()--\n");
      fprintf(stderr, "NDalign::realignForward()--  Save better alignment\n");
      fprintf(stderr, "NDalign::realignForward()--\n");

      sprintf(origScore, "NDalign::realignForward()--  OLD length %u erate %f score %u (%d-%d %d-%d)\n",
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    }

    _bestResult.save(aLo, aHi, bLo, bHi, _editDist->score(), olapType, _editDist->Left_Delta_Len, _editDist->Left_Delta);

    display("NDalign::realignForward()-- ", displayAlign);

    _bestResult.setErate(1.0 - (double)(_matches + _gapmatches) / (length() - _freegaps));

    if (displayAlgorithm) {
      fprintf(stderr, "%sNDalign::realignForward()--  NEW length %u erate %f score %u (%d-%d %d-%d)\n",
              origScore,
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    }
  }

  else if (displayAlgorithm) {
    fprintf(stderr, "NDalign::realignForward()-- Alignment no better   - OLD length %u erate %f score %u (%d-%d %d-%d)\n",
            length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    fprintf(stderr, "NDalign::realignForward()-- Alignment no better   - NEW length %u erate %f score %u (%d-%d %d-%d)\n",
            ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2, 0.0, _editDist->score(), aLo, aHi, bLo, bHi);
    //display("NDalign::realignForward(NB)--", aLo, aHi, bLo, bHi, _editDist->Left_Delta, _editDist->Left_Delta_Len, true, false);
  }
}



void
NDalign::realignBackward(bool displayAlgorithm, bool displayAlign) {
  Match_Node_t  match;

  match.Start  = aend();    //  Begin position in a
  match.Offset = bend();    //  Begin position in b
  match.Len    = 0;
  match.Next   = 0;         //  Not used here

  int32  aLo=0, aHi=0;
  int32  bLo=0, bHi=0;

  char   origScore[1024];

  if (displayAlgorithm)
    fprintf(stderr, "NDalign::realignBackward()--\n");

  pedOverlapType  olapType = _editDist->Extend_Alignment(&match,        //  Initial exact match, relative to start of string
                                                         _aStr, _aLen,
                                                         _bStr, _bLen,
                                                         aLo,   aHi,    //  Output: Regions which the match extends
                                                         bLo,   bHi);

  aHi++;  //  Add one to the end point because Extend_Alignment returns the base-based coordinate.
  bHi++;

  //  Is this a better overlap than what we have?

  if (((score() <  _editDist->score())) ||
      ((score() <= _editDist->score()) && (length() > ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2))) {

    if (displayAlign) {
      fprintf(stderr, "NDalign::realignBackward()--\n");
      fprintf(stderr, "NDalign::realignBackward()--  Save better alignment\n");
      fprintf(stderr, "NDalign::realignBackward()--\n");

      sprintf(origScore, "NDalign::realignBackward()--  OLD length %u erate %f score %u (%d-%d %d-%d)\n",
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    }

    _bestResult.save(aLo, aHi, bLo, bHi, _editDist->score(), olapType, _editDist->Left_Delta_Len, _editDist->Left_Delta);

    display("NDalign::realignBackward()-- ", displayAlign);

    _bestResult.setErate(1.0 - (double)(_matches + _gapmatches) / (length() - _freegaps));

    if (displayAlign) {
      fprintf(stderr, "%sNDalign::realignBackward()--  NEW length %u erate %f score %u (%d-%d %d-%d)\n",
              origScore,
              length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    }
  }

  else if (displayAlign) {
    fprintf(stderr, "NDalign::realignBackward()-- Alignment no better   - OLD length %u erate %f score %u (%d-%d %d-%d)\n",
            length(), erate(), score(), abgn(), aend(), bbgn(), bend());
    fprintf(stderr, "NDalign::realignBackward()-- Alignment no better   - NEW length %u erate %f score %u (%d-%d %d-%d)\n",
            ((aHi - aLo) + (bHi - bLo) + _editDist->Left_Delta_Len) / 2, 0.0, _editDist->score(), aLo, aHi, bLo, bHi);
    //display("NDalign::realignBackward(NB)--", aLo, aHi, bLo, bHi, _editDist->Left_Delta, _editDist->Left_Delta_Len, true, false);
  }
}






void
NDalign::display(char    *prefix,
                 int32    aLo,   int32   aHi,
                 int32    bLo,   int32   bHi,
                 int32   *delta, int32   deltaLen,
                 bool     displayIt,
                 bool     saveStats) {

  uint32 matches    = 0;
  uint32 errors     = 0;
  uint32 gapmatches = 0;
  uint32 freegaps   = 0;

  if (_topDisplay == NULL) {
    _topDisplay = new char [2 * AS_MAX_READLEN + 1];
    _botDisplay = new char [2 * AS_MAX_READLEN + 1];
    _resDisplay = new char [2 * AS_MAX_READLEN + 1];
  }

  char *a      = astr() + aLo;
  int32 a_len = aHi - aLo;

  char *b      = bstr() + bLo;
  int32 b_len = bHi - bLo;

  int32 top_len = 0;

  {
    int32 i = 0;
    int32 j = 0;

    for (int32 k = 0;  k < deltaLen;  k++) {
      for (int32 m = 1;  m < abs (delta[k]);  m++) {
        _topDisplay[top_len++] = a[i++];
        j++;
      }

      if (delta[k] < 0) {
        _topDisplay[top_len++] = '-';
        j++;

      } else {
        _topDisplay[top_len++] = a[i++];
      }
    }

    while (i < a_len && j < b_len) {
      _topDisplay[top_len++] = a[i++];
      j++;
    }

    _topDisplay[top_len] = 0;
  }


  int32 bot_len = 0;

  {
    int32 i = 0;
    int32 j = 0;

    for (int32 k = 0;  k < deltaLen;  k++) {
      for (int32 m = 1;  m < abs (delta[k]);  m++) {
        _botDisplay[bot_len++] = b[j++];
        i++;
      }

      if (delta[k] > 0) {
        _botDisplay[bot_len++] = '-';
        i++;

      } else {
        _botDisplay[bot_len++] = b[j++];
      }
    }

    while (j < b_len && i < a_len) {
      _botDisplay[bot_len++] = b[j++];
      i++;
    }

    _botDisplay[bot_len] = 0;
  }

  //  Compute stats, build the result display.

  {
    for (int32 i=0; (i < top_len) || (i < bot_len); i++) {
      char T = _topDisplay[i];
      char B = _botDisplay[i];

      char t = tolower(T);
      char b = tolower(B);

      //  A minor-allele match if the lower case matches, but upper case doesn't
      if      ((t == b) && (T != B)) {
        _resDisplay[i] = '\'';
        gapmatches++;
      }

      //  A free-gap if either is lowercase
      else if (((t == T) && (b == '-')) ||  //  T is lowercase, B is a gap
               ((b == B) && (t == '-'))) {  //  B is lowercase, T is a gap
        _resDisplay[i] = '-';
        freegaps++;
      }

      //  A match if the originals match
      else if ((t == b) && (T == B)) {
        _resDisplay[i] = ' ';
        matches++;
      }

      //  Otherwise, an error
      else {
        _resDisplay[i] = '^';
        errors++;
      }
    }

    _resDisplay[ max(top_len, bot_len) ] = 0;
  }

  //  Really display it?

  if (displayIt == true) {
    for (int32 i=0; (i < top_len) || (i < bot_len); i += DISPLAY_WIDTH) {
      fprintf(stderr, "%s%d\n", prefix, i);

      fprintf(stderr, "%sA: ", prefix);
      for (int32 j=0;  (j < DISPLAY_WIDTH) && (i+j < top_len);  j++)
        putc(_topDisplay[i+j], stderr);
      fprintf(stderr, "\n");

      fprintf(stderr, "%sB: ", prefix);
      for (int32 j=0; (j < DISPLAY_WIDTH) && (i+j < bot_len); j++)
        putc(_botDisplay[i+j], stderr);
      fprintf(stderr, "\n");

      fprintf(stderr, "%s   ", prefix);
      for (int32 j=0; (j<DISPLAY_WIDTH) && (i+j < bot_len) && (i+j < top_len); j++)
        putc(_resDisplay[i+j], stderr);
      fprintf(stderr, "\n");
    }
  }

  //  Update statistics?  Avoid a whole ton of if statements by always counting then changing only here.

  if (saveStats) {
    _matches    = matches;
    _errors     = errors;
    _gapmatches = gapmatches;
    _freegaps   = freegaps;
  }
}



void
NDalign::display(char *prefix, bool displayIt) {
  display(prefix,
          _bestResult._aLo,   _bestResult._aHi,
          _bestResult._bLo,   _bestResult._bHi,
          _bestResult._delta, _bestResult._deltaLen,
          displayIt,
          true);
}
