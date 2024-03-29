
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "sqStore.H"
#include "tgStore.H"

#include <map>
#include <algorithm>


//  Replace the children list in tig with one that has fewer contains.
//  The original list is saved in the tig.
//  When contains are 'unstashed', positions are updated.



struct readInfo {
  uint32    idx;   //  index of the read in the original unsorted _children
  int32     len;   //  length of the read
  bool      use;   //  true if the read is not contained in some other read
  bool      ign;   //  true if the read should not be used at all in consensus
};



void
tgTig::filterContains(double  maxCov, bool enableStash) {

  clearStash();

  if (_childrenLen == 1)  //  And do nothing if only a single read.
    return;

  //  Flag the read for stashing if it doesn't extend the tig or is to be ignored.

  std::sort(_children, _children + _childrenLen);

  readInfo   *posLen = new readInfo [_childrenLen];   //  Sorting by length of child
  int32       hiEnd  = -1;

  for (uint32 ci=0; ci<_childrenLen; ci++) {
    int32  lo = _children[ci].min();
    int32  hi = _children[ci].max();
    bool   ign = _children[ci].skipConsensus();    //  Completely ignore the read if told to.
    bool   use = (ign == false) && (hiEnd < hi);   //  Use the read if it extends the tig.

    if (ign)  assert(use == false);
    if (use)  assert(hiEnd < hi);

    if (use)                      //  Extend the useful length of the tig.
      hiEnd = hi;

    posLen[ci].idx = ci;          //  Save the read info.
    posLen[ci].len = hi - lo;
    posLen[ci].use = use;
    posLen[ci].ign = ign;

    if      (ign == true) { _stashIgnr[0] += 1;  _stashIgnr[1] += posLen[ci].len; }   //  Count where
    else if (use == true) { _stashBack[0] += 1;  _stashBack[1] += posLen[ci].len; }   //  reads go.
    else                  { _stashStsh[0] += 1;  _stashStsh[1] += posLen[ci].len; }   //
  }

  //  Save the longer of the contained reads, until we reach the maximum
  //  coverage desired.  

  auto longestFirst = [](readInfo const &A, readInfo const &B) {
                        if (A.len != B.len)  return A.len > B.len;
                        else                 return A.idx < B.idx;
                      };

  std::sort(posLen, posLen + _childrenLen, longestFirst);   //  Sort by length.

  for (uint32 ci=1; ci<_childrenLen; ci++)                  //  and verify.
    assert(posLen[ci-1].len >= posLen[ci].len);

  uint64  bLimit = (uint64)floor(maxCov * hiEnd);

  for (uint32 ci=0; ci<_childrenLen; ci++) {
    if ((posLen[ci].use == true) ||              //  Already a backbone read or
        (posLen[ci].ign == true))                //  something we want to ignore;
      continue;                                  //  skip this read.

    if (_stashCont[1] + _stashBack[1] > bLimit)  //  Exceeded coverage limit.
      break;                                     //  Bail.
    
    posLen[ci].use = true;                       //  Do not stash the read.

    _stashStsh[0] -= 1;                          //  Do not stash this read.
    _stashStsh[1] -= posLen[ci].len;             //

    _stashCont[0] += 1;                          //  Do use it for consensus.
    _stashCont[1] += posLen[ci].len;             //
  }

  //  Transfer the 'use' flag to the reads in the tig.

  for (uint32 ci=0; ci<_childrenLen; ci++)
    _children[ posLen[ci].idx ].skipConsensus( !posLen[ci].use );

  //  Optionally remove the unused reads from the tig.  Read correction (with
  //  falconsense) uses this mode; utgcns does not.

  if ((enableStash) && (_stashStsh[0] > 0)) {
    _stashedLen = _childrenLen;                      //  Save the original list
    _stashedMax = _childrenMax;                      //  so we can restore it later.
    _stashed    = _children;

    uint32 c = 0;
    for (uint32 ci=0; ci<_stashedLen; ci++)          //  Count how many reads we
      if (posLen[ci].use == true)                    //  want to save.
        c++;

    _childrenLen  = 0;                               //  Allocate a new list for
    _childrenMax  = c;                               //  exactly what we need to save.
    _children     = new tgPosition [_childrenMax];

    for (uint32 ci=0; ci<_stashedLen; ci++)          //  Copy used reads to the new list.
      if (posLen[ci].use == true)
        _children[_childrenLen++] = _stashed[posLen[ci].idx];

    std::sort(_children, _children + _childrenLen);  //  Sort by position.
  }

  //  Cleanup.

  delete [] posLen;
}


#if 0
void
tgTig::unstashContains(void) {

  if (_stashed == NULL)   //  If no saved list, nothing to unstash.
    return;

  //  For the reads we 'stashed', we need to compute a position in the new
  //  consensus sequence.  We'll just linearly scale the positions.

  int32   oldMax = 0;
  int32   newMax = 0;
  double  sf     = 1.0;

  for (uint32 ci=0; ci<_stashedLen; ci++)
    oldMax = std::max(oldMax, _stashed[ci].max());

  for (uint32 ci=0; ci<_childrenLen; ci++)
    newMax = std::max(newMax, _children[ci].max());

  if (oldMax > 0)
    sf = (double)newMax / oldMax;

  //  Build a map from child ID to it's current position.

  std::map<uint32, uint32>   idmap;

  for (uint32 ci=0; ci < _childrenLen; ci++)
    idmap[_children[ci].ident()] = ci;

  //  Over all the reads in the original list (currently held in _stashed),
  //  update the position, either from its current position out of consensus
  //  or by extrapolating the original position.

  for (uint32 ci=0; ci<_stashedLen; ci++) {
    uint32  id = _stashed[ci].ident();

    if (idmap.find(id) != idmap.end()) {      //  If we find the ID in the current list,
      _stashed[ci] = _children[idmap[id]];    //  simply copy the new position data.
      idmap.erase(id);
    }

    else {
      int32  nmin = sf * _stashed[ci].min();   //  If not, scale the positions based
      int32  nmax = sf * _stashed[ci].max();   //  on new vs old tig length.

      if (nmin > newMax)  nmin = newMax;       //  But don't let them exceed the limit!
      if (nmax > newMax)  nmax = newMax;

      _stashed[ci].setMinMax(nmin, nmax);
    }
  }

  //  Make sure that we updated all the children.

  if (idmap.empty() == false)
    fprintf(stderr, "Failed to unstash the contained reads.  Still have " F_SIZE_T " reads unplaced.\n",
            idmap.size());
  assert(idmap.empty() == true);

  //  Throw out the stashed list, and restore the original.

  delete [] _children;

  _childrenLen  = _stashedLen;
  _childrenMax  = _stashedMax;
  _children     = _stashed;

  _stashedLen = 0;
  _stashedMax = 0;
  _stashed    = NULL;
}
#endif



