
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

#include "runtime.H"
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
};


void
tgTig::stashContains(double  maxCov, tgTigStashed &S) {

  //  Initialize.
  //    Declare that we have no stashed reads.
  //    Clear the return statistics.

  _stashed    = nullptr;
  _stashedLen = 0;
  _stashedMax = 0;

  S.clear();

  if (_childrenLen == 1)
    return;

  //  Sort the original children by position.

  std::sort(_children, _children + _childrenLen);

  //  Decide which children to save.

  readInfo        *posLen  = new readInfo [_childrenLen];   //  Sorting by length of child

  //  Flag the read for stashing if it doesn't extend the tig.

  int32         hiEnd = -1;

  for (uint32 ci=0; ci<_childrenLen; ci++) {
    int32  lo = _children[ci].min();
    int32  hi = _children[ci].max();

    if (hi <= hiEnd) {
      posLen[ci].idx = ci;
      posLen[ci].len = hi - lo;
      posLen[ci].use = false;

      S.nStsh += 1;
      S.bStsh += posLen[ci].len;
    }

    else {
      posLen[ci].idx = ci;
      posLen[ci].len = hi - lo;
      posLen[ci].use = true;

      S.nBack += 1;
      S.bBack += posLen[ci].len;
    }

    hiEnd = max(hi, hiEnd);
  }

  //  Sort by length, longest first, then verify we're sorted.

  auto longestFirst = [](readInfo const &A, readInfo const &B) {
                        return(A.len > B.len);
                      };

  std::sort(posLen, posLen + _childrenLen, longestFirst);

  for (uint32 ci=1; ci<_childrenLen; ci++)
    assert(posLen[ci-1].len >= posLen[ci].len);

  //  Save the longer of the contained reads, until we reach the maximum
  //  coverage desired.  

  uint64  bLimit = (uint64)floor(maxCov * hiEnd);

  for (uint32 ci=0; ci<_childrenLen; ci++) {
    if (posLen[ci].use == true)            //  Already a backbone read.
      continue;                            //  Skip this read.

    if (S.bCont + S.bBack > bLimit)        //  Exceeded coverage limit.
      break;                               //  Bail.

    posLen[ci].use = true;                 //  Do not stash the read.

    S.nStsh -= 1;                          //  Do not stash this read.
    S.bStsh -= posLen[ci].len;             //

    S.nCont += 1;                          //  And do use it for consensus.
    S.bCont += posLen[ci].len;             //
  }

  //  Make a new list of reads, while saving the original, if there reads
  //  we want to stash.

  if (S.nStsh > 0) {
    _stashedLen = _childrenLen;                      //  Save the original list
    _stashedMax = _childrenMax;                      //  so we can restore it later.
    _stashed    = _children;

    _childrenLen  = 0;                               //  Allocate a new list for
    _childrenMax  = S.nBack + S.nCont;               //  exactly what we need to save.
    _children     = new tgPosition [_childrenMax];

    for (uint32 ci=0; ci<_stashedLen; ci++)
      if (posLen[ci].use == true)                                //  If used, we want to keep the
        _children[_childrenLen++] = _stashed[posLen[ci].idx];    //  read, so copy it to the new list.
  }

  // since we added the reads using length sorted order, re-sort them by position to make everyone downstream  happy
  std::sort(_children, _children + _childrenLen);

  //  Cleanup and return the statistics.

  delete [] posLen;
}



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
    oldMax = max(oldMax, _stashed[ci].max());

  for (uint32 ci=0; ci<_childrenLen; ci++)
    newMax = max(newMax, _children[ci].max());

  if (oldMax > 0)
    sf = (double)newMax / oldMax;

  //  Build a map from child ID to it's current position.

  map<uint32, uint32>   idmap;

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

