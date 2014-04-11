// This file is part of A2Amapper.
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Brian Walenz
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received (LICENSE.txt) a copy of the GNU General Public 
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef ATAC_MATCHORDER_H
#define ATAC_MATCHORDER_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"


class atacMatchOrder {
private:
  void           initialize(atacMatchList *ML) {
    _matchesLen    = ML->numberOfMatches();
    _matchesMax    = ML->numberOfMatches();
    _matches       = new atacMatch * [_matchesLen];
    _matchIIDtoIdx = new uint32      [_matchesLen];
    for (uint32 i=0; i<_matchesLen; i++) {
      _matches[i] = ML->getMatch(i);
      _matchIIDtoIdx[i] = ~uint32ZERO;
    }
  };
public:
  atacMatchOrder(atacMatchList *ML) {
    initialize(ML);
  };
  atacMatchOrder(atacMatchList &ML) {
    initialize(&ML);
  };
  ~atacMatchOrder() {
    delete [] _matches;
  };

  uint32            numMatches(void)      { return(_matchesLen);  };
  uint32            numberOfMatches(void) { return(_matchesLen);  };

  //  Return match i in our list.  Common usage would be
  //    getMatch(index(iid) + 1)  //  return the match after the one we have
  //    getMatch(index(iid) - 1)
  //
  atacMatch        *getMatch(uint32 i)    {
    return(_matches[i]);
  };
  atacMatch        *operator[](uint32 i)  {
    return(_matches[i]);
  };

  //  Return the index, in our sorted list, of the matchiid supplied.
  //
  uint32            index(uint32 matchiid) {
    return(_matchIIDtoIdx[matchiid]);
  };

  //  Merge the r match into the l match.  The l match gets the result.
  void              mergeMatches(atacMatch *l, atacMatch *r, uint32 mergeuid);

  void              sortA(uint32 first=0, uint32 len=0);
  void              sortB(uint32 first=0, uint32 len=0);
  void              sortDiagonal(uint32 first=0, uint32 len=0);
  void              sortMatchUID(uint32 first=0, uint32 len=0);
  void              sortParentUID(uint32 first=0, uint32 len=0);

private:
  uint32           _matchesLen;
  uint32           _matchesMax;
  atacMatch      **_matches;

  atacMatchList   *_ML;

  void              updateIndex(void) {
    for (uint32 i=0; i<_matchesLen; i++)
      _matchIIDtoIdx[i] = ~uint32ZERO;
    for (uint32 i=0; i<_matchesLen; i++)
      _matchIIDtoIdx[_matches[i]->matchiid] = i;
  };

  uint32          *_matchIIDtoIdx;
};


#endif  //  ATAC_MATCHORDER_H
