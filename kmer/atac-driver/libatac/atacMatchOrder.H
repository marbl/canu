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
    _matchIIDtoIdx = new u32bit      [_matchesLen];
    for (u32bit i=0; i<_matchesLen; i++) {
      _matches[i] = ML->getMatch(i);
      _matchIIDtoIdx[i] = ~u32bitZERO;
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

  u32bit            numMatches(void)      { return(_matchesLen);  };
  u32bit            numberOfMatches(void) { return(_matchesLen);  };

  //  Return match i in our list.  Common usage would be
  //    getMatch(index(iid) + 1)  //  return the match after the one we have
  //    getMatch(index(iid) - 1)
  //
  atacMatch        *getMatch(u32bit i)    {
    return(_matches[i]);
  };
  atacMatch        *operator[](u32bit i)  {
    return(_matches[i]);
  };

  //  Return the index, in our sorted list, of the matchiid supplied.
  //
  u32bit            index(u32bit matchiid) {
    return(_matchIIDtoIdx[matchiid]);
  };

  //  Merge the r match into the l match.  The l match gets the result.
  void              mergeMatches(atacMatch *l, atacMatch *r, u32bit mergeuid);

  void              sortA(u32bit first=0, u32bit len=0);
  void              sortB(u32bit first=0, u32bit len=0);
  void              sortDiagonal(u32bit first=0, u32bit len=0);
  void              sortMatchUID(u32bit first=0, u32bit len=0);
  void              sortParentUID(u32bit first=0, u32bit len=0);

private:
  u32bit           _matchesLen;
  u32bit           _matchesMax;
  atacMatch      **_matches;

  atacMatchList   *_ML;

  void              updateIndex(void) {
    for (u32bit i=0; i<_matchesLen; i++)
      _matchIIDtoIdx[i] = ~u32bitZERO;
    for (u32bit i=0; i<_matchesLen; i++)
      _matchIIDtoIdx[_matches[i]->matchiid] = i;
  };

  u32bit          *_matchIIDtoIdx;
};


#endif  //  ATAC_MATCHORDER_H
