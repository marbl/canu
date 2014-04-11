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

#ifndef ATAC_MATCHLIST_H
#define ATAC_MATCHLIST_H

#include <string>
#include <map>
using namespace std;

#include "bio++.H"


class atacMatchList {
public:
  atacMatchList();
  ~atacMatchList();

  void           add(atacMatch &m);
  //void           delete();

  uint32         numMatches(void)      { return(_matchesLen);  };
  uint32         numberOfMatches(void) { return(_matchesLen);  };

  atacMatch     *getMatch(uint32 i)    { return(_matches + i); };
  atacMatch     *operator[](uint32 i)  { return(_matches + i); };

private:
  uint32        _matchesLen;
  uint32        _matchesMax;
  atacMatch    *_matches;
};

#endif  //  ATAC_MATCHLIST_H
