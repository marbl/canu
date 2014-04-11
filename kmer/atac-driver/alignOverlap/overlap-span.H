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

#ifndef OVERLAP_SPAN_H
#define OVERLAP_SPAN_H

#include "overlap.H"

#define COLORSHIFT 24
#define COLORMASK  0x00ffffff

class span_t {
public:
  uint32  _iid;
  uint32  _beg;
  uint32  _end;
  uint32  _matchesLen;
  uint32  _matchesMax;
  uint32 *_matches;

  span_t(uint32 iid, uint32 beg, uint32 end) {
    _iid = iid;
    _beg = beg;
    _end = end;
    _matchesLen = 0;
    _matchesMax = 0;
    _matches    = 0L;
  };

  ~span_t() {
    delete [] _matches;
  };

  //  The top X bits of the _matches is for storing the color.  This
  //  does cut down the number of matches we can store.  Human-Human
  //  is ~1 million matches.

  void   addMatch(uint32 matchiid, uint32 color) {
    if (_matchesLen >= _matchesMax) {
      if (_matchesMax == 0)
        _matchesMax = 2;
      _matchesMax *= 2;
      uint32 *X = new uint32 [_matchesMax];
      memcpy(X, _matches, sizeof(uint32) * _matchesLen);
      delete [] _matches;
      _matches = X;
    }

    if (matchiid >> COLORSHIFT)
      fprintf(stderr, "ERROR!  span_t::addMatch()-- match id too big, decrease the color space.\n"), exit(1);

    _matches[_matchesLen++] = (color << COLORSHIFT) | (matchiid);
  };

  //  Split this span at position, return two new spans
  //
  void   split(uint32 position, span_t* &l, span_t* &r) {

    if ((position < _beg) || (_end < position)) {
      fprintf(stderr, "span_t::split()-- _beg="uint32FMT" _end="uint32FMT" postition="uint32FMT"?\n", _beg, _end, position);
      exit(1);
    }

    l = new span_t(_iid, _beg, position);
    r = new span_t(_iid, position, _end);

    l->_matchesLen = _matchesLen;
    l->_matchesMax = _matchesMax;
    l->_matches    = new uint32 [_matchesMax];
    memcpy(l->_matches, _matches, sizeof(uint32) * _matchesLen);

    r->_matchesLen = _matchesLen;
    r->_matchesMax = _matchesMax;
    r->_matches    = new uint32 [_matchesMax];
    memcpy(r->_matches, _matches, sizeof(uint32) * _matchesLen);
  };
};

#endif  //  OPERLAP_SPAN_H
