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
  u32bit  _iid;
  u32bit  _beg;
  u32bit  _end;
  u32bit  _matchesLen;
  u32bit  _matchesMax;
  u32bit *_matches;

  span_t(u32bit iid, u32bit beg, u32bit end) {
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

  void   addMatch(u32bit matchiid, u32bit color) {
    if (_matchesLen >= _matchesMax) {
      if (_matchesMax == 0)
        _matchesMax = 2;
      _matchesMax *= 2;
      u32bit *X = new u32bit [_matchesMax];
      memcpy(X, _matches, sizeof(u32bit) * _matchesLen);
      delete [] _matches;
      _matches = X;
    }

    if (matchiid >> COLORSHIFT)
      fprintf(stderr, "ERROR!  span_t::addMatch()-- match id too big, decrease the color space.\n"), exit(1);

    _matches[_matchesLen++] = (color << COLORSHIFT) | (matchiid);
  };

  //  Split this span at position, return two new spans
  //
  void   split(u32bit position, span_t* &l, span_t* &r) {

    if ((position < _beg) || (_end < position)) {
      fprintf(stderr, "span_t::split()-- _beg="u32bitFMT" _end="u32bitFMT" postition="u32bitFMT"?\n", _beg, _end, position);
      exit(1);
    }

    l = new span_t(_iid, _beg, position);
    r = new span_t(_iid, position, _end);

    l->_matchesLen = _matchesLen;
    l->_matchesMax = _matchesMax;
    l->_matches    = new u32bit [_matchesMax];
    memcpy(l->_matches, _matches, sizeof(u32bit) * _matchesLen);

    r->_matchesLen = _matchesLen;
    r->_matchesMax = _matchesMax;
    r->_matches    = new u32bit [_matchesMax];
    memcpy(r->_matches, _matches, sizeof(u32bit) * _matchesLen);
  };
};

#endif  //  OPERLAP_SPAN_H
