// This file is part of A2Amapper.
// Copyright (c) 2004 Applera Corporation
// Author: Dan Fasulo
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

#include <new>
#include <assert.h>
#include "MEMatch.H"

using namespace std;

MEMatch * MEMatch::__freeList = NULL;



MEMatch *
MEMatch::getNewMatch(const string& id,
                     FastASequenceInCore *s1,
                     FastASequenceInCore *s2, 
		     SeqOff_t start1,
                     SeqOff_t start2,
                     SeqOff_t len, 
		     bool s2_match_rev) {
  MEMatch *result;

  if (__freeList) {
    result = __freeList;
    __freeList = result->__next;
    result = new (result) MEMatch(id, s1, s2, start1, start2, len, 
				  s2_match_rev);
  }
  else {
    result = new MEMatch(id, s1, s2, start1, start2, len, s2_match_rev);
  }

  return result;
}


void 
MEMatch::freeMatch(MEMatch *m) {

  // Get rid of references to chromosomal sequences before putting
  // these in the free list!
  //
  m->_s1   = 0L;
  m->_s2   = 0L;

  m->__next = __freeList;
  __freeList = m;
}


MEMatch::MEMatch(const string& id, FastASequenceInCore *s1, FastASequenceInCore *s2, 
		 SeqOff_t start1, 
		 SeqOff_t start2, SeqOff_t len,
		 bool s2_match_rev) : _s1(s1), _s2(s2) {
  _id = id;
  _start1 = start1;
  _start2 = start2;
  _len = len;
  _rc = s2_match_rev;
  _del = false;
  __next = NULL;
}


SeqOff_t 
MEMatch::rcPos2(void) const {
  return _s2->sequenceLength() - _start2 - _len;
}


u64bit 
MEMatch::diagID(void) const {
  if (!_rc) {
    return _start2 + _s1->sequenceLength() - _start1 - 1;
  }
  else {			// reverse complement case
    return _start1 + _start2 + _len - 1;
  }
}


void 
MEMatch::extendForward(SeqOff_t num_chars) {
  _len += num_chars;
  if (_rc)
    _start2 -= num_chars;
}


void 
MEMatch::extendBackward(SeqOff_t num_chars) {
  _start1 -= num_chars;
  _len += num_chars;
  if (!_rc)
    _start2 -= num_chars;
}


void 
MEMatch::trimFront(SeqOff_t num_chars) {
  _start1 += num_chars;
  _len -= num_chars;
  if (!_rc)
    _start2 += num_chars;
}


void 
MEMatch::trimEnd(SeqOff_t num_chars) {
  _len -= num_chars;
  if (_rc)
    _start2 += num_chars;
}


bool 
MEMatch::canMergeWith(const MEMatch& m2) const {
  if (diagID() != m2.diagID())
    return false;

  return _start1 + _len >= m2._start1;
}


void 
MEMatch::consume(const MEMatch& m2) {
  SeqOff_t len;

  len = m2._start1 + m2._len - _start1;
  if (len > _len) {
    _len = len;
    if (_rc) 
      _start2 = m2._start2;
  }
}


void 
MEMatch::textDump(ostream& dest, bool show_seq, unsigned int margin) {
  u64bit p1 = pos1();
  u64bit p2 = pos2();
  SeqOff_t num_print = len() + 2 * margin;
  SeqOff_t i;

  FastAAccessor   A(s1(), false);
  FastAAccessor   B(s2(), isReversed());

  B.setReverseComplementRange(pos2(), len());

  dest << "ID : " << id() << "  Pos1: " << pos1() << "  Pos2: " << pos2()
       << "  Len: " << len() << '\n'
       << "Rev: ";
  if (isReversed()) 
    dest << "Y (" << rcPos2() << ")";
  else
    dest << "N";
  dest << "  Diag: " << diagID() << "\n";
      
  if (show_seq) {
    dest << ">>> ";

    //  Print the margin, or move the sequence back
    //
    if (p1 < margin) {
      for (i = 0; i < (margin - p1); ++i)
	dest << ' ';
      p1 = 0;
    } else {
      p1 = p1 - margin;
      i = 0;
    }
    A.setPosition(p1);

    for (; i < num_print; ++i, p1++, ++A) {
      if ((i == margin) || (i == (num_print - margin)))
	dest << '|';
      if (A.isValid() == false)
	dest << ' ';
      else
	dest << *A;
    }
    dest << '\n';



    //  Print the margin, or move the sequence back
    //
    if (isReversed()) {
      dest << "<<< ";

      if (p2 + len() + margin > _s2->sequenceLength()) {
        for (i = 0; i < (p2 + len() + margin - _s2->sequenceLength()); ++i)
          dest << ' ';
        //p2 = 0;
      } else {
        p2 = p2 - margin;
        i = 0;
      }
    } else {
      dest << ">>> ";

      if (p2 < margin) {
        for (i = 0; i < (margin - p2); ++i)
          dest << ' ';
        p2 = 0;
      } else {
        p2 = p2 - margin;
        i = 0;
      }
    }
    B.setPosition(p2);

    for (; i < num_print; ++i, p2++, ++B) {
      if ((i == margin) || (i == (num_print - margin)))
	dest << '|';
      if (B.isValid() == false)
	dest << ' ';
      else
	dest << *B;
    }
    dest << '\n';
  }
}


bool 
MEMatch::operator<(const MEMatch& m2) const {
  u64bit d = diagID(), d2 = m2.diagID();

  if (d < d2)
    return true;
  else if (d == d2)
    return (_start1 < m2._start1);
  else
    return false;
}


bool 
MEMatch::operator<=(const MEMatch& m2) const {
  u64bit d = diagID(), d2 = m2.diagID();

  if (d < d2)
    return true;
  else if (d == d2)
    return (_start1 <= m2._start1);
  else
    return false;
}
