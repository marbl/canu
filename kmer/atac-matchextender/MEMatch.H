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

#ifndef _ME_MATCH_H_
#define _ME_MATCH_H_

#include <string>
#include <iostream>
#include "../libbio/bio++.H"

using namespace std;

class MEMatch {

public:

  //
  // CLASS CONSTANTS AND METHODS
  //

  /* Class does its own memory management via a free list for efficiency */
  static MEMatch *getNewMatch(const string   &id,
                              FastASequenceInCore   *s1,
                              FastASequenceInCore   *s2, 
			      SeqOff_t        start1,
                              SeqOff_t        start2,
                              SeqOff_t        len, 
			      bool            s2_match_rev);
  static void freeMatch(MEMatch *m);
  
  //
  // CONSTRUCTORS
  //
  
  /* Note that s2 is always the sequence of the forward strand, even
     if the match is a forward-reverse match. */
  MEMatch(const string& id,
          FastASequenceInCore *s1,
          FastASequenceInCore *s2,
          SeqOff_t start1, 
	  SeqOff_t start2,
          SeqOff_t len,
          bool s2_match_rev);
  
  //
  // METHODS
  //

  /* Basic accessors */

  inline string id(void) const { return _id; }

  /* Note that s2 is always the sequence of the forward strand, even
     if the match is a forward-reverse match. */
  inline FastASequenceInCore *s1(void)   const { return _s1; }
  inline FastASequenceInCore *s2(void)   const { return _s2; }

  inline SeqOff_t pos1(void) const { return _start1; }
  inline SeqOff_t pos2(void) const { return _start2; }

  inline SeqOff_t len(void) const { return _len; }

  inline bool isDeleted(void) const { return _del; }
  inline void setDeleted(bool d = true) { _del = d; }

  /* Other operations */

  /* Sorts primarily on diagonal, then on pos1 */
  bool operator<(const MEMatch& m2) const;

  /* Sorts primarily on diagonal, then on pos1 */
  bool operator<=(const MEMatch& m2) const;

  /* Returns true if and only if the hit is to the reverse complement
     of Seq2. */
  inline bool isReversed(void) const { return _rc; }

  /* Returns the index of the first character of the match in the reverse
     complement of s2(); only useful if isReversed() is true. */
  SeqOff_t rcPos2(void) const;

  /* Calculates the diagonal ID of the match.  Its absolute values is
     in the range 0 to s1.size() + s2.size().  A negative value indicates
     that this is a forward-to-reverse match; otherwise it is a forward-forward
     match. */
  SeqOff_t diagID(void) const;

  /* Extending and merging ops */

  /* Simply updates the match lengths and starts.  Note that "forward" is
     relative to s1. */
  void extendForward(SeqOff_t num_chars);

  /* Simply updates the match lengths and starts.  Note that "backward" is
     relative to s1. */
  void extendBackward(SeqOff_t num_chars);

  /* Removes k positions from the beginning of the match.  Assumes that
     length will not be reduced below 1. */
  void trimFront(SeqOff_t num_chars);

  /* Removes k positions from the end of the match.  Assumes that
     length will not be reduced below 1. */
  void trimEnd(SeqOff_t num_chars);

  /* Returns true if and only if this match hits m2.  The matches must
     be on the same diagonal between the same two sequences. */
  bool canMergeWith(const MEMatch& m2) const;

  /* Modifies this match to consume the other match.  It must be
     the case that this < m2 and this.canMergeWith(m2). */
  void consume(const MEMatch& m2);

  /* Mainly for debugging */

  /* Print a text/debugging version to dest.  If show_seq is true,
     show "alignment" of the two matches.  If margin is non-zero,
     include up to margin characters on either side of the match
     (if present). */
  void textDump(ostream& dest, bool show_seq, unsigned int margin = 0);

protected:

  string                _id;
  FastASequenceInCore  *_s1;
  FastASequenceInCore  *_s2;
  SeqOff_t              _start1;
  SeqOff_t              _start2;
  SeqOff_t              _len;
  bool                  _rc;
  bool                  _del;

private:

  // Next pointer for pooled memory allocation
  MEMatch *__next;

  // Pool of free structures for allocation
  static MEMatch * __freeList;
};

#endif

