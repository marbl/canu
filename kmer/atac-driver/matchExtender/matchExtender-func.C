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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <vector>
#include <algorithm>

#include "bio++.H"
#include "atac-common.H"
#include "match.H"

using namespace std;


u32bit   DEF_MIN_END_RUN_LEN   = 10;
u32bit   DEF_MAX_MM_BLOCK      = 3;
u32bit   DEF_MIN_BLOCK_SEP     = 20;
double   DEF_MIN_IDENTITY      = 0.95;
u32bit   DEF_MAX_NBR_SEP       = 100;
u32bit   DEF_MAX_NBR_PATH_MM   = 5;




bool
trim_to_pct(vector<match_s *>& matches, u32bit midx, double pct) {
  fprintf(stderr, "trim_to_pct()\n");

  u32bit start        = 0;
  u32bit len          = 0;
  u32bit best_start   = 0;
  u32bit best_len     = 0;

  match_s *m = matches[midx];

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

  A.setPosition(m->pos1());
  B.setPosition(m->pos2());

  for (start = 0; start < m->len(); ++start) {
    u32bit best_run_len = 0;
    u32bit sum          = 0;

    A.setPosition(m->pos1() + start);
    B.setPosition(m->pos2() + start);

    for (len = 1; len <= (m->len() - start); ++len) {
      char c1 = *A;
      char c2 = *B;

      if (validSymbol[c1] &&
          validSymbol[c2] &&
          IUPACidentity[c1][c2])
	sum++;

      if (sum >= pct * len)
	best_run_len = len;

      ++A;
      ++B;
    }

    // Special case : if the whole string is okay, don't check any
    // subranges
    //
    if ((start == 0) && (best_run_len == m->len()))
      return(false);

    if (best_run_len > best_len) {
      best_len = best_run_len;
      best_start = start;   
    }
  }

  if (best_len < m->len()) {
    fprintf(stderr, "Trimming to substring with start="u32bitFMT" and len="u32bitFMT" for percent identity\n",
            best_start, best_len);

    m->extendLeft(-(s32bit)best_start);
    m->extendRight(-(s32bit)(m->len() - best_len));

    return(true);
  }

  return(false);
}


void
extend_match_backward(vector<match_s *>& matches,
                      u32bit midx,
		      u32bit min_start_pos) {
  fprintf(stderr, "extend_match_backward()-- min_start_pos="u32bitFMT"\n", min_start_pos);

  // Assumes when traveling backwards that we will never run into
  // another match (otherwise, that match would have been forward
  // extended previously).

  u32bit     num_recent_mismatches = 0;
  match_s   *m = matches[midx];
  u32bit     good_run_len = (int) m->len();
  u32bit     num_pending = 0;

  m->dump(stderr, "extend_back()-1- ");

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

  A.setPosition(m->_acc1->getRangeBegin());
  B.setPosition(m->_acc2->getRangeBegin());

  //  Decrement, instead of subtract one from the position above, to
  //  avoid any issues with overflow (0 - 1).
  //
  fprintf(stderr, "extend_back(outside)-- %9d/%c  --  %9d/%c\n", A.getPosition(), *A, B.getPosition(), *B);
  --A;
  --B;
  fprintf(stderr, "extend_back(outside)-- %9d/%c  --  %9d/%c\n", A.getPosition(), *A, B.getPosition(), *B);

  while ((A.getPosition() > min_start_pos) && A.isValid() && B.isValid()) {
    char c1 = *A;
    char c2 = *B;

    fprintf(stderr, "extend_back(before)-- %9d/%c  --  %9d/%c\n", A.getPosition(), c1, B.getPosition(), c2);
    m->dump(stderr, "extend_back(before)-- ");

    if (validSymbol[c1] &&
        validSymbol[c2] &&
        IUPACidentity[c1][c2]) {
      good_run_len++; 

      //  If we've gone long enough, erase our mismatch record (BLOCK_SEP=20)
      if (good_run_len == DEF_MIN_BLOCK_SEP)
	num_recent_mismatches = 0;

      //  If we're in the middle of a long good run, add the character to
      //  the match (END_RUN_LEN=10)
      //
      if (good_run_len > DEF_MIN_END_RUN_LEN) {
	m->extendLeft(1);
      } else if (good_run_len == DEF_MIN_END_RUN_LEN) {
        //  else if we just made the minimum extension length, add
        //  all of the pending characters.
        //
	m->extendLeft(num_pending + 1);
	num_pending = 0;
      } else {
        //  Otherwise, this character is pending.  However, still do output
        //  if we're run out of sequence.
        //
	num_pending++;
	if ((A.getPosition() == min_start_pos) || (B.getPosition() == 0))
	  m->extendLeft(num_pending);
      }
    } else {
      //  Else this character is a mismatch.
      //
      num_pending++;
      good_run_len = 0;
      num_recent_mismatches++;
      //  MM_BLOCK=3
      if (num_recent_mismatches > DEF_MAX_MM_BLOCK)
	break;
    }

    m->dump(stderr, "extend_back(after)-- ");

    --A;
    --B;
  }

  fprintf(stderr, "extend_match_backward()-- pos1="u32bitFMT",valid=%d pos2="u32bitFMT",valid=%d\n",
          A.getPosition(), A.isValid(), B.getPosition(), B.isValid());
  m->dump(stderr, "extend_match_backward()--");
}


bool
can_reach_nearby_match(match_s *src, match_s *dest) {
  fprintf(stderr, "can_reach_nearby_match()\n");

  u32bit   num_mismatch = 0;

  if (dest->pos1() - (src->pos1() + src->len()) > (u32bit) DEF_MAX_NBR_SEP)
    return false;

  FastAAccessor  &A = *src->_acc1;
  FastAAccessor  &B = *src->_acc2;

  A.setPosition(src->_acc1->getRangeBegin());
  B.setPosition(src->_acc2->getRangeBegin());

  while ((num_mismatch    <= DEF_MAX_NBR_PATH_MM) &&
         (A.getPosition() <  dest->pos1()) &&
         (A.isValid()) &&
         (B.isValid())) {
    char c1 = *A;
    char c2 = *B;

    if (validSymbol[c1] &&
        validSymbol[c2] &&
        !IUPACidentity[c1][c2])
      num_mismatch++;

    ++A;
    ++B;
  }

  return (num_mismatch <= DEF_MAX_NBR_PATH_MM);
}



  //  Stops and returns true if we hit the next match
bool
extend_match_forward(vector<match_s *>& matches, u32bit midx, match_s *target) {
  fprintf(stderr, "extend_match_forward()\n");

  match_s     *m = matches[midx];
  u32bit       num_recent_mismatches = 0;
  u32bit       num_pending = 0;

  u32bit       good_run_len = (int) m->len();

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

  A.setPosition(m->_acc1->getRangeBegin());
  B.setPosition(m->_acc2->getRangeBegin());

  while (A.isValid() && B.isValid()) {
    char c1 = *A;
    char c2 = *B;

    // If there's a match ...
    //
    if (validSymbol[c1] &&
        validSymbol[c2] &&
        IUPACidentity[c1][c2]) {
      good_run_len++;

      // Pass Go and collect $200
      //
      if (good_run_len == DEF_MIN_BLOCK_SEP)
	num_recent_mismatches = 0;

      //  Not enough good characters yet, so output is pending ...
      //
      if (good_run_len < DEF_MIN_END_RUN_LEN) {
	num_pending++;	

	//  If we've got a short good run but have hit the end of
	//  a sequence, do extension.
        //

	//if (!A.isValid(1) || !B.isValid(1))

        if (((A.getPosition() + 1) == A.getLength()) || ((B.getPosition() + 1) == B.getLength()))
	  m->extendRight(num_pending);

      } else if (good_run_len == DEF_MIN_END_RUN_LEN) {
        //  Else if we've just made the minimum good end run length,
        //  commit to extending.
        //
	m->extendRight(num_pending + 1);
	num_pending = 0;
      } else if (good_run_len > DEF_MIN_END_RUN_LEN) {
        //  Otherwise, if we're above the minimum good run length,
	//  extend by one character.
	m->extendRight(1);
      }

      //  If we've run into (and possibly over) another seed match,
      //  consume it and return so the main loop can restart
      //
      if (target && (m->canMergeWith(target))) {
	m->consume(target);
	return true;
      }
    } else {
      // Otherwise, process mismatch ...
      num_pending++;
      good_run_len = 0;
      num_recent_mismatches++;
      if (num_recent_mismatches > DEF_MAX_MM_BLOCK)
	break;
    }

    ++A;
    ++B;
  }
  
  return false;
}







u32bit
extend_matches_on_diagonal(vector<match_s *>& matches, u32bit diag_start) {
  fprintf(stderr, "extend_matches_on_diagonal()\n");

  u32bit     diag_id = matches[diag_start]->_diagonal;
  u32bit     idx;
  u32bit     prev_end = 0;
  match_s   *m;
  match_s   *next_m = NULL;

  //  Back extend each match as far as possible (but never over the
  //  preceding match
  //
  for (idx = diag_start; 
       (idx < matches.size()) && (matches[idx]->_diagonal == diag_id);
       ++idx) {

    m = matches[idx];

    m->dump(stderr, "Before back extension:");
    extend_match_backward(matches, idx, prev_end);
    m->dump(stderr, "After back extension:");

    prev_end = m->pos1() + m->len();

    if ((m->pos1() > m->s1()->sequenceLength()) || (m->pos2() > m->s2()->sequenceLength()))
      m->dump(stderr, "NEGATIVE after back extend!\n"), abort();
  }


  //  Now forward extend each match


  idx = diag_start; 
  while ((idx < matches.size()) && 
	 (matches[idx]->_diagonal == diag_id)) {

    if (matches[idx]->isDeleted()) {
      idx++;
      continue;
    }
    
    m        = matches[idx];
    next_m   = 0L;

    for (u32bit next_idx=idx+1; ((next_idx < matches.size()) && 
                                 (matches[next_idx]->_diagonal == diag_id) &&
                                 (next_m == 0L)); next_idx++)
      if (matches[next_idx]->isDeleted() == false)
	next_m = matches[next_idx];

    //  First, try to reach the next match with the simple "maximum of
    //  k mismatches" rule. If we made it, consume the next match and
    //  start the loop again with the same match (now extended)
    //
    if (next_m && can_reach_nearby_match(m, next_m)) {
      m->consume(next_m);
      next_m->setDeleted();
      m->dump(stderr, "Extended through next match via neighbor search:");
      continue;
    }

    //  Otherwise, try to make it to the next match with the
    //  character-at- a-time extension rules.  If we make it, restart
    //  the loop with the same match (now extended).  Otherwise, trim
    //  the extended match as necessary and move on to the next
    //  match.
    //
    if (extend_match_forward(matches, idx, next_m)) {
      next_m->setDeleted();
      m->dump(stderr, "Extended through next match via forward extension:");
      continue;
    } else {
      m->dump(stderr, "Failed to make next match.  Final extended version:");
    }
    
    //  Didn't make it, so trim and move on
    //
    if (trim_to_pct(matches, idx, DEF_MIN_IDENTITY)) {
      m->dump(stderr, "After trimming:");
    } else {
      fprintf(stderr, "No trimming done.\n");
    }

    if ((m->pos1() > m->s1()->sequenceLength()) || (m->pos2() > m->s2()->sequenceLength()))
      m->dump(stderr, "NEGATIVE after forward extend!"), abort();

    fprintf(stderr, "\n==============\n\n");

    ++idx;
  }
  
  return idx;
}



