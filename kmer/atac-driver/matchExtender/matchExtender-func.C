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
#include "atac.H"
#include "match.H"

using namespace std;


extern u32bit  minEndRunLen;
extern u32bit  maxMMBlock;
extern u32bit  minBlockSep;
extern double  minIdentity;
extern u32bit  maxNbrSep;
extern u32bit  maxNbrPathMM;


//#define DEBUG_TRACE
//#define DEBUG_TRIMTOPERCENT
//#define DEBUG_EXTEND
//#define DEBUG_EXTEND_CONSUME
//#define DEBUG_EXTEND_BACK
//#define DEBUG_EXTEND_FORWARD



//  Return true if c1 and c2 are identities, false otherwise.
//
bool
isIdentity(char c1, char c2) {
  return((letterToBits[(int)c1] != 0xff) &&
         (letterToBits[(int)c2] != 0xff) &&
         IUPACidentity[(int)c1][(int)c2]);
}



//  Finds the largest block >= 'pct' (95%) identity.
//
bool
trim_to_pct(vector<match_s *>& matches, u32bit midx, double pct) {
#ifdef DEBUG_TRACE
  fprintf(stderr, "trim_to_pct()\n");
#endif

  u32bit     best_start   = 0;
  u32bit     best_len     = 0;
  match_s   *m = matches[midx];

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

  A.setPosition(m->pos1());
  B.setPosition(m->pos2());

#ifdef DEBUG_TRIMTOPERCENT
  //m->dump(stderr, "TrimToPercent", false);
#endif

  //  For all starting positions:
  //
  //  We could short-circuit here - once (m->len() - start) becomes
  //  shorter than our best_len, we have no hope in finding a better
  //  one.
  //
  for (u32bit start=0;
       (start< m->len()) && (m->len() - start > best_len);
       ++start) {
    u32bit best_run_len = 0;
    u32bit sum          = 0;

    A.setPosition(m->pos1() + start);
    B.setPosition(m->pos2() + start);

    //  And all ending positions:
    //
    //  Compute the number of identities we've seen, and remember the
    //  length of the highest identity.
    //
    for (u32bit len = 1; start + len <= m->len(); ++len) {
      char c1 = *A;
      char c2 = *B;

      //  We just extend the last result by one, rather than recompute
      //  the whole value for our new range (start, len).
      //
      if (isIdentity(c1, c2))
	sum++;

      //  If the sum is more than 'pct' identities, we are by
      //  construction the longest run at this starting point, so
      //  remember it.
      //
      if (sum >= pct * len)
	best_run_len = len;

      ++A;
      ++B;
    }

    //  Special case: if the whole string is okay, don't check any
    //  subranges
    //
    if ((start == 0) && (best_run_len == m->len()))
      return(false);

    //  If we've just found a longer subrange, remember it.
    //
    if (best_run_len > best_len) {
      best_start = start;
      best_len   = best_run_len;
    }
  }

  if (best_len < m->len()) {
#ifdef DEBUG_TRIMTOPERCENT
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "Trimming to substring with start="u32bitFMT" and len="u32bitFMT" for percent identity\n",
            best_start, best_len);
    m->dump(stderr, "BEFORE", true);
#endif

    m->extendLeft(-(s32bit)best_start);
    m->extendRight(-(s32bit)(m->len() - best_len));

#ifdef DEBUG_TRIMTOPERCENT
    m->dump(stderr, "AFTER", true);
    fprintf(stderr, "============================================================\n");
#endif

    return(true);
  }

  return(false);
}


void
extend_match_backward(vector<match_s *>& matches,
                      u32bit midx,
		      u32bit min_start_pos) {
#ifdef DEBUG_TRACE
  fprintf(stderr, "extend_match_backward()-- min_start_pos="u32bitFMT"\n", min_start_pos);
#endif

  // Assumes when traveling backwards that we will never run into
  // another match (otherwise, that match would have been forward
  // extended previously).

  u32bit     num_recent_mismatches = 0;
  match_s   *m = matches[midx];
  u32bit     good_run_len = (int) m->len();
  u32bit     num_pending = 0;

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

  A.setPosition(m->_acc1->getRangeBegin());
  B.setPosition(m->_acc2->getRangeBegin());

  //  Decrement, instead of subtract one from the position above, to
  //  avoid any issues with overflow (e.g., 0 - 1).
  //
  --A;
  --B;

  while ((A.getPosition() > min_start_pos) && A.isValid() && B.isValid()) {
    char c1 = *A;
    char c2 = *B;

    if (isIdentity(c1, c2)) {
      good_run_len++; 

      //  If we've gone long enough, erase our mismatch record
      //
      if (good_run_len == minBlockSep)  //  20 by default
	num_recent_mismatches = 0;

      //  If we're in the middle of a long good run, add the character
      //  to the match (END_RUN_LEN=10)
      //
      //  Otherwise, if we just made the minimum extension length, add
      //  all of the pending characters.
      //
      //  Otherwise, this character is pending.  However, still do
      //  output if we're run out of sequence.
      //
      if (good_run_len > minEndRunLen) {
	m->extendLeft(1);
      } else if (good_run_len == minEndRunLen) {
	m->extendLeft(num_pending + 1);
	num_pending = 0;
      } else {
	num_pending++;
      }
    } else {
      good_run_len = 0;
      num_pending++;
      num_recent_mismatches++;
      if (num_recent_mismatches > maxMMBlock)  //  3 by default
	break;
    }

    --A;
    --B;
  }

  //  If we hit the end of the sequence, and are good, do extension
  //
  //  if ((A.getPosition() == min_start_pos) || (B.getPosition() == 0))
  //
  if (!A.isValid() || !B.isValid() || (A.getPosition() <= min_start_pos))
    m->extendLeft(num_pending);


#ifdef DEBUG_EXTEND_BACK
  fprintf(stderr, "extend_back()-- M u %s . %s %d %d 1 %s %d %d 1\n",
          m->_matchId,
          m->_id1, m->_acc1->getRangeBegin(), m->_acc1->getRangeLength(),
          m->_id2, m->_acc2->getRangeBegin(), m->_acc2->getRangeLength());
#endif
}


bool
can_reach_nearby_match(match_s *src, match_s *dest) {
#ifdef DEBUG_TRACE
  fprintf(stderr, "can_reach_nearby_match()\n");
#endif

  if (dest->pos1() - (src->pos1() + src->len()) > (u32bit) maxNbrSep)  // 100
    return false;

#if 0
  src->dump(stderr, "src:");
  dest->dump(stderr, "dst:");
#endif

  FastAAccessor  &A = *src->_acc1;
  FastAAccessor  &B = *src->_acc2;

  A.setPosition(A.getRangeEnd() - 1);
  B.setPosition(B.getRangeEnd() - 1);

  ++A;
  ++B;

  u32bit  num_mismatch = 0;

  while ((num_mismatch    <= maxNbrPathMM) &&  // 5
         (A.getPosition() <  dest->pos1()) &&
         (A.isValid()) &&
         (B.isValid())) {
    if (!isIdentity(*A, *B))
      num_mismatch++;

    ++A;
    ++B;
  }

#if 0
  fprintf(stderr, "num_mismatch=%d  pos: %d %d   valid: A:%d B:%d\n",
          num_mismatch, A.getPosition(), dest->pos1(), A.isValid(), B.isValid());
#endif

  return(num_mismatch <= maxNbrPathMM);  // 5
}



//  Stops and returns true if we hit the next match
//
bool
extend_match_forward(vector<match_s *>& matches, u32bit midx, match_s *target) {
#ifdef DEBUG_TRACE
  fprintf(stderr, "extend_match_forward()\n");
#endif

  match_s     *m = matches[midx];
  u32bit       num_recent_mismatches = 0;
  u32bit       num_pending = 0;

  u32bit       good_run_len = (int) m->len();

  FastAAccessor  &A = *m->_acc1;
  FastAAccessor  &B = *m->_acc2;

#ifdef DEBUG_EXTEND_FORWARD
  fprintf(stderr, "extend_match_forward()-- A:%4d-%4d B:%4d-%4d\n",
          A.getRangeBegin(), A.getRangeLength(),
          B.getRangeBegin(), B.getRangeLength());
#endif

  //  Set our position to the last valid base in the range, then move
  //  to the next one.
  //
  A.setPosition(A.getRangeEnd() - 1);
  B.setPosition(B.getRangeEnd() - 1);

  ++A;
  ++B;

  while (A.isValid() && B.isValid()) {
    char c1 = *A;
    char c2 = *B;

    if (isIdentity(c1, c2)) {
      good_run_len++;

      //fprintf(stderr, "extend-forward %c %c\n", c1, c2);

      // Pass Go and collect $200
      //
      if (good_run_len == minBlockSep)
	num_recent_mismatches = 0;

      //  If not enough good characters yet, increase the length
      //  pending.  We used to check for the hitting the end of the
      //  sequence here.
      //
      //  Otherwise, if we have just made the minumum good run length,
      //  do the extension.
      //
      //  Otherwise, if we're above the minimum good length, extend by
      //  another character.
      //
      if (good_run_len < minEndRunLen) {
	num_pending++;
      } else if (good_run_len == minEndRunLen) {
	m->extendRight(num_pending + 1);
	num_pending = 0;
      } else if (good_run_len > minEndRunLen) {
	m->extendRight(1);
      }

      //  If we've run into (and possibly over) another seed match,
      //  return so the main loop can consume and restart.
      //
      if (m->canMergeWith(target))
	return(true);
    } else {
      good_run_len = 0;
      num_pending++;
      num_recent_mismatches++;

      if (num_recent_mismatches > maxMMBlock)
	return(false);
    }

    ++A;
    ++B;
  }

  //  If we've got a short good run but have hit the end of
  //  a sequence, do extension.
  //
  if ((!A.isValid() || !B.isValid()) && (good_run_len < minEndRunLen))
    m->extendRight(num_pending);

#ifdef DEBUG_EXTEND_FORWARD
  fprintf(stderr, "extend_match_forward(finish)-- A:%4d-%4d B:%4d-%4d\n",
          A.getRangeBegin(), A.getRangeLength(),
          B.getRangeBegin(), B.getRangeLength());
#endif
  
  return(false);
}







u32bit
extend_matches_on_diagonal(vector<match_s *>& matches, u32bit diag_start) {
#ifdef DEBUG_TRACE
  fprintf(stderr, "extend_matches_on_diagonal()\n");
#endif

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

#ifdef DEBUG_EXTEND_BACK
    m->dump(stderr, "Before back extension:", true);
#endif

    extend_match_backward(matches, idx, prev_end);

#ifdef DEBUG_EXTEND_BACK
    m->dump(stderr, "After back extension:", true);
#endif

#ifdef DEBUG_EXTEND
    fprintf(stderr, "1M u %s . %s %d %d 1 %s %d %d 1\n",
            matches[idx]->_matchId,
            matches[idx]->_id1, matches[idx]->_acc1->getRangeBegin(), matches[idx]->_acc1->getRangeLength(),
            matches[idx]->_id2, matches[idx]->_acc2->getRangeBegin(), matches[idx]->_acc2->getRangeLength());
#endif

    prev_end = m->pos1() + m->len();

    if ((m->pos1() > m->seq1()->sequenceLength()) || (m->pos2() > m->seq2()->sequenceLength()))
      m->dump(stderr, "NEGATIVE after back extend!\n", true), abort();
  }


  //  Now forward extend each match


  idx = diag_start; 
  while ((idx < matches.size()) && 
	 (matches[idx]->_diagonal == diag_id)) {

    if (matches[idx]->isDeleted()) {
      idx++;
      continue;
    }
    
#ifdef DEBUG_EXTEND
    fprintf(stderr, "2M u %s . %s %d %d 1 %s %d %d 1\n",
            matches[idx]->_matchId,
            matches[idx]->_id1, matches[idx]->_acc1->getRangeBegin(), matches[idx]->_acc1->getRangeLength(),
            matches[idx]->_id2, matches[idx]->_acc2->getRangeBegin(), matches[idx]->_acc2->getRangeLength());
#endif

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
#ifdef DEBUG_EXTEND_CONSUME
      m->dump(stderr, "I can_reach_nearby_match and extend this", true);
      next_m->dump(stderr, "with this", true);
#endif
      m->consume(next_m);
      next_m->setDeleted();
#ifdef DEBUG_EXTEND_CONSUME
      m->dump(stderr, "Extended through next match via neighbor search:", true);
#endif
      continue;
    }

    //  Otherwise, try to make it to the next match with the
    //  character-at- a-time extension rules.  If we make it, restart
    //  the loop with the same match (now extended).  Otherwise, trim
    //  the extended match as necessary and move on to the next
    //  match.
    //
    if (extend_match_forward(matches, idx, next_m)) {
#ifdef DEBUG_EXTEND_CONSUME
        m->dump(stderr, "I extend_match_forward and extend this", true);
        next_m->dump(stderr, "with this", true);
#endif
      m->consume(next_m);
      next_m->setDeleted();
#ifdef DEBUG_EXTEND_CONSUME
      m->dump(stderr, "Extended through next match via forward extension:", true);
#endif
      continue;
    }

#ifdef DEBUG_EXTEND
    //m->dump(stderr, "Failed to make next match.  Final extended version:", true);
#endif

#ifdef DEBUG_EXTEND
    fprintf(stderr, "3M u %s . %s %d %d 1 %s %d %d 1\n",
            matches[idx]->_matchId,
            matches[idx]->_id1, matches[idx]->_acc1->getRangeBegin(), matches[idx]->_acc1->getRangeLength(),
            matches[idx]->_id2, matches[idx]->_acc2->getRangeBegin(), matches[idx]->_acc2->getRangeLength());
#endif

    //  Didn't make it, so trim and move on
    //
    if (trim_to_pct(matches, idx, minIdentity)) {
#ifdef DEBUG_EXTEND_TRIMMING
      m->dump(stderr, "After trimming:", true);
#endif
    } else {
#ifdef DEBUG_EXTEND_TRIMMING
      fprintf(stderr, "No trimming done.\n");
#endif
    }

#ifdef DEBUG_EXTEND
    fprintf(stderr, "4M u %s . %s %d %d 1 %s %d %d 1\n",
            matches[idx]->_matchId,
            matches[idx]->_id1, matches[idx]->_acc1->getRangeBegin(), matches[idx]->_acc1->getRangeLength(),
            matches[idx]->_id2, matches[idx]->_acc2->getRangeBegin(), matches[idx]->_acc2->getRangeLength());
#endif

#ifdef DEBUG_EXTEND
    if ((m->pos1() > m->seq1()->sequenceLength()) || (m->pos2() > m->seq2()->sequenceLength()))
      m->dump(stderr, "NEGATIVE after forward extend!", true), abort();

    fprintf(stderr, "\n==============\n\n");
#endif

    ++idx;
  }
  
  return idx;
}



