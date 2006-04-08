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

#include "overlap.H"

//  I really wanted this to be parameterized with two macros, but the
//  preprocessor merges, then replaces:
//    #define INDEXA 1
//    #define INDEXB 2
//    #define NODE node ## INDEXA
//  results in 'nodeINDEXA' not 'node1'

void
NAME(FILE             *outfile,
     spanTree         *S,
     atacMatchList    *M1,
     atacMatchList    *M2,
     overlapStats     &stats,
     annoList         *AL,
     u32bit           &ALlen,
     u32bit           &ALmax) {

  dnode_t  *node = dict_first(S->_tree);

  while (node) {
    span_t  *span = (span_t *)dnode_getkey(node);
    u32bit   spanLen = span->_end - span->_beg;

    if (span->_matchesLen == 0) {
      stats.unmapped += spanLen;
      printAnno(outfile, AL, ALlen, 'U', INDEX, span);
    } else if (span->_matchesLen == 1) {
      u32bit    match = span->_matches[0];
      atacMatch  *m;

      if (match >> COLORSHIFT) {
        m = &M2->_matches[match & COLORMASK];
        stats.map2unique += spanLen;
      } else {
        m = &M1->_matches[match & COLORMASK];
        stats.map1unique += spanLen;
      }

      printAnno(outfile, AL, ALlen, '1', INDEX, span, match, m);
    } else if ((span->_matchesLen == 2) &&
               ((span->_matches[0] >> COLORSHIFT) == (span->_matches[1] >> COLORSHIFT))) {
      stats.inconsistent += spanLen;
      printAnno(outfile, AL, ALlen, '?', INDEX, span);
    } else if (span->_matchesLen == 2) {
      u32bit match1 = span->_matches[0];
      u32bit match2 = span->_matches[1];

      if (match1 >> COLORSHIFT) {
        match1 = span->_matches[1];
        match2 = span->_matches[0];
      }

      atacMatch  *m1 = &M1->_matches[match1 & COLORMASK];
      atacMatch  *m2 = &M2->_matches[match2 & COLORMASK];

      if (m1->iid2 == m2->iid2) {
        u32bit off1  = span->_beg - m1->POS1;
        u32bit pos1l = m1->POS2 + off1;
        u32bit pos1r = m1->POS2 + m1->LEN2 - off1;

        u32bit off2  = span->_beg - m2->POS1;
        u32bit pos2l = m2->POS2 + off2;
        u32bit pos2r = m2->POS2 + m2->LEN2 - off2;

        if ((pos1l == pos2l) || (pos1r == pos2r)) {
          stats.same += spanLen;
          printAnno(outfile, AL, ALlen, 'Y', INDEX, span, match1, m1, match2, m2);
        } else {
          stats.different += spanLen;
          printAnno(outfile, AL, ALlen, 'N', INDEX, span, match1, m1, match2, m2);
        }
      } else {
        //  Wildly different matches!  Mapped to different scaffolds!
        stats.wilddiff += spanLen;
        printAnno(outfile, AL, ALlen, '!', INDEX, span, match1, m1, match2, m2);
      }
    } else {
      stats.inconsistent += spanLen;
      printAnno(outfile, AL, ALlen, '?', INDEX, span);
    }

    node = dict_next(S->_tree, node);
  }
}
  
