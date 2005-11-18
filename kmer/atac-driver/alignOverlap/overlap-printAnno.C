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

void
printAnno(FILE *F, annoList *AL, u32bit &ALlen,
          char label,
          u32bit axis,
          span_t *span,
          u32bit match1, match_t *m1,
          u32bit match2, match_t *m2) {

  //  If we're just given match1, make it match2 if it is the second mapping
  //
  if ((match1 >> COLORSHIFT) && (match2 == u32bitZERO)) {
    match2 = match1; m2 = m1;
    match1 = 0;      m1 = 0;
  }

  u32bit len  = span->_end - span->_beg;

  //  axis is 1 or 2; if we're the first axis (B35 centric) make a
  //  list of the matches for later processing

  if (axis == 1)
    AL[ALlen++].add(label, span->_iid, span->_beg, len,
                    match1 & COLORMASK, m1,
                    match2 & COLORMASK, m2);

  fprintf(F, "%c "u32bitFMTW(4)":"u32bitFMTW(09)"-"u32bitFMTW(09)"["u32bitFMTW(6)"] ",
          label,
          span->_iid, span->_beg, span->_end, len);

  if (m1) {
    fprintf(F, u32bitFMTW(07)" ", m1->matchuid);
    u32bit off1 = span->_beg - m1->pos1;

    if (axis == 1) {
      u32bit  sta = m1->pos2 + off1;
      u32bit  end = m1->pos2 + off1 + len;

      if (m1->fwd2 == 0) {
        sta = m1->pos2 + m1->len2 - off1;
        end = m1->pos2 + m1->len2 - off1 - len;
      }

      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m1->iid2, sta, end);
    } else {
      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m1->iid1, m1->pos1 + off1, m1->pos1 + off1 + len);
    }
  } else {
    fprintf(F, u32bitFMTW(07)" ", u32bitZERO);
    fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", u32bitZERO, u32bitZERO, u32bitZERO);
  }

  if (m2) {
    fprintf(F, u32bitFMTW(07)" ", m2->matchuid);
    u32bit off2 = span->_beg - m2->pos1;

    if (axis == 1) {
      u32bit  sta = m2->pos2 + off2;
      u32bit  end = m2->pos2 + off2 + len;
      
      if (m2->fwd2 == 0) {
        sta = m2->pos2 + m2->len2 - off2;
        end = m2->pos2 + m2->len2 - off2 - len;
      }

      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m2->iid2, sta, end);
    } else {
      fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", m2->iid1, m2->pos1 + off2, m2->pos1 + off2 + len);
    }
  } else {
    fprintf(F, u32bitFMTW(07)" ", u32bitZERO);
    fprintf(F, "("u32bitFMTW(8)": "u32bitFMTW(9)"-"u32bitFMTW(9)") ", u32bitZERO, u32bitZERO, u32bitZERO);
  }

  fprintf(F, "\n");
}
