// This file is part of A2Amapper.
// Copyright (c) 2004 Applied Biosystems
// Copyright (c) 2005 J. Craig Venter Institute
// Author: Dan Fasulo
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

#include "bio++.H"
#include "match.H"


void
match_s::dump(FILE *out, const char *descr, bool showSeq) {
  fprintf(out, "%s: ID:%s range1:"u32bitFMT","u32bitFMT" _pos="u32bitFMT" (seqlen="u32bitFMT")\n",
          descr, _matchId,
          _acc1->getRangeBegin(), _acc1->getRangeLength(), _acc1->_pos, _seq1->sequenceLength());
  fprintf(out, "%s  ID:%s range2:"u32bitFMT","u32bitFMT" _pos="u32bitFMT" (seqlen="u32bitFMT")  diag:"u32bitFMT" %s\n",
          descr, _matchId,
          _acc2->getRangeBegin(), _acc2->getRangeLength(), _acc2->_pos, _seq2->sequenceLength(),
          _diagonal, (_ori1 != _ori2) ? "reversed" : "");

  if (showSeq) {
    FastAAccessor  &A = *_acc1;
    FastAAccessor  &B = *_acc2;

    //  Save the position of the accessors
    //
    u32bit  acc1pos = A._pos;
    u32bit  acc2pos = B._pos;

    A.setPosition(A.getRangeBegin());
    B.setPosition(B.getRangeBegin());

    u32bit    margin = 5;
    u32bit    i = 0;
    char     *seq = new char [_acc1->getRangeEnd() - _acc1->getRangeBegin() + margin + margin + 32];
    char     *las = seq;

    strcpy(seq, ">>> ");
    while (*las) las++;

    for (i=0; i<margin; i++)
      --A;
    for (i=0; i<margin; i++, ++A)
      if (A.isValid())
        *las++ = *A;
      else
        *las++ = ' ';
    *las++ = ':';
    for (i=0; i<_acc1->getRangeEnd() - _acc1->getRangeBegin(); i++, ++A)
      if (A.isValid())
        *las++ = *A;
      else
        *las++ = ' ';
    *las++ = ':';
    for (i=0; i<margin; i++, ++A)
      if (A.isValid())
        *las++ = *A;
      else
        *las++ = ' ';

    *las++ = 0;
    fprintf(out, "%s\n", seq);

    las = seq;

    strcpy(seq, (_ori1 == _ori2) ? ">>> " : "<<< ");
    while (*las) las++;

    for (i=0; i<margin; i++)
      --B;
    for (i=0; i<margin; i++, ++B)
      if (B.isValid())
        *las++ = *B;
      else
        *las++ = ' ';
    *las++ = ':';
    for (i=0; i<_acc1->getRangeEnd() - _acc1->getRangeBegin(); i++, ++B)
      if (B.isValid())
        *las++ = *B;
      else
        *las++ = ' ';
    *las++ = ':';
    for (i=0; i<margin; i++, ++B)
      if (B.isValid())
        *las++ = *B;
      else
        *las++ = ' ';

    *las++ = 0;
    fprintf(out, "%s\n", seq);

    delete [] seq;

    //  Restore positions
    _acc1->_pos = acc1pos;
    _acc2->_pos = acc2pos;
  }
}









