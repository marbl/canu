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
  fprintf(stderr, "%s\n", descr);

  showSeq = true;

  fprintf(out, "ID: %s  pos1:"u32bitFMT","u32bitFMT" _pos="u32bitFMT" (seqlen="u32bitFMT") pos2:"u32bitFMT","u32bitFMT" _pos="u32bitFMT" (seqlen="u32bitFMT")  diag:"u32bitFMT" %s\n",
          _matchId,
          _acc1->getRangeBegin(), _acc1->getRangeEnd() - _acc1->getRangeBegin(), _acc1->_pos, _seq1->sequenceLength(),
          _acc2->getRangeBegin(), _acc2->getRangeEnd() - _acc2->getRangeBegin(), _acc2->_pos, _seq2->sequenceLength(),
          _diagonal,
          (_ori1 != _ori2) ? "reversed" : "");

  if (showSeq) {
    FastAAccessor  &A = *_acc1;
    FastAAccessor  &B = *_acc2;

    //  Save the position of the accessors
    //
    u32bit  acc1pos = _acc1->getPosition();
    u32bit  acc2pos = _acc2->getPosition();

    A.setPosition(_acc1->getRangeBegin());
    B.setPosition(_acc2->getRangeBegin());

    u32bit    margin = 15;
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
    *las++ = '|';
    for (i=0; i<_acc1->getRangeEnd() - _acc1->getRangeBegin(); i++, ++A)
      if (A.isValid())
        *las++ = *A;
      else
        *las++ = ' ';
    *las++ = '|';
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
    *las++ = '|';
    for (i=0; i<_acc1->getRangeEnd() - _acc1->getRangeBegin(); i++, ++B)
      if (B.isValid())
        *las++ = *B;
      else
        *las++ = ' ';
    *las++ = '|';
    for (i=0; i<margin; i++, ++B)
      if (B.isValid())
        *las++ = *B;
      else
        *las++ = ' ';

    *las++ = 0;
    fprintf(out, "%s\n", seq);

    delete [] seq;

    //  Restore positions
    _acc1->setPosition(acc1pos);
    _acc2->setPosition(acc2pos);
  }

#if 0
  return;



  if (showSeq) {
    u32bit    p1 = _acc1->getRangeBegin();
    u32bit    p2 = _acc2->getRangeBegin();
    u32bit    margin = 25;
    u32bit    num_print = _acc1->getRangeEnd() - _acc1->getRangeBegin() + 2 * margin;
    u32bit    i;
 
    char  *seq = new char [_acc1->getRangeEnd() - _acc1->getRangeBegin() + 32];
    char  *las = seq;

    strcpy(seq, ">>> ");
    while (*las) las++;

    //  Print the margin, or move the sequence back
    //
    if (p1 < margin) {
      for (i = 0; i < (margin - p1); ++i)
	*las++ = ' ';
      p1 = 0;
    } else {
      p1 = p1 - margin;
      i = 0;
    }
    A.setPosition(p1);

    for (; i < num_print; ++i, p1++, ++A) {
      if ((i == margin) || (i == (num_print - margin)))
	*las++ = '|';
      if (A.isValid() == false)
	*las++ = ' ';
      else
	*las++ = *A;
    }
    *las++ = 0;
    fprintf(out, "%s\n", seq);


    //  Print the margin, or move the sequence back
    //
    las = seq;

    if (_ori1 != _ori2) {
      strcpy(seq, "<<< ");
      while (*las) las++;

      if (p2 + _acc1->getRangeEnd() - _acc1->getRangeBegin() + margin > _seq2->sequenceLength()) {
        for (i = 0; i < (p2 + _acc1->getRangeEnd() - _acc1->getRangeBegin() + margin - _seq2->sequenceLength()); ++i)
          *las++ = ' ';
        //p2 = 0;
      } else {
        p2 = p2 - margin;
        i = 0;
      }
    } else {
      strcpy(seq, ">>> ");
      while (*las) las++;

      if (p2 < margin) {
        for (i = 0; i < (margin - p2); ++i)
          *las++ = ' ';
        p2 = 0;
      } else {
        p2 = p2 - margin;
        i = 0;
      }
    }
    B.setPosition(p2);

    for (; i < num_print; ++i, p2++, ++B) {
      if ((i == margin) || (i == (num_print - margin)))
	*las++ = '|';
      if (B.isValid() == false)
	*las++ = ' ';
      else
	*las++ = *B;
    }
    *las++ = 0;
    fprintf(out, "%s\n", seq);

    delete [] seq;
  }
#endif
}
