
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2005-2007, J. Craig Venter Institute.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

//static const char *rcsid = "$Id$";

#include "finalTrim.H"


//  sort the position values on the 5' end -- this sorts increasingly
int
position_compare5(const void *a, const void *b) {
  OVSoverlap  *A = (OVSoverlap *)a;
  OVSoverlap  *B = (OVSoverlap *)b;

  if (A->dat.obt.a_beg < B->dat.obt.a_beg)  return(-1);
  if (A->dat.obt.a_beg > B->dat.obt.a_beg)  return(1);
  return(0);
}


//  sort the position values on the 3' end -- this sorts decreasingly
int
position_compare3(const void *a, const void *b) {
  OVSoverlap  *A = (OVSoverlap *)a;
  OVSoverlap  *B = (OVSoverlap *)b;

  if (A->dat.obt.a_end < B->dat.obt.a_end)  return(1);
  if (A->dat.obt.a_end > B->dat.obt.a_end)  return(-1);
  return(0);
}


class obtConsolidate {
public:
  obtConsolidate(OVSoverlap *ovl, uint32 ovlLen, uint32 ibgn);
  ~obtConsolidate() {
  };

  uint32 _min5;    //  minimum value we've ever seen
  uint32 _minm5;   //  minimum value we've ever seen more than once
  uint32 _minm5c;  //  number of times we've seen minm5
  uint32 _mode5;   //  mode
  uint32 _mode5c;  //  number of time we've seen the mode
  uint32 _max3;
  uint32 _maxm3;
  uint32 _maxm3c;
  uint32 _mode3;
  uint32 _mode3c;
};


obtConsolidate::obtConsolidate(OVSoverlap *ovl, uint32 ovlLen, uint32 ibgn) {
  _min5   = UINT32_MAX;
  _minm5  = UINT32_MAX;
  _minm5c = 0;
  _mode5  = UINT32_MAX;
  _mode5c = 1;
  _max3   = UINT32_MAX;
  _maxm3  = UINT32_MAX;
  _maxm3c = 0;
  _mode3  = UINT32_MAX;
  _mode3c = 1;

  if ((ovl == NULL) || (ovlLen == 0))
    return;

  uint32  mtmp;  //  Mode we are computing, value
  uint32  mcnt;  //  Mode we are computing, number of times we've seen it

  qsort(ovl, ovlLen, sizeof(OVSoverlap), position_compare5);

  _min5   = ovl[0].dat.obt.a_beg;
  _mode5  = mtmp = ovl[0].dat.obt.a_beg;
  _mode5c = mcnt = 1;

  for (uint32 i=1; i<ovlLen; i++) {

    //  5' end.  Scan the list, remembering the best mode we've seen so far.  When a better one
    //  arrives, we copy it to the saved one -- and keep copying it as it gets better.

    if (mtmp == ovl[i].dat.obt.a_beg) {  //  Same mode?  Count.
      mcnt++;
    } else {
      mtmp = ovl[i].dat.obt.a_beg;  //  Different mode, restart.
      mcnt = 1;
    }

    if (mcnt > _mode5c) {  //  Bigger mode?  Save it.
      _mode5  = mtmp;
      _mode5c = mcnt;
    }

    //  If our mode is more than one and we've not seen a multiple hit before
    //  save this position.
    //
    if ((_mode5c > 1) && (_minm5 == UINT32_MAX))
      _minm5  = _mode5;

    if (_minm5 == _mode5)
      _minm5c = _mode5c;
  }

  //  Do it all again for the 3' -- remember that we've sorted this decreasingly.
  //
  qsort(ovl, ovlLen, sizeof(OVSoverlap), position_compare3);

  _max3   = ovl[0].dat.obt.a_end;
  _mode3  = mtmp = ovl[0].dat.obt.a_end;
  _mode3c = mcnt = 1;

  for (uint32 i=1; i<ovlLen; i++) {
    if (mtmp == ovl[i].dat.obt.a_end) {
      mcnt++;
    } else {
      mtmp = ovl[i].dat.obt.a_end;
      mcnt = 1;
    }

    if (mcnt > _mode3c) {
      _mode3  = mtmp;
      _mode3c = mcnt;
    }

    if ((_mode3c > 1) && (_maxm3 == UINT32_MAX))
      _maxm3  = _mode3;

    if (_maxm3 == _mode3)
      _maxm3c = _mode3c;
  }

  //  Update to be relative to the start of the read, not the start of the clear range.
  //  The count variables (*c) are left in for 'documentation'.

  _min5   += ibgn;
  _minm5  += ibgn;
  _minm5c += 0;
  _mode5  += ibgn;
  _mode5c += 0;
  _max3   += ibgn;
  _maxm3  += ibgn;
  _maxm3c += 0;
  _mode3  += ibgn;
  _mode3c += 0;
}
