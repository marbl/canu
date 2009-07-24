
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute.
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

static const char *rcsid = "$Id: AS_MER_gkpStoreChain.C,v 1.1 2009-07-24 12:09:56 brianwalenz Exp $";

#include "AS_MER_gkpStoreChain.H"

gkpStoreChain::gkpStoreChain(char const *gkpName,
                             uint32 clr) {
  clear();

  _gkp    = new gkStore(gkpName, FALSE, FALSE);

  strcpy(_filename, _gkp->gkStore_path());

  //  maxChainLen must be > 1024*1024.  Set to 128*1024 for valgrind fun.

  _maxChains           = 1024;
  _maxChainLen         = 512 * 1024 * 1024;

  _numberOfSequences = 0;

  _chainIID    = ~u32bitZERO;
  _chainLen    = new u32bit [_maxChains];
  _chainBgnFrg = new u32bit [_maxChains];

  for (uint32 i=0; i<_maxChains; i++) {
    _chainLen[i]    = 0;
    _chainBgnFrg[i] = 0;
  }

  _chainBgnFrg[0] = 1;

  gkStream *stm = new gkStream(_gkp, 0, 0, GKFRAGMENT_INF);

  while (stm->next(&_frg)) {
    if (_frg.gkFragment_getIsDeleted() == false) {
      u32bit  iid = _frg.gkFragment_getReadIID();
#ifdef FULLREAD
      u32bit  len = _frg.gkFragment_getSequenceLength();
#else
      u32bit  len = _frg.gkFragment_getClearRegionLength();
#endif

      if (_chainLen[_numberOfSequences] + len > _maxChainLen) {
        _numberOfSequences++;
        _chainBgnFrg[_numberOfSequences] = iid;
      }

      //  This'll take 1TB of sequence to trip
      assert(_numberOfSequences < _maxChains);

      _chainLen[_numberOfSequences] += len + 1;
    }
  }

  delete stm;

  if (_chainLen[_numberOfSequences] != 0)
    _numberOfSequences++;

  _chainBgnFrg[_numberOfSequences] = _gkp->gkStore_getNumFragments();

  //for (uint32 i=0; i<=_numberOfSequences; i++)
  //  fprintf(stderr, "gkpStoreChain::gkpStoreChain()-- bgn=%u len=%u\n", _chainBgnFrg[i], _chainLen[i]);
}



gkpStoreChain::gkpStoreChain() {
  clear();
}



gkpStoreChain::~gkpStoreChain() {
  delete    _gkp;
  delete [] _chainLen;
  delete [] _chainBgnFrg;
  delete [] _frgLengths;
}



seqFile*
gkpStoreChain::openFile(const char *name) {
  char  *p;
  char   n[FILENAME_MAX + 32];

  //  We need to make a copy of the name so we can munge out the IID's
  //  and clear range specification.
  //
  strcpy(n, name);

  //  Parse the name, look for :'s.  If we find a ':', terminate the
  //  gkpName.  Then decode the clear range spec.
  //
  p = strchr(n, ':');
  if (p) {
    *p = 0;

    if (strcmp(p+1, "chain") != 0)
      p = 0L;
  }

  if (p == 0L)
    return(0L);

  return(new gkpStoreChain(n, AS_READ_CLEAR_LATEST));
};



bool
gkpStoreChain::getSequence(u32bit iid,
                           char *&h, u32bit &hLen, u32bit &hMax,
                           char *&s, u32bit &sLen, u32bit &sMax) {

  if (_chainLen[iid] == 0)
    return(false);

  if (sMax < _chainLen[iid]) {
    delete [] s;
    sMax = _chainLen[iid];
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 1024;
    h    = new char [hMax];
  }

  sprintf(h, "iid:"u32bitFMT"-"u32bitFMT, _chainBgnFrg[iid], _chainBgnFrg[iid+1]-1);

  hLen = strlen(h);
  sLen = 0;

  for (uint32 frg=_chainBgnFrg[iid]; frg<_chainBgnFrg[iid+1]; frg++) {
    _gkp->gkStore_getFragment(frg, &_frg, GKFRAGMENT_SEQ);

#if FULLREAD
    strcpy(s + sLen, _frg.gkFragment_getSequence());
    sLen += _frg.gkFragment_getSequenceLength();
    s[sLen++] = '.';
#else
    strncpy(s + sLen, _frg.gkFragment_getSequence() + _frg.gkFragment_getClearRegionBegin(), _frg.gkFragment_getClearRegionLength());
    sLen += _frg.gkFragment_getClearRegionLength();
    s[sLen++] = '.';
#endif
  }

  s[sLen] = 0;

  return(true);
}



bool
gkpStoreChain::getSequence(u32bit iid,
                           u32bit bgn, u32bit end, char *s) {

  char   *t = s;

  *s = 0;

  //fprintf(stderr, "gkpStoreChain::getSequence()-- iid=%u len=%u bgn=%u end=%u\n", iid, _chainLen[iid], bgn, end);

  if (_chainLen[iid] == 0)
    return(false);

  if ((_chainIID != iid) || (_frgLengths == 0L)) {
    delete [] _frgLengths;
    _frgLengths = new u32bit [_chainBgnFrg[iid+1] - _chainBgnFrg[iid]];
    _chainIID   = iid;

    for (u32bit i=0; i<_chainBgnFrg[iid+1] - _chainBgnFrg[iid]; i++)
      _frgLengths[i] = ~u32bitZERO;

    _lastFrg = _chainBgnFrg[iid];
    _lastPos = 0;
  }


  //  Jump to where were last time.  Usually, this function is called to stream through a sequence,
  //  so we can generally skip the "skip fragments that come before bgn" loop.  We can't do this if
  //  the bgn is before the _pastPos saved below.
  //
  uint32 bas = _chainBgnFrg[iid];
  uint32 frg = _lastFrg;
  uint32 pos = _lastPos;

  if (bgn < _lastPos) {
    bas = _chainBgnFrg[iid];
    frg = _chainBgnFrg[iid];
    pos = 0;
  }

  //  Skip fragments that come before bgn
  //
  while ((frg < _chainBgnFrg[iid+1]) &&
         (_frgLengths[frg-bas] != ~u32bitZERO) &&
         (pos + _frgLengths[frg-bas] < bgn)) {
    pos += _frgLengths[frg-bas];
    frg++;
  }

  //fprintf(stderr, "i=%d pos=%d\n", i, pos);

  //  We're either ready to start copying fragments, or we're before bgn, but we don't know the
  //  fragment length.  Grab the fragment.
  //
  while ((frg < _chainBgnFrg[iid+1]) &&
         (pos < end)) {

    _gkp->gkStore_getFragment(frg, &_frg, GKFRAGMENT_SEQ);

#ifdef FULLREAD
    u32bit  clr = 0;
    u32bit  len = _frg.gkFragment_getSequenceLength();
#else
    u32bit  clr = _frg.gkFragment_getClearRegionBegin();
    u32bit  len = _frg.gkFragment_getClearRegionLength();
#endif

    _frgLengths[frg-bas] = len + 1;

    if        (_frg.gkFragment_getIsDeleted() == true) {
      //  0)  Fragment is deleted
      frg++;
      continue;

    } else if (pos + len + 1 <= bgn) {
      //  1)  Fragment completely before bgn

    } else if (pos < bgn) {
      //  2) Fragment partially at the start
      //
      //  Nothing to copy if pos + len + 1 == bgn, just add the '.', otherwise, add pos + len - bgn
      //  bases, then add the '.'.
      //  
      if (pos + len + 1 > bgn) {
        strncpy(t, _frg.gkFragment_getSequence() + clr + (bgn - pos), len - (bgn - pos));
        t += len - (bgn - pos);
      }
      *t++ = '.';
      *t = 0;
#if 0
      if (bgn == 175112190)
        fprintf(stderr, "2 pos=%d len=%d bgn=%d %s\n", pos, len, bgn, t - len + (bgn - pos) - 1);
#endif
      //assert(strlen(s) == pos - bgn + len + 1);

    } else if ((bgn <= pos) && (pos + len + 1 <= end)) {
      //  3)  Fragment and '.' completely contained
      strncpy(t, _frg.gkFragment_getSequence() + clr, len);
      t += len;
      *t++ = '.';
      *t = 0;
#if 0
      if (bgn == 175112190)
        fprintf(stderr, "3 pos=%d len=%d bgn=%d %s\n", pos, len, bgn, t - len - 1);
#endif
      //assert(strlen(s) == pos - bgn + len + 1);

    } else if (pos < end) {
      //  4)  Fragment partially at the end
      strncpy(t, _frg.gkFragment_getSequence() + clr, (end - pos));
      t += (end - pos);
      *t = 0;
#if 0
      if (bgn == 175112190)
        fprintf(stderr, "4 pos=%d len=%d bgn=%d %s\n", pos, len, bgn, t - (end - pos));
#endif
      //assert(strlen(s) == pos - bgn + (end - pos));

    } else if (end <= pos) {
      //  5)  Fragment completely after end

    } else {
      assert(0);
    }

    //  Remember which fragment we're at, and where it was located
    _lastFrg = frg;
    _lastPos = pos;

    //  Advance to the next fragment (or past the end, in which case we don't care about the values
    //  anymore.
    pos += len + 1;
    frg++;
  }

  //  Terminate the sequence.
  *t = 0;

  //fprintf(stderr, "i=%d pos=%d len=%d\n", i, pos, strlen(s));

  return(true);
}



void
gkpStoreChain::clear(void) {
  _gkp          = 0L;
  _maxChains    = 0;
  _maxChainLen  = 0;

  _chainIID     = 0;
  _chainLen     = 0L;
  _chainBgnFrg  = 0L;

  _frgLengths   = 0L;

  strcpy(_filename, "");
  strcpy(_typename, "gkpStoreChain");

  _numberOfSequences = 0;
}
