
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

#include "AS_MER_gkpStore_to_FastABase.H"


gkpStoreSequence::gkpStoreSequence(char const *gkpName, uint32 clearRangeSpec) {
  _gkp = openGateKeeperStore(gkpName, FALSE);
  _frg = new_fragRecord();
  _clr = clearRangeSpec;
  _iid = 1;
  _eof = false;

  assert(_clr < AS_READ_CLEAR_NUM);

  //  even though, yes, we don't strictly need to be loading the
  //  sequence lengths, and we could get them from the store, this
  //  code is much simpler, and a whole lot faster, if we always load
  //  them.
  //
  _seqLengths = new uint32 [getNumberOfSequences() + 1];

  FragStream *stm = openFragStream(_gkp, FRAG_S_INF);
  while (nextFragStream(stm, _frg)) {
    uint32  beg = getFragRecordClearRegionBegin(_frg, _clr);
    uint32  end = getFragRecordClearRegionEnd  (_frg, _clr);

    if (getFragRecordIsDeleted(_frg))
      end = beg;

    _seqLengths[getFragRecordIID(_frg)] = end - beg;
  }

  closeFragStream(stm);
}


gkpStoreSequence::~gkpStoreSequence() {
  del_fragRecord(_frg);
  closeGateKeeperStore(_gkp);
  delete [] _seqLengths;
}


void
gkpStoreSequence::getSequence(uint32 iid,
                              uint32 &hLen, char *&h,
                              uint32 &sLen, char *&s) {

  _iid = iid + 1;

  if (_iid != getFragRecordIID(_frg))
    find(_iid - 1);

  h    = new char [65];
  sprintf(h, F_UID","F_IID,
          getFragRecordUID(_frg),
          getFragRecordIID(_frg));
  hLen = strlen(h);

  sLen = _seqLengths[_iid];
  s    = new char [sLen + 1];

  if (sLen > 0)
    strncpy(s, getFragRecordSequence(_frg) + getFragRecordClearRegionBegin(_frg, _clr), sLen);

  s[sLen] = 0;
}


FastASequenceInCore *
gkpStoreSequence::getSequence(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  getSequence(_iid, hLen, h, sLen, s);

  return(new FastASequenceInCore(_iid++ - 1, h, hLen, s, sLen));
}


FastASequenceOnDisk *
gkpStoreSequence::getSequenceOnDisk(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  getSequence(_iid, hLen, h, sLen, s);

  return(new FastASequenceOnDisk(_iid++ - 1, h, hLen, s, sLen));
}


u32bit
gkpStoreSequence::sequenceLength(IID_t iid) {
  return(_seqLengths[iid+1]);
}


u32bit
gkpStoreSequence::headerLength(IID_t iid) {
  //  Not critical, just an upper bound, I think.
  return(64);
}


bool
gkpStoreSequence::find(IID_t  iid) {
  iid++;
  if (iid > getLastElemFragStore(_gkp)) {
    _eof = true;
    return(false);
  }
  getFrag(_gkp, iid, _frg, FRAG_S_SEQ);
  _iid = iid;
  return(true);
}


bool
gkpStoreSequence::find(char  *id) {
  return(false);
}
