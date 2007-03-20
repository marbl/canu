
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


gkpStoreSequence::gkpStoreSequence(char const *gkpName) {
  _gkp = openGateKeeperStore(gkpName, FALSE);
  _stm = openFragStream(_gkp, FRAG_S_SEQ);

  _frg = new_fragRecord();

  _iid = 1;
  _eof = false;

  _seqLengths = 0L;
}


gkpStoreSequence::~gkpStoreSequence() {
  del_fragRecord(_frg);
  closeFragStream(_stm);
  closeGateKeeperStore(_gkp);
  delete [] _seqLengths;
}


FastASequenceInCore *
gkpStoreSequence::getSequence(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  if (_iid != getFragRecordIID(_frg))
    find(_iid-1);

  h    = new char [65];
  sprintf(h, F_UID","F_IID,
          getFragRecordUID(_frg),
          getFragRecordIID(_frg));
  hLen = strlen(h);

  sLen = getFragRecordSequenceLength(_frg);
  s    = new char [sLen + 1];
  strcpy(s, getFragRecordSequence(_frg));

  //  The same stuff above is used for InCore and OnDisk.

  FastASequenceInCore *f = new FastASequenceInCore(_iid - 1,
                                                   h, hLen,
                                                   s, sLen);
  _iid++;
  return(f);
}


FastASequenceOnDisk *
gkpStoreSequence::getSequenceOnDisk(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  if (_iid != getFragRecordIID(_frg))
    find(_iid-1);

  h    = new char [65];
  sprintf(h, F_UID","F_IID,
          getFragRecordUID(_frg),
          getFragRecordIID(_frg));
  hLen = strlen(h);

  sLen = getFragRecordSequenceLength(_frg);
  s    = new char [sLen + 1];
  strcpy(s, getFragRecordSequence(_frg));

  //  The same stuff above is used for InCore and OnDisk.

  FastASequenceOnDisk *f = new FastASequenceOnDisk(_iid - 1,
                                                   h, hLen,
                                                   s, sLen);
  _iid++;
  return(f);
}


u32bit
gkpStoreSequence::sequenceLength(IID_t iid) {

  if (_seqLengths) {
    return(_seqLengths[iid+1]);
  } else {
    assert(0);
    if (iid+1 != getFragRecordIID(_frg))
      find(iid);
    return(getFragRecordSequenceLength(_frg));
  }
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



void
gkpStoreSequence::openIndex(u32bit indextypetoload) {
  FragStream *stm = openFragStream(_gkp, FRAG_S_INF);
  fragRecord *frg = new_fragRecord();

  fprintf(stderr, "gkpStoreSequence::openIndex()--  opening.\n");

  _seqLengths = new uint32 [getNumberOfSequences() + 1];

  while (nextFragStream(stm, frg))
    _seqLengths[getFragRecordIID(frg)] = getFragRecordSequenceLength(frg);

  del_fragRecord(frg);
  closeFragStream(stm);
}
