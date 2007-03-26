
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

gkpStoreSequence::gkpStoreSequence() {
  _gkp = 0L;
  _frg = 0L;
  _bgn = 0;
  _end = 0;
  _clr = AS_READ_CLEAR_LATEST;
  _iid = 0;
  _eof = false;
}

gkpStoreSequence::gkpStoreSequence(char const *gkpName,
                                   uint32 bgn,
                                   uint32 end,
                                   uint32 clr) {
  _gkp = openGateKeeperStore(gkpName, FALSE);
  _frg = new_fragRecord();
  _bgn = bgn;
  _end = end;
  _clr = clr;
  _iid = _bgn;
  _eof = false;

  assert(_clr < AS_READ_CLEAR_NUM);

  if (_end < _bgn)
    fprintf(stderr, "gkpStoreSequence()--  ERROR:  begin IID = %u > end IID = %u\n",
            _bgn, _end);
  assert(_bgn <= _end);

  //  even though, yes, we don't strictly need to be loading the
  //  sequence lengths, and we could get them from the store, this
  //  code is much simpler, and a whole lot faster, if we always load
  //  them.
  //
  _seqLengths = new uint32 [getNumberOfSequences() + 1];

  FragStream *stm = openFragStream(_gkp, FRAG_S_INF);
  while (nextFragStream(stm, _frg)) {
    bgn = getFragRecordClearRegionBegin(_frg, _clr);
    end = getFragRecordClearRegionEnd  (_frg, _clr);

    if (getFragRecordIsDeleted(_frg))
      end = bgn;
    if ((getFragRecordIID(_frg) < _bgn) || (_end < getFragRecordIID(_frg)))
      end = bgn;

    _seqLengths[getFragRecordIID(_frg)] = end - bgn;
  }

  closeFragStream(stm);
}


gkpStoreSequence::~gkpStoreSequence() {
  del_fragRecord(_frg);
  closeGateKeeperStore(_gkp);
  delete [] _seqLengths;
}


seqFile*
gkpStoreSequence::openFile(const char *name) {
  char  *p;
  char  *q;
  char   n[FILENAME_MAX + 32];

  uint32 bgn = 1;
  uint32 end = ~(uint32)0;
  uint32 clr = AS_READ_CLEAR_LATEST;

  //fprintf(stderr, "openFile gkpStoreSequence '%s'\n", name);

  //  We need to make a copy of the name so we can munge out the IID's
  //  and clear range specification.
  //
  strcpy(n, name);

  //  Parse the name, look for :'s.  If we find a ':', terminate the
  //  gkpName.  Then check if there is a '-'.  If so, decode the two
  //  IID's, and move to the next ':'.  Finally, decode the clear
  //  range spec.
  //
  p = strchr(n, ':');
  if (p) {
    *p = 0;
    q = strchr(p+1, '-');
    if (q) {
      bgn = atoi(p+1);
      end = atoi(q+1);
      p = strrchr(q, ':');
    }

    if (p)
      clr = AS_PER_decodeClearRangeLabel(p+1);
  }


  if (testOpenGateKeeperStore(n, FALSE))
    return(new gkpStoreSequence(n, bgn, end, clr));
  else
    return(0L);
};


bool
gkpStoreSequence::getSequence(uint32 &hLen, char *&h,
                              uint32 &sLen, char *&s) {

  if (_iid > _end)
    return(false);

  if (_iid != getFragRecordIID(_frg))
    if (find(_iid - 1) == false)
      return(false);

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

  return(true);
}


seqInCore *
gkpStoreSequence::getSequenceInCore(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  if (getSequence(hLen, h, sLen, s) == false)
    return(0L);

  return(new seqInCore(_iid++ - 1, h, hLen, s, sLen));
}


seqOnDisk *
gkpStoreSequence::getSequenceOnDisk(void) {
  char   *h = NULL, *s = NULL;
  uint32  hLen=0,    sLen=0;

  if (getSequence(hLen, h, sLen, s) == false)
    return(0L);

  return(new seqOnDisk(_iid++ - 1, h, hLen, s, sLen));
}



bool
gkpStoreSequence::find(seqIID  iid) {
  iid++;
  if (iid > getLastElemFragStore(_gkp)) {
    _eof = true;
    _iid = getLastElemFragStore(_gkp) + 1;
    return(false);
  }
  if ((iid < _bgn) || (_end < iid)) {
    //  Not sure this is a useful warning.  Certainly, if we're in a
    //  loop from iid=0 to iid=last, we'll get this LOTS.
    //  
    //fprintf(stderr, "gkpStoreSequence::find()--  WARNING: looked for gkp iid=%u, not in range [%u,%u]\n",
    //        iid, _bgn, _end);
    _iid = getLastElemFragStore(_gkp) + 1;
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
