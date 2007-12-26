
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
  _bgn = 0;
  _end = 0;
  _clr = AS_READ_CLEAR_LATEST;
  _eof = false;
  _tst = 0;
}

gkpStoreSequence::gkpStoreSequence(char const *gkpName,
                                   uint32 bgn,
                                   uint32 end,
                                   uint32 clr) {

  _gkp = openGateKeeperStore(gkpName, FALSE);
  _bgn = bgn;
  _end = end;
  _clr = clr;
  _curIID = _bgn;
  _eof = false;
  _tst = 0;

  {
    char    pt[FILENAME_MAX];
    stat_s  st;

    sprintf(pt, "%s/frg", gkpName);

    if (stat(pt, &st) != 0) {
      fprintf(stderr, "Couldn't stat() '%s'\n%s\n", pt, strerror(errno));
      exit(1);
    }

    _tst = st.st_mtime;
  }

  if (_end < _bgn)
    fprintf(stderr, "gkpStoreSequence()--  ERROR:  begin IID = %u > end IID = %u\n",
            _bgn, _end);
  assert(_bgn <= _end);

  uint32    max = getNumberOfSequences() + 1;

  _seqLen = new uint16 [max];
  _clrBeg = new uint16 [max];
  _clrEnd = new uint16 [max];
  
  for (uint32 i=0; i<max; i++) {
    _seqLen[i] = 0;
    _clrBeg[i] = 0;
    _clrEnd[i] = 0;
  }

  FragStream *stm = openFragStream(_gkp, FRAG_S_INF);

  if ((_bgn < max) && (end < max))
    resetFragStream(stm, _bgn, _end);

  while (nextFragStream(stm, &_frg)) {
    if (!getFragRecordIsDeleted(&_frg)) {
      uint32  iid = getFragRecordIID(&_frg);
      _clrBeg[iid] = getFragRecordClearRegionBegin(&_frg, _clr);
      _clrEnd[iid] = getFragRecordClearRegionEnd  (&_frg, _clr);
      _seqLen[iid] = getFragRecordSequenceLength  (&_frg);
    }
  }

  closeFragStream(stm);
}


gkpStoreSequence::~gkpStoreSequence() {
  closeGateKeeperStore(_gkp);
  delete [] _seqLen;
  delete [] _clrBeg;
  delete [] _clrEnd;
}


seqFile*
gkpStoreSequence::openFile(const char *name) {
  char  *p;
  char  *q;
  char   n[FILENAME_MAX + 32];

  uint32 bgn = 0;
  uint32 end = ~(uint32)0;
  uint32 clr = AS_READ_CLEAR_LATEST;

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
gkpStoreSequence::getSequence(u32bit &hLen, char *&h,
                              u32bit &sLen, char *&s) {

  if (_curIID > _end)
    return(false);

  if ((_curIID > 0) && (_curIID != getFragRecordIID(&_frg)))
    if (find(_curIID) == false)
      return(false);

  h    = new char [65];
  sprintf(h, "%s,"F_IID, AS_UID_toString(getFragRecordUID(&_frg)), getFragRecordIID(&_frg));
  hLen = strlen(h);

  sLen = _clrEnd[_curIID] - _clrBeg[_curIID];
  s    = new char [sLen + 1];

  if (sLen > 0)
    strncpy(s, getFragRecordSequence(&_frg) + _clrBeg[_curIID], sLen);

  s[sLen] = 0;

  return(true);
}


bool
gkpStoreSequence::find(seqIID  iid) {

  if (iid == 0) {
    return(false);
  }
  if (iid > getLastElemFragStore(_gkp) + 1) {
    _eof = true;
    _curIID = getLastElemFragStore(_gkp) + 2;
    //fprintf(stderr, "gkpStoreSequence()-- find(%d) failed, too big\n", iid);
    return(false);
  }
  if ((iid < _bgn) || (_end < iid)) {
    _curIID = getLastElemFragStore(_gkp) + 2;
    //fprintf(stderr, "gkpStoreSequence()-- find(%d) failed, out of range\n", iid);
    return(false);
  }
  getFrag(_gkp, iid, &_frg, FRAG_S_SEQ);
  _curIID = iid;
  return(true);
}


bool
gkpStoreSequence::find(char  *id) {
  fprintf(stderr, "gkpStoreSequence::find()-- %s (NOT IMPLEMENTED!)\n", id);
  return(false);
}
