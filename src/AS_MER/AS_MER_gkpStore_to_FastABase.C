
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

static const char *rcsid = "$Id: AS_MER_gkpStore_to_FastABase.C,v 1.14 2008-10-08 22:02:57 brianwalenz Exp $";

#include "AS_MER_gkpStore_to_FastABase.H"

gkpStoreFile::gkpStoreFile(char const *gkpName,
                           uint32 bgn,
                           uint32 end,
                           uint32 clr) {
  clear();

  _gkp    = openGateKeeperStore(gkpName, FALSE);
  _bgn    = bgn;
  _end    = end;

  strcpy(_filename, _gkp->storePath);

  _numberOfSequences = getLastElemFragStore(_gkp) + 1;

  if (_end > _numberOfSequences)
    _end = _numberOfSequences;

  if (_end < _bgn)
    fprintf(stderr, "gkpStoreFile()--  ERROR:  begin IID = %u > end IID = %u\n",
            _bgn, _end);

  assert(_bgn <= _end);

  _clrBeg = new uint16 [_numberOfSequences];
  _clrEnd = new uint16 [_numberOfSequences];

  for (uint32 i=0; i<_numberOfSequences; i++) {
    _clrBeg[i] = 0;
    _clrEnd[i] = 0;
  }

  FragStream *stm = openFragStream(_gkp, FRAG_S_INF);

  if ((_bgn < _numberOfSequences) && (end < _numberOfSequences))
    resetFragStream(stm, _bgn, _end);

  while (nextFragStream(stm, &_frg)) {
    if (!getFragRecordIsDeleted(&_frg)) {
      uint32  iid = getFragRecordIID(&_frg);
      _clrBeg[iid] = getFragRecordClearRegionBegin(&_frg, clr);
      _clrEnd[iid] = getFragRecordClearRegionEnd  (&_frg, clr);
    }
  }

  closeFragStream(stm);
}



gkpStoreFile::gkpStoreFile() {
  clear();
}



gkpStoreFile::~gkpStoreFile() {
  closeGateKeeperStore(_gkp);
  delete [] _clrBeg;
  delete [] _clrEnd;
}



seqFile*
gkpStoreFile::openFile(const char *name) {
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
    return(new gkpStoreFile(n, bgn, end, clr));
  else
    return(0L);
};



bool
gkpStoreFile::getSequence(u32bit iid,
                          char *&h, u32bit &hLen, u32bit &hMax,
                          char *&s, u32bit &sLen, u32bit &sMax) {

  if (sMax == 0) {
    sMax = AS_FRAG_MAX_LEN;
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 2048;
    h    = new char [hMax];
  }

  if (iid == 0) {
    h[0] = 0;
    s[0] = 0;

    hLen = 0;
    sLen = 0;

    return(true);
  }

  if ((iid < _bgn) || (_end < iid))
    return(false);

  getFrag(_gkp, iid, &_frg, FRAG_S_SEQ);

  sprintf(h, "%s,"F_IID, AS_UID_toString(getFragRecordUID(&_frg)), getFragRecordIID(&_frg));

  hLen = strlen(h);
  sLen = _clrEnd[iid] - _clrBeg[iid];

  if (sLen > 0)
    strncpy(s, getFragRecordSequence(&_frg) + _clrBeg[iid], sLen);

  s[sLen] = 0;

  return(true);
}



bool
gkpStoreFile::getSequence(u32bit iid,
                          u32bit bgn, u32bit end, char *s) {

  if (iid == 0)
    fprintf(stderr, "gkpStoreFile::getSequence(part)-- someone requested iid==0?\n"), exit(1);

  if ((iid < _bgn) || (_end < iid))
    return(false);

  getFrag(_gkp, iid, &_frg, FRAG_S_SEQ);

  strncpy(s, getFragRecordSequence(&_frg) + _clrBeg[iid] + bgn, end - bgn);

  s[end - bgn] = 0;

  return(true);
}



void
gkpStoreFile::clear(void) {
  _gkp = 0L;
  _bgn = 0;
  _end = 0;

  strcpy(_filename, "");
  strcpy(_typename, "gkpStore");

  _numberOfSequences = 0;
}
