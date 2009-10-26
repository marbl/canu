
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

static const char *rcsid = "$Id: AS_MER_gkpStore_to_FastABase.C,v 1.18 2009-10-26 13:20:26 brianwalenz Exp $";

#include "AS_MER_gkpStore_to_FastABase.H"

gkpStoreFile::gkpStoreFile(char const *gkpName,
                           uint32 bgn,
                           uint32 end,
                           uint32 clr) {
  clear();

  _gkp    = new gkStore(gkpName, FALSE, FALSE);
  _bgn    = bgn;
  _end    = end;

  strcpy(_filename, _gkp->gkStore_path());

  _numberOfSequences = _gkp->gkStore_getNumFragments();

  if (_end > _numberOfSequences)
    _end = _numberOfSequences;

  if (_end < _bgn)
    fprintf(stderr, "gkpStoreFile()--  ERROR:  begin IID = %u > end IID = %u\n", _bgn, _end);
  assert(_bgn <= _end);

  _clrBeg = new uint16 [_numberOfSequences + 1];
  _clrEnd = new uint16 [_numberOfSequences + 1];

  for (uint32 i=0; i<=_numberOfSequences; i++) {
    _clrBeg[i] = 0;
    _clrEnd[i] = 0;
  }

  gkStream *stm = new gkStream(_gkp, _bgn, _end, GKFRAGMENT_INF);

  while (stm->next(&_frg)) {
    if (!_frg.gkFragment_getIsDeleted()) {
      uint32  iid = _frg.gkFragment_getReadIID();
      _clrBeg[iid] = _frg.gkFragment_getClearRegionBegin(clr);
      _clrEnd[iid] = _frg.gkFragment_getClearRegionEnd  (clr);
    }
  }

  delete stm;
}



gkpStoreFile::gkpStoreFile() {
  clear();
}



gkpStoreFile::~gkpStoreFile() {
  delete    _gkp;
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
      clr = gkStore_decodeClearRegionLabel(p+1);
  }

  if (clr == AS_READ_CLEAR_ERROR)
    return(0L);

  //  if can't open store, or file not there, return null

  return(new gkpStoreFile(n, bgn, end, clr));
};



bool
gkpStoreFile::getSequence(u32bit iid,
                          char *&h, u32bit &hLen, u32bit &hMax,
                          char *&s, u32bit &sLen, u32bit &sMax) {

  if (sMax == 0) {
    sMax = AS_READ_MAX_NORMAL_LEN+1;
    s    = new char [sMax];
  }

  if (hMax == 0) {
    hMax = 1024;
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

  _gkp->gkStore_getFragment(iid, &_frg, GKFRAGMENT_SEQ);

  sprintf(h, "%s,"F_IID, AS_UID_toString(_frg.gkFragment_getReadUID()), _frg.gkFragment_getReadIID());

  hLen = strlen(h);
  sLen = _clrEnd[iid] - _clrBeg[iid];

  if (sLen > 0)
    strncpy(s, _frg.gkFragment_getSequence() + _clrBeg[iid], sLen);

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

  _gkp->gkStore_getFragment(iid, &_frg, GKFRAGMENT_SEQ);

  strncpy(s, _frg.gkFragment_getSequence() + _clrBeg[iid] + bgn, end - bgn);

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
