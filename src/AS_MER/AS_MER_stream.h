
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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
#ifndef MERSTREAM_H
#define MERSTREAM_H

extern "C" {
#include "AS_PER_gkpStore.h"
#include "AS_PER_gkpStore.h"
}

extern unsigned char   compressSymbol[256];
extern unsigned char   validSymbol[256];
extern unsigned char   decompressSymbol[256];
extern unsigned char   complementSymbol[256];
extern unsigned char   validCompressedSymbol[256];

class merStream {
public:
  merStream(cds_uint32 merSize, char *gkpstore);
  merStream(cds_uint32 merSize, char *gkpstore, int skipNum);
  ~merStream();

  cds_uint64         theFMer(void)        { return(_theFMer); };
  cds_uint64         theRMer(void)        { return(_theRMer); };

  bool               nextMer(void);
private:
  void               loadMer(cds_uint32 s);
  unsigned char      nextSymbol(void) {

    //  If the sequence in valid, just return the next letter.
    //
    if (_thePos < _endPos)
      return((unsigned char)_theSeq[_thePos++]);

    //  Otherwise, we need to read another fragment and return an 'N'
    //  to break the merstream.

    if (_iid >= _max)
      return(0);

    getFrag(_fs, _iid, _fr, FRAG_S_SEQ);
    _iid += _skipNum;

    _thePos = getFragRecordClearRegionBegin(_fr, AS_READ_CLEAR_OBT);
    _endPos = getFragRecordClearRegionEnd  (_fr, AS_READ_CLEAR_OBT);

    _theSeq = getFragRecordSequence(_fr);

    return('N');
  };

  char                     *_theSeq;
  cds_uint32                _theLen;
  unsigned int              _thePos;
  unsigned int              _endPos;

  cds_uint32                _skipNum;
  cds_uint32                _merSize;
  cds_int32                 _timeUntilValid;
  cds_uint64                _theMerMask;

  cds_uint64                _theFMer;
  cds_uint64                _theRMer;

  cds_uint32                _theRMerShift;

  char                      _theMerString[33];

  GateKeeperStore      *_fs;
  fragRecord           *_fr;

  cds_uint64            _iid;
  cds_uint64            _max;
};




//  We have a problem; at the start of a stream, we want to
//  initialize the mer with merSize-1 bases.  In the middle of the
//  stream, we need to load the mer with merSize bases.
//
//  loadMer will push s bases onto the mer, restarting if it hits a
//  mer break.
//
//  No masking of the mer is performed.
//
inline
void
merStream::loadMer(cds_uint32 s) {
  cds_uint64   ch = 255;
  cds_uint64   cf = 0;
  cds_uint64   cr = 0;

  _timeUntilValid = s;

  //  While we are invalid, and still in the sequence
  //  push characters onto the mer.  The valid
  //  count is updated if we hit an invalid base.
  //
  while ((_timeUntilValid > 0) && (ch != 0)) {
    ch = nextSymbol();
    cf = validCompressedSymbol[ch];

    //  Rather than take the chance of generating a cache miss accessing
    //  another array, we just mask out the upper bits of validCompressedSymbol[];
    //  this is exactly the same as using compressSymbol[].
    //
    //  We need to mask the upper bits for reverse (but not for forward)
    //  because, in reverse, we shift to the right.  If we don't mask these
    //  out, we will have extra bits set in the mer.
    //
    //  Example: Consider placing 255 (== 11111111, the invalid symbol returned
    //  from validCompressedSymbol) at the fourth base in a four mer:
    //    000000**xxxxxx -- the ** are is the fourth base.
    //
    //  Without masking, we'd set all the 0's and all the *'s to one, in effect,
    //  preloading the next three bases of the mer.
    //
    cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

    _timeUntilValid--;

    if (cf & 0xfc)
      _timeUntilValid = s;

    //  If the ch is valid, we can obviously add it to the mer.
    //  If it is invalid, we don't care; by the time the mer is valid
    //  all bits of any invalid mers are removed.  Sure, we could
    //  put this into the else, but I suspect that it will be faster outside.
    //
    _theFMer <<= 2;
    _theRMer >>= 2;

    _theFMer |= cf;
    _theRMer |= cr << _theRMerShift;
  }
}

inline
bool
merStream::nextMer(void) {
  cds_uint64  ch = nextSymbol();
  cds_uint64  cf = validCompressedSymbol[ch];
  cds_uint64  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

  //  EOF?
  if (ch == 0)
    return(false);

  //  Push the ch onto the mer.
  //
  _theFMer <<= 2;
  _theRMer >>= 2;

  _theFMer |= cf;
  _theRMer |= cr << _theRMerShift;

  //  If the ch is invalid, we need to make a whole new mer.
  //
  if (cf & (unsigned char)0xfc)
    loadMer(_merSize);

  _theFMer &= _theMerMask;
  _theRMer &= _theMerMask;

  //  Still need to check if we are valid -- we could
  //  have run off the end of the sequence before a valid
  //  mer was created.
  //
  return(_timeUntilValid <= 0);
}

#endif  //  MERSTREAM_H
