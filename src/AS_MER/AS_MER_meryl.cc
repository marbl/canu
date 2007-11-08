
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_gkpStore.h"
}


//  Parameters for the mer counter, and a nice utility function.
//
class mcDescription {
public:
  uint32      _merSizeInBases;
  uint32      _merSizeInBits;

  uint32      _tableSizeInBits;
  uint64      _tableSizeInEntries;

  uint32      _chckBits;
  uint64      _chckMask;

  uint32      _hashWidth;
  uint64      _hashMask;

  uint64      _actualNumberOfMers;

  uint64 HASH(uint64 a) {
    return((a >> _chckBits) & _hashMask);
  }
};

mcDescription   mcd;

uint64     *_chck;
uint64     *_hash;

uint64       mersToCount;

unsigned char   compressSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   validSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   decompressSymbol[256] = { 65, 67, 71, 84, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   complementSymbol[256] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84, 86, 71, 72, 0, 0, 67, 68, 0, 0, 77, 0, 75, 78, 0, 0, 0, 89, 87, 65, 65, 66, 83, 0, 82, 0, 0, 0, 0, 0, 0, 0, 116, 118, 103, 104, 0, 0, 99, 100, 0, 0, 109, 0, 107, 110, 0, 0, 0, 121, 119, 97, 97, 98, 115, 0, 114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
unsigned char   validCompressedSymbol[256] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 1, 255, 255, 255, 2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

#ifdef TRUE64BIT
#define UINT64_ZERO   0x0000000000000000UL
#define UINT64_ONE    0x0000000000000001UL
#else
#define UINT64_ZERO   0x0000000000000000ULL
#define UINT64_ONE    0x0000000000000001ULL
#endif

#define UINT64_MASK(X)   ((~UINT64_ZERO) >> (64 - (X)))


class merStream {
public:
  merStream(uint32 merSize, char *gkpstore, int skipNum);
  ~merStream();

  uint64         theFMer(void)        { return(_theFMer); };
  uint64         theRMer(void)        { return(_theRMer); };

  bool               nextMer(void);
private:
  void               loadMer(uint32 s);
  unsigned char      nextSymbol(void) {

    //  If the sequence in valid, just return the next letter.
    //
    if (_thePos < _endPos)
      return((unsigned char)_theSeq[_thePos++]);

    //  Otherwise, we need to read another fragment and return an 'N'
    //  to break the merstream.

    if (_iid >= _max)
      return(0);

    getFrag(_fs, _iid, &_fr, FRAG_S_SEQ);
    _iid += _skipNum;

    _thePos = getFragRecordClearRegionBegin(&_fr, AS_READ_CLEAR_OBT);
    _endPos = getFragRecordClearRegionEnd  (&_fr, AS_READ_CLEAR_OBT);

    _theSeq = getFragRecordSequence(&_fr);

    return('N');
  };

  char                 *_theSeq;
  uint32                _theLen;
  uint32                _thePos;
  uint32                _endPos;

  uint32                _skipNum;
  uint32                _merSize;
  int32                 _timeUntilValid;
  uint64                _theMerMask;

  uint64                _theFMer;
  uint64                _theRMer;

  uint32                _theRMerShift;

  GateKeeperStore      *_fs;
  fragRecord            _fr;

  uint64                _iid;
  uint64                _max;
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
merStream::loadMer(uint32 s) {
  uint64   ch = 255;
  uint64   cf = 0;
  uint64   cr = 0;

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
  uint64  ch = nextSymbol();
  uint64  cf = validCompressedSymbol[ch];
  uint64  cr = validCompressedSymbol[complementSymbol[ch]] & 0x03;

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


merStream::merStream(uint32 merSize, char *gkpstore, int skipNum) {
  _theSeq          = 0L;
  _theLen          = 0;
  _thePos          = 0;
  _endPos          = 0;

  _merSize         = merSize;
  _timeUntilValid  = 0;
  _theMerMask      = UINT64_MASK(_merSize << 1);

  _theFMer         = 0;
  _theRMer         = 0;

  _theRMerShift    = (_merSize << 1) - 2;

  _fs = openGateKeeperStore(gkpstore, FALSE);
  if (_fs == NULL) {
    fprintf(stderr, "ERROR:  couldn't open the gatekeeper store '%s'\n", gkpstore);
    exit(1);
  }

  _iid = 1;
  _max = getLastElemFragStore(_fs);

  _skipNum=skipNum;
  loadMer(_merSize - 1);
}


merStream::~merStream() {
  closeGateKeeperStore(_fs);
}


inline
uint64
getDecodedValue(uint64 *ptr,
                uint64  pos,
                uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 b2  = siz - b1;  //  Only used if siz > b1
  uint64 ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);
  } else {
    ret  = (ptr[wrd] & UINT64_MASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & UINT64_MASK(b2);
  }

  ret &= UINT64_MASK(siz);

  return(ret);
}



inline
void
setDecodedValue(uint64 *ptr,
                uint64  pos,
                uint64  siz,
                uint64  val) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 b2  = siz - b1;  //  Only used if siz > b1

  val &= UINT64_MASK(siz);

  if (b1 >= siz) {
    ptr[wrd] &= ~( UINT64_MASK(siz) << (b1-siz) );
    ptr[wrd] |= val << (b1-siz);
  } else {
    ptr[wrd] &= ~UINT64_MASK(b1);
    ptr[wrd] |= (val & (UINT64_MASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(UINT64_MASK(b2) << (64-b2));
    ptr[wrd+1] |= (val & (UINT64_MASK(b2))) << (64-b2);
  }
}



inline
uint64
preDecrementDecodedValue(uint64 *ptr,
                         uint64  pos,
                         uint64  siz) {
  uint64 wrd = (pos >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (pos     ) & 0x000000000000003fllu;
  uint64 b1  = 64 - bit;
  uint64 b2  = siz - b1;  //  Only used if siz > b1
  uint64 ret = 0;

  if (b1 >= siz) {
    ret = ptr[wrd] >> (b1 - siz);

    ret--;
    ret &= UINT64_MASK(siz);

    ptr[wrd] &= ~( UINT64_MASK(siz) << (b1-siz) );
    ptr[wrd] |= ret << (b1-siz);
  } else {
    ret  = (ptr[wrd] & UINT64_MASK(b1)) << b2;
    ret |= (ptr[wrd+1] >> (64 - b2)) & UINT64_MASK(b2);

    ret--;
    ret &= UINT64_MASK(siz);

    ptr[wrd] &= ~UINT64_MASK(b1);
    ptr[wrd] |= (ret & (UINT64_MASK(b1) << (b2))) >> (b2);

    ptr[wrd+1] &= ~(UINT64_MASK(b2) << (64-b2));
    ptr[wrd+1] |= (ret & (UINT64_MASK(b2))) << (64-b2);
  }

  return(ret);
}



////////////////////////////////////////////////////////////////////////////////
//
//  BUILD
//
//

void
createHashTable(char *inputFile, uint32 skipNum) {
  uint64  mer;

  uint32 *_ctbl = new uint32 [ mcd._tableSizeInEntries ];
  for (uint32 i=mcd._tableSizeInEntries; i--; )
    _ctbl[i] = 0;

  //  I think we can get by with 24 bits of count here; the fragstore
  //  seems to have only 9 million mers in the first bucket.  Should
  //  probably check for overflow, though.
  //  
  merStream          M(mcd._merSizeInBases, inputFile, skipNum);

  while (M.nextMer() && (mcd._actualNumberOfMers < mersToCount)) {
    mer = M.theFMer();
    if (mer > M.theRMer())
      mer = M.theRMer();

    _ctbl[ mcd.HASH(mer) ]++;

    mcd._actualNumberOfMers++;
  }


  //
  //  Allocate a PACKED array for the hash table.  This needs to be
  //  packed only for mcd._actualNumberOfMers > 4 billion, really.
  //

  //  Determine how many bits we need to hold the value
  //  mcd._actualNumberOfMers.....then....
  //
  //  This is mcd._actualNumberOfMers+1 because we need to store the
  //  first position after the last mer.  That is, if there are two
  //  mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  mcd._hashWidth  = 1;
  while ((mcd._actualNumberOfMers+1) > (UINT64_ONE << mcd._hashWidth))
    mcd._hashWidth++;

  //  ....allocate a hash table that is that many bits wide.
  //
  _hash = new uint64 [(mcd._tableSizeInEntries+1) * mcd._hashWidth / 64 + 2];

  //
  //  Create the hash index using the counts.  The hash points
  //  to the end of the bucket; when we add a word, we move the
  //  hash bucket pointer down one.
  //
  //  When done, we can deallocate the counting table.
  //
  uint64 i=0;
  uint64 j=0;
  uint64 c=0;

  while (i < mcd._tableSizeInEntries) {
    c += _ctbl[i++];
    setDecodedValue(_hash, j, mcd._hashWidth, c);
    j += mcd._hashWidth;
  }

  //  Add the location of the end of the table.  This is not
  //  modified when adding words, but is used to determine
  //  the size of the last bucket.
  //
  setDecodedValue(_hash, j, mcd._hashWidth, c);

  delete [] _ctbl;
}


void
verifyHashTable(void) {
  uint64 i=0, j=0, c=0, d=0;
  uint64 fail=0;

  while (i <= mcd._tableSizeInEntries) {
    c = getDecodedValue(_hash, j, mcd._hashWidth);

    if (c < d) {
      fprintf(stderr, "ERROR:  Table[%lu] out of order.\n", i);
      fail++;
    }

    d = c;

    j += mcd._hashWidth;
    i++;
  }

  assert(fail == 0);
}


void
fillCheckTable(char *inputFile, uint32 skipNum) {
  uint64  mer, b, c;

  //  Allocate space for mcd._actualNumberOfMers mers in the _chck array.
  //  This doesn't need to be cleared.
  //
  _chck = new uint64 [mcd._actualNumberOfMers * mcd._chckBits / 64 + 1];

  merStream   M(mcd._merSizeInBases, inputFile, skipNum);
  uint64      nummers=0;

  while (M.nextMer() && (nummers++ < mersToCount)) {
    mer = M.theFMer();
    if (mer > M.theRMer())
      mer = M.theRMer();

    b = mcd.HASH(mer) * mcd._hashWidth;
    c = preDecrementDecodedValue(_hash, b, mcd._hashWidth) * mcd._chckBits;

    setDecodedValue(_chck, c, mcd._chckBits, mer & mcd._chckMask);
  }
}



////////////////////////////////////////////////////////////////////////////////
//
//  OUTPUT
//
void
adjustHeap(uint64 *M, int64 i, int64 n) {
  uint64   m = M[i];
  int64    j = (i << 1) + 1;  //  let j be the left child

  while (j < n) {
    if (j<n-1 && M[j] < M[j+1])
      j++;                   //  j is the larger child

    if (m >= M[j])           //  a position for M[i] has been found
      break;

    M[(j-1)/2] = M[j];       //  Move larger child up a level

    j = (j << 1) + 1;
  }

  M[(j-1)/2] = m;
}



void
sortAndOutput(char   *outfilename,
              uint32  targetCount) {
  uint64 m     = UINT64_ONE << mcd._tableSizeInBits;
  uint32 count = 0;
  uint32 items = 0;

  uint64  *_sortedList    = 0L;
  uint32   _sortedListMax = 0;
  uint32   _sortedListLen = 0;

  unsigned char  theMerString[33];
  uint64         mer;

  //  Open the output file
  //
  FILE *outFile = fopen(outfilename, "w");


  //  For each bucket, sort it.  The output is done
  //  in the sort.
  //
  for (uint64 B=0, b=0; b<m; b++) {
    uint64 st = getDecodedValue(_hash, B, mcd._hashWidth);
    B        += mcd._hashWidth;
    uint64 ed = getDecodedValue(_hash, B, mcd._hashWidth);

    assert(ed >= st);

    _sortedListLen = ed - st;

    count = 0;
    items = 0;

    if (_sortedListLen > 0) {

      //  Allocate more space, if we need to.
      //
      if (_sortedListLen > _sortedListMax) {
        delete [] _sortedList;
        _sortedList    = new uint64 [_sortedListLen + 1];
        _sortedListMax = _sortedListLen;
      }

      //  Unpack the check values
      //
      for (uint64 i=st, J=st*mcd._chckBits; i<ed; i++, J += mcd._chckBits)
        _sortedList[i-st] = getDecodedValue(_chck, J, mcd._chckBits);

      //  Sort if there is more than one item
      //
      if (_sortedListLen > 1) {

        //  Create the heap of lines.
        //
        for (int64 t=(_sortedListLen-2)/2; t>=0; t--)
          adjustHeap(_sortedList, t, _sortedListLen);

        //  Interchange the new maximum with the element at the end of the tree
        //
        for (int64 t=_sortedListLen-1; t>0; t--) {
          uint64           tv = _sortedList[t];
          _sortedList[t]      = _sortedList[0];
          _sortedList[0]      = tv;

          adjustHeap(_sortedList, 0, t);
        }
      }


      //  Scan the list of sorted mers, counting them.  Whenever we 
      //  know the count, output it.
      //
      count = 1;
      if (_sortedListLen > 0) {
        uint32 t;
        for (t=1; t<_sortedListLen; t++) {
          if (_sortedList[t] != _sortedList[t-1]) {
            if (targetCount <= count) {
              mer = (b << mcd._chckBits) | _sortedList[t-1];
              for (uint32 i=0; i<mcd._merSizeInBases; i++)
                theMerString[mcd._merSizeInBases-i-1] = decompressSymbol[(mer >> (2*i)) & 0x03];
              theMerString[mcd._merSizeInBases] = 0;
              fprintf(outFile, ">%d\n%s\n", count, theMerString);
            }
            count = 0;
          }

          count++;
        }

        //  Dump the last mer
        //
        if (targetCount <= count) {
          mer = (b << mcd._chckBits) | _sortedList[t-1];
          for (uint32 i=0; i<mcd._merSizeInBases; i++)
            theMerString[mcd._merSizeInBases-i-1] = decompressSymbol[(mer >> (2*i)) & 0x03];
          theMerString[mcd._merSizeInBases] = 0;
          fprintf(outFile, ">%d\n%s\n", count, theMerString);
        }
      }
    }
  }

  fclose(outFile);
}



int
main(int argc, char **argv) {
  uint32            merSize          = 0;
  char             *fragStore        = 0L;
  char             *outputFile       = 0L;
  uint64            minimumCount     = 0;
  uint32            skipNum = 1;

  if (argc == 1) {
    fprintf(stderr, "usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "Given a sequence file (-s) and lots of parameters, compute\n");
    fprintf(stderr, "the mer-count tables.  By default, both strands are processed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "        -m #                    (size of a mer)\n");
    fprintf(stderr, "        -s /path/to/fragstore   (a fragstore)\n");
    fprintf(stderr, "        -n #                    (threshold to output)\n");
    fprintf(stderr, "        -N #                    (number of mers to examine)\n");
    fprintf(stderr, "        -K #                    (read only every Kth ffragment)\n");
    fprintf(stderr, "        -o tblprefix            (output table prefix)\n");
    exit(1);
  }

  mersToCount = ~UINT64_ZERO;

  for (int arg=1; arg < argc; arg++) {
    if (argv[arg][0] != '-') {
      fprintf(stderr, "Not an option: '%s'.\n", argv[arg]);
      exit(1);
    } else {
      switch (argv[arg][1]) {
        case 'm':
          arg++;
          merSize = atoi(argv[arg]);
          break;
        case 's':
          arg++;
          fragStore = argv[arg];
          break;
        case 'n':
          arg++;
          minimumCount = strtoull(argv[arg], NULL, 10);
          break;
        case 'K':
	  arg++;
	  skipNum = atoi(argv[arg]);
	  break;
        case 'N':
	  arg++;
          mersToCount = strtoull(argv[arg], NULL, 10);
 	  if(mersToCount <= 0){
 	    fprintf(stderr,"Trouble getting number of mers to count from %s (option -N)\n",
 		    argv[arg]);
 	    exit(1);
 	  }
	  break;
        case 'o':
          arg++;
          outputFile = argv[arg];
          break;
        case 'V':
          fprintf(stdout, "version: CA $Id: AS_MER_meryl.cc,v 1.11 2007-11-08 12:38:12 brianwalenz Exp $\n");
          exit(0);
          break;
        default:
          fprintf(stderr, "Unknown option '%s'.\n", argv[arg]);
          break;
      }
    }
  }

  if (fragStore == 0L) {
    fprintf(stderr, "ERROR - no fragstore loaction specified.\n");
    exit(1);
  }

  if (outputFile == 0L) {
    fprintf(stderr, "ERROR - no output file specified.\n");
    exit(1);
  }

  if (merSize == 0) {
    fprintf(stderr, "ERROR - no merSize specified.\n");
    exit(1);
  }

  mcd._merSizeInBases      = merSize;
  mcd._merSizeInBits       = mcd._merSizeInBases << 1;

  //  XXX:  Hardcoded; need to estimate based on the fragstore
  //  estimateTableSize(fragStore, merSize),
  //
  mcd._tableSizeInBits     = 26;
  mcd._tableSizeInEntries  = 1 << mcd._tableSizeInBits;

  mcd._chckBits            = mcd._merSizeInBits - mcd._tableSizeInBits;
  _chck                    = 0L;
  mcd._chckMask            = UINT64_MASK(mcd._chckBits);

  mcd._hashWidth           = 0;
  _hash                    = 0L;
  mcd._hashMask            = UINT64_MASK(mcd._tableSizeInBits);  //  unused?

  mcd._actualNumberOfMers  = 0;

  createHashTable(fragStore,skipNum);
  verifyHashTable();
  fillCheckTable(fragStore,skipNum);
  sortAndOutput(outputFile, minimumCount);

  delete [] _chck;
  delete [] _hash;

  return(0);
}
