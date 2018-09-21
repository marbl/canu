
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#ifndef MERYLCOUNTARRAY_H
#define MERYLCOUNTARRAY_H

#include "meryl.H"
#include "bits.H"
#include "system.H"

class merylCountArray {
public:
  merylCountArray();
  ~merylCountArray();

  uint64    initialize(uint64 prefix, uint32 width, uint32 segsize);

private:
  void      removeSegments(void);
  void      addSegment(uint32 seg);

  //  Add a value to the table.
  //
  //  wordPos is 0 for the high bits and 63 for the bit that represents integer 1.
public:
  uint64    add(uint64 value) {
    uint64  seg       = _nBits / _segSize;   //  Which segment are we in?
    uint64  segPos    = _nBits % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;         //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _width;

    //  Increment the position.

    _nBits += _width;

    //  If the first word and the first position, we need to allocate a segment.
    //  This catches both the case when _nBits=0 (we've added nothing) and when
    //  _nBits=_segSize (we've added exactly one segment worth of kmers).

    if ((word    == 0) &&
        (wordBgn == 0))
      addSegment(seg);

    //  If there is enough space in the current word, stash the bits there.
    //  If it's a new word, we need to special case an initialization of the word.

    if      (wordBgn == 0) {
      _segments[seg][word]  = (value << (64 - wordEnd));
    }

    else if (wordEnd <= 64) {
      _segments[seg][word] |= (value << (64 - wordEnd));
    }

    //  Otherwise, the value spans two words.  If these can be in the same block,
    //  stash the bits there.

    else if (segPos + _width <= _segSize) {
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      _segments[seg][word+0] |= (value >>        extraBits);
      _segments[seg][word+1]  = (value << (64 -  extraBits));
    }

    //  Otherwise, the word spans two segments.  At least we know it
    //  hits the last word in the first segment and the first word in
    //  the second segment.  And that the second segment hasn't been
    //  allocated yet.  And that it's just the same math as the previous
    //  case, just in different segments instead of words.

    else {
      uint32 W         = word;  //  Just to keep things pretty.  I love my optimizer!
      uint32 extraBits = wordEnd - 64;

      addSegment(seg+1);

      _segments[seg+0][W] |= (value >>        extraBits);
      _segments[seg+1][0]  = (value << (64 -  extraBits));
    }

    return(usedSizeDelta());
  };


private:
  //
  //  Return the kkth kmer suffix stored in the array.  This is only used in sort(),
  //  and only to convert the bit-packed _seg data into unpacked words, so could
  //  be optimized for that case.  I don't expect much of a performance gain.
  //
  uint64    get(uint64 kk) {
    uint64  bitPos    = kk * _width;

    uint64  seg       = bitPos / _segSize;   //  Which segment are we in?
    uint64  segPos    = bitPos % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;         //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _width;

    uint64  bits      = 0;

    //  If the bits are entirely in a single word, copy them out.

    if      (wordEnd <= 64) {
      bits = (_segments[seg][word] >> (64 - wordEnd)) & uint64MASK(_width);
    }

    //  Otherwise, the value spans two words.  If these are in the same block,
    //  grab them.

    else if (segPos + _width <= _segSize) {
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      bits  = (_segments[seg][word+0] & uint64MASK(_width - extraBits)) << extraBits;
      bits |= (_segments[seg][word+1] >> (64 - extraBits) & uint64MASK(extraBits));
    }

    //  Otherwise, the word spans two segments.  At least we know it
    //  hits the last word in the first segment and the first word in
    //  the second segment.  And that the second segment hasn't been
    //  allocated yet.  And that it's just the same math as the previous
    //  case, just in different segments instead of words.

    else {
      uint32 W         = word;  //  Just to keep things pretty.  I love my optimizer!
      uint32 extraBits = wordEnd - 64;

      bits  = (_segments[seg+0][W] & uint64MASK(_width - extraBits)) << extraBits;
      bits |= (_segments[seg+1][0] >> (64 - extraBits) & uint64MASK(extraBits));
    }

    return(bits);
  };


public:
  uint64           numBits(void)        {  return(_nBits);  };

  //  Returns the number of bytes in pages touched by data in this object.
  //  It's a pretty tight bound.  The extra 1 was added to make it
  //  be an overestimate of what 'top' is reporting.  Without it,
  //  it underestimates by a significant amount.

  //  This is underestimating the actual resident memory usage.  The constant
  //  +5 tries to adjust, but still underestimates on large data sets.

  uint64           usedSize(void) {
    uint64  fullSegs  = (_nBits / _segSize);        //  Number of fully filled segments
    uint64  partSeg   = (_nBits % _segSize) + 64;   //  Number of bites (rounded to next word) in the last (partially filled) segment

    uint64  pagesUsed = 0;

    pagesUsed += fullSegs * (_segSize / _bitsPerPage) + fullSegs * (((_segSize % _bitsPerPage) == 0) ? 0 : 1);
    pagesUsed +=            (partSeg  / _bitsPerPage) +            (((partSeg  % _bitsPerPage) == 0) ? 0 : 1);
    pagesUsed += 5;
    pagesUsed += _nReAlloc;

    return(pagesUsed * _bitsPerPage / 8 + _segAlloc * sizeof(uint64 **) + sizeof(merylCountArray));
  };

  uint64           usedSizeDelta(void) {

    if (_nBits < _nBitsTrigger)
      return(0);

    uint64  newSize   = usedSize();
    uint64  sizeDelta = newSize - _nBitsOldSize;

    _nBitsTrigger += _bitsPerPage / 16;
    _nBitsOldSize  = newSize;

    return(sizeDelta);
  };


  void             countKmers(void);
  void             dumpCountedKmers(kmerCountBlockWriter *out);
  void             removeCountedKmers(void);


private:
  uint32           _width;        //  Size of the element we're storing

  uint64           _prefix;       //  The kmer prefix we're storing data for
  uint64          *_suffix;       //  After sorting, the suffix of each kmer
  uint32          *_counts;       //  After sorting, the number of times we've seen this kmer

  uint64           _nKmers;       //  Number of kmers.

  uint64           _bitsPerPage;
  uint64           _nReAlloc;

  uint32           _segSize;      //  Number of bits in each segment.
  uint32           _segAlloc;     //  Number of segments we're allowed to allocate  (size of the array below).
  uint64         **_segments;     //  An array of blocks of data.

  uint64           _nBits;        //  Number of bits stored.
  uint64           _nBitsTrigger; //  Number of bits we need to store for a size recalculation to occur.
  uint64           _nBitsOldSize; //  Last computed size.
};



#endif  //  MERYLCOUNTARRAY_H
