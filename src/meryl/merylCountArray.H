
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef MERYLCOUNTARRAY_H
#define MERYLCOUNTARRAY_H

#include "meryl.H"
#include "bits.H"
#include "system.H"




template<typename VALUE>
class swv {   //  That's suffix-with-value
public:
  void      set(uint64 suffix, VALUE value) {
    _s1 = suffix >> 32;
    _s2 = suffix  & uint64MASK(32);
    _v  = value;
  };

  uint64    getSuffix(void) {
    uint64   s;

    s   = _s1;
    s <<= 32;
    s  |= _s2;

    return(s);
  };

  VALUE     getValue(void) {
    return(_v);
  };

  bool      operator<(swv<VALUE> const that) const {
    return(((_s1  < that._s1)) ||
           ((_s1 == that._s1) && (_s2  < that._s2)) ||
           ((_s1 == that._s1) && (_s2 == that._s2) && (_v < that._v)));
  };

private:
  uint32   _s1;  //  This bit of ugly is, splitting the suffix into
  uint32   _s2;  //  two words, allows an swv to be aligned to
  VALUE    _v;   //  4-byte boundaries (if VALUE is uint32).
};





template<typename VALUE>
class merylCountArray {
public:
  merylCountArray();
  ~merylCountArray();

  uint64    initialize(uint64 prefix, uint32 width, uint32 segsize=64);

  uint64    initializeValues(VALUE maxValue=0);

  void      enableMultiSet(bool enable)    {  _multiSet = enable;  };

private:
  void      removeSegments(void);
  void      removeValues(void);
  void      addSegment(uint32 seg);

  //  Add a suffix to the table.
  //
  //  wordPos is 0 for the high bits and 63 for the bit that represents integer 1.
public:
  uint64    add(uint64 suffix) {
    uint64  seg       = _nBits / _segSize;   //  Which segment are we in?
    uint64  segPos    = _nBits % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;         //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _sWidth;

    //  Increment the position.

    _nBits += _sWidth;

    //  If the first word and the first position, we need to allocate a segment.
    //  This catches both the case when _nBits=0 (we've added nothing) and when
    //  _nBits=_segSize (we've added exactly one segment worth of kmers).

    if ((word    == 0) &&
        (wordBgn == 0))
      addSegment(seg);

    //  If there is enough space in the current word, stash the bits there.
    //  If it's a new word, we need to special case an initialization of the word.

    if      (wordBgn == 0) {
      _segments[seg][word]  = (suffix << (64 - wordEnd));
    }

    else if (wordEnd <= 64) {
      _segments[seg][word] |= (suffix << (64 - wordEnd));
    }

    //  Otherwise, the suffix spans two words.  If these can be in the same block,
    //  stash the bits there.

    else if (segPos + _sWidth <= _segSize) {
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      _segments[seg][word+0] |= (suffix >>        extraBits);
      _segments[seg][word+1]  = (suffix << (64 -  extraBits));
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

      _segments[seg+0][W] |= (suffix >>        extraBits);
      _segments[seg+1][0]  = (suffix << (64 -  extraBits));
    }

    return(usedSizeDelta());
  };

  uint64    addValue(VALUE value) {

    if (_vals == NULL)
      return(0);

    if (_vWidth == 0)
      return(_vals->setEliasDelta(value));
    else
      return(_vals->setBinary(_vWidth, value));
  };

private:
  uint64      *unpackSuffixes(uint64 nSuffixes);
  swv<VALUE>  *unpackSuffixesAndValues(uint64 nSuffixes);

private:
  //
  //  Return the kkth kmer suffix stored in the array.  This is only used in sort(),
  //  and only to convert the bit-packed _seg data into unpacked words, so could
  //  be optimized for that case.  I don't expect much of a performance gain.
  //
  uint64    get(uint64 kk) {
    uint64  bitPos    = kk * _sWidth;

    uint64  seg       = bitPos / _segSize;   //  Which segment are we in?
    uint64  segPos    = bitPos % _segSize;   //  Bit position in that segment.

    uint32  word      = segPos / 64;         //  Which word are we in=?
    uint32  wordBgn   = segPos % 64;         //  Bit position in that word.
    uint32  wordEnd   = wordBgn + _sWidth;

    uint64  bits      = 0;

    //  If the bits are entirely in a single word, copy them out.

    if      (wordEnd <= 64) {
      bits = (_segments[seg][word] >> (64 - wordEnd)) & uint64MASK(_sWidth);
    }

    //  Otherwise, the suffix spans two words.  If these are in the same block,
    //  grab them.

    else if (segPos + _sWidth <= _segSize) {
      uint32   extraBits = wordEnd - 64;

      assert(wordEnd > 64);

      bits  = (_segments[seg][word+0] & uint64MASK(_sWidth - extraBits)) << extraBits;
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

      bits  = (_segments[seg+0][W] & uint64MASK(_sWidth - extraBits)) << extraBits;
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


private:
  void             countSingleKmers(void);
  void             countSingleKmersWithValues(void);
  void             countMultiSetKmers(void);
public:
  void             countKmers(void);
  void             dumpCountedKmers(kmerCountBlockWriter *out);
  void             removeCountedKmers(void);


private:
  uint32           _sWidth;       //  Size of the suffix we're storing
  uint32           _vWidth;       //  Size of the values we're storing

  uint64           _prefix;       //  The kmer prefix we're storing data for
  uint64          *_suffix;       //  After sorting, the suffix of each kmer
  VALUE           *_counts;       //  After sorting, the number of times we've seen this kmer

  uint64           _nKmers;       //  Number of kmers.

  uint64           _bitsPerPage;
  uint64           _nReAlloc;

  uint32           _segSize;      //  Number of bits in each segment.
  uint32           _segAlloc;     //  Number of segments we're allowed to allocate  (size of the array below).
  uint64         **_segments;     //  An array of blocks of data.

  stuffedBits     *_vals;         //  An array of compressed values.

  uint64           _nBits;        //  Number of bits stored.
  uint64           _nBitsTrigger; //  Number of bits we need to store for a size recalculation to occur.
  uint64           _nBitsOldSize; //  Last computed size.

  bool             _multiSet;     //  Treat the input kmers as a multiset.
};



#endif  //  MERYLCOUNTARRAY_H
