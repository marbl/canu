
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



class swv;   //  Defined and only used in merylCountArray.C


//  The basic meryl counting object.  It collects a list of kmers from the
//  input, sorts them, the returns a list of distinct kmers with the number
//  of times each occurs in the list.
//
//  A kmer is split into a prefix and a suffix.  All kmers in one
//  merylCountArray object will share the same prefix (the client needs to do
//  this for us).  The merylCountArray will save the lower _sWidth bits of
//  each input kmer.
//
//  If the kmer comes with a count or color, those are saved too.  Otherwise,
//  only the kmer suffix bits are saved.
//
//  Since we don't know how many kmers will be input to us, and don't want to
//  keep reallocating data, our list is grown in segments of size _segSize
//  (segsize in initialize()) kilo bytes.

#undef ADD_VERBOSE
#define ADD_INSTRUMENT

class merylCountArray {
public:
  merylCountArray();
  ~merylCountArray();

  uint64    initialize(uint64 prefix, uint32 width, uint32 segsize=64);

  uint64    initializeValues(kmvalu maxValue=0);

  void      enableMultiSet(bool enable)    {  _multiSet = enable;  };

public:
  void      initializeForTesting(uint32 width, uint32 nwords);
  void      dumpData(void);

private:
  void      removeSegments(void);
  void      removeValues(void);
  void      addSegment(uint32 seg);

  //  Count the number of times we hit each case, and where the word starts in that case
public:
  void      clearStats(void);
  void      dumpStats(void);

#ifdef ADD_INSTRUMENT
  uint64    nTests[6];
  uint64    nStart[6][64];
#endif

  //  Add a suffix to the table.
  //
public:
  uint64    add(kmdata suffix);
  uint64    addValue(kmvalu value);

private:
  kmdata      *unpackSuffixes(uint64 nSuffixes);
  swv         *unpackSuffixesAndValues(uint64 nSuffixes);

  //
  //  Return the kkth kmer suffix stored in the array.  This is only used in sort(),
  //  and only to convert the bit-packed _seg data into unpacked words, so could
  //  be optimized for that case.  I don't expect much of a performance gain.
  //
public:
  kmdata    get(uint64 kk);
  kmdata    getSimple(uint64 kk);


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
  void             dumpCountedKmers(merylBlockWriter *out);
  void             removeCountedKmers(void);


private:
  uint32           _sWidth;       //  Size of the suffix we're storing
  uint32           _vWidth;       //  Size of the values we're storing

  uint64           _prefix;       //  The kmer prefix we're storing data for
  kmdata          *_suffix;       //  After sorting, the suffix of each kmer
  kmvalu          *_counts;       //  After sorting, the number of times we've seen this kmer

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
