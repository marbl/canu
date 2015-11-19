
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
 *  This file is derived from:
 *
 *    kmer/libutil/bitPackedArray.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2005-JUL-12 to 2014-APR-11
 *      are Copyright 2005-2006,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef BITPACKEDARRAY_H
#define BITPACKEDARRAY_H

#include "AS_global.H"

#undef  DEBUG_BPH_ADD
#undef  DEBUG_BPH_GET

////////////////////////////////////////
//
//  bitPackedArray
//
//  implements an integer array using bit-widths less than word-sizes,
//  e.g., a memory efficient way to store 23 bit numbers.  Numbers may
//  be up to 64 bits wide.
//
//  The array is variable length, and it is implemented as an array,
//  not a list or tree -- accessing element 1,000,000 will allocate
//  elements 0 through 999,999.
//
class bitPackedArray {
public:

  //  Create a bitpacked array with elements of width 'width' using
  //  'segmentSize' KB per segment.  If you know your array is going
  //  to be much bigger or smaller, crank this value.
  //
  bitPackedArray(uint32 valueWidth, uint32 segmentSize = 1024);
  ~bitPackedArray();

  //  No array operator is provided, because we cannot return a
  //  reference to a value that is split across two words (or even a
  //  reference to a value that is not bit aligned in the word).
  //
  uint64   get(uint64 idx);
  void     set(uint64 idx, uint64 val);

  //  Clear the array.  Since the array is variable sized, you must add
  //  things to a new array before clearing it.
  void     clear(void);

private:
  uint32   _valueWidth;
  uint32   _segmentSize;
  uint64   _nextElement;  //  the first invalid element
  uint64   _valuesPerSegment;

  uint64   _numSegments;
  uint64   _maxSegments;
  uint64 **_segments;
};


//  An array of bits.  Exactly the same as the bitPackedArray, but
//  optimized for width=1.
//
class bitArray {
public:

  bitArray(uint32 segmentSize = 1024);
  ~bitArray();

  uint64   get(uint64 idx);

  uint64   getAndSet(uint64 idx);

  void     set(uint64 idx);
  void     clr(uint64 idx);

  void     clear(void);

private:
  void     resize(uint64 s);

  uint32   _segmentSize;
  uint64   _valuesPerSegment;

  uint64   _numSegments;
  uint64   _maxSegments;
  uint64 **_segments;
};


//  Uses the bitPackedArray to implement a heap.  The bitPackedArray is dynamically sized,
//  so this can be too.
//
class bitPackedHeap {
public:
  bitPackedHeap(uint32 width, uint64 size=16) {
    _array    = new bitPackedArray(width, size);
    _array->set(0, 0);
    _lastVal  = 0;
  };

  ~bitPackedHeap() {
    delete _array;
  };

  uint64    get(void) {
    uint64  biggestVal = ~uint64ZERO;

    if (_lastVal == 0)
      return(biggestVal);

    biggestVal = _array->get(0);
    _lastVal--;

    if (_lastVal == 0)
      return(biggestVal);

    uint64  t    = _array->get(_lastVal);

    _array->set(0, t);

    uint64  pidx = 0;
    uint64  pval = t;
    uint64  cidx = 1;
    uint64  cval = 0;  //  set below

    while (cidx < _lastVal) {
      //  Set cval here, so we can first test if cidx is in range.
      cval = _array->get(cidx);

      //  Pick the smallest of the two kids
      if (cidx+1 < _lastVal) {
        t = _array->get(cidx+1);
        if (cval > t) {
          cidx++;
          cval = t;
        }
      }

#ifdef DEBUG_BPH_GET
      fprintf(stderr, "test c="F_U64" and p="F_U64" lastVal="F_U64"\n",
              cidx, pidx, _lastVal);
      fprintf(stderr, "test c="F_U64"="F_U64"\n",
              cidx, cval);
      fprintf(stderr, "test p="F_U64"="F_U64"\n",
              pidx, pval);
      fprintf(stderr, "test c="F_U64"="F_U64" and p="F_U64"="F_U64"\n",
              cidx, cval, pidx, pval);
#endif

      if (cval < pval) {

#ifdef DEBUG_BPH_GET
        fprintf(stderr, "swap c="F_U64"="F_U64" and p="F_U64"="F_U64"\n",
                cidx, cval, pidx, pval);
#endif

        //  Swap p and c
        _array->set(pidx, cval);
        _array->set(cidx, pval);

        //  Move down the tree -- pval doesn't change, we moved it into cidx!
        pidx = cidx;
        cidx = cidx * 2 + 1;
      } else {
        cidx = _lastVal;
      }
    }

    return(biggestVal);
  };

  void      add(uint64 value) {
    uint64  cidx = _lastVal;
    uint64  cval = value;
    uint64  pidx = 0;
    uint64  pval = 0;
    bool    more = false;

#ifdef DEBUG_BPH_ADD
    fprintf(stderr, "add  c="F_U64"="F_U64" -- lastVal="F_U64"\n",
            cidx, cval, _lastVal);
#endif

    _array->set(cidx, cval);

    if (cidx > 0)
      more = true;

    while (more) {
      pidx = (cidx-1) / 2;

#ifdef DEBUG_BPH_ADD
      fprintf(stderr, "more c="F_U64" and p="F_U64"\n", cidx, pidx);
#endif

      pval = _array->get(pidx);

#ifdef DEBUG_BPH_ADD
      fprintf(stderr, "test c="F_U64"="F_U64" and p="F_U64"="F_U64"\n",
              cidx, cval, pidx, pval);
#endif

      if (pval > cval) {

#ifdef DEBUG_BPH_ADD
        fprintf(stderr, "swap c="F_U64"="F_U64" and p="F_U64"="F_U64"\n",
                cidx, cval, pidx, pval);
#endif

        //  Swap p and c
        _array->set(cidx, pval);
        _array->set(pidx, cval);

        //  Move up the tree -- cval doesn't change, we moved it into pidx!
        cidx = pidx;
      } else {
        more = false;
      }
      if (cidx == 0)
        more = false;
    }

    _lastVal++;

    //dump();
  };

  void      dump(void) {
    for (uint32 i=0; i<_lastVal; i++)
      fprintf(stderr, "HEAP["F_U32"]="F_U64"\n", i, _array->get(i));
  }

  void      clear(void) {
    _array->clear();
    _lastVal = 0;
  };

private:
  bitPackedArray   *_array;
  uint64            _lastVal;
};



inline
uint64
bitArray::get(uint64 idx) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = idx % _valuesPerSegment;

  uint64 wrd = (p >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (p     ) & 0x000000000000003fllu;

  return((_segments[s][wrd] >> bit) & 0x0000000000000001llu);
}


inline
void
bitArray::resize(uint64 s) {

  if (s < _numSegments)
    return;

  if (s > _maxSegments) {
    _maxSegments = s + 16;
    uint64 **S = new uint64 * [_maxSegments];
    for (uint32 i=0; i<_numSegments; i++)
      S[i] = _segments[i];
    delete [] _segments;
    _segments = S;
  }

  while (_numSegments <= s)
    _segments[_numSegments++] = new uint64 [_segmentSize * 1024 / 8];
}


inline
uint64
bitArray::getAndSet(uint64 idx) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = idx % _valuesPerSegment;

  uint64 wrd = (p >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (p     ) & 0x000000000000003fllu;

  uint64 ret = (_segments[s][wrd] >> bit) & 0x0000000000000001llu;

  _segments[s][wrd] |= uint64ONE << bit;

  return(ret);
}


inline
void
bitArray::set(uint64 idx) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = idx % _valuesPerSegment;

  resize(s);

  uint64 wrd = (p >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] |= uint64ONE << bit;
}


inline
void
bitArray::clr(uint64 idx) {
  uint64 s = idx / _valuesPerSegment;
  uint64 p = idx % _valuesPerSegment;

  resize(s);

  uint64 wrd = (p >> 6) & 0x0000cfffffffffffllu;
  uint64 bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] &= ~(0x0000000000000001llu << bit);
}


#endif  // BITPACKEDARRAY_H
