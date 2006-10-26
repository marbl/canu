#ifndef BITPACKEDARRAY_H
#define BITPACKEDARRAY_H

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
  bitPackedArray(u32bit valueWidth, u32bit segmentSize = 1024);
  ~bitPackedArray();

  //  No array operator is provided, because we cannot return a
  //  reference to a value that is split across two words (or even a
  //  reference to a value that is not bit aligned in the word).
  //
  u64bit   get(u64bit idx);
  void     set(u64bit idx, u64bit val);

  //  Clear the array.  Since the array is variable sized, you must add
  //  things to a new array before clearing it.
  void     clear(void);

private:
  u32bit   _valueWidth;
  u32bit   _segmentSize;
  u64bit   _nextElement;  //  the first invalid element
  u64bit   _valuesPerSegment;

  u64bit   _numSegments;
  u64bit   _maxSegments;
  u64bit **_segments;
};


//  An array of bits.  Exactly the same as the bitPackedArray, but
//  optimized for width=1.
//
class bitArray {
public:

  bitArray(u32bit segmentSize = 1024);
  ~bitArray();

  u64bit   get(u64bit idx);

  u64bit   getAndSet(u64bit idx);

  void     set(u64bit idx);
  void     clr(u64bit idx);

  void     clear(void);

private:
  void     resize(u64bit s);

  u32bit   _segmentSize;
  u64bit   _valuesPerSegment;

  u64bit   _numSegments;
  u64bit   _maxSegments;
  u64bit **_segments;
};


//  Uses the bitPackedArray to implement a heap.  The bitPackedArray is dynamically sized,
//  so this can be too.
//
class bitPackedHeap {
public:
  bitPackedHeap(u32bit width, u64bit size=16) {
    _array    = new bitPackedArray(width, size);
    _array->set(0, 0);
    _lastVal  = 0;
  };

  ~bitPackedHeap() {
    delete _array;
  };

  u64bit    get(void) {
    u64bit  biggestVal = ~u64bitZERO;

    if (_lastVal == 0)
      return(biggestVal);

    biggestVal = _array->get(0);
    _lastVal--;

    if (_lastVal == 0)
      return(biggestVal);

    u64bit  t    = _array->get(_lastVal);

    _array->set(0, t);

    u64bit  pidx = 0;
    u64bit  pval = t;
    u64bit  cidx = 1;
    u64bit  cval = 0;  //  set below

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
      fprintf(stderr, "test c="u64bitFMT" and p="u64bitFMT" lastVal="u64bitFMT"\n",
              cidx, pidx, _lastVal);
      fprintf(stderr, "test c="u64bitFMT"="u64bitFMT"\n",
              cidx, cval);
      fprintf(stderr, "test p="u64bitFMT"="u64bitFMT"\n",
              pidx, pval);
      fprintf(stderr, "test c="u64bitFMT"="u64bitFMT" and p="u64bitFMT"="u64bitFMT"\n",
              cidx, cval, pidx, pval);
#endif

      if (cval < pval) {

#ifdef DEBUG_BPH_GET
        fprintf(stderr, "swap c="u64bitFMT"="u64bitFMT" and p="u64bitFMT"="u64bitFMT"\n",
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

  void      add(u64bit value) {
    u64bit  cidx = _lastVal;
    u64bit  cval = value;
    u64bit  pidx = 0;
    u64bit  pval = 0;
    bool    more = false;

#ifdef DEBUG_BPH_ADD
    fprintf(stderr, "add  c="u64bitFMT"="u64bitFMT" -- lastVal="u64bitFMT"\n",
            cidx, cval, _lastVal);
#endif

    _array->set(cidx, cval);

    if (cidx > 0)
      more = true;

    while (more) {
      pidx = (cidx-1) / 2;

#ifdef DEBUG_BPH_ADD
      fprintf(stderr, "more c="u64bitFMT" and p="u64bitFMT"\n", cidx, pidx);
#endif

      pval = _array->get(pidx);

#ifdef DEBUG_BPH_ADD
      fprintf(stderr, "test c="u64bitFMT"="u64bitFMT" and p="u64bitFMT"="u64bitFMT"\n",
              cidx, cval, pidx, pval);
#endif

      if (pval > cval) {

#ifdef DEBUG_BPH_ADD
        fprintf(stderr, "swap c="u64bitFMT"="u64bitFMT" and p="u64bitFMT"="u64bitFMT"\n",
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
    for (u32bit i=0; i<_lastVal; i++)
      fprintf(stderr, "HEAP["u32bitFMT"]="u64bitFMT"\n", i, _array->get(i));
  }

  void      clear(void) {
    _array->clear();
    _lastVal = 0;
  };

private:
  bitPackedArray   *_array;
  u64bit            _lastVal;
};



inline
u64bit
bitArray::get(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  return((_segments[s][wrd] >> bit) & 0x0000000000000001llu);
}


inline
void
bitArray::resize(u64bit s) {

  if (s < _numSegments)
    return;

  if (s > _maxSegments) {
    _maxSegments = s + 16;
    u64bit **S = new u64bit * [_maxSegments];
    for (u32bit i=0; i<_numSegments; i++)
      S[i] = _segments[i];
    delete [] _segments;
    _segments = S;
  }

  while (_numSegments <= s)
    _segments[_numSegments++] = new u64bit [_segmentSize * 1024 / 8];
}


inline
u64bit
bitArray::getAndSet(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  u64bit ret = (_segments[s][wrd] >> bit) & 0x0000000000000001llu;
  
  _segments[s][wrd] |= u64bitONE << bit;

  return(ret);
}


inline
void
bitArray::set(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  resize(s);

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] |= u64bitONE << bit;
}


inline
void
bitArray::clr(u64bit idx) {
  u64bit s = idx / _valuesPerSegment;
  u64bit p = idx % _valuesPerSegment;

  resize(s);

  u64bit wrd = (p >> 6) & 0x0000cfffffffffffllu;
  u64bit bit = (p     ) & 0x000000000000003fllu;

  _segments[s][wrd] &= ~(0x0000000000000001llu << bit);
}


#endif  // BITPACKEDARRAY_H
