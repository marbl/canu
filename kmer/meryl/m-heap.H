#ifndef M_HEAP_H
#define M_HEAP_H

#include "util++.H"
#include "bio++.H"

//
//  This is a bit packed heap, derived from bitPackedHeap.  It is
//  customized to decode a kmer from a merStream, given the location
//  of the kmer in the stream.  This kmer is used for the value of the
//  data item in the heap, instead of the value stored in the heap.
//

class bitPackedMerHeap {
public:
  bitPackedMerHeap(seqStore *SS, u32bit width, u64bit size=16) {
    _array    = new bitPackedArray(width, size);
    _array->set(0, 0);
    _lastVal  = 0;
    _mers     = SS;
  };

  ~bitPackedMerHeap() {
    delete _array;
  };

  //  Get the mer with index idx in the merStream
  //
  kMer const  &getMer(u64bit idx) {
    _mers->setIterationStart(idx);
    _mers->nextMer();
    if (_mers->theRMer() < _mers->theFMer())
      return(_mers->theRMer());
    return(_mers->theFMer());
  }

  u64bit       get(kMer &mer) {
    u64bit  pos = ~u64bitZERO;

    if (_lastVal == 0)
      return(pos);

    pos = _array->get(0);
    mer = getMer(pos);

    if (--_lastVal == 0)
      return(pos);

    //  Rebalance the heap

    u64bit  tval = _array->get(_lastVal);
    kMer    tmer;

    _array->set(0, tval);

    u64bit  pidx = 0;
    u64bit  pval = tval;
    kMer    pmer = getMer(pval);
    u64bit  cidx = 1;
    u64bit  cval = 0;  //  set below
    kMer    cmer;

    while (cidx < _lastVal) {
      //  Set cval here, so we can first test if cidx is in range.
      cval = _array->get(cidx);
      cmer = getMer(cval);

      //  Pick the smallest of the two kids
      if (cidx+1 < _lastVal) {
        tval = _array->get(cidx+1);
        tmer = getMer(tval);

        if (cmer > tmer) {
          cidx++;
          cval = tval;
          cmer = tmer;
        }
      }

      if (cmer < pmer) {

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

    return(pos);
  };

  void      add(u64bit value) {
    u64bit  cidx = _lastVal++;
    u64bit  cval = value;
    kMer    cmer;
    u64bit  pidx = 0;
    u64bit  pval = 0;
    kMer    pmer;
    bool    more = true;

    _array->set(cidx, cval);

    if (cidx == 0)
      return;

    cmer = getMer(cval);

    while (more) {
      pidx = (cidx-1) / 2;
      pval = _array->get(pidx);
      pmer = getMer(pval);

      if (pmer > cmer) {

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
  bitPackedArray       *_array;
  u64bit                _lastVal;
  seqStore             *_mers;
};


#endif  //  M_HEAP_H
