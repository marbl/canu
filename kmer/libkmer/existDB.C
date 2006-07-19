#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "positionDB.H"
#include "bio++.H"


existDB::existDB(char const  *filename,
                 bool         loadData) {
  clear();

  //  XXX  we should extend this function to take flags like the other constructor.

  _compressedHash   = false;
  _compressedBucket = false;
  _beVerbose        = false;

  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "existDB::existDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}


existDB::existDB(char const    *filename,
                 u32bit         merSize,
                 u32bit         tblBits,
                 u32bit         lo,
                 u32bit         hi,
                 existDBflags   flags) {
  clear();

  _compressedHash   = flags & existDBcompressHash;
  _compressedBucket = flags & existDBcompressBuckets;
  _beVerbose        = flags & existDBverbose;


  //  Try to read state from the filename.  If successful, make sure that
  //  the merSize and tblBits are correct.
  //
  if (loadState(filename)) {
    bool fail = false;

    if (_merSizeInBases != merSize) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different mer sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", _merSizeInBases, merSize);
      fail = true;
    }
    if (_mask1 != u64bitMASK(tblBits)) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different table sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", 2 * _merSizeInBases - _shift1, tblBits);
      fail = true;
    }

    if (fail)
      exit(1);

    return;
  }

  //  If no direction flags are set, set the default direction of
  //  forward.  Stupid precedence rules.
  //
  if ((flags & (existDBcanonical | existDBforward | existDBreverse)) == u32bitZERO)
    flags |= existDBforward;

  //  If we can open 'filename' for reading, then we assume the file is a multi-fasta, and
  //  we build an existDB using merSize = p1 and tblBits = p2
  //
  //  Otherwise, we assume that 'filename' is really the prefix for a meryl database.


  if (fileExists(filename))
    createFromFastA(filename, merSize, tblBits, flags);
  else
    createFromMeryl(filename, lo, hi, tblBits, flags);
}


existDB::~existDB() {
  delete [] _hashTable;
  delete [] _buckets;
}





bool
existDB::exists(u64bit mer) {
  u64bit h, st, ed;

  if (_compressedHash) {
    h  = HASH(mer) * _hashWidth;
    st = getDecodedValue(_hashTable, h,              _hashWidth);
    ed = getDecodedValue(_hashTable, h + _hashWidth, _hashWidth);
  } else {
    h  = HASH(mer);
    st = _hashTable[h];
    ed = _hashTable[h+1];
  }

  if (st == ed)
    return(false);

  u64bit  c = CHECK(mer);

  if (_compressedBucket) {
    st *= _chckWidth;
    ed *= _chckWidth;

    for (; st<ed; st += _chckWidth) {
      if (getDecodedValue(_buckets, st, _chckWidth) == c)
        return(true);
    }
  } else {
    for (; st<ed; st++) {
      if (_buckets[st] == c)
        return(true);
    }
  }

  return(false);
}


//  Special case used from inside the posDB.  Errrr, it _was_ used
//  inside the posDB.
//
//  If:
//    The existDB and posDB are built using the same hashWidth.
//    We know the bucket and check we want to find.
//    The checks are in increasing order.
//  Then we can bypass a lot of the overhead of checking for existence.
//
#if 0
bool
existDB::exists(u64bit b, u64bit c) {

  //  Are we in a new bucket?  Reset!
  //
  if (b != _es_bucket) {
    _es_bucket = b;

    if (_compressedHash) {
      b      *= _hashWidth;
      _es_st  = getDecodedValue(_hashTable, b,              _hashWidth);
      _es_ed  = getDecodedValue(_hashTable, b + _hashWidth, _hashWidth);
    } else {
      _es_st  = _hashTable[b];
      _es_ed  = _hashTable[b+1];
    }

    if (_compressedBucket) {
      _es_st *= _chckWidth;
      _es_ed *= _chckWidth;
    }
  }

  //  If we're all done, return early.  Probably doesn't do much, just
  //  skips the setup of a for loop.
  //
  if (_es_st == _es_ed)
    return(false);

#ifdef ES_SORTED
  //  This would work great, except existDB isn't sorted.
  if (_compressedBucket) {
    for (; _es_st<_es_ed; _es_st += _chckWidth) {
      if (getDecodedValue(_buckets, _es_st, _chckWidth) == c)
        return(true);
    }
  } else {
    for (; _es_st<_es_ed; _es_st++) {
      if (_buckets[_es_st] == c)
        return(true);
    }
  }
#else
  if (_compressedBucket) {
    for (u64bit st=_es_st; st<_es_ed; st += _chckWidth) {
      if (getDecodedValue(_buckets, st, _chckWidth) == c)
        return(true);
    }
  } else {
    for (u64bit st=_es_st; st<_es_ed; st++) {
      if (_buckets[st] == c)
        return(true);
    }
  }
#endif

  return(false);
}
#endif
