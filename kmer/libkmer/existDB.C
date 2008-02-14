#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "bio++.H"


existDB::existDB(char const  *filename,
                 bool         loadData) {
  clear();

  _compressedHash   = false;
  _compressedBucket = false;

  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "existDB::existDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}


existDB::existDB(char const    *filename,
                 u32bit         merSize,
                 existDBflags   flags,
                 u32bit         lo,
                 u32bit         hi) {
  clear();

  _compressedHash   = flags & existDBcompressHash;
  _compressedBucket = flags & existDBcompressBuckets;
  _compressedCounts = flags & existDBcompressCounts;

  //  Try to read state from the filename.  If successful, make sure
  //  that the merSize is correct.
  //
  if (loadState(filename)) {
    bool fail = false;

    if (_merSizeInBases != merSize) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different mer sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", _merSizeInBases, merSize);
      fail = true;
    }

    if (fail)
      exit(1);

    return;
  }

  //  If no direction flags are set, set the default direction of
  //  forward.  Stupid precedence rules.
  //
  if ((flags & (existDBcanonical | existDBforward)) == u32bitZERO)
    flags |= existDBforward;

  //  If we can open 'filename' for reading, then we assume the file
  //  is a multi-fasta, and we build an existDB/
  //
  //  Otherwise, we assume that 'filename' is really the prefix for a
  //  meryl database.


  if (fileExists(filename))
    createFromFastA(filename, merSize, flags);
  else
    createFromMeryl(filename, lo, hi, flags);
}


existDB::~existDB() {
  delete [] _hashTable;
  delete [] _buckets;
  delete [] _counts;
}





bool
existDB::exists(u64bit mer) {
  u64bit c, h, st, ed;

  if (_compressedHash) {
    h  = HASH(mer) * _hshWidth;
    st = getDecodedValue(_hashTable, h,             _hshWidth);
    ed = getDecodedValue(_hashTable, h + _hshWidth, _hshWidth);
  } else {
    h  = HASH(mer);
    st = _hashTable[h];
    ed = _hashTable[h+1];
  }

  if (st == ed)
    return(false);

  c = CHECK(mer);

  if (_compressedBucket) {
    st *= _chkWidth;
    ed *= _chkWidth;

    for (; st<ed; st += _chkWidth) {
      if (getDecodedValue(_buckets, st, _chkWidth) == c)
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


u64bit
existDB::count(u64bit mer) {
  u64bit c, h, st, ed, ct;

  if (_counts == 0L)
    return(0);

  if (_compressedHash) {
    h  = HASH(mer) * _hshWidth;
    st = getDecodedValue(_hashTable, h,             _hshWidth);
    ed = getDecodedValue(_hashTable, h + _hshWidth, _hshWidth);
    ct = st;
  } else {
    h  = HASH(mer);
    st = _hashTable[h];
    ed = _hashTable[h+1];
    ct = st;
  }

  if (st == ed)
    return(false);

  c = CHECK(mer);

  if (_compressedBucket) {
    st *= _chkWidth;
    ed *= _chkWidth;

    for (; st<ed; st += _chkWidth) {
      if (getDecodedValue(_buckets, st, _chkWidth) == c)
        goto returncount;
    }
  } else {
    for (; st<ed; st++) {
      if (_buckets[st] == c)
        goto returncount;
    }
  }

  return(0);

 returncount:
  if (_compressedCounts)
    return(getDecodedValue(_counts, ct * _cntWidth, _cntWidth));
  else
    return(_counts[ct]);
}




#if 0
//  Special case used from inside the posDB.  Errrr, it _was_ used
//  inside the posDB.
//
//  If:
//    The existDB and posDB are built using the same hshWidth.
//    We know the bucket and check we want to find.
//    The checks are in increasing order.
//  Then we can bypass a lot of the overhead of checking for existence.
//
bool
existDB::exists(u64bit b, u64bit c) {

  //  Are we in a new bucket?  Reset!
  //
  if (b != _es_bucket) {
    _es_bucket = b;

    if (_compressedHash) {
      b      *= _hshWidth;
      _es_st  = getDecodedValue(_hashTable, b,             _hshWidth);
      _es_ed  = getDecodedValue(_hashTable, b + _hshWidth, _hshWidth);
    } else {
      _es_st  = _hashTable[b];
      _es_ed  = _hashTable[b+1];
    }

    if (_compressedBucket) {
      _es_st *= _chkWidth;
      _es_ed *= _chkWidth;
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
    for (; _es_st<_es_ed; _es_st += _chkWidth) {
      if (getDecodedValue(_buckets, _es_st, _chkWidth) == c)
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
    for (u64bit st=_es_st; st<_es_ed; st += _chkWidth) {
      if (getDecodedValue(_buckets, st, _chkWidth) == c)
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
