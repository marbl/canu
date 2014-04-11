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
                 uint32         merSize,
                 existDBflags   flags,
                 uint32         lo,
                 uint32         hi) {
  clear();

  _compressedHash   = flags & existDBcompressHash;
  _compressedBucket = flags & existDBcompressBuckets;
  _compressedCounts = flags & existDBcompressCounts;

  _searchForDupe = false;

  //  Try to read state from the filename.  If successful, make sure
  //  that the merSize is correct.
  //
  if (loadState(filename)) {
    bool fail = false;

    if (_merSizeInBases != merSize) {
      fprintf(stderr, "existDB::existDB()-- Read state from '%s', but got different mer sizes\n", filename);
      fprintf(stderr, "existDB::existDB()-- Got "uint32FMT", expected "uint32FMT"\n", _merSizeInBases, merSize);
      fail = true;
    }

    if (fail)
      exit(1);

    return;
  }

  //  If no direction flags are set, set the default direction of
  //  forward.  Stupid precedence rules.
  //
  if ((flags & (existDBcanonical | existDBforward)) == uint32ZERO)
    flags |= existDBforward;

  //  If we can open 'filename' for reading, then we assume the file
  //  is a multi-fasta, and we build an existDB/
  //
  //  Otherwise, we assume that 'filename' is really the prefix for a
  //  meryl database.


  if (fileExists(filename))
    createFromFastA(filename, merSize, flags);
  else
    createFromMeryl(filename, merSize, lo, hi, flags);
}


existDB::existDB(char const    *sequence,
                 uint32         merSize,
                 existDBflags   flags) {
  clear();

  _compressedHash   = flags & existDBcompressHash;
  _compressedBucket = flags & existDBcompressBuckets;
  _compressedCounts = flags & existDBcompressCounts;

  if ((flags & (existDBcanonical | existDBforward)) == uint32ZERO)
    flags |= existDBforward;

  createFromSequence(sequence, merSize, flags);
}


existDB::~existDB() {
  delete [] _hashTable;
  delete [] _buckets;
  delete [] _counts;
}





bool
existDB::exists(uint64 mer) {
  uint64 c, h, st, ed;

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


uint64
existDB::count(uint64 mer) {
  uint64 c, h, st, ed;

  if (_counts == 0L)
    return(0);

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
    return(0);

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
    return(getDecodedValue(_counts, st * _cntWidth, _cntWidth));
  else
    return(_counts[st]);
}
