#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "existDB.H"
#include "positionDB.H"
#include "bio++.H"


existDB::existDB(char const  *filename,
                 bool         loadData) {
  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "existDB::existDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}


existDB::existDB(char const  *filename,
                 u32bit       merSize,
                 u32bit       tblBits,
                 u32bit       lo,
                 u32bit       hi,
                 positionDB  *posDB,
                 bool         compressHash,
                 bool         compressBuckets,
                 bool         beVerbose) {

  _compressedHash   = compressHash;
  _compressedBucket = compressBuckets;
  _beVerbose        = beVerbose;

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

  //
  //  If we can open 'filename' for reading, then we assume the file is a multi-fasta, and
  //  we build an existDB using merSize = p1 and tblBits = p2
  //
  //  Otherwise, we assume that 'filename' is really the prefix for a meryl database.
  //  

  if (fileExists(filename)) {
    //fprintf(stderr, "Loading mers from fasta '%s'\n", filename);
    createFromFastA(filename, merSize, tblBits, posDB);
  } else {
    //fprintf(stderr, "Loading mers from meryl '%s'\n", filename);
    createFromMeryl(filename, lo, hi, tblBits, posDB);
  }
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

