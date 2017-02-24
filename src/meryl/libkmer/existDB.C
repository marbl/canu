
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
 *    Brian P. Walenz from 2003-JAN-02 to 2003-OCT-20
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-12 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-20 to 2014-APR-11
 *      are Copyright 2005-2008,2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

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
