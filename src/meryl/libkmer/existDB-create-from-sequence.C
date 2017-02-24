
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
 *    Brian P. Walenz from 2012-MAY-08 to 2014-APR-11
 *      are Copyright 2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-AUG-31
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
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
#include "seqCache.H"
#include "seqStream.H"
#include "merStream.H"

bool
existDB::createFromSequence(char const  *sequence,
                            uint32       merSize,
                            uint32       flags) {

  bool               beVerbose  = false;
  bool               rebuilding = false;

  _hashTable = 0L;
  _buckets   = 0L;
  _counts    = 0L;

  _merSizeInBases        = merSize;

  _searchForDupe = true;

  if ((flags & existDBcompressHash) ||
      (flags & existDBcompressBuckets) ||
      (flags & existDBcompressCounts))
    fprintf(stderr, "existDB::createFromSequence: compression not supported.\n"), exit(1);

  //  This (at =22) eats up 16MB, and should allow a lot of mers at big sizes.  Unfortunately, we
  //  know nothing about how man mers are going to be in the input.
  //
  //  Setting this too high drastically reduces performance, suspected because of cache misses.
  //  Setting this too low will also reduce performance, by increasing the search time in a bucket.
  //
  uint32 tblBits = logBaseTwo64(strlen(sequence));

 rebuild:
  _shift1                = 2 * _merSizeInBases - tblBits;
  _shift2                = _shift1 / 2;
  _mask1                 = uint64MASK(tblBits);
  _mask2                 = uint64MASK(_shift1);

  _hshWidth              = uint32ZERO;
  _chkWidth              = 2 * merSize - tblBits;
  _cntWidth              = 16;

  uint64  tableSizeInEntries = uint64ONE << tblBits;
  uint64  numberOfMers       = uint64ZERO;
  uint64 *countingTable      = new uint64 [tableSizeInEntries + 1];

  for (uint64 i=tableSizeInEntries+1; i--; )
    countingTable[i] = 0;

  _isCanonical = flags & existDBcanonical;
  _isForward   = flags & existDBforward;

  assert(_isCanonical + _isForward == 1);

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //
  merStream    *M = new merStream(new kMerBuilder(_merSizeInBases),
                                  new seqStream(sequence, strlen(sequence)),
                                  true, true);

  while (M->nextMer()) {
    if (_isForward) {
      countingTable[ HASH(M->theFMer()) ]++;
      numberOfMers++;
    }

    if (_isCanonical) {
      countingTable[ HASH(M->theCMer()) ]++;
      numberOfMers++;
    }
  }

  delete M;

#ifdef STATS
  uint64  dist[32] = {0};
  uint64  maxcnt = 0;
  for (uint64 i=tableSizeInEntries+1; i--; ) {
    if (countingTable[i] > maxcnt)
      maxcnt = countingTable[i];

    if (countingTable[i] < 32)
      dist[countingTable[i]]++;
  }

  for(uint64 i=0; i<32; i++)
    fprintf(stderr, "existDB::usage[%2d] = %d\n", i, dist[i]);
  fprintf(stderr, "existDB::maxcnt    = %d\n", maxcnt);
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Determine how many bits we need to hold the value
  //  numberOfMers.....then....
  //
  //  This is numberOfMers+1 because we need to store the
  //  first position after the last mer.  That is, if there are two
  //  mers, we will store that the first mer is at position 0, the
  //  second mer is at position 1, and the end of the second mer is at
  //  position 2.
  //
  if (_compressedHash) {
    _hshWidth = 1;
    while ((numberOfMers+1) > (uint64ONE << _hshWidth))
      _hshWidth++;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate a hash table and some mer storage buckets.
  //
  _hashTableWords = tableSizeInEntries + 2;
  if (_compressedHash)
    _hashTableWords = _hashTableWords * _hshWidth / 64 + 1;

  _bucketsWords = numberOfMers + 2;
  if (_compressedBucket)
    _bucketsWords = _bucketsWords * _chkWidth / 64 + 1;

  _countsWords = numberOfMers + 2;
  if (_compressedCounts)
    _countsWords = _countsWords * _cntWidth / 64 + 1;

  if (beVerbose) {
    fprintf(stderr, "existDB::createFromSequence()-- hashTable is "uint64FMT"MB\n", _hashTableWords >> 17);
    fprintf(stderr, "existDB::createFromSequence()-- buckets is "uint64FMT"MB\n", _bucketsWords >> 17);
    if (flags & existDBcounts)
      fprintf(stderr, "existDB::createFromSequence()-- counts is "uint64FMT"MB\n", _countsWords >> 17);
  }

  _hashTable   = new uint64 [_hashTableWords];
  _buckets     = new uint64 [_bucketsWords];
  _countsWords = (flags & existDBcounts) ?             _countsWords  : 0;
  _counts      = (flags & existDBcounts) ? new uint64 [_countsWords] : 0L;

  //  These aren't strictly needed.  _buckets is cleared as it is initialied.  _hashTable
  //  is also cleared as it is initialized, but in the _compressedHash case, the last
  //  few words might be uninitialized.  They're unused.

  //memset(_hashTable, 0, sizeof(uint64) * _hashTableWords);
  //memset(_buckets,   0, sizeof(uint64) * _bucketsWords);  //  buckets is cleared as it is built
  //memset(_counts,    0, sizeof(uint64) * _countsWords);

  _hashTable[_hashTableWords-1] = 0;
  _hashTable[_hashTableWords-2] = 0;
  _hashTable[_hashTableWords-3] = 0;
  _hashTable[_hashTableWords-4] = 0;

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Make the hash table point to the start of the bucket, and reset
  //  the counting table -- we're going to use it to fill the buckets.
  //
  uint64  tmpPosition = 0;
  uint64  begPosition = 0;
  uint64  ptr         = 0;

  if (_compressedHash) {
    for (uint64 i=0; i<tableSizeInEntries; i++) {
      tmpPosition    = countingTable[i];
      countingTable[i] = begPosition;

      setDecodedValue(_hashTable, ptr, _hshWidth, begPosition);
      ptr         += _hshWidth;

      begPosition += tmpPosition;
    }

    setDecodedValue(_hashTable, ptr, _hshWidth, begPosition);
  } else {
    for (uint64 i=0; i<tableSizeInEntries; i++) {
      tmpPosition    = countingTable[i];
      countingTable[i] = begPosition;

      _hashTable[i] = begPosition;

      begPosition += tmpPosition;
    }

    //  Set the last position in the hash, but we don't care about
    //  the temporary counting table.
    //
    _hashTable[tableSizeInEntries] = begPosition;
  }




  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers, placed into buckets
  //
  M  = new merStream(new kMerBuilder(_merSizeInBases),
                     new seqStream(sequence, strlen(sequence)),
                     true, true);

  while (M->nextMer()) {
    if (_isForward)
      insertMer(HASH(M->theFMer()), CHECK(M->theFMer()), 1, countingTable);

    if (_isCanonical)
      insertMer(HASH(M->theCMer()), CHECK(M->theCMer()), 1, countingTable);
  }

  delete M;

  //  Compress out the gaps we have from redundant kmers.

  uint64  pos = 0;
  uint64  frm = 0;
  uint64  len = 0;

  for (uint64 i=0; i<tableSizeInEntries; i++) {
    frm = _hashTable[i];
    len = countingTable[i] - _hashTable[i];

    _hashTable[i] = pos;

    for (uint64 j=0; j<len; j++) {
      if (_counts)
        _counts[pos]  = _counts[frm];

      _buckets[pos++] = _buckets[frm++];
    }
  }

  if (beVerbose)
    fprintf(stderr, "Compressed from "uint64FMT" to "uint64FMT" ("uint64FMT" bits)\n",
            _hashTable[tableSizeInEntries], pos, logBaseTwo64(pos));

  while (pos < _bucketsWords)
    _buckets[pos++] = 0;

  _hashTable[tableSizeInEntries] = pos;

  //  All done.  Delete temporary stuff
  //
  delete [] countingTable;

  //  But if we horribly screwed up the estimate of tblBits, reset and recompute

  if ((logBaseTwo64(pos) < tblBits) &&
      (rebuilding == false)) {
    rebuilding = true;

    delete [] _hashTable;
    delete [] _buckets;
    delete [] _counts;

    _hashTable = 0L;
    _buckets   = 0L;
    _counts    = 0L;

    tblBits = logBaseTwo64(pos);

    goto rebuild;
  }

  return(true);
}
