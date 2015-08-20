
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
 *    Brian P. Walenz from 2003-JAN-02 to 2004-APR-05
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-JAN-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Clark Mobarry on 2004-FEB-12
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-19 to 2014-APR-11
 *      are Copyright 2005-2008,2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <new>

#include "bio++.H"
#include "positionDB.H"
#include "existDB.H"
#include "libmeryl.H"

#undef ERROR_CHECK_COUNTING
#undef ERROR_CHECK_COUNTING_ENCODING
#undef ERROR_CHECK_EMPTY_BUCKETS

//  This tests Chunlin Xiao's discovered bug -- if there are a small
//  number of unique mers, compared to distinct mers (2 * #unique_mers
//  < #distinct_mers, we would overflow the position pointer in
//  buckets.  This enables a check that it doesn't occur.
//
//  This has a fixed allocation size, and crashes on larger inputs.
//
#undef  TEST_NASTY_BUGS

//  Tests that mers are masked out properly.  Doesn't handle canonical
//  mers though.
//
#undef  MER_REMOVAL_TEST


positionDB::positionDB(char const        *filename,
                       uint32             merSize,
                       uint32             merSkip,
                       uint32             maxMismatch,
                       bool               loadData) {
  memset(this, 0, sizeof(positionDB));

  //  loadData == false only for driver-posDB.C, and only so it can
  //  dump stats on a posDB file.

  if (loadState(filename, true, false) == false) {
    fprintf(stderr, "positionDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }

  if ((loadData) && (merSize != _merSizeInBases)) {
    fprintf(stderr, "positionDB()-- Tried to read state from '%s', but mer size is wrong (found "uint32FMT", wanted "uint32FMT").\n",
            filename, _merSizeInBases, merSize);
    exit(1);
  }

  if ((loadData) && (merSkip != _merSkipInBases)) {
    fprintf(stderr, "positionDB()-- Tried to read state from '%s', but mer skip is wrong (found "uint32FMT", wanted "uint32FMT").\n",
            filename, _merSkipInBases, merSkip);
    exit(1);
  }

if ((loadData) && (maxMismatch != _nErrorsAllowed)) {
    fprintf(stderr, "positionDB()-- Tried to read state from '%s', but max number of mismatches is wrong (found "uint32FMT", wanted "uint32FMT").\n",
            filename, _nErrorsAllowed, maxMismatch);
    exit(1);
  }

  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "positionDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}


positionDB::positionDB(merStream          *MS,
                       uint32              merSize,
                       uint32              merSkip,
                       existDB            *mask,
                       existDB            *only,
                       merylStreamReader  *counts,
                       uint32              minCount,
                       uint32              maxCount,
                       uint32              maxMismatch,
                       uint32              maxMemory,
                       bool                beVerbose) {

  memset(this, 0, sizeof(positionDB));

  //  Guesstimate a nice table size based on the number of input mers
  //  and the mersize, unless the user gave us a table size.
  //
  //  We need to ensure that
  //    2 * merSize + posnWidth + 1 - 64 <= tblBits <= 2 * merSize - 4
  //
  //  The catch is that we don't exactly know posnWidth right now.  We
  //  can overestimate it, though, based on the size of the sequence
  //  that is backing the merStream.
  //
  //  The second catch is that we don't want to make tblBits too big
  //  or too small.  If too big, we waste a lot of memory in the hash
  //  table pointers, and if too small, we waste even more memory in
  //  the data table (not to mention the algorithm dies because it
  //  assumed buckets in the data table are small).
  //
  //  The memory size is (roughly):
  //
  //    2^tblBits * log(numDistinctMers) +
  //    numDistinctMers * (2*merSize - tblBits + 1 + log(numMers) +
  //    (numMers - numUniqieMers) * log(numMers)
  //
  //  this is approximately proportional to:
  //
  //    2^tblBits * posnWidth +
  //    approxMers * (2*merSize - tblBits + 1 + posnWidth)
  //
  uint64  approxMers = MS->approximateNumberOfMers();
  uint64  posnWidth  = logBaseTwo64(approxMers + 1);

  //  Find the smallest and largest tblBits we could possibly use.
  //
  uint64  sm = 2 * merSize + posnWidth + 1 - 64;
  uint64  lg = 2 * merSize - 4;

  if (2 * merSize + posnWidth + 1 < 64)
    sm = 2;

  if (sm < 16)
    sm = 16;

  if (sm > lg) {
    fprintf(stderr, "ERROR:  too many mers for this mersize!\n");
    fprintf(stderr, "        sm         = "uint64FMT"\n", sm);
    fprintf(stderr, "        lg         = "uint64FMT"\n", lg);
    fprintf(stderr, "        merSize    = "uint32FMT" bits\n", 2 * merSize);
    fprintf(stderr, "        approxMers = "uint64FMT" mers\n", approxMers);
    fprintf(stderr, "        posnWidth  = "uint64FMT" bits\n", posnWidth);
    exit(1);
  }


  //  Iterate through all the choices, picking the one with the
  //  smallest expected footprint.
  //
  {

    if (beVerbose) {
      fprintf(stderr, "potential configurations for approximately "uint64FMT" "uint32FMT"-mers (posnW="uint64FMT").\n",
              approxMers, merSize, posnWidth);
    }

    uint64  mini = 0;      //  tblSize of the smallest found
    uint64  minm = ~mini;  //  memory size of the smallest found
    double  minw = 0.0;    //  work of the smallest found

    uint64  memory    = 0;
    double  effort    = 0;

    if (maxMemory == 0)
      maxMemory = ~uint32ZERO;

    for (uint64 i=sm; i<=lg; i++) {

      //  These are only needed if maxMismatch is set, but it's
      //  simpler to always set.
      //
      _merSizeInBases        = merSize;
      _merSizeInBits         = 2 * _merSizeInBases;
      _merSkipInBases        = merSkip;
      _tableSizeInBits       = i;
      _tableSizeInEntries    = uint64ONE << _tableSizeInBits;
      _hashWidth             = uint32ZERO;
      _hashMask              = uint64MASK(_tableSizeInBits);
      _chckWidth             = _merSizeInBits - _tableSizeInBits;
      _posnWidth             = uint64ZERO;
      _sizeWidth             = 0;

      _shift1                = _merSizeInBits - _tableSizeInBits;
      _shift2                = _shift1 / 2;
      _mask1                 = uint64MASK(_tableSizeInBits);
      _mask2                 = uint64MASK(_shift1);

      //  Everyone wants to know the memory size (in MB).
      //
      memory = ((uint64ONE << i) * posnWidth + approxMers * (2*merSize - i + 1 + posnWidth)) >> 23;

      //  If we know we're looking for mismatches, we compute the amount
      //  of work needed per lookup, and use that, instead of strict
      //  memory sizing, to deicde the table size.
      //
      if (maxMismatch > 0)
        effort = setUpMismatchMatcher(maxMismatch, approxMers);

      //  If our memory size is smaller than allowed, AND it's the
      //  smallest, or the work is smaller, save the table size.
      //
      if ((memory < maxMemory) &&
          ((memory < minm) ||
           (effort < minw))) {
        mini = i;
        minm = memory;
        minw = effort;
      }

      if (beVerbose) {
        fprintf(stderr, "tblBits="uint64FMTW(2)" shifts="uint32FMTW(02)","uint32FMTW(02)" -- size %8.3fGB -- work %8.3f%s\n",
                i, _shift1, _shift2, memory / 1024.0, effort, (mini == i) ? " ***" : "");
      }
    }

    _tableSizeInBits = mini;
  }


  if (_tableSizeInBits == 0) {
    fprintf(stderr, "ERROR:  No positionDB parameters within allowed memory limit.\n");
    exit(1);
  }


  if (beVerbose) {
    uint32 s1 = 2*merSize-_tableSizeInBits;
    fprintf(stderr, "tblBits="uint32FMT" s1="uint32FMT" s2="uint32FMT" -- merSize="uint32FMT" bits + posnWidth="uint64FMT" bits (est "uint64FMT" mers) FINAL\n",
            _tableSizeInBits, s1, s1/2, merSize, posnWidth, approxMers);
  }


  _merSizeInBases        = merSize;
  _merSizeInBits         = 2 * _merSizeInBases;
  _merSkipInBases        = merSkip;
  _tableSizeInEntries    = uint64ONE << _tableSizeInBits;
  _hashWidth             = uint32ZERO;
  _hashMask              = uint64MASK(_tableSizeInBits);
  _chckWidth             = _merSizeInBits - _tableSizeInBits;
  _posnWidth             = uint64ZERO;
  _sizeWidth             = 0;

  if (maxCount == 0)
    maxCount = ~uint32ZERO;

  if (counts)
    _sizeWidth = (maxCount < ~uint32ZERO) ? logBaseTwo64(maxCount+1) : 32;

  _shift1                = _merSizeInBits - _tableSizeInBits;
  _shift2                = _shift1 / 2;
  _mask1                 = uint64MASK(_tableSizeInBits);
  _mask2                 = uint64MASK(_shift1);

#if 0
  fprintf(stderr, "merSizeInBits   "uint32FMT"\n", _merSizeInBits);
  fprintf(stderr, "hashWidth       "uint32FMT"\n", _hashWidth);
  fprintf(stderr, "chckWidth       "uint32FMT"\n", _chckWidth);
  fprintf(stderr, "shift1          "uint32FMT"\n", _shift1);
  fprintf(stderr, "shift2          "uint32FMT"\n", _shift2);
#endif

  if (maxMismatch > 0)
    setUpMismatchMatcher(maxMismatch, approxMers);

  build(MS, mask, only, counts, minCount, maxCount, beVerbose);
}



void
positionDB::build(merStream          *MS,
                  existDB            *mask,
                  existDB            *only,
                  merylStreamReader  *counts,
                  uint32              minCount,
                  uint32              maxCount,
                  bool                beVerbose) {

  _bucketSizes           = 0L;
  _countingBuckets       = 0L;
  _hashTable_BP          = 0L;
  _hashTable_FW          = 0L;
  _buckets               = 0L;
  _positions             = 0L;

  _wCnt                  = 0;
  _wFin                  = 0;

  //  For get/setDecodedValues().
  uint64  lensC[4] = {~uint64ZERO, ~uint64ZERO, ~uint64ZERO, ~uint64ZERO};
  uint64  lensF[4] = {~uint64ZERO, ~uint64ZERO, ~uint64ZERO, ~uint64ZERO};
  uint64  vals[4]  = {0};
  uint64  nval     = (_sizeWidth == 0) ? 3 : 4;

  _numberOfMers          = uint64ZERO;
  _numberOfPositions     = uint64ZERO;
  _numberOfDistinct      = uint64ZERO;
  _numberOfUnique        = uint64ZERO;
  _numberOfEntries       = uint64ZERO;
  _maximumEntries        = uint64ZERO;

  //  We assume later that these are already allocated.
  _sortedMax             = 16384;
  _sortedChck            = new uint64 [_sortedMax];
  _sortedPosn            = new uint64 [_sortedMax];

  if (MS == 0L) {
    fprintf(stderr, "positionDB()-- ERROR: No merStream?  Nothing to build a table with!\n");
    exit(1);
  }

  MS->rewind();


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //

  //  We'll later want to reuse the _bucketSizes space for storing the
  //  hash table.  To make it somewhat safe, we allocate the space as
  //  uint64, then cast it to be uint32.
  //
  //  bktAllocIsJunk tells us if we should release this memory (if we
  //  need to allocate separate space for the hash table).  We'd need
  //  to do this if the hashWidth is more than 32 bits, but we won't
  //  know that for a little bit.
  //
  //  The _bucketSizes is offset by one from bktAlloc so that we don't
  //  overwrite _bucketSizes when we are constructing hash table.
  //
  uint64 *bktAlloc;
  try {
    bktAlloc = new uint64 [_tableSizeInEntries / 2 + 4];
  } catch (std::bad_alloc) {
    fprintf(stderr, "positionDB()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "positionDB()-- bktAlloc = new uint64 ["uint64FMT"]\n", _tableSizeInEntries / 2 + 4);
    exit(1);
  }
  bool     bktAllocIsJunk = false;

  bzero(bktAlloc, sizeof(uint64) * (_tableSizeInEntries / 2 + 4));

  //  Why +2?  We try to reuse the bktAlloc space for the hash table,
  //  which is constructed from the bucketSizes.  The hashTable is
  //  built from the bucketSizes.  It definitely needs to be +1, and
  //  so we use +2 just in case the human is being stupid again.
  //
  _bucketSizes = (uint32 *)(bktAlloc + 2);

#ifdef ERROR_CHECK_COUNTING
  fprintf(stdout, "ERROR_CHECK_COUNTING is defined.\n");
  uint32 *_errbucketSizes = new uint32 [_tableSizeInEntries + 2];
  for (uint64 i=0; i<_tableSizeInEntries + 2; i++)
    _errbucketSizes[i] = uint32ZERO;
#endif

  if (beVerbose)
    fprintf(stderr, "    Allocated bucket size counting space with total size "uint64FMT" KB\n", _tableSizeInEntries >> 8);


  speedCounter  *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

  //  Two choices here
  //
  //  1)  No masking or onlying is done.  Stream the mers and just
  //      count the positions.  This is the original behavior.
  //
  //  2)  Masking or onlying is done.  Open the output stream file,
  //      stream the mers by, checking for mask/only of both
  //      forward and reverse mers.  If either is found, push
  //      the (forward) mer and position onto the stream.
  //      close the output stream.
  //
  //      Save the mer if it doesn't exist in the mask (both f and r),
  //      or does exist in the only (either f or r), add it.
  //
  //      The input databases for mask and only are (currently) made
  //      using canonical mers.  We halve the number of exists() by
  //      also using canonical mers here.
  //

  MS->rewind();

  while (MS->nextMer(_merSkipInBases)) {
    _bucketSizes[ HASH(MS->theFMer()) ]++;

#ifdef ERROR_CHECK_COUNTING
    _errbucketSizes[ HASH(MS->theFMer()) ]++;
#endif

    _numberOfMers++;
    _numberOfPositions = MS->thePositionInStream();
    assert((_numberOfPositions >> 60) == 0);
    C->tick();
  }


  delete C;
  C = 0L;

  if (beVerbose)
    fprintf(stderr, "    Found "uint64FMT" mers (max position = "uint64FMT")\n", _numberOfMers, _numberOfPositions);

  //  This caught a nasty bug in merStream rewind(), and it's pretty
  //  cheap, so I left it in.  Search for the other DEBUGnumPositions.
  //
  uint64 DEBUGnumPositions = _numberOfPositions + 1;

  //  This is _numberOfMers+1 because we need to store the first
  //  position after the last mer.  That is, if there are two mers, we
  //  will store that the first mer is at position 0, the second mer
  //  is at position 1, and the end of the second mer is at position
  //  2.
  //
  //  In reality, it should be the number of distinct mers, not the
  //  total number of mers, but we don't know that yet.  And so
  //  occasionally we'll make things too big and waste a bit of
  //  memory.
  //
  _hashWidth = logBaseTwo64(_numberOfMers+1);
  _posnWidth = logBaseTwo64(_numberOfPositions+1);



  ///////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate buckets and make bucketSizes be a pointer into them
  //
  _wCnt          = _chckWidth + _posnWidth + 1 + _sizeWidth;

  lensC[0] = _chckWidth;
  lensC[1] = _posnWidth;
  lensC[2] = 1;
  lensC[3] = _sizeWidth;

  uint64   bucketsSpace  = (_numberOfMers+1) * _wCnt / 64 + 1;
  uint32   endPosition   = 0;

  if (beVerbose)
    fprintf(stderr, "    Allocated "uint64FMT"KB for buckets ("uint64FMT" 64-bit words)\n", bucketsSpace >> 7, bucketsSpace);
  try {
    _countingBuckets = new uint64 [bucketsSpace];
  } catch (std::bad_alloc) {
    fprintf(stderr, "positionDB()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "positionDB()-- _countingBuckets = new uint64 ["uint64FMT"]\n", bucketsSpace);
    exit(1);
  }

  for (uint64 i=0; i<bucketsSpace; i++)
    _countingBuckets[i] = ~uint64ZERO;

  for (uint64 i=0; i<_tableSizeInEntries; i++) {
    endPosition     += _bucketSizes[i];
    _bucketSizes[i]  = endPosition;
  }
  _bucketSizes[_tableSizeInEntries] = endPosition;

#ifdef ERROR_CHECK_COUNTING
  if (endPosition != _numberOfMers)
    fprintf(stdout, "ERROR_CHECK_COUNTING: BUCKETSIZE COUNTING PROBLEM -- endPos="uint32FMT" != numMers="uint64FMT"\n",
            endPosition, _numberOfMers);
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers with positions
  //
  if (beVerbose)
    fprintf(stderr, "    Building lists with positions.\n");

  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);

#ifdef ERROR_CHECK_COUNTING_ENCODING
  fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING is defined!\n");
#endif


  MS->rewind();

  while (MS->nextMer(_merSkipInBases)) {
    uint64 h = HASH(MS->theFMer());

#ifdef ERROR_CHECK_COUNTING
    if (_bucketSizes[h] == 0) {
      char  str[33];
      fprintf(stderr, "positionDB()-- ERROR_CHECK_COUNTING: Bucket "uint64FMT" ran out of things!  '%s'\n", h, MS->theFMer().merToString(str));
      fprintf(stderr, "positionDB()-- ERROR_CHECK_COUNTING: Stream is at "uint64FMT"\n", MS->thePositionInStream());
    }
#endif

    _bucketSizes[h]--;

#ifdef ERROR_CHECK_COUNTING
    _errbucketSizes[h]--;
#endif


#ifdef ERROR_CHECK_EMPTY_BUCKETS
    //  Check that everything is empty.  Empty is defined as set to all 1's.
    getDecodedValues(_countingBuckets, (uint64)_bucketSizes[h] * (uint64)_wCnt, nval, lensC, vals);

    if (((~vals[0]) & uint64MASK(lensC[0])) ||
        ((~vals[1]) & uint64MASK(lensC[1])) ||
        ((~vals[2]) & uint64MASK(lensC[2])) ||
        ((lensC[3] > 0) && ((~vals[3]) & uint64MASK(lensC[3]))))
      fprintf(stdout, "ERROR_CHECK_EMPTY_BUCKETS: countingBucket not empty!  pos=%lu 0x%016lx 0x%016lx 0x%016lx 0x%016lx\n",
              _bucketSizes[h] * _wCnt,
              (~vals[0]) & uint64MASK(lensC[0]),
              (~vals[1]) & uint64MASK(lensC[1]),
              (~vals[2]) & uint64MASK(lensC[2]),
              (~vals[3]) & uint64MASK(lensC[3]));
#endif

    vals[0] = CHECK(MS->theFMer());
    vals[1] = MS->thePositionInStream();
    vals[2] = 0;
    vals[3] = 0;

    setDecodedValues(_countingBuckets, (uint64)_bucketSizes[h] * (uint64)_wCnt, nval, lensC, vals);

#ifdef ERROR_CHECK_COUNTING_ENCODING
    getDecodedValues(_countingBuckets, (uint64)_bucketSizes[h] * (uint64)_wCnt, nval, lensC, vals);

    if (vals[0] != CHECK(MS->theFMer()))
      fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  CHCK corrupted!  Wanted "uint64HEX" got "uint64HEX"\n",
              CHECK(MS->theFMer()), vals[0]);
    if (vals[1] != MS->thePositionInStream())
      fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  POSN corrupted!  Wanted "uint64HEX" got "uint64HEX"\n",
              MS->thePositionInStream(), vals[1]);
    if (vals[2] != 0)
      fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  UNIQ corrupted.\n");
    if (vals[3] != 0)
      fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  SIZE corrupted.\n");
#endif

    C->tick();
  }


  delete C;
  C = 0L;

#ifdef ERROR_CHECK_COUNTING
  for (uint64 i=0; i<_tableSizeInEntries; i++)
    if (_errbucketSizes[i] != 0)
      fprintf(stdout, "ERROR_CHECK_COUNTING: Bucket "uint32FMT" wasn't filled fully?  "uint32FMT" left over.\n", i, _errbucketSizes[i]);

  delete [] _errbucketSizes;
  _errbucketSizes = 0L;
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  4)  Sort each bucket -- count:
  //        1) number of distinct mers
  //        2) number of unique mers
  //        3) number of entries in position table ( sum mercount+1 for all mercounts > 1)
  //      also need to repack the sorted things
  //
  if (beVerbose)
    fprintf(stderr, "    Sorting and repacking buckets ("uint64FMT" buckets).\n", _tableSizeInEntries);

  C = new speedCounter("    %7.2f Mbuckets -- %5.2f Mbuckets/second\r", 1000000.0, 0x1ffffff, beVerbose);
  for (uint64 i=0; i<_tableSizeInEntries; i++) {
    sortAndRepackBucket(i);
    C->tick();
  }
  delete C;
  C = 0L;

  if (beVerbose)
    fprintf(stderr,
            "    Found "uint64FMTW(12)" total mers\n"
            "    Found "uint64FMTW(12)" distinct mers\n"
            "    Found "uint64FMTW(12)" unique mers\n"
            "    Need "uint64FMT" non-unique position list entries ("uint64FMT" maximum count)\n",
            _numberOfMers, _numberOfDistinct, _numberOfUnique, _numberOfEntries, _maximumEntries);



  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Compute the size of the final bucket position entry.  It's
  //  either a position into the sequence, or a pointer into a list of
  //  positions.  In rare cases, the pointer is larger than the
  //  sequence position, and we need to do extra work.
  //
  //  The width of position pointers (in buckets) is the max of
  //  _posnWidth (a pointer to the sequence position) and
  //  _pptrWidth (a pointer to an entry in the positions table).
  //
  _pptrWidth = logBaseTwo64(_numberOfEntries+1);
  if (_pptrWidth < _posnWidth)
    _pptrWidth = _posnWidth;

  _wFin = _chckWidth + _pptrWidth + 1 + _sizeWidth;

  lensF[0] = _chckWidth;
  lensF[1] = _pptrWidth;
  lensF[2] = 1;
  lensF[3] = _sizeWidth;

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  5)  Allocate: real hash table, buckets and position table.
  //

  //  XXXX how do we count the number of buckets/positions we never
  //  use because they are masked out??
  //
  //  If we are just thresholding (ignore things with count > 100)
  //  it's easy, a simple loop over something.
  //
  //  If we have an exist/only db....are they in the same order?  Can
  //  we loop over both at the same time and count that way?  That'd
  //  be cool!  Mersize is the same, why can the table size be the
  //  same too -- OK, if the existDB has a small number of mers in it,
  //  then we don't need a large table.

  uint64  hs = _tableSizeInEntries * _hashWidth / 64 + 1;
  uint64  bs = _numberOfDistinct   * _wFin      / 64 + 1;
  uint64  ps = _numberOfEntries    * _posnWidth / 64 + 1;

  if (_hashWidth <= 32) {
    if (beVerbose)
      fprintf(stderr, "    Reusing bucket counting space for hash table.\n");

#ifdef UNCOMPRESS_HASH_TABLE
    _hashTable_BP  = 0L;
    _hashTable_FW  = (uint32 *)bktAlloc;
#else
    _hashTable_BP  = bktAlloc;
    _hashTable_FW  = 0L;
#endif

    bktAllocIsJunk = false;
  } else {

    //  Can't use the full-width hash table, since the data size is >
    //  32 bits -- we'd need to allocate 64-bit ints for it, and
    //  that'll likely be too big...and we'd need to have
    //  _hashTable_FW64 or something.

    if (beVerbose)
      fprintf(stderr, "    Allocated "uint64FMTW(10)"KB for hash table ("uint64FMT" 64-bit words)\n", hs >> 7, hs);
    try {
      _hashTable_BP = new uint64 [hs];
      _hashTable_FW = 0L;
    } catch (std::bad_alloc) {
      fprintf(stderr, "positionDB()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "positionDB()-- _hashTable_BP = new uint64 ["uint64FMT"]\n", hs);
      exit(1);
    }
    bktAllocIsJunk = true;
  }


  //  If we have enough space to reuse the counting space, reuse it.
  //  Else, allocate more space.
  //
  //  We need to ensure that there are enough bits and that the size
  //  of a bucket didn't increase.  If the bucket size did increase,
  //  and we see more unique buckets than total mers (up to some
  //  point) we overwrite data.
  //
  //  Recall that bucketSpace ~= numberOfMers * wCnt
  //
  if ((bs < bucketsSpace) && (_wFin <= _wCnt)) {
    if (beVerbose)
      fprintf(stderr, "    Reusing bucket space; Have: "uint64FMT"  Need: "uint64FMT" (64-bit words)\n", bucketsSpace, bs);

    _buckets = _countingBuckets;

    bs = bucketsSpace;  // for output at the end
  } else {
    if (beVerbose)
      fprintf(stderr, "    Allocated "uint64FMTW(10)"KB for buckets    ("uint64FMT" 64-bit words)\n", bs >> 7, bs);
    try {
      _buckets   = new uint64 [bs];
    } catch (std::bad_alloc) {
      fprintf(stderr, "positionDB()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "positionDB()-- _buckets = new uint64 ["uint64FMT"]\n", bs);
      exit(1);
    }
  }

  if (beVerbose)
    fprintf(stderr, "    Allocated "uint64FMTW(10)"KB for positions  ("uint64FMT" 64-bit words)\n", ps >> 7, ps);
  try {
    _positions = new uint64 [ps];
  } catch (std::bad_alloc) {
    fprintf(stderr, "positionDB()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "positionDB()-- _positions = new uint64 ["uint64FMT"\n", ps);
    exit(1);
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  6)  Transfer from the sorted buckets to the hash table.
  //
  if (beVerbose)
    fprintf(stderr, "    Transferring to final structure ("uint64FMT" buckets).\n", _tableSizeInEntries);

  uint64   bucketStartPosition = 0;

  //  Current positions and bit positions in the buckets and position list.
  //
  uint64  currentBbit = uint64ZERO;  //  Bit position into bucket
  uint64  currentPbit = uint64ZERO;  //  Bit position into positions
  uint64  currentPpos = uint64ZERO;  //  Value position into positions

#ifdef TEST_NASTY_BUGS
  //  Save the position array pointer of each bucket for debugging.
  //
  uint64  currentBpos = uint64ZERO;  //  Value position into bucket
  uint32 *posPtrCheck = new uint32 [65826038];
#endif

  //  We also take this opportunity to reset some statistics that are
  //  wrong.
  //
  _numberOfMers      = 0;
  _numberOfPositions = 0;
  _numberOfDistinct  = 0;
  _numberOfUnique    = 0;
  _numberOfEntries   = 0;
  _maximumEntries    = 0;

  C = new speedCounter("    %7.2f Mbuckets -- %5.2f Mbuckets/second\r", 1000000.0, 0x1ffffff, beVerbose);

  //  We need b outside the loop!
  //
  uint64  b;
  for (b=0; b<_tableSizeInEntries; b++) {
    C->tick();

    //  Set the start of the bucket -- we took pains to ensure that
    //  we don't overwrite _bucketSizes[b], if we are reusing that
    //  space for the hash table.
    //
    if (_hashTable_BP)
      setDecodedValue(_hashTable_BP, (uint64)b * (uint64)_hashWidth, _hashWidth, bucketStartPosition);
    else
      _hashTable_FW[b] = bucketStartPosition;

    //  Get the number of mers in the counting bucket.  The error
    //  checking and sizing of _sortedChck and _sortedPosn was already
    //  done in the sort.
    //
    uint64 st = _bucketSizes[b];
    uint64 ed = _bucketSizes[b+1];
    uint32 le = ed - st;

    //  Unpack the check values
    //
    for (uint64 i=st, J=st * _wCnt; i<ed; i++, J += _wCnt) {
      getDecodedValues(_countingBuckets, J, 2, lensC, vals);
      _sortedChck[i-st] = vals[0];
      _sortedPosn[i-st] = vals[1];
    }


    //  Walk through the counting bucket, adding things to the real
    //  bucket as we see them.  Mers with more than one position are
    //  inserted into the bucket, and the positions inserted into the
    //  position list.

    //  start and end locations of the mer.  For mers with only
    //  one occurrance (unique mers), stM+1 == edM.
    //
    uint32  stM = uint32ZERO;
    uint32  edM = uint32ZERO;

    while (stM < le) {

      //  Move to the next mer.
      //
      edM++;

      //  Keep moving while the two mers are the same.
      //
      while ((edM < le) && (_sortedChck[stM] == _sortedChck[edM]))
        edM++;

      //  edM is now the mer after the last.  Write all mers from stM
      //  up to edM to the final structure.  If there is one mer, put
      //  it in the bucket.  If not, put a pointer to the position
      //  array there.

      //  We're in bucket b, looking at mer _sortedChck[stM].  Ask the
      //  only/mask if that exists, if so do/do not include the mer.
      //
      bool    useMer = true;

      if (edM - stM < minCount)
        useMer = false;

      if (edM - stM > maxCount)
        useMer = false;

      if ((useMer == true) && (mask || only)) {

        //  MER_REMOVAL_DURING_XFER.  Great.  The existDB has
        //  (usually) the canonical mer.  We have the forward mer.
        //  Well, no, we have the forward mers' hash and check.  So,
        //  we reconstruct the mer, reverse complement it, and then
        //  throw the mer out if either the forward or reverse exists
        //  (or doesn't exist).

        uint64 m = REBUILD(b, _sortedChck[stM]);
        uint64 r;

        if (mask) {
          if (mask->isCanonical()) {
            r = reverseComplementMer(_merSizeInBases, m);
            if (r < m)
              m = r;
          }
          if (mask->exists(m))
            useMer = false;
        }

        if (only) {
          if (only->isCanonical()) {
            r = reverseComplementMer(_merSizeInBases, m);
            if (r < m)
              m = r;
          }
          if (only->exists(m) == false)
            useMer = false;
        }
      }

      if (useMer) {
        _numberOfMers      += edM - stM;
        _numberOfPositions += edM - stM;
        _numberOfDistinct++;

        if (stM+1 == edM) {
          _numberOfUnique++;

#ifdef TEST_NASTY_BUGS
          posPtrCheck[currentBpos++] = _sortedPosn[stM];
#endif

          vals[0] = _sortedChck[stM];
          vals[1] = _sortedPosn[stM];
          vals[2] = 1;
          vals[3] = 0;

          currentBbit = setDecodedValues(_buckets, currentBbit, nval, lensF, vals);
          bucketStartPosition++;
        } else {
          _numberOfEntries  += edM - stM;
          if (_maximumEntries < edM - stM)
            _maximumEntries = edM - stM;

#ifdef TEST_NASTY_BUGS
          posPtrCheck[currentBpos++] = currentPpos;
#endif

          vals[0] = _sortedChck[stM];
          vals[1] = currentPpos;
          vals[2] = 0;
          vals[3] = 0;

          currentBbit = setDecodedValues(_buckets, currentBbit, nval, lensF, vals);
          bucketStartPosition++;

          //  Store the positions.  Store the number of positions
          //  here, then store all positions.
          //
          //  The positions are in the proper place in _sortedPosn,
          //  and setDecodedValue masks out the extra crap, so no
          //  temporary needed.  Probably should be done with
          //  setDecodedValues, but then we need another array telling
          //  the sizes of each piece.
          //
          setDecodedValue(_positions, currentPbit, _posnWidth, edM - stM);
          currentPbit += _posnWidth;
          currentPpos++;

          for (; stM < edM; stM++) {
            if (_sortedPosn[stM] >= DEBUGnumPositions) {
              fprintf(stderr, "positionDB()-- ERROR:  Got position "uint64FMT", but only "uint64FMT" available!\n",
                      _sortedPosn[stM], DEBUGnumPositions);
              abort();
            }
            setDecodedValue(_positions, currentPbit, _posnWidth, _sortedPosn[stM]);
            currentPbit += _posnWidth;
            currentPpos++;
          }
        }
      }  //  useMer

      //  All done with this mer.
      //
      stM = edM;
    }  //  while (stM < le)
  }  //  for each bucket

  //  Set the end of the last bucket
  //
  if (_hashTable_BP)
    setDecodedValue(_hashTable_BP, b * _hashWidth, _hashWidth, bucketStartPosition);
  else
    _hashTable_FW[b] = bucketStartPosition;

  delete C;

  //  Clear out the end of the arrays -- this is only so that we can
  //  checksum the result.
  //
  if (_hashTable_BP) {
    b = b * _hashWidth + _hashWidth;
    setDecodedValue(_hashTable_BP, b,           64 - (b % 64),           uint64ZERO);
  }
  setDecodedValue(_buckets,   currentBbit, 64 - (currentBbit % 64), uint64ZERO);
  setDecodedValue(_positions, currentPbit, 64 - (currentPbit % 64), uint64ZERO);


  if (beVerbose) {
    fprintf(stderr, "    Avail: Bucket "uint64FMTW(12)"    Position "uint64FMTW(12)" (64-bit words)\n", bs, ps);
    fprintf(stderr, "    Avail: Bucket "uint64FMTW(12)"    Position "uint64FMTW(12)" (entries)\n", _numberOfDistinct, _numberOfEntries);
    fprintf(stderr, "    Used:  Bucket "uint64FMTW(12)"    Position "uint64FMTW(12)" (64-bit words)\n", currentBbit / 64, currentPbit / 64);
  }

  //  Reset the sizes to what we actually found.  If we then
  //  dump/reload, we shrink our footprint.
  //
  _numberOfDistinct = currentBbit / _wFin;
  _numberOfEntries  = currentPbit / _posnWidth;

  if (beVerbose) {
    fprintf(stderr, "    Used:  Bucket "uint64FMTW(12)"    Position "uint64FMTW(12)" (entries)\n", _numberOfDistinct, _numberOfEntries);
    fprintf(stderr,
            "    Found "uint64FMTW(12)" total mers\n"
            "    Found "uint64FMTW(12)" distinct mers\n"
            "    Found "uint64FMTW(12)" unique mers\n"
            "    Need "uint64FMT" non-unique position list entries ("uint64FMT" maximum count)\n",
            _numberOfMers, _numberOfDistinct, _numberOfUnique, _numberOfEntries, _maximumEntries);
  }


  //  If we removed mers, there is a small chance that our hash table
  //  is too big -- we might have removed enoough mers to make the
  //  width smaller.  If so, rebuild the hash table.
  //
  //  Also, hooray, we finally know the number of distinct mers, so we
  //  can make this nice and tight
  //
  if (_hashTable_BP) {
    uint32 newHashWidth = 1;
    while ((_numberOfDistinct+1) > (uint64ONE << newHashWidth))
      newHashWidth++;

    if (newHashWidth != _hashWidth) {
      uint64 npos = 0;
      uint64 opos = 0;

      if (beVerbose)
        fprintf(stderr, "    Rebuilding the hash table, from "uint32FMT" bits wide to "uint32FMT" bits wide.\n",
                _hashWidth, newHashWidth);

      for (uint64 z=0; z<_tableSizeInEntries+1; z++) {
        setDecodedValue(_hashTable_BP,
                        npos,
                        newHashWidth,
                        getDecodedValue(_hashTable_BP, opos, _hashWidth));
        npos += newHashWidth;
        opos += _hashWidth;
      }

      //  Clear the end again.
      setDecodedValue(_hashTable_BP, npos, 64 - (npos % 64), uint64ZERO);
    }

    _hashWidth = newHashWidth;
  }


  //  If supplied, add in any counts.  The meryl table is, sadly, in
  //  the wrong order, and we must hash and search.
  //
  //  Meryl _should_ be storing only forward mers, but we have no way
  //  of checking.
  //
  //  After all counts are loaded, check if we can compress the counts
  //  space any.  Check if the largestMerylCount is much smaller than
  //  the space it is stored in.  If so, we can compress the table.
  //
  uint64  largestMerylCount = 0;
  uint64  countsLoaded      = 0;

  if (counts) {
    if (beVerbose)
      fprintf(stderr, "    Loading "uint64FMT" mercounts.\n", counts->numberOfDistinctMers());

    C = new speedCounter("    %7.2f Mmercounts -- %5.2f Mmercounts/second\r", 1000000.0, 0x1fffff, beVerbose);

    while (counts->nextMer()) {
      kMer    k = counts->theFMer();
      uint64  c = counts->theCount();
      uint64  f = setCount(k, c);
      k.reverseComplement();
      uint64  r = setCount(k, c);

      if (f + r > 0) {
        countsLoaded++;
        if (largestMerylCount < c)
          largestMerylCount = c;
      }

      C->tick();
    }

    delete C;

    if (beVerbose)
      fprintf(stderr, "    Loaded "uint64FMT" mercounts; largest is "uint64FMT".\n", countsLoaded, largestMerylCount);

    if (logBaseTwo64(largestMerylCount + 1) < _sizeWidth) {
      if (beVerbose)
        fprintf(stderr, "    Compress sizes from "uint32FMT" bits to "uint32FMT" bits.\n",
                _sizeWidth,
                (uint32)logBaseTwo64(largestMerylCount + 1));

      uint64 oSiz[4] = { _chckWidth, _pptrWidth, 1, _sizeWidth };
      uint64 nSiz[4] = { _chckWidth, _pptrWidth, 1, logBaseTwo64(largestMerylCount + 1) };
      uint64 tVal[4] = { 0, 0, 0, 0 };

      uint64  oP = 0, oS = oSiz[0] + oSiz[1] + oSiz[2] + oSiz[3];
      uint64  nP = 0, nS = nSiz[0] + nSiz[1] + nSiz[2] + nSiz[3];

      assert(nS < oS);

      C = new speedCounter("    %7.2f Mmercounts -- %5.2f Mmercounts/second\r", 1000000.0, 0x1fffff, beVerbose);

      for (uint64 bu=0; bu<_numberOfDistinct; bu++) {
        getDecodedValues(_buckets, oP, 4, oSiz, tVal);
        setDecodedValues(_buckets, nP, 4, nSiz, tVal);

        oP += oS;
        nP += nS;

        C->tick();
      }

      delete C;

      _sizeWidth = nSiz[3];
      _wFin      = _chckWidth + _pptrWidth + 1 + _sizeWidth;
    }
  }


#ifdef TEST_NASTY_BUGS
  //  Unpack the bucket positions and check.  Report the first one
  //  that is broken.
  //
  for(uint64 bb=0; bb<currentBpos; bb++)
    if (posPtrCheck[bb] != getDecodedValue(_buckets, bb * _wFin + _chckWidth, _pptrWidth))
      fprintf(stderr, "Bucket %lu (at bitpos %lu) failed position check (wanted %lu got %lu)\n",
              bb,
              bb * _wFin,
              posPtrCheck[bb],
              getDecodedValue(_buckets, bb * _wFin + _chckWidth, _pptrWidth));
  delete [] posPtrCheck;
#endif


#ifdef MER_REMOVAL_TEST
#warning MER_REMOVAL_TEST was not updated to deal with canonical mers
  if (beVerbose)
    fprintf(stderr, "positionDB()--     TESTING MER REMOVAL\n");

  MS->rewind();
  if (mask) {
    C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);
    uint32  extraMer   = 0;
    while (MS->nextMer(_merSkipInBases)) {
      uint64  mer = MS->theFMer();
      if (mask->exists(mer) && exists(mer))
        extraMer++;
      C->tick();
    }
    delete C;
    fprintf(stderr, "positionDB()-- mask: "uint32FMT" mers extra!\n", extraMer);
  } else if (only) {
    C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);
    uint32  missingMer = 0;
    while (MS->nextMer(_merSkipInBases)) {
      uint64  mer = MS->theFMer();
      if (only->exists(mer) && !exists(mer))
        missingMer++;
      C->tick();
    }
    delete C;
    fprintf(stderr, "positionDB()-- only: "uint32FMT" mers missing!\n", missingMer);
  }
#endif

  //  Free the counting buckets if we aren't using the space for
  //  something else.
  //
  if (_buckets != _countingBuckets)
    delete [] _countingBuckets;

  //  In theory, we could move these to be immediately after the data
  //  is useless.
  //
  _bucketSizes     = 0L;
  _countingBuckets = 0L;

  delete [] _sortedChck;
  delete [] _sortedPosn;

  _sortedMax  = 0;
  _sortedChck = 0L;
  _sortedPosn = 0L;

  if (bktAllocIsJunk)
    delete [] bktAlloc;
}

positionDB::~positionDB() {
  delete [] _hashTable_BP;
  delete [] _hashTable_FW;
  delete [] _buckets;
  delete [] _positions;
  delete [] _hashedErrors;
}
