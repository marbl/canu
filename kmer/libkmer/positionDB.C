#include <stdio.h>
#include <stdlib.h>
#include <new>

#include "positionDB.H"
#include "existDB.H"
#include "bri++.H"

#define REUSE_BUCKETS

//#define ERROR_CHECK_COUNTING
//#define ERROR_CHECK_COUNTING_ENCODING
//#define ERROR_CHECK_EMPTY_BUCKETS


//  Progress messages disabled if we're silent.
//
#ifndef SILENTPOSITIONDB

#define MSG_OUTPUT  stderr

#ifdef TRUE64BIT
const char *buckCountMsg = "    Allocated bucket size counting space with total size %lu KB\n";
const char *foundMersMsg = "\n    Found %lu mers\n";
const char *cbktAllocMsg = "    Allocated %12lu KB for buckets (%lu 64-bit words)\n";
const char *mersStatsMsg = "    Found %10lu total mers\n"
                           "    Found %10lu distinct mers\n"
                           "    Found %10lu unique mers\n"
                           "    Need  %10lu non-unique position list entries (%lu maximum count)\n";
const char *hashAllocMsg = "    Allocated %12luKB for hash table (%lu 64-bit words)\n";
const char *buckAllocMsg = "    Allocated %12luKB for buckets    (%lu 64-bit words)\n";
const char *posnAllocMsg = "    Allocated %12luKB for positions  (%lu 64-bit words)\n";
const char *buckReuseMsg = "    Reusing bucket space; Have: %lu  Need: %lu (64-bit words)\n";
const char *spceAvalMsg  = "    Avail: Bucket %12lu    Position %12lu (64-bit words)\n";
const char *spceUsedMsg  = "    Used:  Bucket %12lu    Position %12lu (64-bit words)\n";
#else
const char *buckCountMsg = "    Allocated bucket size counting space with total size %llu KB\n";
const char *foundMersMsg = "\n    Found %llu mers\n";
const char *cbktAllocMsg = "    Allocated %12llu KB for buckets (%llu 64-bit words)\n";
const char *mersStatsMsg = "    Found %10llu total mers\n"
                           "    Found %10llu distinct mers\n"
                           "    Found %10llu unique mers\n"
                           "    Need  %10llu non-unique position list entries (%llu maximum count)\n";
const char *hashAllocMsg = "    Allocated %12lluKB for hash table (%llu 64-bit words)\n";
const char *buckAllocMsg = "    Allocated %12lluKB for buckets    (%llu 64-bit words)\n";
const char *posnAllocMsg = "    Allocated %12lluKB for positions  (%llu 64-bit words)\n";
const char *buckReuseMsg = "    Reusing bucket space; Have: %llu  Need: %llu (64-bit words)\n";
const char *spceAvalMsg  = "    Avail: Bucket %12llu    Position %12llu (64-bit words)\n";
const char *spceUsedMsg  = "    Used:  Bucket %12llu    Position %12llu (64-bit words)\n";
#endif

#endif  //  SILENTPOSITIONDB


//  Error messages are always present.
//
#ifdef TRUE64BIT
const char *cntSizeMsg   = "ERROR: wCnt=%u and wFin=%u (should be 64 or less).\n";
const char *buckReuseErr = "ERROR:  I didn't allocate enough space to reuse buckets!\n"
                           "        Have: %lu  Need: %lu (64-bit words)\n";
#else
const char *cntSizeMsg   = "ERROR: data sizes to big: wCnt=%lu and wFin=%lu.\n";
const char *buckReuseErr = "ERROR:  I didn't allocate enough space to reuse buckets!\n"
                           "        Have: %llu  Need: %llu (64-bit words)\n";
#endif  //  TRUE64BIT



positionDB::positionDB(char const    *filename,
                       bool           loadData) {
  if (loadState(filename, true, loadData) == false) {
    fprintf(stderr, "positionDB::positionDB()-- Tried to read state from '%s', but failed.\n", filename);
    exit(1);
  }
}


positionDB::positionDB(char const  *seq,
                       char const  *filename,
                       u32bit       merSize,
                       u32bit       merSkip,
                       u32bit       tblBits,
                       bool         beVerbose,
                       existDB       *mask,
                       existDB       *only,
                       char const    *streamFileName) {
  _bucketSizes           = 0L;
  _countingBuckets       = 0L;
  _hashTable             = 0L;
  _buckets               = 0L;
  _positions             = 0L;

  _merSizeInBases        = merSize;
  _merSizeInBits         = 2 * _merSizeInBases;
  _merSkipInBases        = merSkip;
  _tableSizeInBits       = tblBits;
  _tableSizeInEntries    = u64bitONE << _tableSizeInBits;
  _hashWidth             = u32bitZERO;
  _hashMask              = u64bitMASK(_tableSizeInBits);
  _chckWidth             = _merSizeInBits - _tableSizeInBits;
  _chckMask              = u64bitMASK(_chckWidth);
  _posnWidth             = u64bitZERO;
  _posnMask              = u64bitZERO;

  _wCnt                  = 0;
  _wFin                  = 0;

  _shift1                = _merSizeInBits - _tableSizeInBits;
  _shift2                = _shift1 / 2;
  _mask1                 = u64bitMASK(_tableSizeInBits);
  _mask2                 = u64bitMASK(_shift1);

  _numberOfMers          = u64bitZERO;
  _numberOfPositions     = u64bitZERO;
  _numberOfDistinct      = u64bitZERO;
  _numberOfUnique        = u64bitZERO;
  _numberOfEntries       = u64bitZERO;
  _maximumEntries        = u64bitZERO;

  _sortedListMax         = u32bitZERO;
  _sortedListLen         = u32bitZERO;
  _sortedList            = 0L;


  if (filename && loadState(filename)) {
    bool fail = false;

    if (_merSizeInBases != merSize) {
      fprintf(stderr, "positionDB::positionDB()-- Read state from '%s', but got different mer sizes\n", filename);
      fprintf(stderr, "positionDB::positionDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", _merSizeInBases, merSize);
      fail = true;
    }
    if (_merSkipInBases != merSkip) {
      fprintf(stderr, "positionDB::positionDB()-- Read state from '%s', but got different mer skips\n", filename);
      fprintf(stderr, "positionDB::positionDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", _merSkipInBases, merSkip);
      fail = true;
    }
    if (_tableSizeInBits != tblBits) {
      fprintf(stderr, "positionDB::positionDB()-- Read state from '%s', but got different table sizes\n", filename);
      fprintf(stderr, "positionDB::positionDB()-- Got "u32bitFMT", expected "u32bitFMT"\n", _tableSizeInBits, tblBits);
      fail = true;
    }

    if (fail)
      exit(1);

    return;
  }

  if (!filename && !seq) {
    fprintf(stderr, "positionDB()-- ERROR: No sequence and no filename?  Nothing to build a table with!\n");
    fprintf(stderr, "positionDB()--        The code is probably broken.\n");
    exit(1);
  }

  if ((mask || only) && !streamFileName) {
    fprintf(stderr, "positionDB()-- WARNING: You supplied a mask/only existDB but no temporary file!  The\n");
    fprintf(stderr, "positionDB()--          temporary file is required for masking to occur!\n");
  }

  if (!(mask || only) && streamFileName) {
    fprintf(stderr, "positionDB()-- WARNING: You supplied a temporary file, but no mask/only existDB!  The\n");
    fprintf(stderr, "positionDB()--          temporary file is required only for masking.  I'll disable it,\n");
    fprintf(stderr, "positionDB()--          but you should fix your command line.\n");
  }


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  1)  Count bucket sizes
  //

  //  We'll later want to reuse the _bucketSizes space for storing
  //  the _hashTable.  To make it somewhat safe, we allocate the
  //  space as u64bit, then cast it to be u32bit.
  //
  //  bktAllocIsJunk tells us if we should release this memory (if
  //  we need to allocate separate space for the _hashTable).
  //
  //  The _bucketSizes is offset by one from bktAlloc so that
  //  we don't overwrite _bucketSizes when we are constructing
  //  _hashTable.
  //
  u64bit *bktAlloc;
  try {
    bktAlloc = new u64bit [_tableSizeInEntries / 2 + 2];
  } catch (std::bad_alloc) {
    fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "hitMatrix::filter()-- bktAlloc = new u64bit ["u64bitFMT"]\n", _tableSizeInEntries / 2 + 2);
    exit(1);
  }
  bool     bktAllocIsJunk = false;

  for (u64bit i=_tableSizeInEntries / 2 + 2; i--; )
    bktAlloc[i] = 0;

  _bucketSizes = (u32bit *)(bktAlloc + 1);

#ifdef ERROR_CHECK_COUNTING
  fprintf(stdout, "ERROR_CHECK_COUNTING is defined.\n");
  u32bit *_errbucketSizes = new u32bit [_tableSizeInEntries + 2];
  for (u64bit i=_tableSizeInEntries + 2; i--; )
    _errbucketSizes[i] = 0;
#endif

#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, buckCountMsg, _tableSizeInEntries >> 8);
#endif

  merStream     *M = 0L;

  if (seq)
    M = new merStream(_merSizeInBases, seq, 0);

  if (filename)
    M = new merStream(_merSizeInBases, filename);

  if (M == 0L) {
    fprintf(stderr, "ERROR:  Nothing to initialize with!\n");
    exit(1);
  }

#ifndef SILENTPOSITIONDB
  speedCounter  *C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);
#endif

  //  Two choices here
  //
  //  1)  No masking or onlying is done.  Stream the mers and just
  //      count the positions.  This is the original behavior.
  //
  //  2)  Masking or onlyint is done.  Open the output stream file,
  //      stream the mers by, checking for mask/only of both
  //      forward and reverse mers.  If either is found, push
  //      the (forward) mer and position onto the stream.
  //      close the output stream.
  //
  if (streamFileName == 0L) {
    while (M->nextMer(_merSkipInBases)) {
      _bucketSizes[ HASH(M->theFMer()) ]++;

#ifdef ERROR_CHECK_COUNTING
      _errbucketSizes[ HASH(M->theFMer()) ]++;
#endif

      _numberOfMers++;
      _numberOfPositions = M->thePosition();
#ifndef SILENTPOSITIONDB
      C->tick();
#endif
    }
  } else {

    //  XXXX:  HACK THIS!

    bitPackedFile *ms = new bitPackedFile(streamFileName);

    while (M->nextMer(_merSkipInBases)) {

      //  Save the mer if it doesn't exist in the mask (both f and r),
      //  or does exist in the only (either f or r), add it.
      //
      //  The input databases for mask and only are (currently) made
      //  using canonical mers.  We halve the number of exists() by
      //  also using canonical mers here.
      //
      bool  save = false;

      u64bit  canonicalmer = M->theFMer();
      if (canonicalmer > M->theRMer())
        canonicalmer = M->theRMer();

      if (mask && !mask->exists(canonicalmer))
        save = true;
      if (only &&  only->exists(canonicalmer))
        save = true;

      if (save) {
        _bucketSizes[ HASH(M->theFMer()) ]++;
        _numberOfMers++;
        _numberOfPositions = M->thePosition();

        //  We probably don't want to use putNumber instead of putBits
        //  because our distribution of number is uniform up to
        //  (usually) 31 bits.

        ms->putBits(M->theFMer(),     _merSizeInBits);
        ms->putBits(M->thePosition(), 32);
      }

#ifndef SILENTPOSITIONDB
      C->tick();
#endif
    }

    delete ms;
  }

  delete M;
#ifndef SILENTPOSITIONDB
  delete C;
#endif

  M = 0L;
  C = 0L;

#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, foundMersMsg, _numberOfMers);
#endif

  //  This is _numberOfMers+1 because we need to store the first
  //  position after the last mer.  That is, if there are two mers, we
  //  will store that the first mer is at position 0, the second mer
  //  is at position 1, and the end of the second mer is at position 2.
  //
  _hashWidth = 1;
  while ((_numberOfMers+1) > (u64bitONE << _hashWidth))
    _hashWidth++;

  _posnWidth = 1;
  while ((_numberOfPositions+1) > (u64bitONE << _posnWidth))
    _posnWidth++;
  _posnMask = u64bitMASK(_posnWidth);


#if 0
  fprintf(stderr, "hashWidth = %d\n", _hashWidth);
  fprintf(stderr, "chckWidth = %d\n", _chckWidth);
  fprintf(stderr, "posnWidth = %d\n", _posnWidth);
#endif


  ///////////////////////////////////////////////////////////////////////////////
  //
  //  2)  Allocate buckets and make bucketSizes be a pointer into them
  //
#ifdef REUSE_BUCKETS
  _wCnt          = _chckWidth + _posnWidth + 1;
  _wFin          = _chckWidth + _posnWidth + 1;
#else
  _wCnt          = _chckWidth + _posnWidth;
  _wFin          = _chckWidth + _posnWidth + 1;
#endif
  u64bit   bucketsSpace  = (_numberOfMers+1) * _wCnt / 64 + 1;
  u32bit   endPosition   = 0;

  if ((_wCnt > 64) || (_wFin > 64)) {
    fprintf(stderr, "ERROR: data size too big.  Reduce mersize, number of mers, or both.\n");
    fprintf(stderr, cntSizeMsg, _wCnt, _wFin);
    exit(1);
  }

#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, cbktAllocMsg, bucketsSpace >> 7, bucketsSpace);
#endif
  try {
    _countingBuckets = new u64bit [bucketsSpace];
  } catch (std::bad_alloc) {
    fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "hitMatrix::filter()-- _countingBuckets = new u64bit ["u64bitFMT"]\n", bucketsSpace);
    exit(1);
  }

  for (u64bit i=bucketsSpace; i--; )
    _countingBuckets[i] = ~u64bitZERO;

  for (u32bit i=0; i<_tableSizeInEntries; i++) {
    endPosition     += _bucketSizes[i];
    _bucketSizes[i]  = endPosition;
  }
  _bucketSizes[_tableSizeInEntries] = endPosition;

#ifdef ERROR_CHECK_COUNTING
  fprintf(stdout, "ERROR_CHECK_COUNTING:  endposition = "u32bitFMT"\n", endPosition);
  fprintf(stdout, "ERROR_CHECK_COUNTING:  nummers     = "u64bitFMT"\n", _numberOfMers);
  if (endPosition != _numberOfMers)
    fprintf(stdout, "ERROR_CHECK_COUNTING: BUCKETSIZE COUNTING PROBLEM -- endPos != numMers\n");
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  3)  Build list of mers with positions
  //
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, "    Building lists with positions.\n");
#endif

#ifndef SILENTPOSITIONDB
  C = new speedCounter("    %7.2f Mmers -- %5.2f Mmers/second\r", 1000000.0, 0x1fffff, beVerbose);
#endif

#ifdef ERROR_CHECK_COUNTING_ENCODING
  fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING is defined!\n");
#endif


  //  XXX:  Two choices again
  //
  //  If no file of pre-encoded mers, read the mers from a merStream.
  //  Otherwise, read the pre-encoded mers.
  //

  if (streamFileName == 0L) {
    if (seq)
      M = new merStream(_merSizeInBases, seq, 0);

    if (filename)
      M = new merStream(_merSizeInBases, filename);


    while (M->nextMer(_merSkipInBases)) {
      u64bit h = HASH(M->theFMer());

      _bucketSizes[h]--;

#ifdef ERROR_CHECK_COUNTING
      _errbucketSizes[h]--;
#endif


#ifdef ERROR_CHECK_EMPTY_BUCKETS
      if ((~getDecodedValue(_countingBuckets, (u64bit)_bucketSizes[h] * (u64bit)_wCnt, _wCnt)) & u64bitMASK(_wCnt))
        fprintf(stdout, "ERROR_CHECK_EMPTY_BUCKETS: countingBucket not empty!  pos=%lu\n", _bucketSizes[h] * _wCnt);
#endif

      setDecodedValue(_countingBuckets, (u64bit)_bucketSizes[h] * (u64bit)_wCnt, _wCnt,
                      (CHECK(M->theFMer()) << _posnWidth) | (M->thePosition() & _posnMask));


#ifdef ERROR_CHECK_COUNTING_ENCODING
      u64bit v = getDecodedValue(_countingBuckets, (u64bit)_bucketSizes[h] * (u64bit)_wCnt, _wCnt);

      //  This test is only valid if we have an extra bit at the start --
      //  if we are planning on reusing the counting space for buckets.
      //
      if ((_wCnt == _wFin) && (0 != (v >> (_wCnt - 1))))
        fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error: HBIT is set!      Wanted "u64bitHEX" got "u64bitHEX"\n",
                (CHECK(M->theFMer()) << _posnWidth) | (M->thePosition() & _posnMask), v);
      if (CHECK(M->theFMer()) != ((v >> _posnWidth) & _chckMask))
        fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  CHCK corrupted!  Wanted "u64bitHEX" got "u64bitHEX"\n",
                (CHECK(M->theFMer()) << _posnWidth) | (M->thePosition() & _posnMask), v);
      if (M->thePosition() != (v & _posnMask))
        fprintf(stdout, "ERROR_CHECK_COUNTING_ENCODING error:  POSN corrupted!  Wanted "u64bitHEX" got "u64bitHEX"\n",
                (CHECK(M->theFMer()) << _posnWidth) | (M->thePosition() & _posnMask), v);
#endif

#ifndef SILENTPOSITIONDB
      C->tick();
#endif
    }

    delete M;
    M = 0L;
  } else {
    bitPackedFile *ms = new bitPackedFile(streamFileName);

    for (u64bit i=0; i<_numberOfMers; i++) {
      u64bit  mer = ms->getBits(_merSizeInBits);
      u64bit  pos = ms->getBits(32);
      
      u64bit  h = HASH(mer);
        
      _bucketSizes[h]--;

      setDecodedValue(_countingBuckets, (u64bit)_bucketSizes[h] * (u64bit)_wCnt, _wCnt,
                      (CHECK(mer) << _posnWidth) | (pos & _posnMask));

#ifndef SILENTPOSITIONDB
      C->tick();
#endif
    }

    delete ms;
  }


#ifndef SILENTPOSITIONDB
  delete C;
#endif
  C = 0L;

#ifdef ERROR_CHECK_COUNTING
  fprintf(stdout, "Checking for unfilled buckets\n");
  for (u32bit i=0; i<_tableSizeInEntries; i++)
    if (_errbucketSizes[i] != 0)
      fprintf(stdout, "ERROR_CHECK_COUNTING: Bucket "u32bitFMT" wasn't filled fully?  "u64bitFMT" left over.\n", i, _errbucketSizes[i]);
  fprintf(stdout, "ERROR_CHECK_COUNTING\n");
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  4)  Sort each bucket -- count:
  //        1) number of distinct mers
  //        2) number of unique mers
  //        3) number of entries in position table ( sum mercount+1 for all mercounts > 1)
  //      also need to repack the sorted things
  //
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, "\n    Sorting and repacking buckets.\n");
#endif

  for (u32bit i=0; i<_tableSizeInEntries; i++)
    sortAndRepackBucket(i);

#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, mersStatsMsg, _numberOfMers, _numberOfDistinct, _numberOfUnique, _numberOfEntries, _maximumEntries);
#endif


  ////////////////////////////////////////////////////////////////////////////////
  //
  //  5)  Allocate: real hash table, buckets and position table.
  //
  u64bit  hs = _tableSizeInEntries * _hashWidth / 64 + 1;
  u64bit  bs = _numberOfDistinct   * _wFin      / 64 + 1;
  u64bit  ps = _numberOfEntries    * _posnWidth / 64 + 1;

  if (_hashWidth < 32) {
#ifndef SILENTPOSITIONDB
    if (beVerbose)
      fprintf(MSG_OUTPUT, "    Reusing bucket counting space for hash table.\n");
#endif

    _hashTable     = bktAlloc;
    bktAllocIsJunk = false;
  } else {
#ifndef SILENTPOSITIONDB
    if (beVerbose)
      fprintf(MSG_OUTPUT, hashAllocMsg, hs >> 7, hs);
#endif

    try {
      _hashTable     = new u64bit [hs];
    } catch (std::bad_alloc) {
      fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
      fprintf(stderr, "hitMatrix::filter()-- _hashTable = new u64bit ["u64bitFMT"]\n", hs);
      exit(1);
    }
    bktAllocIsJunk = true;
  }

#ifdef REUSE_BUCKETS
  if (bucketsSpace < bs) {
    fprintf(stderr, buckReuseErr, bucketsSpace, bs);
    exit(1);
  }
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, buckReuseMsg, bucketsSpace, bs);
#endif

  _buckets   = _countingBuckets;

  bs = bucketsSpace; // for output at the end
#else
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, buckAllocMsg, bs >> 7, bs);
#endif
  try {
    _buckets   = new u64bit [bs];
  } catch (std::bad_alloc) {
    fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "hitMatrix::filter()-- _buckets = new u64bit ["u64bitFMT"]\n", bs);
    exit(1);
  }
#endif

#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, posnAllocMsg, ps >> 7, ps);
#endif
  try {
    _positions = new u64bit [ps];
  } catch (std::bad_alloc) {
    fprintf(stderr, "hitMatrix::filter()-- caught std::bad_alloc in %s at line %d\n", __FILE__, __LINE__);
    fprintf(stderr, "hitMatrix::filter()-- _positions = new u64bit ["u64bitFMT"\n", ps);
    exit(1);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  6)  Transfer from the sorted buckets to the hash table.
  //
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, "    Transferring to final structure.\n");
#endif

  u64bit   bucketStartPosition = 0;
  u64bit   checkMask = ~_posnMask;

  //  Current positions and bit positions in the buckets and position list.
  //
  u64bit  currentBbit = u64bitZERO;  //  Bit position into bucket
  u64bit  currentPbit = u64bitZERO;  //  Bit position into positions
  u64bit  currentPpos = u64bitZERO;  //  Value position into positions

  //  We need b outside the loop!
  //
  u32bit  b;
  for (b=0; b<_tableSizeInEntries; b++) {

    //  Set the start of the bucket -- we took pains to ensure that
    //  we don't overwrite _bucketSizes[b], if we are reusing that
    //  space for _hashTable.
    //
    setDecodedValue(_hashTable, (u64bit)b * (u64bit)_hashWidth, _hashWidth, bucketStartPosition);

    //  Get the number of mers in the counting bucket.  The error checking
    //  and sizing of _sortedList was already done in the sort.
    //
    u64bit st = _bucketSizes[b];
    u64bit ed = _bucketSizes[b+1];
    u64bit le = ed - st;

    //  NOTE:  If all the mers are unique, _sortedList will not be allocated
    //  in sortAndRepackBucket().  In that case, we need to allocate some space.
    //
    if (_sortedList == 0L) {
      _sortedListMax = 1024;
      _sortedListLen = u32bitZERO;
      _sortedList    = new heapbit [_sortedListMax];
    }

    //  Unpack the check values
    //
    for (u64bit i=st, J=st * _wCnt; i<ed; i++, J += _wCnt)
      _sortedList[i-st] = getDecodedValue(_countingBuckets, J, _wCnt);

    //
    //  Walk through the counting bucket, adding things to the real
    //  bucket as we see them.  Mers with more than one position are
    //  inserted into the bucket, and the positions inserted into the
    //  position list.
    //

    //  start and end locations of the mer.  For mers with only
    //  one occurrance (unique mers), stM+1 == edM.
    //
    u64bit  stM = u64bitZERO;
    u64bit  edM = u64bitZERO;
    u64bit  v;

    while (stM < le) {

      //  Move to the next mer.
      //
      edM++;

      //  Keep moving while the two mers are the same.
      //
      while ((edM < le) &&
             ((_sortedList[stM] & checkMask) == (_sortedList[edM] & checkMask)))
        edM++;

      //  edM is now the mer after the last.  Write all mers from stM up to edM
      //  to the final structure.  If there is one mer, put it in the bucket.
      //  If not, put a pointer to the position array there.
      //
      if (stM+1 == edM) {

        //  Rearrange the check and position, insert a bit telling us that
        //  this is a position and not an offset.
        //
        v  = u64bitONE << (_wFin - 1);
        v |= (_sortedList[stM] & _posnMask) << _chckWidth;
        v |= (_sortedList[stM] & checkMask) >> _posnWidth;

        setDecodedValue(_buckets, currentBbit, _wFin, v);
        currentBbit += _wFin;
        bucketStartPosition++;
      } else {

        //  Create a pointer to the position list.
        //
        v  = u64bitZERO;
        v |= (currentPpos      & _posnMask) << _chckWidth;
        v |= (_sortedList[stM] & checkMask) >> _posnWidth;

        setDecodedValue(_buckets, currentBbit, _wFin, v);
        currentBbit += _wFin;
        bucketStartPosition++;

        //  Store the positions.
        //
        setDecodedValue(_positions, currentPbit, _posnWidth, edM - stM);
        currentPbit += _posnWidth;
        currentPpos++;

        //  The positions are in the proper place in _sortedList, and setDecodedValue
        //  masks out the extra crap, so no temporary needed.
        //
        for (; stM < edM; stM++) {
          setDecodedValue(_positions, currentPbit, _posnWidth, _sortedList[stM]);
          currentPbit += _posnWidth;
          currentPpos++;
        }
      }

      //  All done with this mer.
      //
      stM = edM;
    }
  }

  //  Set the end of the last bucket
  //
  setDecodedValue(_hashTable, b * _hashWidth, _hashWidth, bucketStartPosition);

#ifndef SILENTPOSITIONDB
  if (beVerbose) {
    fprintf(MSG_OUTPUT, spceAvalMsg, bs, ps);
    fprintf(MSG_OUTPUT, spceUsedMsg, currentBbit / 64, currentPbit / 64);
  }
#endif

  //  Clean up our temporary tables
  //
  //  Make sure we don't delete the _hashTable if we are reusing the
  //  space allocated by _bucketSizes!
  //
#ifndef REUSE_BUCKETS
#ifndef SILENTPOSITIONDB
  if (beVerbose)
    fprintf(MSG_OUTPUT, "    Deleting counting buckets.\n");
#endif
  delete [] _countingBuckets;
#endif

  //  These aren't immediately after the data is useless, but we know
  //  the data is useless here.
  //
  _bucketSizes     = 0L;
  _countingBuckets = 0L;

  delete [] _sortedList;

  _sortedListMax = 0;
  _sortedListLen = 0;
  _sortedList    = 0L;

  if (bktAllocIsJunk)
    delete [] bktAlloc;
}

positionDB::~positionDB() {
  delete [] _hashTable;
  delete [] _buckets;
  delete [] _positions;
}

