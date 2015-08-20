
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
 *    Brian P. Walenz from 2003-JAN-02 to 2003-OCT-21
 *      are Copyright 2003 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-SEP-13 to 2014-APR-11
 *      are Copyright 2005-2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef POSITIONDB_H
#define POSITIONDB_H

#include <stdio.h>
#include "bio++.H"
#include "merStream.H"

//  The two existDB inputs can be either forward or canonical.  If
//  canonical, we are smart enough to search exist/only with the
//  canonical mer.

//  Returns position in posn, resizing it if needed.  Space is
//  allocated if none supplied.  The following is valid:
//
//    uint64  *posn    = 0L;
//    uint64   posnMax = 0;
//    uint64   posnLen = 0;
//    if (get(somemer, posn, posnMax, posnLen)) {
//      do something with the positions
//    }
//
//  exists() returns T/F if mer exists or not
//  count() returns the number of times that mer is present

//  Define this to use an uncompressed hash table when the width is 32
//  bits or less.  Doing so is A LOT faster in mismatch lookups, but
//  does use more memory.
#undef UNCOMPRESS_HASH_TABLE

//  Define this to leave out references to getTime(), speedCounter()
//  and make the positionDB build very quietly.
#undef SILENTPOSITIONDB

//  Define these to enable some debugging methods
#undef DEBUGPOSDB
#undef DEBUGREBUILD

class existDB;
class merylStreamReader;

class positionDB {
public:
  positionDB(char const        *filename,
             uint32             merSize,
             uint32             merSkip,
             uint32             maxMismatch,
             bool               loadData=true);

  positionDB(merStream         *MS,
             uint32             merSize,
             uint32             merSkip,
             existDB           *mask,
             existDB           *only,
             merylStreamReader *counts,
             uint32             minCount,
             uint32             maxCount,
             uint32             maxMismatch,
             uint32             maxMemory,
             bool               beVerbose);

  ~positionDB();

private:
  void  build(merStream         *MS,
              existDB           *mask,
              existDB           *only,
              merylStreamReader *counts,
              uint32             minCount,
              uint32             maxCount,
              bool               beVerbose);

private:
  void        reallocateSpace(uint64*&    posn,
                              uint64&     posnMax,
                              uint64&     posnLen,
                              uint64      len);

  void        loadPositions(uint64      v,
                            uint64*&    posn,
                            uint64&     posnMax,
                            uint64&     posnLen,
                            uint64&     count);

public:
  bool        getExact(uint64      mer,
                       uint64*&    posn,
                       uint64&     posnMax,
                       uint64&     posnLen,
                       uint64&     count);
  bool        existsExact(uint64   mer);
  uint64      countExact(uint64    mer);

public:
  void        filter(uint64 lo, uint64 hi);

private:
  double      setUpMismatchMatcher(uint32 nErrorsAllowed, uint64 approxMers);
public:
  bool        getUpToNMismatches(uint64      mer,
                                 uint32      maxMismatches,
                                 uint64*&    posn,
                                 uint64&     posnMax,
                                 uint64&     posnLen);
private:
  uint64      setCount(uint64 mer, uint64 count);

  //  Save or load a built table
  //
public:
  void        saveState(char const *filename);
  bool        loadState(char const *filename, bool beNoisy=false, bool loadData=true);

  void        printState(FILE *stream);

  //  Only really useful for debugging.  Don't use.
  //
  void        dump(char *name);


  bool         checkREBUILD(uint64 m) {
#define DEBUGREBUILD
#ifdef DEBUGREBUILD
    uint64 h = HASH(m);
    uint64 c = CHECK(m);
    uint64 r = REBUILD(h, c);
    if (r != m) {
      fprintf(stderr, "shift1 = "uint32FMT"\n", _shift1);
      fprintf(stderr, "shift2 = "uint32FMT"\n", _shift2);
      fprintf(stderr, "M = "uint64HEX"\n", m);
      fprintf(stderr, "H = "uint64HEX"\n", h);
      fprintf(stderr, "C = "uint64HEX"\n", c);
      fprintf(stderr, "R = "uint64HEX"\n", r);
      return(false);
    }
    return(true);
#else
    return(REBUILD(HASH(m), CHECK(m)) == m);
#endif
  };

private:

  uint64       HASH(uint64 k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  uint64       CHECK(uint64 k) {
    return(k & _mask2);
  };

  uint64       REBUILD(uint64 h, uint64 c) {
    //  Decode a HASH and a CHECK to get back the mer.  You'd better
    //  bloody PRAY you don't break this (test/test-rebuild.C).  It
    //  was a headache++ to write.

    uint64 sha = _shift1 - _shift2;
    uint64 msk = uint64MASK(sha);

    //  The check is exactly the mer....just not all there.
    uint64 mer = c;

    uint64 shf = sha - (_tableSizeInBits % 2);
    uint64 shg = 0;
    uint64 shh = _shift1;

    //  Unrolling this is troublesome - we still need the tests,
    //  bizarre merSize, tblSize combinations use lots of iterations
    //  (when the merSize and tblSize are about the same, the CHECK is
    //  small, and so we need to do lots of iterations).

    //fprintf(stderr, "shf="uint64FMTW(2)" shg="uint64FMTW(2)" shh="uint64FMTW(2)" mer="uint64HEX"\n", shf, shg, shh, mer);

    do {
      mer |= (((h >> shg) ^ (mer >> shg) ^ (mer >> shf)) & msk) << shh;
      //fprintf(stderr, "shf="uint64FMTW(2)" shg="uint64FMTW(2)" shh="uint64FMTW(2)" mer="uint64HEX"\n", shf, shg, shh, mer);

      shf += sha;
      shg += sha;
      shh += sha;
    } while ((shf < _merSizeInBits) && (shh < 64));

    mer &= uint64MASK(_merSizeInBits);

    return(mer);
  };

  void         sortAndRepackBucket(uint64 b);

  uint32     *_bucketSizes;
  uint64     *_countingBuckets;
  uint64     *_hashTable_BP;  //  Bit packed
  uint32     *_hashTable_FW;  //  Full width
  uint64     *_buckets;

  uint64     *_positions;

  uint32      _merSizeInBases;
  uint32      _merSizeInBits;

  uint32      _merSkipInBases;

  uint64      _tableSizeInEntries;
  uint32      _tableSizeInBits;

  uint32      _hashWidth;  // Hash bith
  uint32      _chckWidth;  // Check bits
  uint32      _posnWidth;  // Positions in the sequence
  uint32      _pptrWidth;  // Pointers to positions
  uint32      _sizeWidth;  // Extra number in the table

  uint64      _hashMask;

  uint32      _wCnt;
  uint32      _wFin;

  uint32      _shift1;
  uint32      _shift2;
  uint64      _mask1;
  uint64      _mask2;

  uint64      _numberOfMers;
  uint64      _numberOfPositions;
  uint64      _numberOfDistinct;
  uint64      _numberOfUnique;
  uint64      _numberOfEntries;
  uint64      _maximumEntries;

  //  For sorting the mers
  //
  uint32      _sortedMax;
  uint64     *_sortedChck;
  uint64     *_sortedPosn;

  //  For the mismatch matcher
  uint32      _nErrorsAllowed;
  uint32      _hashedErrorsLen;
  uint32      _hashedErrorsMax;
  uint64     *_hashedErrors;
};

#endif  //  POSITIONDB_H
