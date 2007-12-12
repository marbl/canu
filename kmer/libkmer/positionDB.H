#ifndef POSITIONDB_H
#define POSITIONDB_H

#include <stdio.h>
#include "bio++.H"

//  The two existDB inputs can be either forward or canonical.  If
//  canonical, we are smart enough to search exist/only with the
//  canonical mer.

//  Returns position in posn, resizing it if needed.  Space is
//  allocated if none supplied.  The following is valid:
//
//    u64bit  *posn    = 0L;
//    u64bit   posnMax = 0;
//    u64bit   posnLen = 0;
//    if (get(somemer, posn, posnMax, posnLen)) {
//      do something with the positions
//    }
//
//  exists() returns T/F if mer exists or not
//  count() returns the number of times that mer is present

//  Define this to leave out references to getTime(), speedCounter()
//  and make the positionDB build very quietly.
#undef SILENTPOSITIONDB

//  Define these to enable some debugging methods
#undef DEBUGPOSDB
#undef DEBUGREBUILD

class existDB;

class positionDB {
public:
  positionDB(char const    *filename,
             bool           loadData=true);

  positionDB(merStream     *MS,
             u32bit         merSize,
             u32bit         merSkip,
             existDB       *mask,
             existDB       *only,
             u32bit         maxCount,
             bool           beVerbose,
             bool           isForMismatches = false);

  ~positionDB();

private:
  void  build(merStream   *MS,
              existDB     *mask,
              existDB     *only,
              u32bit       maxCount,
              bool         beVerbose);

private:
  void        loadPositions(u64bit      v,
                            u64bit*&    posn,
                            u64bit&     posnMax,
                            u64bit&     posnLen);

public:
  bool        get(u64bit      mer,
                  u64bit*&    posn,
                  u64bit&     posnMax,
                  u64bit&     posnLen);
  bool        getMismatch(u64bit      mer,
                          u64bit      numMismatches,
                          u64bit*&    posn,
                          u64bit&     posnMax,
                          u64bit&     posnLen);
  bool        exists(u64bit   mer);
  u64bit      count(u64bit    mer);

  //  Save or load a built table
  //
  void        saveState(char const *filename);
  bool        loadState(char const *filename, bool beNoisy=false, bool loadData=true);

  void        printState(FILE *stream);

  //  Only really useful for debugging.  Don't use.
  //
#ifdef DEBUGPOSDB
  void        dump(u64bit mer);
  void        dumpTable(void);
#endif


  bool         checkREBUILD(u64bit m) {
#ifdef DEBUGREBUILD
    u64bit h = HASH(m);
    u64bit c = CHECK(m);
    u64bit r = REBUILD(h, c);
    if (r != m) {
      fprintf(stderr, "shift1 = "u32bitFMT"\n", _shift1);
      fprintf(stderr, "shift2 = "u32bitFMT"\n", _shift2);
      fprintf(stderr, "M = "u64bitHEX"\n", m);
      fprintf(stderr, "H = "u64bitHEX"\n", h);
      fprintf(stderr, "C = "u64bitHEX"\n", c);
      fprintf(stderr, "R = "u64bitHEX"\n", r);
      return(false);
    }
    return(true);
#else
    return(REBUILD(HASH(m), CHECK(m)) == m);
#endif
  };

private:

  u64bit       HASH(u64bit k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  u64bit       CHECK(u64bit k) {
    return(k & _mask2);
  };

  u64bit       REBUILD(u64bit h, u64bit c) {
    //  Decode a HASH and a CHECK to get back the mer.  You'd better
    //  bloody PRAY you don't break this (test/test-rebuild.C).  It
    //  was a headache++ to write.

    u64bit sha = _shift1 - _shift2;
    u64bit msk = u64bitMASK(sha);

    //  The check is exactly the mer....just not all there.
    u64bit mer = c;

    u64bit shf = sha - (_tableSizeInBits % 2);
    u64bit shg = 0;
    u64bit shh = _shift1;

    //  Unrolling this is troublesome - we still need the tests,
    //  bizarre merSize, tblSize combinations use lots of iterations
    //  (when the merSize and tblSize are about the same, the CHECK is
    //  small, and so we need to do lots of iterations).

    //fprintf(stderr, "shf="u64bitFMTW(2)" shg="u64bitFMTW(2)" shh="u64bitFMTW(2)" mer="u64bitHEX"\n", shf, shg, shh, mer);

    do {
      mer |= (((h >> shg) ^ (mer >> shg) ^ (mer >> shf)) & msk) << shh;
      //fprintf(stderr, "shf="u64bitFMTW(2)" shg="u64bitFMTW(2)" shh="u64bitFMTW(2)" mer="u64bitHEX"\n", shf, shg, shh, mer);

      shf += sha;
      shg += sha;
      shh += sha;
    } while ((shf < _merSizeInBits) && (shh < 64));

    mer &= u64bitMASK(_merSizeInBits);

    return(mer);
  };

  void         sortAndRepackBucket(u64bit b);

  u32bit     *_bucketSizes;
  u64bit     *_countingBuckets;
  u64bit     *_hashTable;
  u64bit     *_buckets;

  u64bit     *_positions;

  u32bit      _merSizeInBases;
  u32bit      _merSizeInBits;

  u32bit      _merSkipInBases;

  u64bit      _tableSizeInEntries;
  u32bit      _tableSizeInBits;

  u32bit      _hashWidth;
  u32bit      _chckWidth;
  u32bit      _posnWidth;
  u32bit      _posnPtrWidth;
  u32bit      _sizeWidth;

  u64bit      _hashMask;

  u32bit      _wCnt;
  u32bit      _wFin;

  u32bit      _shift1;
  u32bit      _shift2;
  u64bit      _mask1;
  u64bit      _mask2;

  u64bit      _numberOfMers;
  u64bit      _numberOfPositions;
  u64bit      _numberOfDistinct;
  u64bit      _numberOfUnique;
  u64bit      _numberOfEntries;
  u64bit      _maximumEntries;

  //  For sorting the mers
  //
  u32bit      _sortedMax;
  u64bit     *_sortedChck;
  u64bit     *_sortedPosn;
};

#endif  //  POSITIONDB_H
