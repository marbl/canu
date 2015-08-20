
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
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-20 to 2014-APR-11
 *      are Copyright 2005,2007,2012-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Liliana Florea on 2010-DEC-06
 *      are Copyright 2010 Liliana Florea, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef EXISTDB_H
#define EXISTDB_H

//  Used by wgs-assembler, to determine if a rather serious bug was patched.
#define EXISTDB_H_VERSION 1960

#include "bio++.H"

//  Takes as input a list of mers (in a file) and builds a searchable
//  structure listing those mers.  Duplicate mers are not removed and
//  will be stored multiple times.
//
//  Using a compressed hash is allowed, but somewhat useless -- it is
//  really slow and doesn't save that much.
//
//  If existDBcanonical is requested, this will store only the
//  canonical mer.  It is up to the client to be sure that is
//  appropriate!  See positionDB.H for more.

//#define STATS

typedef uint32 existDBflags;
const existDBflags  existDBnoFlags         = 0x0000;
const existDBflags  existDBcompressHash    = 0x0001;
const existDBflags  existDBcompressBuckets = 0x0002;
const existDBflags  existDBcompressCounts  = 0x0004;
const existDBflags  existDBcanonical       = 0x0008;
const existDBflags  existDBforward         = 0x0010;
const existDBflags  existDBcounts          = 0x0020;

class existDB {
public:

  //  Read state from an existDB file
  existDB(char const  *filename,
          bool         loadData=true);

  //  Load mers from an existing existDB file, a fastafile, or a meryl database
  existDB(char const    *filename,
          uint32         merSize,
          existDBflags   flags,
          uint32         lo,
          uint32         hi);

  //  Load mers from a character string
  existDB(char const    *sequence,
          uint32         merSize,
          existDBflags   flags);

  ~existDB();

  void        saveState(char const *filename);

  void        printState(FILE *stream);

  bool        isForward(void)    { return(_isForward);   };
  bool        isCanonical(void)  { return(_isCanonical); };

  bool        exists(uint64 mer);
  uint64      count(uint64 mer);

private:
  bool        loadState(char const *filename, bool beNoisy=false, bool loadData=true);
  bool        createFromFastA(char const  *filename,
                              uint32       merSize,
                              uint32       flags);
  bool        createFromMeryl(char const  *filename,
                              uint32       merSize,
                              uint32       lo,
                              uint32       hi,
                              uint32       flags);
  bool        createFromSequence(char const  *sequence,
                                 uint32       merSize,
                                 uint32       flags);

  uint64       HASH(uint64 k) {
    return(((k >> _shift1) ^ (k >> _shift2) ^ k) & _mask1);
  };

  uint64       CHECK(uint64 k) {
    return(k & _mask2);
  };

  void         insertMer(uint64 hsh, uint64 chk, uint64 cnt, uint64 *countingTable) {

    //  If the mer is already here, just update the count.  This only
    //  works if not _compressedBucket, and only makes sense for loading from
    //  fasta or sequence.

    if ((_compressedBucket == false) &&
        (_searchForDupe)) {
      uint64 st = _hashTable[hsh];
      uint64 ed = countingTable[hsh];

      for (; st<ed; st++) {
        if (_buckets[st] == chk) {
          if (_counts)
            _counts[st] += cnt;
          return;
        }
      }
    }

    if (_compressedBucket)
      setDecodedValue(_buckets, countingTable[hsh] * _chkWidth, _chkWidth, chk);
    else
      _buckets[countingTable[hsh]] = chk;

    if (_counts) {
      if (_compressedCounts) {
        setDecodedValue(_counts, countingTable[hsh] * _cntWidth, _cntWidth, cnt);
      } else {
        _counts[countingTable[hsh]] = cnt;
      }
    }

    countingTable[hsh]++;
  };

  bool        _compressedHash;
  bool        _compressedBucket;
  bool        _compressedCounts;
  bool        _isForward;
  bool        _isCanonical;

  bool        _searchForDupe;

  uint32      _merSizeInBases;

  uint32      _shift1;
  uint32      _shift2;
  uint64      _mask1;
  uint64      _mask2;

  uint32      _hshWidth;  //  Only for the compressed hash
  uint32      _chkWidth;  //  Only for the compressed bucket
  uint32      _cntWidth;  //  Only for the compressed counts

  uint64      _hashTableWords;
  uint64      _bucketsWords;
  uint64      _countsWords;

  uint64     *_hashTable;
  uint64     *_buckets;
  uint64     *_counts;

  void clear(void) {
  };
};

#endif  //  EXISTDB_H
