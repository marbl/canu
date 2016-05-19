
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
 *  This file is derived from:
 *
 *    kmer/libmeryl/libmeryl.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-SEP-08 to 2004-APR-08
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAY-23 to 2014-APR-11
 *      are Copyright 2005,2007-2008,2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-05 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef LIBMERYL_H
#define LIBMERYL_H

#include "kMer.H"

//  A merStream reader/writer for meryl mercount data.
//
//  merSize is used to check that the meryl file is the correct size.
//  If it isn't the code fails.
//
//  The reader returns mers in lexicographic order.  No random access.
//  The writer assumes that mers come in sorted increasingly.
//
//  numUnique    the total number of mers with count of one
//  numDistinct  the total number of distinct mers in this file
//  numTotal     the total number of mers in this file


class merylStreamReader {
public:
  merylStreamReader(const char *fn, uint32 ms=0);
  ~merylStreamReader();

  kMer           &theFMer(void)      { return(_thisMer);          };
  uint64          theCount(void)     { return(_thisMerCount);     };

  bool            hasPositions(void)    { return(_POS != 0L); };
  uint32         *thePositions(void)    { return(_thisMerPositions); };
  uint32          getPosition(uint32 i) { return(((_POS) && (i < _thisMerCount)) ? _thisMerPositions[i] : ~uint32ZERO); };

  uint32          merSize(void)         { return(_merSizeInBits >> 1); };
  uint32          merCompression(void)  { return(_merCompression); };

  uint32          prefixSize(void) { return(_prefixSize); };

  uint64          numberOfUniqueMers(void)   { return(_numUnique); };
  uint64          numberOfDistinctMers(void) { return(_numDistinct); };
  uint64          numberOfTotalMers(void)    { return(_numTotal); };

  uint64          histogram(uint32 i)         { return((i < _histogramLen) ? _histogram[i] : ~uint64ZERO); };
  uint64          histogramLength(void)       { return(_histogramLen); };
  uint64          histogramMaximumCount(void) { return(_histogramMaxValue); };

  bool            nextMer(void);
  bool            validMer(void) { return(_validMer); };
private:
  char                   _filename[FILENAME_MAX];

  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;
  bitPackedFile         *_POS;

  uint64                  getIDXnumber(void) {
    uint64 n = 1;

    if (_idxIsPacked)
      n = _IDX->getNumber();
    else
      n = _IDX->getBits(32);

    return(n);
  };
  uint64                  getDATnumber(void) {
    uint64 n = 1;

    if (_datIsPacked) {
      if (_DAT->getBits(1))
        n = _DAT->getNumber() + 2;
    } else {
      n = _DAT->getBits(32);
    }

    return(n);
  };

  //  Why not bool?  Seems like the bitPackedFile is incompatible
  //  with bools.
  uint32                 _idxIsPacked;
  uint32                 _datIsPacked;
  uint32                 _posIsPacked;

  uint32                 _merSizeInBits;
  uint32                 _merCompression;
  uint32                 _prefixSize;
  uint32                 _merDataSize;
  uint64                 _thisBucket;
  uint64                 _thisBucketSize;
  uint64                 _numBuckets;

  kMer                   _thisMer;
  uint64                 _thisMerCount;

  uint32                 _thisMerPositionsMax;
  uint32                *_thisMerPositions;

  uint64                 _numUnique;
  uint64                 _numDistinct;
  uint64                 _numTotal;

  uint64                 _histogramPos;        // position of the histogram data in IDX
  uint64                 _histogramLen;        // number of entries in the histo
  uint64                 _histogramMaxValue;   // highest count ever seen
  uint64                *_histogram;

  bool                   _validMer;
};


class merylStreamWriter {
public:
  merylStreamWriter(const char *filePrefix,
                    uint32 merSize,          //  In bases
                    uint32 merComp,          //  A length, bases
                    uint32 prefixSize,       //  In bits
                    bool   positionsEnabled);
  ~merylStreamWriter();

  void                    addMer(kMer &mer, uint32 count=1, uint32 *positions=0L);
  void                    addMer(uint64 prefix, uint32 prefixBits,
                                 uint64 mer,    uint32 merBits,
                                 uint32 count=1,
                                 uint32 *positions=0L);

private:
  void                    writeMer(void);

  void                    setIDXnumber(uint64 n) {
    if (_idxIsPacked)
      _IDX->putNumber(n);
    else
      _IDX->putBits(n, 32);
  };
  void                    setDATnumber(uint64 n) {
    if (_datIsPacked) {
      if (n == 1) {
        _DAT->putBits(uint64ZERO, 1);
      } else {
        _DAT->putBits(uint64ONE, 1);
        _DAT->putNumber(n-2);
      }
    } else {
      _DAT->putBits(n, 32);
    }
  };

  char                   _filename[FILENAME_MAX];

  bitPackedFile         *_IDX;
  bitPackedFile         *_DAT;
  bitPackedFile         *_POS;

  uint32                 _idxIsPacked;
  uint32                 _datIsPacked;
  uint32                 _posIsPacked;

  uint32                 _merSizeInBits;
  uint32                 _merCompression;
  uint32                 _prefixSize;
  uint32                 _merDataSize;
  uint64                 _thisBucket;
  uint64                 _thisBucketSize;
  uint64                 _numBuckets;

  uint64                 _numUnique;
  uint64                 _numDistinct;
  uint64                 _numTotal;

  uint64                 _histogramPos;        // position of the histogram data in IDX
  uint64                 _histogramLen;        // number of entries in the histogram
  uint64                 _histogramMaxValue;   // highest count ever seen
  uint64                *_histogram;

  bool                   _thisMerIsBits;
  bool                   _thisMerIskMer;

  kMer                   _thisMer;

  uint64                 _thisMerPre;
  uint64                 _thisMerMer;

  uint32                 _thisMerPreSize;
  uint32                 _thisMerMerSize;

  uint64                 _thisMerCount;
};

#endif  //  LIBMERYL_H
