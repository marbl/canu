
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef OVERLAPINCORE_KMERREF_H
#define OVERLAPINCORE_KMERREF_H

#include "runtime.H"
#include "types.H"

//
//  The object we store in the hash table, a fancy pointer to a position
//  in an input string (stored in _reads (hashReadInfo) and _bases (char).
//

//
//  Initialize:
//   - _stringPos is set to the equivalent of -1 so the keyIterator can
//     track position properly.
//   - _ignore (previously "Empty") is true if this ref is to a frequent
//     kmer and we should ignore it.  It must be initialized to false.
//

class kmerRef {
public:
  kmerRef() {
    _stringID  = buildLowBitMask<uint64>(64 - 2 - AS_MAX_READLEN_BITS);
    _stringPos = buildLowBitMask<uint64>(         AS_MAX_READLEN_BITS);
    _isPointer = 0x0llu;
    _ignore    = 0x0llu;
    _isLast    = 0x1llu;
  };

  kmerRef(uint32 sID) {
    _stringID  = sID;
    _stringPos = buildLowBitMask<uint64>(         AS_MAX_READLEN_BITS);
    _isPointer = 0x0llu;
    _ignore    = 0x0llu;
    _isLast    = 0x1llu;
  };

  void      setLinearPointer(uint64 pointer) {
    _stringID  = (pointer >> AS_MAX_READLEN_BITS) & buildLowBitMask<uint64>(64 - 3 - AS_MAX_READLEN_BITS);
    _stringPos = (pointer                       ) & buildLowBitMask<uint64>(         AS_MAX_READLEN_BITS);
    _isPointer = 1;
  };
  uint64    getLinearPointer(void) {
    return(_stringID << AS_MAX_READLEN_BITS | _stringPos);
  };

  uint64   _stringID  : 64 - 3 - AS_MAX_READLEN_BITS;
  uint64   _stringPos :          AS_MAX_READLEN_BITS;
  uint64   _isPointer : 1;
  uint64   _ignore    : 1;
  uint64   _isLast    : 1;
};


#endif
