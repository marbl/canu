
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

#ifndef OVERLAPINCORE_HASHREADINFO_H
#define OVERLAPINCORE_HASHREADINFO_H

#include "runtime.H"
#include "types.H"

//
//  Data for each sequence - either a read or a frequent kmer - in the hash
//  table.
//
//    _length    - the length of the read in bases
//    _bases     - pointer to the sequence stored in hashTable::_bases
//    _hOffset   - position of the sequence in hashTable::_bases
//
//    _lScreened - true if a frequent kmer occurs near the end
//    _rScreened   of the read
//
//    _isKmer    - true if the sequence is from a frequent kmer string
//

class hashReadInfo {
public:
  uint32      _length       = 0;
  char       *_bases        = nullptr;
  uint64      _hOffset      = uint64max;

  bool        _lScreened    = false;
  bool        _rScreened    = false;
  bool        _isKmer       = false;
};

#endif
