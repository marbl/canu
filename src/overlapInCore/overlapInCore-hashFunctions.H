
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

#ifndef OVERLAPINCORE_HASHFUNCTIONS_H
#define OVERLAPINCORE_HASHFUNCTIONS_H

#include "runtime.H"
#include "types.H"

class hashFunctions {
public:
  hashFunctions() {
  };

  void   initialize(uint32 hashBits, uint32 kmerLen) {
    MASK  = buildLowBitMask<uint64>(hashBits);

    HSF1  =     kmerLen - hashBits / 2;
    HSF2  = 2 * kmerLen - hashBits;

    SV1   = HSF1 + 2;
    SV2   = (HSF1 + HSF2) / 2;   //  == 3/2 kmerLen - 3/4 hashBits, not that it helps to understand.
    SV3   = HSF2 - 2;
  };

  static uint64 hashFunction(uint64 k)   { return(             (((k) ^ (k >> HSF1) ^ (k >> HSF2)) & MASK) ); }
  static uint32 hashProbe(uint64 k)      { return(uint32one  | (((k) ^ (k >>  SV2) ^ (k >>  SV3)) & 0x3e) ); }
  static uint32 keyCheck(uint64 k)       { return(             (((k) ^ (k >>  SV1) ^ (k >>  SV3)) & 0xff) ); }
  static uint32 isPresentBit(uint64 k)   { return(uint32one << (((k) ^ (k >>  SV1) ^ (k >>  SV2)) & 0x1f) ); }

private:
  static uint64  MASK;   //  These are instiantiated in overlapInCore.C
  static uint32  HSF1;
  static uint32  HSF2;
  static uint32  SV1;
  static uint32  SV2;
  static uint32  SV3;
};

#endif
