
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

#ifndef INCLUDE_AS_BAT_READ_INFO
#define INCLUDE_AS_BAT_READ_INFO

#include "runtime.H"
#include "ovStore.H"
#include "sqStore.H"



class ReadInfo {
public:
  ReadInfo(const char *seqStorePath, const char *prefix, uint32 minReadLen, uint32 maxReadLen);
  ~ReadInfo();

  //  Needed for writing overlaps in BestOverlapGraph.  Do not use for anything else.
  sqStore *seqStore(void) { return(_seqStore); };

  uint64  memoryUsage(void) {
    return(sizeof(uint64) + sizeof(uint32) + sizeof(uint32) + sizeof(ReadStatus) * (_numReads + 1));
  };

  uint64  numBases(void)     { return(_numBases); };
  uint32  numReads(void)     { return(_numReads); };
  uint32  numLibraries(void) { return(_numLibraries); };

  uint32  readLength(uint32 iid)     { return(_readStatus[iid].readLength); };
  uint32  libraryIID(uint32 iid)     { return(_readStatus[iid].libraryID);  };

  bool    isValid(uint32 iid)        { return((_readStatus[iid].isPresent == 1) &&
                                              (_readStatus[iid].isIgnored == 0)); };


  //  DO NOT USE THIS.  By setting a read length to zero, the read
  //  can be effectively removed from the assembly.  There are probably
  //  Bad Things waiting to happen.  This was used to remove reads
  //  from BestOverlapGraph, to output a 'cleaner' version.
  void    deleteRead(uint32 iid, uint32 len=0) { _readStatus[iid].readLength = len; };



private:
  void    overlapLengths  (uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang, uint32 &aovl, uint32 &bovl);
public:
  uint32  overlapLength   (uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang);

private:
  struct ReadStatus {
    uint64  readLength   : AS_MAX_READLEN_BITS;
    uint64  libraryID    : AS_MAX_LIBRARIES_BITS;
    uint64  isPresent    : 1;
    uint64  isIgnored    : 1;
    uint64  unused       : (64 - AS_MAX_READLEN_BITS - AS_MAX_LIBRARIES_BITS - 2);
  };

#if AS_MAX_READLEN_BITS + AS_MAX_LIBRARIES_BITS + 1 > 64
#error Not enough space in struct ReadStatus.
#endif

  sqStore     *_seqStore;       //  Needed for dumping overlaps.

  uint64       _numBases;       //  The actual number of bases we can assemble.
  uint32       _numReads;       //  The total number of reads in the store.
  uint32       _numLibraries;   //  The total number of libraries in the store.

  ReadStatus  *_readStatus;
};


extern ReadInfo     *RI;


inline
void
ReadInfo::overlapLengths(uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang, uint32 &aovl_, uint32 &bovl_) {
  int32  alen = readLength(a_iid);
  int32  blen = readLength(b_iid);
  int32  aovl = 0;
  int32  bovl = 0;

  aovl_ = 0;
  bovl_ = 0;

  if ((isValid(a_iid) == false) ||   //  If either read is 'deleted', then
      (isValid(b_iid) == false))     //  the overlap is of length zero.
    return;

  //  b_hang < 0      ?     ----------  :     ----
  //                  ?  ----------     :  ----------
  //
  if (a_hang < 0) {
    aovl = (b_hang < 0) ? (alen + b_hang) : (alen);
    bovl = (b_hang < 0) ? (blen + a_hang) : (blen + a_hang - b_hang);
  }

  //  b_hang < 0      ?  ----------              :  ----------
  //                  ?     ----                 :     ----------
  //
  else {
    aovl = (b_hang < 0) ? (alen - a_hang + b_hang) : (alen - a_hang);
    bovl = (b_hang < 0) ? (blen)                   : (blen - b_hang);
  }

  if ((aovl <= 0) || (bovl <= 0) || (aovl > alen) || (bovl > blen)) {
    fprintf(stderr, "WARNING: bogus overlap found for A=" F_U32 " B=" F_U32 "\n", a_iid, b_iid);
    fprintf(stderr, "WARNING:                     A len=" F_S32 " hang=" F_S32 " ovl=" F_S32 "\n", alen, a_hang, aovl);
    fprintf(stderr, "WARNING:                     B len=" F_S32 " hang=" F_S32 " ovl=" F_S32 "\n", blen, b_hang, bovl);
  }

  if (aovl < 0)     aovl = 0;
  if (bovl < 0)     bovl = 0;

  if (aovl > alen)  aovl = alen;
  if (bovl > blen)  bovl = blen;

  assert(aovl > 0);
  assert(bovl > 0);
  assert(aovl <= alen);
  assert(bovl <= blen);

  aovl_ = aovl;
  bovl_ = bovl;
}



inline
uint32
ReadInfo::overlapLength(uint32 a_iid, uint32 b_iid, int32 a_hang, int32 b_hang) {
  uint32  aovl;
  uint32  bovl;

  overlapLengths(a_iid, b_iid, a_hang, b_hang, aovl, bovl);

  return(max(aovl, bovl));
}


#endif  //  INCLUDE_AS_BAT_READ_INFO
