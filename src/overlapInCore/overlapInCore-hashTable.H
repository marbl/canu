
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

#ifndef OVERLAPINCORE_HASHTABLE_H
#define OVERLAPINCORE_HASHTABLE_H

#include "runtime.H"
#include "types.H"

#include "overlapInCore-hashReadInfo.H"


//
//  Each hashBucket in the hash table has space for ENTRIES_PER_BUCKET kmers.
//  Each entry has an associated 'check' value, which allows us to skip an
//  expensive string comparison to decide if the hash hit is true of false.
//


//  Recommended values are 21, 31 or 42, depending on cache line size.
#define ENTRIES_PER_BUCKET 21



class hashBucket {
public:
  hashBucket() {
    for (uint32 ii=0; ii<ENTRIES_PER_BUCKET; ii++) {
      _entry [ii] = kmerRef();
      _echeck[ii] = 0;
    }

    _entryLen = 0;
    _present  = 0;
  }

public:
  void   setIsPresent(uint64 k)   {         _present |= hashFunctions::isPresentBit(k);        };
  bool   isNotPresent(uint64 k)   { return((_present  & hashFunctions::isPresentBit(k)) == 0); };

  bool   isFull(void) {
    return(_entryLen == ENTRIES_PER_BUCKET);
  };
  
public:
  kmerRef        _entry [ENTRIES_PER_BUCKET];
  uint8          _echeck[ENTRIES_PER_BUCKET];
  uint32         _entryLen;
  uint32         _present;
};



class hashTable {
public:
  hashTable(uint32 kmerSize,
            uint32 hashBits,
            double maxLoad,
            uint32 minOverlapLen) {
    _kmerSize       = kmerSize;
    _hashBits       = hashBits;
    _maxLoad        = maxLoad;
    _minOverlapLen  = minOverlapLen;

    _bucketsLen     = uint64one << _hashBits;
    _buckets        = new hashBucket [_bucketsLen];

    _hashReadsMax   = uint32max;
    _hashBasesMax   = uint64max;
    _hashEntriesMax = _maxLoad * _bucketsLen * ENTRIES_PER_BUCKET;
  };

  ~hashTable() {
    delete [] _bases;
    delete [] _reads;
    delete [] _buckets;
    delete [] _nextKmer;
    delete [] _kmerList;
  };

public:
  uint32   kmerSize(void) const     { return(_kmerSize); };

private:
  void     insertKmer(uint64 key, uint32 keyid, uint32 keypos, bool isFreq);

private:
  void     computeKmerSize(char const *skipFile);
  void     computeReadSize(sqStore *seqStore);

  void     allocateStorage(void);

  //  Save a chunk of sequence 
private:
  void     loadKmerSequence(char *seq, uint32 seqLen);
  bool     loadReadSequence(sqStore *seqStore, uint32 id);

  void     loadReads(sqStore *seqStore);
  void      loadKmers(char const *skipFile);

public:
  uint32   loadTable(sqStore *seqStore,
                     uint32 bgnID, uint32 endID,
                     char *frequentMersPath);

  //  Returns the position of the kmerRef 'kmer' in the _bases, _nextKmer and
  //  _kmerList arrays.
private:
  uint64   kmerPosition(kmerRef kmer) {
    return(_reads[kmer._stringID]._hOffset + kmer._stringPos);
  };

  //  The range of reads we have loaded into the table.  Usual C semantics.
private:
  uint32           _bgnID          = 0;
  uint32           _endID          = 0;

  //  Parameters and limits on how much data we can load into the table.
private:
  uint32           _kmerSize       = 0;
  uint32           _hashBits       = 0;
  double           _maxLoad        = 0.0;
  uint32           _minOverlapLen  = 0;

  uint32           _hashReadsMax   = 0;
  uint64           _hashBasesMax   = 0;
  uint64           _hashEntriesMax = 0;

  //  Initial count of how much stuff we will load into the table.  These are
  //  used just to size the various arrays; the table can have less data
  //  loaded if the hash table load factor is met early.
  //
  //  _read is a temporary used for loading reads from sqStore.  It's vastly
  //  easier to have it here than to pass it around in the function.
private:
  uint64           _basesFound     = 0;   //  Amount sequence in reads and skip kmers we need to save.
  uint64           _readsFound     = 0;   //  Number of reads;          metadata 0          .. readsFound
  uint64           _kmersFound     = 0;   //  Number of kmer sequences; metadata readsFound .. readsFound+kmersFound

  sqRead           _read;   //  A temporary used when loading sequence into the hash table.

  //  Data!  The sequence data, the read meta data and the hash table itself.
private:
  uint64           _basesLen       = 0;
  char            *_bases          = nullptr;

  uint32           _readsLen       = 0;
  uint32           _kmersLen       = 0;
  hashReadInfo    *_reads          = nullptr;   //  String_Info

  uint64           _bucketsLen     = 0;
  hashBucket      *_buckets        = nullptr;   //  Hash_Table

  uint64           _hashEntries    = 0;

  kmerRef         *_nextKmer       = nullptr;   //  nextRef
  kmerRef         *_kmerList       = nullptr;   //  Extra_Data_Len
};


#endif
