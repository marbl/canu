
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

#include "overlapInCore.H"

#include "sequence.H"
#include "strings.H"


//  Insert a kmer 'key' at position 'keyid,keypos' into the hash table.

void
hashTable::insertKmer(uint64    key,
                      uint32    keyid,      //  ID of sequence in _reads
                      uint32    keypos,     //  position of the key in that sequence
                      bool      isFreq) {
  uint64  bucket    = hashFunctions::hashFunction(key);
  uint64  probe     = hashFunctions::hashProbe(key);
  uint8   keycheck  = hashFunctions::keyCheck(key);
  uint64  keyoffset = _reads[keyid]._hOffset + keypos;

  _buckets[bucket].setIsPresent(key);  //  'key' is present in this bucket, or a later probe

  //  Iterate over buckets until we find the kmer in a bucket, or find space
  //  to add a new kmer in a bucket, or we've visited every bucket in the
  //  table and decide the table is full.

  for (int64 Ct=0; Ct < _bucketsLen; Ct++) {
    for (int32 ee=0; ee < _buckets[bucket]._entryLen; ee++) {
      kmerRef   hashref = _buckets[bucket]._entry[ee];
      char     *s       = _reads[keyid            ]._bases + keypos;
      char     *t       = _reads[hashref._stringID]._bases;

      if ((_buckets[bucket]._echeck[ee] != keycheck) ||   //  If check fails, kmer can't be here,
          (strncmp(s, t, _kmerSize) != 0))                //  otherwise verify it's the same kmer,
        continue;                                         //  and if not, keep looking.

      _nextKmer[keyoffset] = hashref;                     //  Move existing ref to storage,

      _buckets[bucket]._entry[ee]._stringID  = keyid;    //  And set the hash table bucket
      _buckets[bucket]._entry[ee]._stringPos = keypos;   //  pointing to that reference.
      _buckets[bucket]._entry[ee]._isPointer = false;    //  (keyoffset is computed using
      _buckets[bucket]._entry[ee]._ignore    = isFreq;   //  keyid and keypos).
      _buckets[bucket]._entry[ee]._isLast    = false;

      return;
    }

    //  Didn't find the kmer in the table.  If no space in the bucket for a new
    //  entry move to the next probe and try again.

    if (_buckets[bucket]._entryLen == ENTRIES_PER_BUCKET) {
      bucket = (bucket + probe) % _bucketsLen;
      continue;
    }

    //  Didn't find the kmer in the table, but there is space to add a new
    //  entry into this bucket!

    int32 ee = _buckets[bucket]._entryLen;

    _buckets[bucket]._entry[ee]._stringID  = keyid;
    _buckets[bucket]._entry[ee]._stringPos = keypos;
    _buckets[bucket]._entry[ee]._isPointer = false;
    _buckets[bucket]._entry[ee]._ignore    = isFreq;
    _buckets[bucket]._entry[ee]._isLast    = true;

    _buckets[bucket]._echeck[ee] = keycheck;

    _buckets[bucket]._entryLen += 1;

    _hashEntries++;

    return;
  }

  //  If we get here, we searched every bucket for the kmer or for space to
  //  add a new one, and failed.  All we can do is blow up.

  fprintf(stderr, "ERROR:  Hash table full\n");
  exit(1);
}




//  Loads reads bgnID <= x < endID into the hash table.
//
//  Returns the id of the first read that isn't loaded into the table - that
//  is, returns endID if all reads are loaded.
//
uint32
hashTable::loadTable(sqStore *seqStore,
                     uint32 bgnID, uint32 endID,
                     char *frequentMersPath) {

  fprintf(stderr, "Build_Hash_Index from " F_U32 " to " F_U32 "\n", bgnID, endID);

  _bgnID = bgnID;
  _endID = bgnID;

  //
  //  Compute an upper limit on how many reads and how much sequence we will
  //  load into the hash table, then allocate space for all this good stuff.
  //
  //  The loading can still stop early if the hash table load factor gets too
  //  high.
  //

  computeReadSize(seqStore);
  computeKmerSize(frequentMersPath);

  allocateStorage();

  loadReads(seqStore);
  loadKmers(frequentMersPath);




  //
  //  Coalesce all the kmerRefs into a single array.  While doing so, set the
  //  _lScreened and _rScreened flags on reads as appropriate.
  //

  //  Count the amount of space we need to copy all the references out of the hash table.

  uint64  kmerListMax = 0;
  uint64  kmerListPos = 0;

  for (uint64 hh=0; hh < _bucketsLen; hh++) { 
    for (int32 ee=0; ee < _buckets[hh]._entryLen; ee++) {
      kmerRef   ref = _buckets[hh]._entry[ee];
      uint32    num = 1;

      while (ref._isLast == false) {    //  Count the length of the chain.
        ref = _nextKmer[ kmerPosition(ref) ];
        num++;
      }

      if (num > 1)                      //  If more than one, we'll copy it
        kmerListMax += num;             //  to the linear list.
    }
  }

  //  Allocate space for the copies.

  _kmerList = new kmerRef [kmerListMax];

  //  Traverse the hash table again, this time copying references to the linear list.

  for (uint64 hh=0; hh < _bucketsLen; hh++) {
    for (int32 ee=0; ee < _buckets[hh]._entryLen; ee++) {
      kmerRef   ref = _buckets[hh]._entry[ee];

      //  Leave as is if this is the only reference in the chain.

      if (ref._isLast == true)
        continue;

      //  Remember if this is a frequent kmer.

      bool isFreq = ref._ignore;

//  Mark  left/right_end_screened in global  String_Info  for
//   ref  and everything in its list, if they occur near
//  enough to the end of the string.

#if 0
static
void
Mark_Screened_Ends_Single(String_Ref_t ref) {
  uint32 sID     = ref._stringID;

  int32  bgnSize =                           ref._stringPos;
  int32  endSize = String_Info[sID].length - ref._stringPos - _kmerSize + 1;

  if (bgnSize < HOPELESS_MATCH)   String_Info[sID].lfrag_end_screened = true;
  if (endSize < HOPELESS_MATCH)   String_Info[sID].rfrag_end_screened = true;
}
#endif



      //  Otherwise:
      //   - convert the ref in the hash table itself to a pointer to the linear list
      //   - copy the ref to the linear list
      //   - copy the rest of the reference chain to the linear list
      //   - set lScreened and rScreened

      _buckets[hh]._entry[ee].setLinearPointer(kmerListPos);


      {
        uint32 sID     =                       ref._stringID;
        int32  bgnSize =                       ref._stringPos;
        int32  endSize = _reads[sID]._length - ref._stringPos - _kmerSize + 1;

        if (bgnSize < HOPELESS_MATCH)   _reads[sID]._lScreened = true;
        if (endSize < HOPELESS_MATCH)   _reads[sID]._rScreened = true;
      }

      _kmerList[kmerListPos++] = ref;

      while (ref._isLast == false) {
        ref = _nextKmer[ kmerPosition(ref) ];
        _kmerList[kmerListPos++] = ref;
      }
    }
  }

  //  Remove the nextKmer chain.

  delete [] _nextKmer;
  _nextKmer = nullptr;

  //  Finished.  Return the readID of the last read loaded.

  return(_bgnID + _readsLen);
}
