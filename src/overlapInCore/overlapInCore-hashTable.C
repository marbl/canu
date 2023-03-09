
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


//  Estimate how many ignore kmers we'll add.  Most of these should already
//  exist in the table; however, some of them will need to be added.  We'll
//  just allocate space for all of them.  We can't be more precise without
//  actually building the table.
void
hashTable::computeKmerSize(char const *skipFile) {
    uint64  basesInSkip = AS_UTL_sizeOfFile(skipFile);

    uint64  nSkipKmers = basesInSkip / (_kmerSize + 1);
    uint64  nSkipBases = basesInSkip - nSkipKmers;

    _kmersFound += nSkipKmers;
    _basesFound += nSkipBases + nSkipKmers;

    fprintf(stderr, "Estimated %lu skip kmers with %lu bases.\n", nSkipKmers, nSkipBases);
}



//  Count how many reads and bases we'll load.
void
hashTable::computeReadSize(sqStore *seqStore) {
  uint32  nShortReads = 0, nLoadReads = 0;
  uint64  nShortBases = 0, nLoadBases = 0;

  for (uint32 curID=_bgnID; ((_readsFound <  _hashReadsMax) &&
                             (_basesFound <  _hashBasesMax) &&
                             (curID       <= _endID)); curID++) {
    uint32  readLen = seqStore->sqStore_getReadLength(curID);

    if (readLen < _minOverlapLen) {
      nShortReads += 1;
      nShortBases += readLen;
    } else {
      nLoadReads  += 1;
      nLoadBases  += readLen;
    }
  }

  _readsFound += nLoadReads;
  _basesFound += nLoadBases + nLoadReads

  fprintf(stderr, "\n");
  fprintf(stderr, "Loading  %u reads with %lu bases into the hash table.\n",   nLoadReads, nLoadBases);
  fprintf(stderr, "Skipping %u reads with %lu bases shorter than the minimum overlap length %u.\n", nShortReads, nShortBases, _minOverlapLen);
}



//  After computeReadSize() and computeKmerSize() figure out how much data to
//  we need to load, allocate space for it.
void
hashTable::allocateStorage(void) {
  _basesLen   = 0;
  _bases      = new char         [_basesFound];

  _readsLen   = 0;
  _kmersLen   = 0;
  _reads      = new hashReadInfo [_readsFound];

  _nextKmer   = new kmerRef      [_basesFound];
  _kmerList   = new kmerRef      [_basesFound];
}



//  Load read readID from seqStore and save it in our storage.

bool
hashTable::loadReadSequence(sqStore *seqStore, uint32 readID) {
  uint32  idx    = _readsLen;
  uint32  seqLen = seqStore->sqStore_getReadLength(readID);

  assert(idx == readID - _bgnID);

  //  Clear read metadata.
  //
  //  For simplicity later, reads that are too short to be loaded have a
  //  valid bases pointer; it points to a single NUL byte - correctly for a
  //  sequence of length zero - but _hOffset is invalid.

  _reads[idx]._length     = 0;
  _reads[idx]._bases      = _bases + _basesLen;
  _reads[idx]._hOffset    = uint64max;

  _reads[idx]._lScreened  = false;
  _reads[idx]._rScreened  = false;
  _reads[idx]._isKmer     = false;

  //  If the read is too short, stop.  Leave _bases and _hOffset as invalid values.

  if (seqLen < _minOverlapLen) {
    _bases[_basesLen++] = 0;       //  Add the empty read to _bases.
    return(false);
  }

  //  A good read.  Set pointers to where we'll put the bases, and to
  //  where we'll save kmer references.

  _reads[idx]._length     = seqLen;
  _reads[idx]._bases      = _bases + _basesLen;
  _reads[idx]._hOffset    =          _basesLen;

  //  Grab the read from seqStore.  Copy the bases to storage.

  seqStore->sqStore_getRead(readID, &_read);

  char   *seq = _read.sqRead_sequence();

  for (uint32 ii=0; ii<seqLen; ii++)
    _bases[_basesLen++] = seq[ii];
  _bases[_basesLen++] = 0;

  //  Blow up if we've exceeded our allocation.

  if (_basesLen > _basesFound)
    fprintf(stderr, "ERROR: loaded too much data.  basesLen=%lu basesFound=%lu\n", _basesLen, _basesFound);
  assert(_basesLen <= _basesFound);

  return(true);
}



  //  Load sequence from disk and insert into the table.
  //
  //  Once loaded, add kmers in the read to the hash table.  keyIterator is
  //  valid for a sequence of length zero - the first next() returns false -
  //  and we can let loadReadSequence filter reads we don't want to load by
  //  just adding an empty sequence.
  //
  //  Note that curID is the ID of the read in seqStore, while _readsLen is
  //  (abused to be) the ID of the read in our metadata.
  //
void
hashTable::loadReads(sqStore *seqStore) {

  for (uint32 curID=_bgnID; ((_readsFound <  _hashReadsMax) &&
                             (_basesFound <  _hashBasesMax) &&
                             //  LOAD FACTOR
                             (curID       <= _endID)); curID++, _readsLen++) {
    loadReadSequence(seqStore, curID);

    keyIterator   it(_reads[_readsLen]._bases, _kmerSize);

    while (it.next() == true)
      if (it.isValid() == true)
        insertKmer(it._key, _readsLen, it._pos, false);

    if ((_readsLen % 100000) == 0)
      fprintf (stderr, "reads:%12" F_U32P "/%12" F_U32P "  bases:%12" F_U64P "/%12" F_U64P "  kmers:%12" F_U64P "/%12" F_U64P "  Load: %.2f%%\n",
               _readsLen,    _hashReadsMax,
               _basesLen,    _hashBasesMax,
               _hashEntries, _hashEntriesMax,  100.0 * _hashEntries / _hashEntriesMax);
  }

  fprintf(stderr, "HASH LOADING STOPPED: read id  %12" F_U32P " out of %12" F_U32P "\n", _readsLen,    _hashReadsMax);
  fprintf(stderr, "HASH LOADING STOPPED: bases    %12" F_U64P " out of %12" F_U64P "\n", _basesLen,    _hashBasesMax);
  fprintf(stderr, "HASH LOADING STOPPED: entries  %12" F_U64P " out of %12" F_U64P "\n", _hashEntries, _hashEntriesMax);
}


//  Load a chunk of sequence to sequence storage and 
void
hashTable::loadKmerSequence(char *seq, uint32 seqLen) {
  uint32  idx = _readsLen + _kmersLen;

  //  Initialize read metadata.

  _reads[idx]._length     = seqLen;
  _reads[idx]._bases      = _bases + _basesLen;
  _reads[idx]._hOffset    =          _basesLen;

  _reads[idx]._lScreened  = false;
  _reads[idx]._rScreened  = false;
  _reads[idx]._isKmer     = true;

  //  Copy the sequence to storage.

  for (uint32 ii=0; ii<seqLen; ii++)
    _bases[_basesLen++] = tolower(seq[ii]);
  _bases[_basesLen++] = 0;

  //  Blow up if we've exceeded our allocation.

  if (_basesLen > _basesFound)
    fprintf(stderr, "ERROR: loaded too much data.  basesLen=%lu basesFound=%lu\n", _basesLen, _basesFound);
  assert(_basesLen <= _basesFound);
}


  //
  //  Load the over occurring kmers from disk, and flag them in the table.
  //

bool
isACGT(char c) {
  return((c == 'a') || (c == 'A') ||
         (c == 'c') || (c == 'C') ||
         (c == 'g') || (c == 'G') ||
         (c == 't') || (c == 'T'));
}

void
hashTable::loadKmers(char const *skipFile) {
  compressedFileReader  *in = new compressedFileReader(frequentMersPath);

  while (in->readLine() == true) {
    char   *line    = in->line();
    uint32  lineLen = 0;

    //  Convert to uppercase, truncating the line at the first non-acgt.
    //  This lets us trim off any garbage at the end of the line, like the
    //  count.  And ignores lines that aren't startig with ACGT.

    for (lineLen=0; lineLen < in->lineLen(); lineLen++) {
      line[lineLen] = toupper(line[lineLen]);

      if (isACGT(line[lineLen]) == false) {
        line[lineLen]   = 0;   //  Terminate the line,
        break;                 //  terminate the loop.
      }
    }

    //  If smaller than a kmer size, ignore it.

    if (lineLen < _kmerSize)
      continue;

    //  Add the line to storage, then add the kmers to the table and flag as
    //  over occurring.  Reverse the sequence and do it all again.

    loadKmerSequence(line, lineLen);

    for (keyIterator  it(line, _kmerSize); it.next(); )
      insertKmer(it._key, _readsLen + _kmersLen, it._pos, true);

    _kmersLen++;

    reverseComplementSequence(line, lineLen);
    loadKmerSequence(line, lineLen);

    for (keyIterator  it(line, _kmerSize); it.next(); )
      insertKmer(it._key, _readsLen + _kmersLen, it._pos, true);

    _kmersLen++;
  }

  delete    in;

  fprintf(stderr, "\n");
  fprintf(stderr, "Marked %u over-occurring kmers.\n", _kmersLen);
  fprintf(stderr, "\n");
}
