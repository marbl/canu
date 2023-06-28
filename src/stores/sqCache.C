
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

#include "sqCache.H"
#include "sequence.H"

#include <set>
#include <vector>
#include <algorithm>


void
sqCache::loadMetadata(void) {

  _readsMax = _seqStore->sqStore_lastReadID() + 1;
  _readsLen = _seqStore->sqStore_lastReadID() + 1;
  _reads    = new sqCacheEntry [_readsMax];

  for (uint32 id=0; id < _readsMax; id++) {
    _reads[id]._nLen           = 0;
    _reads[id]._sLen           = 0;
    _reads[id]._bgn            = 0;
    _reads[id]._end            = 0;
    _reads[id]._nData          = nullptr;
    _reads[id]._sData          = nullptr;
  }

  for (uint32 id=1; id < _readsLen; id++) {
    if (_seqStore->sqStore_isIgnoredRead(id, _which) == true)
      continue;

    //  Set the length of the data we're going to store.

    if (_which & sqRead_raw)
      _reads[id]._sLen = _seqStore->sqStore_getReadLength(id, sqRead_raw);

    if (_which & sqRead_corrected)
      _reads[id]._sLen = _seqStore->sqStore_getReadLength(id, sqRead_corrected);

    //  Set the portion of the read we should be returning to the user.

    _reads[id]._bgn = (_trimmed == true) ? _seqStore->sqStore_getClearBgn(id, _which) : 0;
    _reads[id]._end = (_trimmed == true) ? _seqStore->sqStore_getClearEnd(id, _which) : _seqStore->sqStore_getReadLength(id, _which);
  }
}


sqCache::sqCache(sqStore       *seqStore,
                 sqRead_which   which) {

  _seqStore        =  seqStore;

  _which           = which;
  _compressed      = ((_which & sqRead_compressed) == sqRead_unset) ? false : true;
  _trimmed         = ((_which & sqRead_trimmed)    == sqRead_unset) ? false : true;

  if (_seqStore)
    loadMetadata();
}



sqCache::~sqCache() {
  if (_data == nullptr)                        //  If no data blocks, we have allocated
    for (uint32 ii=0; ii < _readsLen; ii++)    //  data for each read.
      delete [] _reads[ii]._sData;

  delete [] _reads;                            //  Delete read metadata.

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)   //  Delete any data blocks
    delete [] _dataBlocks[ii];

  delete [] _dataBlocks;                       //  And pointers to data blocks.
}



void
sqCache::loadRead(uint32 id) {

  if ((_reads[id]._sData != nullptr) ||        //  If already loaded, we're done.
      (_reads[id]._sLen == 0))                 //  If the read doesn't exist, we're done.
    return;

  //  Load the encoded blob, without decoding it.

  _read.sqRead_fetchBlob(_seqStore->sqStore_getReadBuffer(id));

  //  Find the encoded read data.  This mirrors sqRead_loadFromBuffer.

  uint32   blobPos  = 0;
  uint8   *rptr     = nullptr;
  uint8   *cptr     = nullptr;

  while (blobPos < _read._blobLen) {
    char   *cName =  (char *)  (_read._blob + blobPos + 0);
    uint32  cLen  = *(uint32 *)(_read._blob + blobPos + 4);

    if (((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
        ((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
        ((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')))
      rptr = _read._blob + blobPos;

    if (((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')) ||
        ((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')) ||
        ((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
      cptr = _read._blob + blobPos;

    blobPos += 8 + cLen;
  }

  //  Save either the raw or corrected sequence.

  uint8   *bptr = (_which & sqRead_raw) ? rptr : cptr;
  uint32   blen = *(uint32 *)(bptr + 4) + 8;

  //  If we have a gigantic storage space for read data, use that, otherwise,
  //  allocate space for this data.

  if (_data == nullptr) {
    _reads[id]._sData = new uint8 [blen];
  }

  else {
    if (_dataLen + blen > _dataMax)
      allocateNewBlock();

    _reads[id]._sData = _data + _dataLen;
  }

  //  Copy the data and release the blob.

  memcpy(_reads[id]._sData, bptr, blen);

  //  Update the pointer to the next free chunk of storage.

  if (_data != nullptr) {
    _dataLen += blen;

    assert(_dataLen <= _dataMax);
  }
}



void
sqCache::loadRead(dnaSeq &seq) {

  //  This avoids a crash in merylutil::encode2bitSequence() when it tries
  //  to check that the sequence is nul-terminated; accessing seq[len-1] is
  //  invalid.  But Canu can't (Oct 2022) update merylutil without huge
  //  changes to use it's new namespace support.
  //
  if (seq.length() == 0)
    return;

  //  Make space for this read.
  //
  //  We need space for both metadata (in _reads) and read sequence (in
  //  _data, etc).  The former we grow by 128k entries whenever needed; the
  //  latter is allocated initially, then whenever the read data has a chance
  //  of not fitting in the current chunk.

  increaseArray(_reads, _readsLen+1, _readsMax, 131072);

  if (_dataBlocksLen == 0)
    allocateNewBlock();

  if (_dataLen + 4 + 4 + seq.length() + 4 > _dataMax)
    allocateNewBlock();

  //  Clear metadata for this read.

  uint32  id = ++_readsLen;
  uint8  *dd = _data + _dataLen + 8;

  _reads[id]._sLen  = seq.length();
  _reads[id]._bgn   = 0;
  _reads[id]._end   = seq.length();
  _reads[id]._sData = _data + _dataLen;

  //  Save the name to id mapping.

  std::string readname = seq.ident();

  _nameToID[ readname ] = id;
  _IDtoName.push_back( _nameToID.find(readname) );

  //  Encode the data as 2-bit, 3-bit, or plain bases, whatever works first.
  //  Note the '+8' is to leave space at the start for the AIFF tag and
  //  length; see sqReadData.C and/or sqReadDataWriter.C for details,

  uint8  tag[4] = { '2', 'S', 'Q', 'R' };

  uint32 el2    =                            encode2bitSequence(dd, seq.bases(), seq.length());
  uint32 el3    = (el2 == 0)               ? encode3bitSequence(dd, seq.bases(), seq.length()) : 0;
  uint32 elu    = (el2 == 0) && (el3 == 0) ? encode8bitSequence(dd, seq.bases(), seq.length()) : 0;
  uint32 ell    = 0;

  //  Store the IFF tag and length at the start of the data block.

  if      (el2 > 0) {
    //fprintf(stderr, "  read %9u name '%s' length %5lu 2-bit encoded in %4u bytes.\n", id, seq.ident(), seq.length(), el2);
    tag[0] = '2';
    ell    = el2;
  }
  else if (el3 > 0) {
    //fprintf(stderr, "  read %9u name '%s' length %5lu 3-bit encoded in %4u bytes.\n", id, seq.ident(), seq.length(), el3);
    tag[0] = '3';
    ell    = el3;
  }
  else {
    //fprintf(stderr, "  read %9u name '%s' length %5lu byte encoded in %4u bytes.\n", id, seq.ident(), seq.length(), elu);
    tag[0] = 'U';
    ell    = elu;
  }

  dd = _data + _dataLen;

  memcpy(dd+0,  tag, sizeof(uint8) * 4);
  memcpy(dd+4, &ell, sizeof(uint32));

  //  Clear any pad bytes at the end of the data.

  uint32 padLen = 4 - (ell % 4);

  for (uint32 ii=0; ii<padLen; ii++)
    dd[ell + ii] = 0;

  //  Advance the storage pointer.

  _dataLen += 4 + 4 + ell + padLen;
}



char *
sqCache::sqCache_getSequence(uint32    id) {
  uint32  seqLen = 0;
  uint32  seqMax = 0;
  char   *seq     = nullptr;

  return(sqCache_getSequence(id, seq, seqLen, seqMax));
}



char *
sqCache::sqCache_getSequence(uint32    id,
                             char    *&seq,
                             uint32   &seqLen,
                             uint32   &seqMax) {

  //  If not loaded, load it.

  if (_reads[id]._sData == nullptr)
    loadRead(id);

  //  Decide how many bases are encoded in the encoding and make space to
  //  decode the entire sequence (that is, the untrimmed sequence).

  resizeArray(seq, 0, seqMax, _reads[id]._sLen + 1, _raAct::doNothing);

  //  Decode it.

  char   *cName =  (char *)  (_reads[id]._sData + 0);
  uint32  cLen  = *(uint32 *)(_reads[id]._sData + 4);
  uint8  *chunk     =        (_reads[id]._sData + 8);

  if      (((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode2bitSequence(chunk, cLen, seq, _reads[id]._sLen);

  else if (((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode3bitSequence(chunk, cLen, seq, _reads[id]._sLen);

  else if (((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode8bitSequence(chunk, cLen, seq, _reads[id]._sLen);

  //  If a compressed read, we need to ... compress it.
  //  If not compressed, the (untrimmed) length is exactly basesLen.

  if (_compressed)
    seqLen = homopolyCompress(seq, _reads[id]._sLen, seq);
  else
    seqLen = _reads[id]._sLen;

  //  If a trimmed read, we need to ... trim it.
  //  If not trimmed, seqLen is already set, as is seq, so we're done.

  if (_trimmed) {
    seqLen = _reads[id]._end - _reads[id]._bgn;

    if (_reads[id]._bgn > 0)
      memmove(seq, seq + _reads[id]._bgn, sizeof(char) * seqLen);

    seq[seqLen] = 0;
  }

  //  Return the sequence.

  return(seq);
}



uint32
sqCache::sqCache_mapNameToID(char const *readName) {
  auto  elt = _nameToID.find(std::string(readName));

  if (elt == _nameToID.end())
    return(0);
  else
    return(elt->second);
}





//  Just load all reads.
void
sqCache::sqCache_loadReads(bool verbose) {
  sqCache_loadReads((uint32)0, _readsLen, verbose);
}



void
sqCache::sqCache_loadReads(uint32 bgnID, uint32 endID, bool verbose) {
  uint32  nReads = 0;
  uint64  nBases = 0;

  for (uint32 id=bgnID; id <= endID; id++) {
    if (_reads[id]._sLen > 0) {
      nReads += 1;
      nBases += _reads[id]._end - _reads[id]._bgn;
    }
  }

  if (verbose)
    fprintf(stderr, "Loading %u reads and %lu bases from range %u-%u inclusive.\n",
            nReads, nBases, bgnID, endID);

  //  Allocate a block.

  allocateNewBlock();

  //

  for (uint32 id=bgnID; id <= endID; id++) {
    loadRead(id);

    if ((verbose) && ((id % 4567) == 0)) {
      double  approxSize = ((_dataBlocksLen-1) * _dataMax + _dataLen) / 1024.0 / 1024.0 / 1024.0;

      fprintf(stderr, "Loading %8u < %8u < %8u - %7.2f%% - %.2f GB\r",
              bgnID, id, endID,
              100.0 * (id - bgnID) / (endID - bgnID), approxSize);
    }
  }

  assert(_dataLen <= _dataMax);

  if (verbose) {
    double  approxSize = ((_dataBlocksLen-1) * _dataMax + _dataLen) / 1024.0 / 1024.0 / 1024.0;

    fprintf(stderr, "Loading %8u < %8u < %8u - %7.2f%% - %.2f GB\n",
            bgnID, endID, endID,
            100.0, approxSize);
  }
}



//  Load all the reads in a set of IDs.
void
sqCache::sqCache_loadReads(std::set<uint32> reads, bool verbose) {
  uint32   nToLoad  = reads.size();
  uint32   nLoaded  = 0;
  uint32   nStep    = nToLoad / 100;

  if (verbose)
    fprintf(stderr, "Loading %u reads.\n", nToLoad);

  for (auto it=reads.begin(); it != reads.end(); ++it) {
    loadRead(*it);
    nLoaded++;

    if ((verbose) && ((nLoaded % nStep) == 0))
      fprintf(stderr, "Loading %u reads - %5.1f%%\r", nToLoad, 100.0 * nLoaded / nToLoad);
  }

  if (verbose)
    fprintf(stderr, "\nLoaded " F_SIZE_T " reads.\n", reads.size());
}



//  Load all the reads in a set of IDs, setting age to the second
//  item in the map.
//
//  falconsense gives us a map of readID -> occurrences, which we _could_
//  use to purge reads from the cache when they're no longer used, or to
//  only load reads used more than once.  But we don't do that.
//
//  There was support for 'expiring' reads and removing their data, which we
//  don't do anymore either.
//
void
sqCache::sqCache_loadReads(std::map<uint32, uint32> reads, bool verbose) {
  uint32   nToLoad  = reads.size();
  uint32   nLoaded  = 0;
  uint32   nSkipped = 0;
  uint32   nStep    = nToLoad / 100;

  if (verbose)
    fprintf(stderr, "Loading %u reads.\n", nToLoad);

  for (auto it=reads.begin(); it != reads.end(); ++it) {
    if (it->second > 0) {
      //loadRead(it->first, it->second);
      loadRead(it->first);
      nLoaded++;

    } else {
      nSkipped++;
    }

    if ((verbose) && (((nLoaded + nSkipped) % nStep) == 0))
      fprintf(stderr, "Loading %u reads - %5.1f%%\r", nToLoad, 100.0 * (nLoaded + nSkipped) / nToLoad);
  }

  if (verbose)
    fprintf(stderr, "\nLoaded %u reads; skipped %u singleton reads.\n", nLoaded, nSkipped);
}



//  For trimming, load all the reads in a set of overlaps.
void
sqCache::sqCache_loadReads(ovOverlap *ovl, uint32 nOvl, bool verbose) {
  std::set<uint32>  reads;

  for (uint32 oo=0; oo<nOvl; oo++) {
    reads.insert(ovl[oo].a_iid);
    reads.insert(ovl[oo].b_iid);
  }

  sqCache_loadReads(reads, verbose);
}



//  For correction, load the read the tig represents, and all evidence reads.
void
sqCache::sqCache_loadReads(tgTig *tig, bool verbose) {
  std::set<uint32>  reads;

  reads.insert(tig->tigID());

  for (uint32 oo=0; oo<tig->numberOfChildren(); oo++)
    if (tig->getChild(oo)->isRead() == true)
      reads.insert(tig->getChild(oo)->ident());

  sqCache_loadReads(reads, verbose);
}



//  Load ALL reads in the (possibly compressed) file.
void
sqCache::sqCache_loadReads(char const *filename) {
  dnaSeqFile *readFile = new dnaSeqFile(filename, false);
  dnaSeq      readSeq;

  fprintf(stderr, "-- Loading reads from '%s': %8lu reads.\r", filename, _readsLen);

  while (readFile->loadSequence(readSeq) == true) {
    if ((_readsLen & 0x1ff) == 0x1ff)
      fprintf(stderr, "-- Loading reads from '%s': %8lu reads.\r", filename, _readsLen);

    if (readSeq.length() > 0)
      loadRead(readSeq);
  }

  fprintf(stderr, "-- Loaded  reads from '%s': %8lu reads.\n", filename, _readsLen);

  delete readFile;
}




void
sqCache::sqCache_saveReadToBuffer(writeBuffer *B, uint32 id, sqRead *rd, sqReadDataWriter *wr) {
  sqReadMeta  readMeta;
  sqReadSeq   rawU;   rawU.sqReadSeq_initialize();
  sqReadSeq   rawC;   rawC.sqReadSeq_initialize();
  sqReadSeq   corU;   corU.sqReadSeq_initialize();
  sqReadSeq   corC;   corC.sqReadSeq_initialize();

  //  Load the sequence data, then write it out.

  char    *seq    = nullptr;
  uint32   seqLen = 0;
  uint32   seqMax = 0;

  sqCache_getSequence(id, seq, seqLen, seqMax);
  assert(seq[seqLen] == 0);

  //  Create read metadata and write it.

  readMeta.sqReadMeta_initialize(id);

  B->write(&readMeta, sizeof(sqReadMeta));

  //  Create sequence metadata and write it.

  rawU.sqReadSeq_setLength(seq, seqLen, false);
  rawC.sqReadSeq_setLength(seq, seqLen, true);
  //corU.sqReadSeq_setLength(seq, seqLen, false);
  //corC.sqReadSeq_setLength(seq, seqLen, false);

  B->write(&rawU, sizeof(sqReadSeq));
  B->write(&rawC, sizeof(sqReadSeq));
  B->write(&corU, sizeof(sqReadSeq));
  B->write(&corC, sizeof(sqReadSeq));

  //  Figure out the name of this read.
  //
  //  When reads are loaded from seqStore, names do NOT exist, and the
  //  _IDtoName vector is empty; we pass an empty name to importData().
  //
  //  When reads are loaded from sequence files, names DO exist.  All id's
  //  will have a name; we just ensure that the id is within range.

  char const *name = "";

  if (id < _IDtoName.size())
    name = _IDtoName[id]->first.c_str();

  //  Import the sequence and metadata to the data writer, then write it.

  wr->sqReadDataWriter_importData(name,
                                  seq,
                                  seqLen, 0, seqLen,
                                  &rawU, &rawC, &corU, &corC);

  wr->sqReadDataWriter_writeBlob(B);

  delete [] seq;
}

