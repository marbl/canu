
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

using namespace std;



sqCache::sqCache(sqStore *seqStore, sqRead_which which,  uint64 memoryLimit) {

  _seqStore        =  seqStore;
  _nReads          = _seqStore->sqStore_lastReadID();

  _trackAge        = true;
  _trackExpiration = false;
  _noMoreLoads     = false;

  _which           = which;
  _compressed      = ((_which & sqRead_compressed) == sqRead_unset) ? false : true;
  _trimmed         = ((_which & sqRead_trimmed)    == sqRead_unset) ? false : true;

  _memoryLimit   = memoryLimit * 1024 * 1024 * 1024;

  if (_memoryLimit == 0) {
    _trackAge    = false;
    _memoryLimit = UINT64_MAX;
  }

  _reads         = new sqCacheEntry [_nReads + 1];

  _dataLen       = 0;
  _dataMax       = 0;
  _data          = NULL;

  _dataBlocksLen = 0;
  _dataBlocksMax = 0;
  _dataBlocks    = NULL;

  uint32  nReads = 0;
  uint64  nBases = 0;

  for (uint32 id=1; id <= _nReads; id++) {
    if (_seqStore->sqStore_isIgnoredRead(id, _which) == true) {
      _reads[id]._basesLength    = 0;
      _reads[id]._bgn            = 0;
      _reads[id]._end            = 0;
      _reads[id]._dataExpiration = UINT32_MAX;
      _reads[id]._data           = NULL;
      continue;
    }

    //  _basesLength is the length of the sequence encoded in _data.  It is NOT
    //  the length of the read we will be returning.

    if (_which & sqRead_raw)
      _reads[id]._basesLength = _seqStore->sqStore_getReadLength(id, sqRead_raw);

    if (_which & sqRead_corrected)
      _reads[id]._basesLength = _seqStore->sqStore_getReadLength(id, sqRead_corrected);

    //  Set bgn and end based on trimming status.

    _reads[id]._bgn = (_trimmed == true) ? _seqStore->sqStore_getClearBgn(id, _which) : 0;
    _reads[id]._end = (_trimmed == true) ? _seqStore->sqStore_getClearEnd(id, _which) : _seqStore->sqStore_getReadLength(id, _which);

    //  Set the age, expiration and clear the data pointer.

    //_reads[id]._dataAge        = 0;
    _reads[id]._dataExpiration = UINT32_MAX;
    _reads[id]._data           = NULL;

    //  Do some accounting.

    if (sqCache_getLength(id) > 0) {
      nReads += 1;
      nBases += _reads[id]._end - _reads[id]._bgn;
    }
  }

  fprintf(stderr, "sqCache: found %u %s reads with %lu bases.\n", nReads, toString(_which), nBases);
}



sqCache::~sqCache() {

  //  If we've got a big block of data allocated, reset all the read
  //  data pointers to NULL so they don't try to delete memory that
  //  can't be deleted.

  if (_data)
    for (uint32 ii=0; ii <= _nReads; ii++)
      _reads[ii]._data = NULL;

  delete [] _reads;

  //  Now just delete!

  for (uint32 ii=0; ii<_dataBlocksLen; ii++)
    delete [] _dataBlocks[ii];

  delete [] _dataBlocks;
}





void
sqCache::loadRead(uint32 id, uint32 expiration) {

  //  Reset the age and/or expiration of this read.

  if (_trackAge)
    _reads[id]._dataExpiration = 0;

  if (_trackExpiration)
    _reads[id]._dataExpiration = expiration;

  //  If already loaded, don't load it again.

  if (_reads[id]._data != NULL)
    return;

  //  If no read to load, don't load it.

  if (_reads[id]._basesLength == 0)
    return;

  //fprintf(stderr, "Loading read %u of length %u with expiration %u\n",
  //        id, _reads[id]._basesLength, expiration);

  assert(_noMoreLoads == false);

  //  Load the encoded blob, without decoding it.

  _read.sqRead_fetchBlob(_seqStore->sqStore_getReadBuffer(id));

  //  Find the encoded read data.  This mirrors sqRead_loadFromBuffer.

  uint32   blobPos  = 0;
  uint8   *rptr     = NULL;
  uint8   *cptr     = NULL;

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

  if (_data == NULL) {
    _reads[id]._data = new uint8 [blen];
  }

  else {
    if (_dataLen + blen > _dataMax)
      allocateNewBlock();

    _reads[id]._data = _data + _dataLen;
  }

  //  Copy the data and release the blob.

  memcpy(_reads[id]._data, bptr, blen);

  //  Update the pointer to the next free chunk of storage.

  if (_data != NULL) {
    _dataLen += blen;

    assert(_dataLen <= _dataMax);
  }
}



void
sqCache::removeRead(uint32 id) {

  if (_data == NULL)
    delete [] _reads[id]._data;

  _reads[id]._data           = NULL;
  //_reads[id]._dataAge        = 0;
  _reads[id]._dataExpiration = 0;
}



char *
sqCache::sqCache_getSequence(uint32    id) {
  uint32  seqLen = 0;
  uint32  seqMax = 0;
  char   *seq     = NULL;

  return(sqCache_getSequence(id, seq, seqLen, seqMax));
}



char *
sqCache::sqCache_getSequence(uint32    id,
                             char    *&seq,
                             uint32   &seqLen,
                             uint32   &seqMax) {

  //  If not loaded, load it.

  if (_reads[id]._data == NULL)
    loadRead(id);

  //  Decide how many bases are encoded in the encoding and make space to
  //  decode the entire sequence (that is, the untrimmed sequence).

  resizeArray(seq, 0, seqMax, _reads[id]._basesLength + 1, resizeArray_doNothing);

  //  Decode it.

  char   *cName =  (char *)  (_reads[id]._data + 0);
  uint32  cLen  = *(uint32 *)(_reads[id]._data + 4);
  uint8  *chunk     =        (_reads[id]._data + 8);

  if      (((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == '2') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode2bitSequence(chunk, cLen, seq, _reads[id]._basesLength);

  else if (((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == '3') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode3bitSequence(chunk, cLen, seq, _reads[id]._basesLength);

  else if (((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'R')) ||
           ((cName[0] == 'U') && (cName[1] == 'S') && (cName[2] == 'Q') && (cName[3] == 'C')))
    decode8bitSequence(chunk, cLen, seq, _reads[id]._basesLength);

  //  If a compressed read, we need to ... compress it.
  //  If not compressed, the (untrimmed) length is exactly basesLen.

  if (_compressed)
    seqLen = homopolyCompress(seq, _reads[id]._basesLength, seq);
  else
    seqLen = _reads[id]._basesLength;

  //  If a trimmed read, we need to ... trim it.
  //  If not trimmed, seqLen is already set, as is seq, so we're done.

  if (_trimmed) {
    seqLen = _reads[id]._end - _reads[id]._bgn;

    if (_reads[id]._bgn > 0)
      memmove(seq, seq + _reads[id]._bgn, sizeof(char) * seqLen);

    seq[seqLen] = 0;
  }

  //  If we're tracking age, reset the age to zero.

  if (_trackAge)
    _reads[id]._dataExpiration = 0;

  //  If we're tracking expiration dates, release the data if we're done.

  if ((_trackExpiration) && (--_reads[id]._dataExpiration == 0)) {
    //fprintf(stderr, "READ %u expired.\n", id);
    removeRead(id);
  }

  //  Return the sequence.

  return(seq);
}



void
sqCache::increaseAge(void) {
  if (_trackAge == false)
    return;

  for (uint32 id=0; id <= _nReads; id++)   //  Add one to the age of
    if (_reads[id]._data)                  //  any read that is loaded.
      _reads[id]._dataExpiration++;
}



//  Just load all reads.
void
sqCache::sqCache_loadReads(bool verbose) {
  sqCache_loadReads((uint32)0, _nReads, verbose);

  _noMoreLoads = true;
}



void
sqCache::sqCache_loadReads(uint32 bgnID, uint32 endID, bool verbose) {
  uint32  nReads = 0;
  uint64  nBases = 0;

  for (uint32 id=bgnID; id <= endID; id++) {
    if (_reads[id]._basesLength > 0) {
      nReads += 1;
      nBases += _reads[id]._end - _reads[id]._bgn;
    }
  }

  if (verbose)
    fprintf(stderr, "Loading %u reads and %lu bases from range %u-%u inclusive.\n",
            nReads, nBases, bgnID, endID);

  //  For 50x human, with N's in the sequence, we need 50 * 3 Gbp / 3 bytes.
  //  We'll allocate that in nice 32 MB chunks, 1490 chunks.
  //
  //  Don't bother pre-allocation dataBlocks; easy enough to do that on the
  //  fly.

  _dataMax       = 32 * 1024 * 1024;
  _dataLen       = 0;
  _data          = NULL;

  _dataBlocksLen = 0;
  _dataBlocksMax = 0;
  _dataBlocks    = NULL;

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

  _noMoreLoads = true;
}



//  Load all the reads in a set of IDs.
void
sqCache::sqCache_loadReads(set<uint32> reads, bool verbose) {
  uint32   nToLoad  = reads.size();
  uint32   nLoaded  = 0;
  uint32   nStep    = nToLoad / 100;

  if (verbose)
    fprintf(stderr, "Loading %u reads.\n", nToLoad);

  for (set<uint32>::iterator it=reads.begin(); it != reads.end(); ++it) {
    loadRead(*it);
    nLoaded++;

    if ((verbose) && ((nLoaded % nStep) == 0))
      fprintf(stderr, "Loading %u reads - %5.1f%%\r", nToLoad, 100.0 * nLoaded / nToLoad);
  }

  if (verbose)
    fprintf(stderr, "\nLoaded " F_SIZE_T " reads.\n", reads.size());

  _noMoreLoads = true;
}



//  Load all the reads in a set of IDs, setting age to the second
//  item in the map.
void
sqCache::sqCache_loadReads(map<uint32, uint32> reads, bool verbose) {
  uint32   nToLoad  = reads.size();
  uint32   nLoaded  = 0;
  uint32   nSkipped = 0;
  uint32   nStep    = nToLoad / 100;

  if (verbose)
    fprintf(stderr, "Loading %u reads.\n", nToLoad);

  _trackExpiration = true;

  for (map<uint32,uint32>::iterator it=reads.begin(); it != reads.end(); ++it) {
    if (it->second > 0) {
      loadRead(it->first, it->second);
      nLoaded++;

    } else {
      nSkipped++;
    }

    if ((verbose) && (((nLoaded + nSkipped) % nStep) == 0))
      fprintf(stderr, "Loading %u reads - %5.1f%%\r", nToLoad, 100.0 * (nLoaded + nSkipped) / nToLoad);
  }

  if (verbose)
    fprintf(stderr, "\nLoaded %u reads; skipped %u singleton reads.\n", nLoaded, nSkipped);

  _noMoreLoads = true;
}



//  For trimming, load all the reads in a set of overlaps.
void
sqCache::sqCache_loadReads(ovOverlap *ovl, uint32 nOvl, bool verbose) {
  set<uint32>     reads;

  increaseAge();

  for (uint32 oo=0; oo<nOvl; oo++) {
    reads.insert(ovl[oo].a_iid);
    reads.insert(ovl[oo].b_iid);
  }

  sqCache_loadReads(reads, verbose);

  _noMoreLoads = true;
}



//  For correction, load the read the tig represents, and all evidence reads.
void
sqCache::sqCache_loadReads(tgTig *tig, bool verbose) {
  set<uint32>     reads;

  increaseAge();

  reads.insert(tig->tigID());

  for (uint32 oo=0; oo<tig->numberOfChildren(); oo++)
    if (tig->getChild(oo)->isRead() == true)
      reads.insert(tig->getChild(oo)->ident());

  sqCache_loadReads(reads, verbose);

  _noMoreLoads = true;
}





#if 0
void
sqCache::sqCache_purgeReads(void) {
  uint32  maxAge     = 0;
  uint64  memoryUsed = 0;

  //  Find maxAge, and sum memory used

  for (uint32 rr=0; rr <= _nReads; rr++) {
    if (maxAge < readAge[rr])
      maxAge = readAge[rr];

    memoryUsed += readLen[rr];
  }

  //  Purge oldest until memory is below watermark

  while ((_memoryLimit < memoryUsed) &&
         (maxAge > 1)) {
    fprintf(stderr, "purgeReads()--  used " F_U64 "MB limit " F_U64 "MB -- purge age " F_U32 "\n", memoryUsed >> 20, _memoryLimit >> 20, maxAge);

    for (uint32 rr=0; rr <= _nReads; rr++) {
      if (maxAge == readAge[rr]) {
        memoryUsed -= readLen[rr];

        delete [] readSeqFwd[rr];  readSeqFwd[rr] = NULL;

        readLen[rr] = 0;
        readAge[rr] = 0;
      }
    }

    maxAge--;
  }
}
#endif
