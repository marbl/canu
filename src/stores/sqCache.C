
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2019-FEB-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqCache.H"

#include <set>
#include <vector>
#include <algorithm>

using namespace std;



sqCache::sqCache(sqStore *seqStore, sqRead_version version,  uint64 memoryLimit) {

  _seqStore        =  seqStore;
  _nReads          = _seqStore->sqStore_getNumReads();

  _trackAge        = true;
  _trackExpiration = false;

  _version         =  version;

  if (_version == sqRead_latest) {
    if (_seqStore->sqStore_getNumRawReads()       > 0)   _version = sqRead_raw;
    if (_seqStore->sqStore_getNumCorrectedReads() > 0)   _version = sqRead_corrected;
    if (_seqStore->sqStore_getNumTrimmedReads()   > 0)   _version = sqRead_trimmed;
  }

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

  for (uint32 id=0; id <= _nReads; id++) {
    sqRead  *read = _seqStore->sqStore_getRead(id);

    if (_version != sqRead_trimmed) {
      _reads[id]._readLength = read->sqRead_sequenceLength(_version);
      _reads[id]._bgn        = 0;
      _reads[id]._end        = _reads[id]._readLength;
    } else {
      _reads[id]._readLength = read->sqRead_sequenceLength(sqRead_corrected);
      _reads[id]._bgn        = read->sqRead_clearBgn();
      _reads[id]._end        = read->sqRead_clearEnd();
    }

    //_reads[id]._dataAge        = 0;
    _reads[id]._dataExpiration = UINT32_MAX;
    _reads[id]._data           = NULL;

    if (sqCache_getLength(id) > 0) {
      nReads += 1;
      nBases += _reads[id]._readLength;
    }
  }

  fprintf(stderr, "sqCache: found %u %s reads with %lu bases.\n", nReads, toString(_version), nBases);
}



sqCache::~sqCache() {

  //  If we've got a big block of data allocated, reset all the read
  //  data pointers to NULL so they don't try to delete memory that
  //  can't be deleted.

  if (_data)
    for (uint32 ii=0; ii <= _nReads; ii++)
      _reads[ii]._data = NULL;

  //  Now just delete!

  delete [] _reads;
  delete [] _data;
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

  if (_reads[id]._readLength == 0)
    return;

  //fprintf(stderr, "Loading read %u of length %u with expiration %u\n",
  //        id, _reads[id]._readLength, expiration);

  //  Load the encoded blob.

  uint8   *blob     = _seqStore->sqStore_loadReadBlob(id);
  uint8   *bptr     = blob + 8;
  uint8   *rptr     = NULL;
  uint8   *cptr     = NULL;

  //  Find the encoded read data.

  while ((bptr[0] != 'S') ||
         (bptr[1] != 'T') ||
         (bptr[2] != 'O') ||
         (bptr[3] != 'P')) {
    uint32  chunkLen = 4 + 4 + *((uint32 *)bptr + 1);

    if (((bptr[0] == '2') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'R')) ||
        ((bptr[0] == 'U') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'R')))
      rptr = bptr;

    if (((bptr[0] == '2') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'C')) ||
        ((bptr[0] == 'U') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'C')))
      cptr = bptr;

    bptr += chunkLen;
  }

  //  Decide which read to save.  Raw?  Corrected?  Trimmed?

  if (_version == sqRead_raw)
    bptr = rptr;
  else
    bptr = cptr;

  //  Decode how much data we need to save.

  uint32  chunkLen = 4 + 4 + *((uint32 *)bptr + 1);

  //  If we have a gigantic storage space for read data, use that, otherwise,
  //  allocate space for this data.

  if (_data == NULL) {
    _reads[id]._data = new uint8 [chunkLen];
  }

  else {
    if (_dataLen + chunkLen > _dataMax)
      allocateNewBlock();

    _reads[id]._data = _data + _dataLen;
  }

  //  Copy the data and release the blob.

  memcpy(_reads[id]._data, bptr, chunkLen);

  delete [] blob;

  //  Update the pointer to the next free chunk of storage.

  if (_data != NULL) {
    _dataLen += chunkLen;

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

  seqLen = _reads[id]._readLength;

  resizeArray(seq, 0, seqMax, seqLen + 1, resizeArray_doNothing);

  //  Decode it.

  uint8  *bptr     = _reads[id]._data;
  uint32  chunkLen = *((uint32 *)bptr + 1);

  if      (((bptr[0] == '2') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'R')) ||
           ((bptr[0] == '2') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'C')))
    _readData.sqReadData_decode2bit(_reads[id]._data + 8, chunkLen, seq, seqLen);

  else if (((bptr[0] == 'U') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'R')) ||
           ((bptr[0] == 'U') && (bptr[1] == 'S') && (bptr[2] == 'Q') && (bptr[3] == 'C'))) {
    memcpy(seq, _reads[id]._data + 8, seqLen);
    seq[seqLen] = 0;
  }

  //  If a trimmed read, we need to ... trim it.

  if (_version == sqRead_trimmed) {
    seqLen = sqCache_getLength(id);

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
  uint32  nReads = 0;
  uint64  nBases = 0;

  for (uint32 id=0; id <= _nReads; id++) {
    if (_reads[id]._readLength > 0) {
      nReads += 1;
      nBases += _reads[id]._readLength;
    }
  }

  if (verbose)
    fprintf(stderr, "Loading %u reads and %lu bases out of %u reads in the store.\n",
            nReads, nBases, _nReads);

  //  For 50x human, with N's in the sequence, we need 50 * 3 Gbp / 3 bytes.
  //  We'll allocate that in nice 32 MB chunks, 1490 chunks.
  //
  //  We expect to need 'nBases / 4 + nReads' bytes (2-bit) or 'nBases / 3 + nReads'
  //  (3-bit) bytes, which will let us, at least, pre-allocate the pointers to blocks.

  _dataMax       = 32 * 1024 * 1024;
  _dataLen       = 0;
  _data          = new uint8 [_dataMax];

  _dataBlocksLen = 0;
  _dataBlocksMax = 0;
  _dataBlocks    = NULL;

  //  Allocate another block.

  allocateNewBlock();

  //

  for (uint32 id=0; id <= _nReads; id++) {
    loadRead(id);

    if ((verbose) && ((id % 4567) == 0)) {
      double  approxSize = ((_dataBlocksLen-1) * _dataMax + _dataLen) / 1024.0 / 1024.0 / 1024.0;

      fprintf(stderr, " %9u - %7.2f%% reads - %.2f GB\r",
              id, 100.0 * id / _nReads, approxSize);
    }
  }

  assert(_dataLen <= _dataMax);

  if (verbose) {
    double  approxSize = ((_dataBlocksLen-1) * _dataMax + _dataLen) / 1024.0 / 1024.0 / 1024.0;

    fprintf(stderr, " %9u - %7.2f%% reads - %.2f GB\r",
            _nReads, 100.0, approxSize);
  }
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
