
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
 *  This file is derived from:
 *
 *    src/stores/gkStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-10
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2015-DEC-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "sqStore.H"

#include "files.H"


sqStore       *sqStore::_instance      = NULL;
uint32          sqStore::_instanceCount = 0;

sqRead_version  sqRead_defaultVersion = sqRead_latest;



static
uint8 *
sqStore_loadBlobFromStream(FILE *file) {
  char    tag[5];
  uint32  size;

  //  Ideally, we'd do one read to get the whole blob.  Without knowing
  //  the length, we're forced to do two.

  loadFromFile(tag,  "sqStore::sqStore_loadDataFromFile::blob", 4, file);
  loadFromFile(size, "sqStore::sqStore_loadDataFromFile::size",    file);

  uint8 *blob = new uint8 [8 + size];

  memcpy(blob,    tag,  sizeof(uint8)  * 4);
  memcpy(blob+4, &size, sizeof(uint32) * 1);

  loadFromFile(blob+8, "sqStore::sqStore_loadDataFromFile::blob", size, file);

  return(blob);
}



void
sqRead::sqRead_loadDataFromStream(sqReadData *readData, FILE *file) {
  uint8 *blob = sqStore_loadBlobFromStream(file);

  readData->sqReadData_loadFromBlob(blob);

  delete [] blob;
}



sqRead *
sqStore::sqStore_getRead(uint32 id) {

  if (sqStore_readInPartition(id) == false)
    return(NULL);

  if (sqStore_readInPartition(id) == false)
    fprintf(stderr, "getRead()--  access to read %u in partition %u is not allowed when partition %u is loaded.\n",
            id, _readIDtoPartitionID[id], _partitionID), assert(0);

  sqRead *read = _reads + (((_readIDtoPartitionID     != NULL) &&
                            (_readIDtoPartitionID[id] == _partitionID)) ? _readIDtoPartitionIdx[id] : id);

  if (sqStore_getNumCorrectedReads() > 0)     //  If there are corrected or trimmed reads in the store,
    read->_cExists = true;                    //  set the flags so the read can return the appropriate data.

  if (sqStore_getNumTrimmedReads() > 0)
    read->_tExists = true;

  return(read);
}



uint8 *
sqStore::sqStore_loadReadBlob(uint32 readID) {

  //  If partitioned data, copy from the already-in-core data.

  assert(_blobsData == NULL);

  //  Otherwise, read from disk.

  uint32   tnum = omp_get_thread_num();

  assert(tnum < _blobsFilesMax);

  sqRead  *read = sqStore_getRead(readID);
  FILE    *file = _blobsFiles[tnum].getFile(_storePath, read);

  return(sqStore_loadBlobFromStream(file));
}



void
sqStore::sqStore_loadReadData(sqRead *read, sqReadData *readData) {

  readData->_read    = read;
  readData->_library = sqStore_getLibrary(read->sqRead_libraryID());

  //  If partitioned data, we can load from the already-in-core data.

  if (_blobsData) {
    readData->sqReadData_loadFromBlob(_blobsData + read->sqRead_mByte());
    return;
  }

  //  Otherwise, we need to read from disk.

  uint32   tnum = omp_get_thread_num();

  assert(tnum < _blobsFilesMax);

  read->sqRead_loadDataFromStream(readData, _blobsFiles[tnum].getFile(_storePath, read));
}



void
sqStore::sqStore_loadReadData(uint32  readID, sqReadData *readData) {

  sqStore_loadReadData(sqStore_getRead(readID), readData);
}



//  Dump a block of encoded data to disk, then update the sqRead to point to it.
//
void
sqStore::sqStore_stashReadData(sqReadData *data) {

  data->sqReadData_encodeBlob();                            //  Encode the data.

  _blobsWriter->writeData(data->_blob, data->_blobLen);     //  Write the data.

  data->_read->_mSegm     = _blobsWriter->writtenIndex();       //  Remember where it was written.
  data->_read->_mByteHigh = _blobsWriter->writtenPosition() >> 32;
  data->_read->_mByteLow  = _blobsWriter->writtenPosition() & 0xffffffffllu;
  data->_read->_mPart     = _partitionID;                       //  (0 if not partitioned)
}



//  Load read metadata and data from a stream.
//
void
sqStore::sqStore_loadReadFromStream(FILE *S, sqRead *read, sqReadData *readData) {
  char    tag[5];
  uint32  size;

  //  Mark this as a read.  Needed for tgTig::loadFromStreamOrLayout(), and loading this stuff in
  //  utgcns.

  loadFromFile(tag, "sqStore::sqStore_loadReadFromStream::tag", 4, S);

  if (strncmp(tag, "READ", 4) != 0)
    fprintf(stderr, "Failed to load sqRead, got tag '%c%c%c%c' (0x%02x 0x%02x 0x%02x 0x%02x), expected 'READ'.\n",
            tag[0], tag[1], tag[2], tag[3],
            tag[0], tag[1], tag[2], tag[3]), exit(1);

  //  Load the read metadata

  loadFromFile(read, "sqStore::sqStore_loadReadFromStream::read", S);

  //  Load the read data.
  //
  //  Sadly, we don't have an actual sqStore here (usually), so we don't have a sqLibrary hanging around.

  readData->_read    = read;
  readData->_library = NULL;  //sqStore_getLibrary(read->sqRead_libraryID());

  read->sqRead_loadDataFromStream(readData, S);
}



//  Dump the read metadata and read data to a stream.
//
void
sqStore::sqStore_saveReadToStream(FILE *S, uint32 id) {
  sqRead  *read   = sqStore_getRead(id);
  uint8   *blob   = NULL;
  uint32  blobLen = 0;

  //  If partitioned -- if _blobsData exists -- we can grab the blob from there.  Otherwise,
  //  we need to load it from dist.

  if (_blobsData) {
    blob = _blobsData + read->sqRead_mByte();
  }

  else {
    uint32  tnum = omp_get_thread_num();

    assert(tnum < _blobsFilesMax);

    blob = sqStore_loadBlobFromStream(_blobsFiles[tnum].getFile(_storePath, read));
  }

  blobLen = 8 + *((uint32 *)blob + 1);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the blob to the stream

  fprintf(S, "READ");
  writeToFile(read, "sqStore::sqStore_saveReadToStream::read",          S);
  writeToFile(blob, "sqStore::sqStore_saveReadToStream::blob", blobLen, S);

  //  And cleanup.

  if (_blobsData == NULL)
    delete [] blob;
}



void
sqReadData::sqReadData_setName(char *H) {
  uint32  Hlen = strlen(H) + 1;

  resizeArray(_name, 0, _nameAlloc, Hlen, resizeArray_doNothing);

  memcpy(_name, H, sizeof(char) * Hlen);
}



//  Based on the library type, and presence of read data, either load the
//  sequence into the 'raw' storage, the 'corrected' storage, or maybe both.
//
void
sqReadData::sqReadData_setBasesQuals(char  *S,
                                     uint8 *Q) {
  uint32      Slen  = strlen(S) + 1;

  //  For PacBio HiFi, our correction amounts to stipping homopolymer runs,
  //  and we can do that here.

  if (_library->sqLibrary_readType() == SQ_READTYPE_PACBIO_HIFI) {
    resizeArray(_rseq, 0, _rseqAlloc, Slen, resizeArray_doNothing);    //  Load the raw version.
    resizeArray(_rqlt, 0, _rqltAlloc, Slen, resizeArray_doNothing);

    memcpy(_rseq, S, sizeof(char)  * Slen);
    memcpy(_rqlt, Q, sizeof(uint8) * Slen);

    uint32  cc = 0;                  //  NOTE:  Also used in utility/sequence-extract.C
    uint32  rr = 1;

    while (rr < Slen) {
      if (S[cc] == S[rr])
        rr++;
      else {
        S[++cc] = S[rr  ];
        Q[  cc] = S[rr++];
      }
    }

    Slen = cc + 1;

    S[Slen] = 0;
    Q[Slen] = 0;

    resizeArray(_cseq, 0, _cseqAlloc, Slen, resizeArray_doNothing);   //  Load the corrected
    resizeArray(_cqlt, 0, _cqltAlloc, Slen, resizeArray_doNothing);   //  version.

    memcpy(_cseq, S, sizeof(char)  * Slen);
    memcpy(_cqlt, Q, sizeof(uint8) * Slen);

    return;
  }

  //  Just gross.  sqLibrary_readType() is incorrect.  It's the type of the
  //  read initially loaded into the store, not the type we're currently
  //  loading (a remnant of having one store for raw, corrected and trimmed
  //  reads).
  //
  //  Instead, we need to check both the original type loaded and the
  //  presence of that data to decide what to do.
  //
  //  So, if the library is a 'raw' type but there is no rseq, load as raw.
  //  Otherwise, load as corrected.

  if ((_rseq == NULL) && ((_library->sqLibrary_readType() == SQ_READTYPE_PACBIO_RAW) ||
                          (_library->sqLibrary_readType() == SQ_READTYPE_NANOPORE_RAW))) {
    resizeArray(_rseq, 0, _rseqAlloc, Slen, resizeArray_doNothing);
    resizeArray(_rqlt, 0, _rqltAlloc, Slen, resizeArray_doNothing);

    memcpy(_rseq, S, sizeof(char)  * Slen);
    memcpy(_rqlt, Q, sizeof(uint8) * Slen);
  }

  else {
    if (_read->_cExists)
      fprintf(stderr, "sqReadData_setBasesQuals()- read %u has existing cseq of length %u, replacing with length %u\n",
              _read->_readID, _read->_cseqLen, (uint32)strlen(S));

    resizeArray(_cseq, 0, _cseqAlloc, Slen, resizeArray_doNothing);
    resizeArray(_cqlt, 0, _cqltAlloc, Slen, resizeArray_doNothing);

    memcpy(_cseq, S, sizeof(char)  * Slen);
    memcpy(_cqlt, Q, sizeof(uint8) * Slen);
  }
}




//  Store the 'len' bytes of data in 'dat' into the class-managed _blob data block.
//  Ensures that the _blob block is appropriately padded to maintain 32-bit alignment.
//
void
sqReadData::sqReadData_encodeBlobChunk(char const *tag,
                                       uint32      len,
                                       void       *dat) {

  //  Allocate an initial blob if we don't have one

  if (_blobMax == 0) {
    _blobLen = 0;
    _blobMax = 1048576;
    _blob    = new uint8 [_blobMax];
  }

  //  Or make it bigger

  while (_blobMax <= _blobLen + 8 + len) {
    _blobMax *= 2;
    uint8 *b  = new uint8 [_blobMax];
    memcpy(b, _blob, sizeof(uint8) * _blobLen);
    delete [] _blob;
    _blob = b;
  }

  //  Figure out how much padding we need to add

  uint32 pad = 4 - (len % 4);

  if (pad == 4)
    pad = 0;

  //  Copy in the chunk id and padded length

  len += pad;

  memcpy(_blob + _blobLen,  tag, sizeof(uint8) * 4);    _blobLen += sizeof(uint8) * 4;
  memcpy(_blob + _blobLen, &len, sizeof(uint32));       _blobLen += sizeof(uint32);

  len -= pad;

  //  Then the unpadded data and any padding.

  memcpy(_blob + _blobLen,  dat, sizeof(uint8) * len);  _blobLen += sizeof(uint8) * len;

  if (pad > 2)  _blob[_blobLen++] = 0;
  if (pad > 1)  _blob[_blobLen++] = 0;
  if (pad > 0)  _blob[_blobLen++] = 0;

  //  Finally, update the total blob length.

  _blobLen -= 8;

  memcpy(_blob + 4, &_blobLen, sizeof(uint32));

  _blobLen += 8;
}



void
sqReadData::sqReadData_encodeBlob(void) {

  _blobLen = 0;

  //  Discover which sequence(s) we're encoding.

  _read->_rseqLen = (_rseq == NULL) ? 0 : strlen(_rseq);
  _read->_cseqLen = (_cseq == NULL) ? 0 : strlen(_cseq);

  //  Encode the data into chunks in the blob.

  sqReadData_encodeBlobChunk("BLOB", 0,  NULL);

  sqReadData_encodeBlobChunk("NAME", strlen(_name), _name);

  //  Compute the preferred encodings.  If either fail, the length is set to zero, and the
  //  non-preferred encoding will be computed.  If this too fails, sequences/qualities will be
  //  stored unencoded.

  if (_read->_rseqLen > 0) {
    uint8   *rseq = NULL;
    uint8   *rqlt = NULL;

    uint32  rseq2Len =                   sqReadData_encode2bit(rseq, _rseq, _read->_rseqLen);
    uint32  rseq3Len = (rseq2Len == 0) ? sqReadData_encode3bit(rseq, _rseq, _read->_rseqLen) : 0;

    uint32  rqv      = sqReadData_encodeConstantQV(_rqlt, _read->_rseqLen, _library->sqLibrary_defaultQV());

    uint32  rqlt4Len = ((rqv == 255))                    ? sqReadData_encode4bit(rqlt, _rqlt, _read->_rseqLen) : 0;
    uint32  rqlt5Len = ((rqv == 255) && (rqlt4Len == 0)) ? sqReadData_encode5bit(rqlt, _rqlt, _read->_rseqLen) : 0;

    if      (rseq2Len > 0)
      sqReadData_encodeBlobChunk("2SQR",         rseq2Len, rseq);    //  Two-bit encoded sequence (ACGT only)
    else if (rseq3Len > 0)
      sqReadData_encodeBlobChunk("3SQR",         rseq3Len, rseq);    //  Three-bit encoded sequence (ACGTN)
    else
      sqReadData_encodeBlobChunk("USQR", _read->_rseqLen, _rseq);    //  Unencoded sequence

    if      (rqv < 255)
      sqReadData_encodeBlobChunk("1QVR",                 4, &rqv);   //  Constant QV for every base
    else if (rqlt4Len > 0)
      sqReadData_encodeBlobChunk("4QVR",         rqlt4Len, rqlt);    //  Four-bit (0-15) encoded QVs
    else if (rqlt5Len > 0)
      sqReadData_encodeBlobChunk("5QVR",         rqlt5Len, rqlt);    //  Five-bit (0-32) encoded QVs
    else
      sqReadData_encodeBlobChunk("UQVR", _read->_rseqLen, _rqlt);    //  Unencoded quality

    delete [] rseq;
    delete [] rqlt;
  }

  if (_read->_cseqLen > 0) {
    uint8   *cseq = NULL;
    uint8   *cqlt = NULL;

    uint32  cseq2Len =                   sqReadData_encode2bit(cseq, _cseq, _read->_cseqLen);
    uint32  cseq3Len = (cseq2Len == 0) ? sqReadData_encode3bit(cseq, _cseq, _read->_cseqLen) : 0;

    uint32  cqv      = sqReadData_encodeConstantQV(_cqlt, _read->_cseqLen, _library->sqLibrary_defaultQV());

    uint32  cqlt4Len = ((cqv == 255))                    ? sqReadData_encode4bit(cqlt, _cqlt, _read->_cseqLen) : 0;
    uint32  cqlt5Len = ((cqv == 255) && (cqlt4Len == 0)) ? sqReadData_encode5bit(cqlt, _cqlt, _read->_cseqLen) : 0;

    if      (cseq2Len > 0)
      sqReadData_encodeBlobChunk("2SQC",         cseq2Len, cseq);    //  Two-bit encoded sequence (ACGT only)
    else if (cseq3Len > 0)
      sqReadData_encodeBlobChunk("3SQC",         cseq3Len, cseq);    //  Three-bit encoded sequence (ACGTN)
    else
      sqReadData_encodeBlobChunk("USQC", _read->_cseqLen, _cseq);    //  Unencoded sequence

    if      (cqv < 255)
      sqReadData_encodeBlobChunk("1QVC",                 4, &cqv);   //  Constant QV for every base
    else if (cqlt4Len > 0)
      sqReadData_encodeBlobChunk("4QVC",         cqlt4Len, cqlt);    //  Four-bit (0-15) encoded QVs
    else if (cqlt5Len > 0)
      sqReadData_encodeBlobChunk("5QVC",         cqlt5Len, cqlt);    //  Five-bit (0-32) encoded QVs
    else
      sqReadData_encodeBlobChunk("UQVC", _read->_cseqLen, _cqlt);    //  Unencoded quality

    delete [] cseq;
    delete [] cqlt;
  }

  sqReadData_encodeBlobChunk("STOP", 0,  NULL);
}




//  Lowest level function to load data into a read.
//
void
sqReadData::sqReadData_loadFromBlob(uint8 *blob) {
  char    chunk[5];
  uint32  chunkLen = 0;

  //  Make sure that our blob is actually a blob.

  if ((blob[0] != 'B') && (blob[1] != 'L') && (blob[2] != 'O') && (blob[3] != 'B'))
    fprintf(stderr, "Index error in read " F_U32 " mSegm " F_U64 " mByte " F_U64 " mPart " F_U64 " expected BLOB, got %02x %02x %02x %02x '%c%c%c%c'\n",
            _read->sqRead_readID(),
            _read->sqRead_mSegm(), _read->sqRead_mByte(), _read->sqRead_mPart(),
            blob[0], blob[1], blob[2], blob[3],
            blob[0], blob[1], blob[2], blob[3]);
  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Skip over the BLOB and blobLen.  We probably should track blobLen.

  blob += 8;

  //  Resize strings.

  resizeArray(_rseq, 0, _rseqAlloc, _read->_rseqLen+1, resizeArray_doNothing);
  resizeArray(_rqlt, 0, _rqltAlloc, _read->_rseqLen+1, resizeArray_doNothing);

  resizeArray(_cseq, 0, _cseqAlloc, _read->_cseqLen+1, resizeArray_doNothing);
  resizeArray(_cqlt, 0, _cqltAlloc, _read->_cseqLen+1, resizeArray_doNothing);

  //  Decode the blob data.

  while ((blob[0] != 'S') ||
         (blob[1] != 'T') ||
         (blob[2] != 'O') ||
         (blob[3] != 'P')) {
    chunk[0] = blob[0];
    chunk[1] = blob[1];
    chunk[2] = blob[2];
    chunk[3] = blob[3];
    chunk[4] = 0;

    chunkLen = *((uint32 *)blob + 1);

    if      (strncmp(chunk, "NAME", 4) == 0) {
      resizeArray(_name, 0, _nameAlloc, chunkLen + 1, resizeArray_doNothing);
      memcpy(_name, blob + 8, chunkLen);
      _name[chunkLen] = 0;
    }

    else if (strncmp(chunk, "2SQR", 4) == 0) {
      sqReadData_decode2bit(blob + 8, chunkLen, _rseq, _read->_rseqLen);
    }
    else if (strncmp(chunk, "3SQR", 4) == 0) {
      sqReadData_decode3bit(blob + 8, chunkLen, _rseq, _read->_rseqLen);
    }
    else if (strncmp(chunk, "USQR", 4) == 0) {
      assert(_read->_rseqLen <= chunkLen);
      assert(_read->_rseqLen <= _rseqAlloc);
      memcpy(_rseq, blob + 8, _read->_rseqLen);
      _rseq[_read->_rseqLen] = 0;
    }

    else if (strncmp(chunk, "4QVR", 4) == 0) {
      sqReadData_decode4bit(blob + 8, chunkLen, _rqlt, _read->_rseqLen);
    }
    else if (strncmp(chunk, "5QVR", 4) == 0) {
      sqReadData_decode5bit(blob + 8, chunkLen, _rqlt, _read->_rseqLen);
    }
    else if (strncmp(chunk, "UQVR", 4) == 0) {
      assert(_read->_rseqLen <= chunkLen);
      assert(_read->_rseqLen <= _rqltAlloc);
      memcpy(_rqlt, blob + 8, _read->_rseqLen);
      _rqlt[_read->_rseqLen] = 0;
    }
    else if (strncmp(chunk, "1QVR", 4) == 0) {
      for (uint32 qval = *((uint32 *)blob + 2), ii=0; ii<_read->_rseqLen; ii++)
        _rqlt[ii] = qval;
    }

    else if (strncmp(chunk, "2SQC", 4) == 0) {
      sqReadData_decode2bit(blob + 8, chunkLen, _cseq, _read->_cseqLen);
    }
    else if (strncmp(chunk, "3SQC", 4) == 0) {
      sqReadData_decode3bit(blob + 8, chunkLen, _cseq, _read->_cseqLen);
    }
    else if (strncmp(chunk, "USQC", 4) == 0) {
      assert(_read->_cseqLen <= chunkLen);
      assert(_read->_cseqLen <= _cseqAlloc);
      memcpy(_cseq, blob + 8, _read->_cseqLen);
      _cseq[_read->_cseqLen] = 0;
    }

    else if (strncmp(chunk, "4QVC", 4) == 0) {
      sqReadData_decode4bit(blob + 8, chunkLen, _cqlt, _read->_cseqLen);
    }
    else if (strncmp(chunk, "5QVC", 4) == 0) {
      sqReadData_decode5bit(blob + 8, chunkLen, _cqlt, _read->_cseqLen);
    }
    else if (strncmp(chunk, "UQVC", 4) == 0) {
      assert(_read->_cseqLen <= chunkLen);
      assert(_read->_cseqLen <= _cqltAlloc);
      memcpy(_cqlt, blob + 8, _read->_cseqLen);
      _cqlt[_read->_cseqLen] = 0;
    }
    else if (strncmp(chunk, "1QVC", 4) == 0) {
      for (uint32 qval = *((uint32 *)blob + 2), ii=0; ii<_read->_cseqLen; ii++)
        _cqlt[ii] = qval;
    }

    else {
      fprintf(stderr, "sqRead::sqRead_loadDataFromBlob()--  unknown chunk type %02x %02x %02x %02x '%c%c%c%c' skipped\n",
              chunk[0], chunk[1], chunk[2], chunk[3],
              chunk[0], chunk[1], chunk[2], chunk[3]);
    }

    blob += 4 + 4 + chunkLen;
  }

  //  Decide what data is active.

  if      (_read->_tExists) {
    _aseq = _tseq = _cseq + _read->_clearBgn;
    _aqlt = _tqlt = _cqlt + _read->_clearBgn;
  }

  else if (_read->_cExists) {
    _aseq = _cseq;
    _aqlt = _cqlt;
  }

  else {
    _aseq = _rseq;
    _aqlt = _rqlt;
  }
}



sqLibrary *
sqStore::sqStore_addEmptyLibrary(char const *name) {

  assert(_info.sqInfo_numLibraries() <= _librariesAlloc);

  //  Just like with reads below, there is no _libraries[0] element.

  _info.sqInfo_addLibrary();

  increaseArray(_libraries, _info.sqInfo_numLibraries(), _librariesAlloc, 128);

  //  Initialize the new library.

  _libraries[_info.sqInfo_numLibraries()] = sqLibrary();
  _libraries[_info.sqInfo_numLibraries()]._libraryID = _info.sqInfo_numLibraries();

  //  Bullet proof the library name - so we can make files with this prefix.

  char   *libname    = _libraries[_info.sqInfo_numLibraries()]._libraryName;
  uint32  libnamepos = 0;

  memset(libname, 0, sizeof(char) * LIBRARY_NAME_SIZE);

  for (char const *orig=name; *orig; orig++) {
    if        (*orig == '/') {
      libname[libnamepos++] = '_';

    } else if (isspace(*orig) == 0) {
      libname[libnamepos++] = *orig;

    } else {
      libname[libnamepos++] = '_';
    }

    if (libnamepos >= LIBRARY_NAME_SIZE) {
      libname[LIBRARY_NAME_SIZE-1] = 0;
      break;
    }
  }

  return(_libraries + _info.sqInfo_numLibraries());
}




sqReadData *
sqStore::sqStore_addEmptyRead(sqLibrary *lib) {

  assert(_info.sqInfo_numReads() < _readsAlloc);
  assert(_mode != sqStore_readOnly);

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.sqInfo_numReads().

  _info.sqInfo_addRead();

  increaseArray(_reads, _info.sqInfo_numReads(), _readsAlloc, _info.sqInfo_numReads()/2);

  //  Initialize the new read.

  _reads[_info.sqInfo_numReads()]            = sqRead();
  _reads[_info.sqInfo_numReads()]._readID    = _info.sqInfo_numReads();
  _reads[_info.sqInfo_numReads()]._libraryID = lib->sqLibrary_libraryID();

  //  With the read set up, set pointers in the readData.  Whatever data is in there can stay.

  sqReadData *readData = new sqReadData;

  readData->_read    = _reads + _info.sqInfo_numReads();
  readData->_library = lib;

  return(readData);
}




void
sqStore::sqStore_setClearRange(uint32 id, uint32 bgn, uint32 end) {
  sqRead  *read = sqStore_getRead(id);

  read->_clearBgn = bgn;
  read->_clearEnd = end;
  read->_tExists  = true;
}


void
sqStore::sqStore_setIgnore(uint32 id) {
  sqRead  *read = sqStore_getRead(id);

  read->_ignore = true;
}
