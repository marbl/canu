
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

#include "gkStore.H"

#include "AS_UTL_fileIO.H"


gkStore *gkStore::_instance      = NULL;
uint32   gkStore::_instanceCount = 0;



void
gkRead::gkRead_loadDataFromStream(gkReadData *readData, FILE *file) {
  char    tag[5];
  uint32  size;

  //  Ideally, we'd do one read to get the whole blob.  Without knowing
  //  the length, we're forced to do two.

  AS_UTL_safeRead(file,  tag,  "gkStore::gkStore_loadDataFromFile::blob", sizeof(int8),   4);
  AS_UTL_safeRead(file, &size, "gkStore::gkStore_loadDataFromFile::size", sizeof(uint32), 1);

  uint8 *blob = new uint8 [8 + size];

  memcpy(blob,    tag,  sizeof(uint8)  * 4);
  memcpy(blob+4, &size, sizeof(uint32) * 1);

  AS_UTL_safeRead(file, blob+8, "gkStore::gkStore_loadDataFromFile::blob", sizeof(char), size);

  readData->gkReadData_loadFromBlob(blob);

  delete [] blob;
}



void
gkRead::gkRead_loadDataFromCore(gkReadData *readData, void *blobs) {
  //fprintf(stderr, "gkRead::gkRead_loadDataFromCore()-- read %lu position %lu\n", _readID, _mPtr);
  readData->gkReadData_loadFromBlob(((uint8 *)blobs) + _mPtr);
}



void
gkRead::gkRead_loadDataFromFile(gkReadData *readData, FILE *file) {
  //fprintf(stderr, "gkRead::gkRead_loadDataFromFile()-- read %lu position %lu\n", _readID, _mPtr);
  AS_UTL_fseek(file, _mPtr, SEEK_SET);
  gkRead_loadDataFromStream(readData, file);
}




gkRead *
gkStore::gkStore_getRead(uint32 id) {

  if (gkStore_readInPartition(id) == false)
    return(NULL);

  if (gkStore_readInPartition(id) == false)
    fprintf(stderr, "getRead()--  access to read %u in partition %u is not allowed when partition %u is loaded.\n",
            id, _readIDtoPartitionID[id], _partitionID), assert(0);

  gkRead *read = _reads + (((_readIDtoPartitionID     != NULL) &&
                            (_readIDtoPartitionID[id] == _partitionID)) ? _readIDtoPartitionIdx[id] : id);

  if (gkStore_getNumCorrectedReads() > 0)     //  If there are corrected or trimmed reads in the store,
    read->_cExists = true;                    //  set the flags so the read can return the appropriate data.

  if (gkStore_getNumTrimmedReads() > 0)
    read->_tExists = true;

  return(read);
}



void
gkStore::gkStore_loadReadData(gkRead *read, gkReadData *readData) {

  readData->_read    = read;
  readData->_library = gkStore_getLibrary(read->gkRead_libraryID());

  if (_blobs)
    read->gkRead_loadDataFromCore(readData, _blobs);

  else if (_blobsFiles)
    read->gkRead_loadDataFromFile(readData, _blobsFiles[omp_get_thread_num()]);

  else
    fprintf(stderr, "No data loaded for read %u: no _blobs or _blobsFiles?\n", read->_readID), assert(0);
}


void
gkStore::gkStore_loadReadData(uint32  readID, gkReadData *readData) {

  gkStore_loadReadData(gkStore_getRead(readID), readData);
}



//  Dump a block of encoded data to disk, then update the gkRead to point to it.
//
void
gkStore::gkStore_stashReadData(gkReadData *data) {

  assert(_blobsWriter != NULL);

  data->gkReadData_encodeBlob();

  data->_read->_mPtr = _blobsWriter->tell();
  data->_read->_pID  = _partitionID;                //  0 if not partitioned

  //fprintf(stderr, "STASH read %u at position " F_U64 " or length " F_U64 "\n", read->gkRead_readID(), read->_mPtr, data->_blobLen);

  _blobsWriter->write(data->_blob, data->_blobLen);
}



//  Load read metadata and data from a stream.
//
void
gkStore::gkStore_loadReadFromStream(FILE *S, gkRead *read, gkReadData *readData) {
  char    tag[5];
  uint32  size;

  //  Mark this as a read.  Needed for tgTig::loadFromStreamOrLayout(), and loading this stuff in
  //  utgcns.

  AS_UTL_safeRead(S, tag, "gkStore::gkStore_loadReadFromStream::tag", sizeof(char), 4);

  if (strncmp(tag, "READ", 4) != 0)
    fprintf(stderr, "Failed to load gkRead, got tag '%c%c%c%c' (0x%02x 0x%02x 0x%02x 0x%02x), expected 'READ'.\n",
            tag[0], tag[1], tag[2], tag[3],
            tag[0], tag[1], tag[2], tag[3]), exit(1);

  //  Load the read metadata

  AS_UTL_safeRead(S, read, "gkStore::gkStore_loadReadFromStream::read", sizeof(gkRead), 1);

  //  Load the read data.

  read->gkRead_loadDataFromStream(readData, S);
}



//  Dump the read metadata and read data to a stream.
//
void
gkStore::gkStore_saveReadToStream(FILE *S, uint32 id) {

  //  Mark this as a read.  Needed for tgTig::loadFromStreamOrLayout(), and loading this stuff in
  //  utgcns.

  fprintf(S, "READ");

  //  Dump the read metadata

  gkRead  *read = gkStore_getRead(id);

  AS_UTL_safeWrite(S, read, "gkStore::gkStore_saveReadToStream::read", sizeof(gkRead), 1);

  //  Figure out where the blob actually is, and make sure that it really is a blob

  uint8  *blob    = (uint8 *)_blobs + read->_mPtr;
  uint32  blobLen = 8 + *((uint32 *)blob + 1);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the blob to the stream

  AS_UTL_safeWrite(S, blob, "gkStore::gkStore_saveReadToStream::blob", sizeof(char), blobLen);
}



void
gkReadData::gkReadData_setName(char *H) {
  uint32  Hlen = strlen(H) + 1;

  resizeArray(_name, 0, _nameAlloc, Hlen, resizeArray_doNothing);

  memcpy(_name, H, sizeof(char) * Hlen);
}



void
gkReadData::gkReadData_setBasesQuals(char  *S,
                                     uint8 *Q) {
  bool        isRaw = ((_library->gkLibrary_readType() == GK_READTYPE_PACBIO_RAW) ||
                       (_library->gkLibrary_readType() == GK_READTYPE_NANOPORE_RAW));

  uint32      Slen  = strlen(S) + 1;

  //  If loading raw reads, and no raw read, save the data there.

  if ((isRaw == true) && (_rseq == NULL)) {
    resizeArray(_rseq, 0, _rseqAlloc, Slen, resizeArray_doNothing);
    resizeArray(_rqlt, 0, _rqltAlloc, Slen, resizeArray_doNothing);

    memcpy(_rseq, S, sizeof(char)  * Slen);
    memcpy(_rqlt, Q, sizeof(uint8) * Slen);
  }

  //  If loading corrected reads, and no corrected read, save the date there.

  else {
    if (_read->_cExists)
      fprintf(stderr, "gkReadData_setBasesQuals()- read %u has existing cseq of length %u, replacing with length %u\n",
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
gkReadData::gkReadData_encodeBlobChunk(char const *tag,
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
gkReadData::gkReadData_encodeBlob(void) {

  _blobLen = 0;

  //  Discover which sequence(s) we're encoding.

  _read->_rseqLen = (_rseq == NULL) ? 0 : strlen(_rseq);
  _read->_cseqLen = (_cseq == NULL) ? 0 : strlen(_cseq);

  //  Compute the preferred encodings.  If either fail, the length is set to zero, and the
  //  non-preferred encoding will be computed.  If this too fails, sequences/qualities will be
  //  stored unencoded.

  uint8   *rseq = NULL, *rqlt = NULL;
  uint8   *cseq = NULL, *cqlt = NULL;

  uint32  rseq2Len =                   gkReadData_encode2bit(rseq, _rseq, _read->_rseqLen);
  uint32  rseq3Len = (rseq2Len == 0) ? gkReadData_encode3bit(rseq, _rseq, _read->_rseqLen) : 0;

  uint32  rqlt4Len =                   gkReadData_encode4bit(rqlt, _rqlt, _read->_rseqLen);
  uint32  rqlt5Len = (rqlt4Len == 0) ? gkReadData_encode5bit(rqlt, _rqlt, _read->_rseqLen) : 0;

  uint32  cseq2Len =                   gkReadData_encode2bit(cseq, _cseq, _read->_cseqLen);
  uint32  cseq3Len = (cseq2Len == 0) ? gkReadData_encode3bit(cseq, _cseq, _read->_cseqLen) : 0;

  uint32  cqlt4Len =                   gkReadData_encode4bit(cqlt, _cqlt, _read->_cseqLen);
  uint32  cqlt5Len = (cqlt4Len == 0) ? gkReadData_encode5bit(cqlt, _cqlt, _read->_cseqLen) : 0;

  uint32  qv       = _library->gkLibrary_defaultQV();

  //  Encode the data into chunks in the blob.

  gkReadData_encodeBlobChunk("BLOB", 0,  NULL);

  gkReadData_encodeBlobChunk("NAME", strlen(_name), _name);

  if      (rseq2Len > 0)
    gkReadData_encodeBlobChunk("2SQR",         rseq2Len, rseq);    //  Two-bit encoded sequence (ACGT only)
  else if (rseq3Len > 0)
    gkReadData_encodeBlobChunk("3SQR",         rseq3Len, rseq);    //  Three-bit encoded sequence (ACGTN)
  else if (_read->_rseqLen > 0)
    gkReadData_encodeBlobChunk("USQR", _read->_rseqLen, _rseq);    //  Unencoded sequence

  if      (rqlt4Len > 0)
    gkReadData_encodeBlobChunk("4QVR",         rqlt4Len, rqlt);    //  Four-bit (0-15) encoded QVs
  else if (rqlt5Len > 0)
    gkReadData_encodeBlobChunk("5QVR",         rqlt5Len, rqlt);    //  Five-bit (0-32) encoded QVs
  else if ((_read->_rseqLen > 0) && (_rqlt[0] < 255))
    gkReadData_encodeBlobChunk("UQVR", _read->_rseqLen, _rqlt);    //  Unencoded quality
  else if (_read->_rseqLen > 0)
    gkReadData_encodeBlobChunk("1QVR",                 4, &qv);    //  Constant QV for every base (if sequence exists)

  if      (cseq2Len > 0)
    gkReadData_encodeBlobChunk("2SQC",         cseq2Len, cseq);    //  Two-bit encoded sequence (ACGT only)
  else if (cseq3Len > 0)
    gkReadData_encodeBlobChunk("3SQC",         cseq3Len, cseq);    //  Three-bit encoded sequence (ACGTN)
  else if (_read->_cseqLen > 0)
    gkReadData_encodeBlobChunk("USQC", _read->_cseqLen, _cseq);    //  Unencoded sequence

  if      (cqlt4Len > 0)
    gkReadData_encodeBlobChunk("4QVC",         cqlt4Len, cqlt);    //  Four-bit (0-15) encoded QVs
  else if (cqlt5Len > 0)
    gkReadData_encodeBlobChunk("5QVC",         cqlt5Len, cqlt);    //  Five-bit (0-32) encoded QVs
  else if ((_read->_cseqLen > 0) && (_cqlt[0] < 255))
    gkReadData_encodeBlobChunk("UQVC", _read->_cseqLen, _cqlt);    //  Unencoded quality
  else if (_read->_cseqLen > 0)
    gkReadData_encodeBlobChunk("1QVC",                 4, &qv);    //  Constant QV for every base (if sequence exists)

  gkReadData_encodeBlobChunk("STOP", 0,  NULL);

  //  Cleanup.

  delete [] rseq;
  delete [] rqlt;

  delete [] cseq;
  delete [] cqlt;
}




//  Lowest level function to load data into a read.
//
void
gkReadData::gkReadData_loadFromBlob(uint8 *blob) {
  char    chunk[5];
  uint32  chunkLen = 0;

  //  Make sure that our blob is actually a blob.

  if ((blob[0] != 'B') && (blob[1] != 'L') && (blob[2] != 'O') && (blob[3] != 'B'))
    fprintf(stderr, "Index error in read " F_U32 " mPtr " F_U64 " pID " F_U64 " expected BLOB, got %02x %02x %02x %02x '%c%c%c%c'\n",
            _read->gkRead_readID(),
            _read->_mPtr, _read->_pID,
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
      gkReadData_decode2bit(blob + 8, chunkLen, _rseq, _read->_rseqLen);
    }
    else if (strncmp(chunk, "3SQR", 4) == 0) {
      gkReadData_decode3bit(blob + 8, chunkLen, _rseq, _read->_rseqLen);
    }
    else if (strncmp(chunk, "USQR", 4) == 0) {
      assert(_read->_rseqLen <= chunkLen);
      assert(_read->_rseqLen <= _rseqAlloc);
      memcpy(_rseq, blob + 8, _read->_rseqLen);
      _rseq[_read->_rseqLen] = 0;
    }

    else if (strncmp(chunk, "4QVR", 4) == 0) {
      gkReadData_decode4bit(blob + 8, chunkLen, _rqlt, _read->_rseqLen);
    }
    else if (strncmp(chunk, "5QVR", 4) == 0) {
      gkReadData_decode5bit(blob + 8, chunkLen, _rqlt, _read->_rseqLen);
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
      gkReadData_decode2bit(blob + 8, chunkLen, _cseq, _read->_cseqLen);
    }
    else if (strncmp(chunk, "3SQC", 4) == 0) {
      gkReadData_decode3bit(blob + 8, chunkLen, _cseq, _read->_cseqLen);
    }
    else if (strncmp(chunk, "USQC", 4) == 0) {
      assert(_read->_cseqLen <= chunkLen);
      assert(_read->_cseqLen <= _cseqAlloc);
      memcpy(_cseq, blob + 8, _read->_cseqLen);
      _cseq[_read->_cseqLen] = 0;
    }

    else if (strncmp(chunk, "4QVC", 4) == 0) {
      gkReadData_decode4bit(blob + 8, chunkLen, _cqlt, _read->_cseqLen);
    }
    else if (strncmp(chunk, "5QVC", 4) == 0) {
      gkReadData_decode5bit(blob + 8, chunkLen, _cqlt, _read->_cseqLen);
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
      fprintf(stderr, "gkRead::gkRead_loadDataFromBlob()--  unknown chunk type %02x %02x %02x %02x '%c%c%c%c' skipped\n",
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



gkLibrary *
gkStore::gkStore_addEmptyLibrary(char const *name) {

  assert(_info.numLibraries <= _librariesAlloc);

  //  Just like with reads below, there is no _libraries[0] element.

  _info.numLibraries++;

  increaseArray(_libraries, _info.numLibraries, _librariesAlloc, 128);

  //  Initialize the new library.

  _libraries[_info.numLibraries] = gkLibrary();
  _libraries[_info.numLibraries]._libraryID = _info.numLibraries;

  //  Bullet proof the library name - so we can make files with this prefix.

  char   *libname    = _libraries[_info.numLibraries]._libraryName;
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

  return(_libraries + _info.numLibraries);
}




gkReadData *
gkStore::gkStore_addEmptyRead(gkLibrary *lib) {

  assert(_info.numReads < _readsAlloc);
  assert(_mode != gkStore_readOnly);

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.numReads.

  _info.numReads++;

  increaseArray(_reads, _info.numReads, _readsAlloc, _info.numReads/2);

  //  Initialize the new read.

  _reads[_info.numReads]            = gkRead();
  _reads[_info.numReads]._readID    = _info.numReads;
  _reads[_info.numReads]._libraryID = lib->gkLibrary_libraryID();

  //  With the read set up, set pointers in the readData.  Whatever data is in there can stay.

  gkReadData *readData = new gkReadData;

  readData->_read    = _reads + _info.numReads;
  readData->_library = lib;

  return(readData);
}




void
gkStore::gkStore_setClearRange(uint32 id, uint32 bgn, uint32 end) {
  gkRead  *read = gkStore_getRead(id);

  read->_clearBgn = bgn;
  read->_clearEnd = end;
}
