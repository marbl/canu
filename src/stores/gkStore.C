
#include "gkStore.H"

#include "AS_UTL_fileIO.H"
#include "AS_UTL_alloc.H"



bool
gkRead::gkRead_loadData(gkReadData *readData, void *blobs) {

  readData->_read = this;

  //  The resize will only increase the space.  if the new is less than the max, it returns immediately.

  resizeArrayPair(readData->_seq, readData->_qlt, readData->_seqAlloc, readData->_seqAlloc, (uint32)_seqLen+1, resizeArray_doNothing);

  //  Where, or where!, is the data?

  uint64  offset = _mPtr;

  //  One might be tempted to set the readData blob to point to the blob data in the mmap,
  //  but doing so will cause it to be written out again.

  readData->_blobLen = 0;
  readData->_blobMax = 0;
  readData->_blob    = NULL;

  //  Instead, we'll use someting horribly similar.

  uint8  *blob    = ((uint8 *)blobs) + offset;
  char    chunk[5];

  if ((blob[0] != 'B') && (blob[1] != 'L') && (blob[2] != 'O') && (blob[3] != 'B'))
    fprintf(stderr, "Index error in read "F_U32" %c mPtr "F_U64" pID "F_U64" expected BLOB, got %c%c%c%c\n",
            gkRead_readID(),
            '?', //(_numberOfPartitions == 0) ? 'm' : 'p',
            _mPtr, _pID,
            blob[0], blob[1], blob[2], blob[3]);
  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  uint32  blobLen = *((uint32 *)blob + 1);

  //fprintf(stderr, "BLOB len %u\n", blobLen);

  blob += 8;

  while ((blob[0] != 'S') ||
         (blob[1] != 'T') ||
         (blob[2] != 'O') ||
         (blob[3] != 'P')) {
    chunk[0] = blob[0];
    chunk[1] = blob[1];
    chunk[2] = blob[2];
    chunk[3] = blob[3];
    chunk[4] = 0;

    uint32   chunkLen = *((uint32 *)blob + 1);

    //fprintf(stderr, "%s len %u\n", chunk, chunkLen);

    if      (strncmp(chunk, "TYPE", 4) == 0) {
      readData->_type = *((uint32 *)blob + 2);
    }

    else if (strncmp(chunk, "VERS", 4) == 0) {
    }

    else if (strncmp(chunk, "QSEQ", 4) == 0) {
      //fprintf(stderr, "QSEQ not supported.\n");
    }

    else if (strncmp(chunk, "USEQ", 4) == 0) {
      assert(_seqLen <= chunkLen);
      assert(_seqLen <= readData->_seqAlloc);
      memcpy(readData->_seq, blob + 8, _seqLen);
      readData->_seq[_seqLen] = 0;
    }

    else if (strncmp(chunk, "UQLT", 4) == 0) {
      assert(_seqLen <= chunkLen);
      assert(_seqLen <= readData->_seqAlloc);
      memcpy(readData->_qlt, blob + 8, _seqLen);
      readData->_qlt[_seqLen] = 0;
    }

    else if (strncmp(chunk, "2SEQ", 4) == 0) {
      //fprintf(stderr, "2SEQ not supported.\n");
    }

    else if (strncmp(chunk, "3SEQ", 4) == 0) {
      //fprintf(stderr, "3SEQ not supported.\n");
    }

    else if (strncmp(chunk, "QVAL", 4) == 0) {
      uint32  qval = *((uint32 *)blob + 2);

      for (uint32 ii=0; ii<_seqLen; ii++)
        readData->_qlt[ii] = qval;
    }

    else if (strncmp(chunk, "IDNT", 4) == 0) {
      //fprintf(stderr, "IDNT = %u\n", *((uint32 *)blob + 2));
    }

    else if (strncmp(chunk, "SLEN", 4) == 0) {
      //fprintf(stderr, "SLEN = %u\n", *((uint32 *)blob + 2));
    }

    else {
      fprintf(stderr, "gkRead::gkRead_loadData()--  unknown chunk type '%s' skipped\n", chunk);
    }

    blob += 4 + 4 + chunkLen;
  }

#if 0
  switch (_type) {
  case GKREAD_ENC_UNKNOWN:
    fprintf(stderr, "gkRead::gkRead_loadData()-- encoding not set.\n");
    assert(0);
    break;
  case GKREAD_ENC_TYPE_QLT:
    return (gkRead_decodeSeqQlt(readData, (char *)blobs + offset));
    break;
  case GKREAD_TYPE_PACBIO:
    return (gkRead_decodePacBio(readData, (char *)blobs + offset));
    break;
  case GKREAD_TYPE_MINION:
    return (gkRead_decodeMinION(readData, (char *)blobs + offset));
    break;
  default:
    fprintf(stderr, "gkRead()-- ERROR: no known technology for type 0x%02lx.\n", _encoding);
    return(false);
    exit(1);
  }
#endif

  return(true);
};


bool
gkRead::gkRead_decodeSeqQlt(gkReadData *readData, void *blob) {
  //assert(_encoding == GKREAD_ENC_SEQ_QLT);

  return(true);
};

bool
gkRead::gkRead_decodePacBio(gkReadData *readData, void *blob) {
  //assert(_encoding == GKREAD_ENC_PACBIO);

  return(true);
};

bool
gkRead::gkRead_decodeMinION(gkReadData *readData, void *blob) {
  //assert(_encoding == GKREAD_ENC_MINION);

  return(true);
};



//  Dump a block of encoded data to disk, then update the gkRead to point to it.
void
gkStore::gkStore_stashReadData(gkRead *read, gkReadData *data) {

  assert(_blobsFile != NULL);

  read->_mPtr = AS_UTL_ftell(_blobsFile);
  read->_pID  = _partitionID;                //  0 if not partitioned

  //fprintf(stderr, "STASH read %u at position "F_SIZE_T"\n", read->gkRead_readID(), AS_UTL_ftell(_blobsFile));

  AS_UTL_safeWrite(_blobsFile,
                   data->_blob,
                   "gkStore_stashReadData::blob",
                   sizeof(char),
                   data->_blobLen);
}


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

  while (_blobMax + 8 + len < _blobMax) {
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

  memcpy(_blob + _blobLen,  tag, sizeof(char) * 4);     _blobLen += sizeof(char) * 4;
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


gkReadData *
gkRead::gkRead_encodeSeqQlt(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  uint32  RID  = _readID;    //  Debugging
  uint32  Slen = strlen(S);
  uint32  Qlen = strlen(Q);

  if (Slen < Qlen) {
    fprintf(stderr, "-- WARNING:  read '%s' sequence length %u != quality length %u; quality bases truncated.\n",
            H, Slen, Qlen);
    Q[Slen] = 0;
  }

  if (Slen > Qlen) {
    fprintf(stderr, "-- WARNING:  read '%s' sequence length %u != quality length %u; quality bases padded.\n",
            H, Slen, Qlen);
    for (uint32 ii=Qlen; ii<Slen; ii++)
      Q[ii] = Q[Qlen-1];
  }

  uint32  blobType = GKREAD_TYPE_SEQ_QLT;
  uint32  blobVers = 0x00000001;

  uint8  *qseq = NULL;

  _seqLen    = Slen;

  rd->gkReadData_encodeBlobChunk("BLOB",    0,  NULL);
  rd->gkReadData_encodeBlobChunk("TYPE",    4, &blobType);
  rd->gkReadData_encodeBlobChunk("VERS",    4, &blobVers);
  rd->gkReadData_encodeBlobChunk("IDNT",    4, &RID);       //  Debugging, the read id we're supposed to be for
  rd->gkReadData_encodeBlobChunk("SLEN",    4, &Slen);      //  Debugging, the length of the read we're supposed to be for
  rd->gkReadData_encodeBlobChunk("QSEQ",    0,  qseq);      //  Encoded sequence and quality
  rd->gkReadData_encodeBlobChunk("USEQ", Slen,  S);         //  Unencoded sequence
  rd->gkReadData_encodeBlobChunk("UQLT", Qlen,  Q);         //  Unencoded quality
  rd->gkReadData_encodeBlobChunk("STOP",    0,  NULL);

  return(rd);
}

gkReadData *
gkRead::gkRead_encodeSeqQlt(char *H, char *S, uint32 qv) {
  gkReadData *rd = new gkReadData;

  uint32  Slen = strlen(S);

  uint32  blobType = GKREAD_TYPE_SEQ_QLT;
  uint32  blobVers = 0x00000002;

  uint8  *qseq = NULL;

  _seqLen    = Slen;

  rd->gkReadData_encodeBlobChunk("BLOB",    0,  NULL);
  rd->gkReadData_encodeBlobChunk("TYPE",    4, &blobType);
  rd->gkReadData_encodeBlobChunk("VERS",    4, &blobVers);
  //rd->gkReadData_encodeBlobChunk("2SEQ",    0,  qseq);      //  Two-bit encoded sequence (ACGT only)
  //rd->gkReadData_encodeBlobChunk("3SEQ",    0,  qseq);      //  Three-bit encoded sequence (ACGTN)
  rd->gkReadData_encodeBlobChunk("USEQ", Slen,  S);         //  Unencoded sequence
  rd->gkReadData_encodeBlobChunk("QVAL",    4, &qv);        //  QV for every base
  rd->gkReadData_encodeBlobChunk("STOP",    0,  NULL);

  return(rd);
}

gkReadData *
gkRead::gkRead_encodePacBio(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  return(rd);
}

gkReadData *
gkRead::gkRead_encodeMinION(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  return(rd);
}


////////////////////////////////////////
//
//  gkLibrary is lightweight, except for three functions that need to parse strings
//

void
gkLibrary::gkLibrary_parsePreset(char *p) {

  //  Corrected PacBio reads.
  if (strcasecmp(p, "pacbio-corrected") == 0) {
    return;
  }

  //  Uncorrected PacBio, for correction.
  if (strcasecmp(p, "pacbio-correct") == 0) {
    _doCorrection = true;
    return;
  }

  if (strcasecmp(p, "pacbio-raw") == 0) {
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "pacbio-ccs") == 0) {
    return;
  }



  //  Specializations of raw PacBio.
  if (strcasecmp(p, "pacbio-raw-p4c2") == 0) {
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "pacbio-raw-p5c3") == 0) {
    _checkForSubReads = true;
    return;
  }

  if (strcasecmp(p, "pacbio-raw-p6c4") == 0) {
    _checkForSubReads = true;
    return;
  }


  if (strcasecmp(p, "minion") == 0) {
    return;
  }


  if (strcasecmp(p, "contigs") == 0) {
    return;
  }

  fprintf(stderr, "gkLibrary::gkLibrary_parsePreset()--  ERROR: unknown preset '%s'\n", p);
  exit(1);
}


void
gkLibrary::gkLibrary_setFinalTrim(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = FINALTRIM_NONE;

  else if (strcasecmp(t, "largest") == 0)
    _finalTrim = FINALTRIM_LARGEST_COVERED;

  //else if (strcasecmp(t, "evidence") == 0)
  //  _finalTrim = FINALTRIM_EVIDENCE_BASED;

  else if (strcasecmp(t, "bestedge") == 0)
    _finalTrim = FINALTRIM_BEST_EDGE;

  else
    fprintf(stderr, "gkLibrary::gkLibrary_setFinalTrim()--  ERROR: unknown final trim '%s'\n",
            t), exit(1);
}



////////////////////////////////////////

//  The N valid modes for a 'new gkpStore' call:
//
//  1)  Add new reads/libraries, modify old ones.  gkStore(path, true, true)
//  2)  No addition, but can modify old ones.      gkStore(path, true)
//  3)  No addition, no modification.              gkStore(path);
//
gkStore::gkStore(char const *path, gkStore_mode mode, uint32 partID) {
  char    name[FILENAME_MAX];

  memset(_storePath, 0, sizeof(char) * FILENAME_MAX);
  memset(_storeName, 0, sizeof(char) * FILENAME_MAX);

  strcpy(_storePath, path);
  strcpy(_storeName, path);  //  Broken.

  sprintf(name, "%s/info", _storePath);

  //  If the info file exists, load it.

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *I = fopen(name, "r");
    AS_UTL_safeRead(I, &_info, "gkStore::_info", sizeof(gkStoreInfo), 1);
    fclose(I);
  }

  //  Check sizes are correct.

  assert(_info.gkLibrarySize == sizeof(gkLibrary));
  assert(_info.gkReadSize    == sizeof(gkRead));

  //  Clear ourself, to make valgrind happier.

  _librariesMMap          = NULL;
  _librariesAlloc         = 0;
  _libraries              = NULL;

  _readsMMap              = NULL;
  _readsAlloc             = 0;
  _reads                  = NULL;

  _blobsMMap              = NULL;
  _blobs                  = NULL;
  _blobsFile              = NULL;

  _mode                   = mode;

  _numberOfPartitions     = 0;
  _partitionID            = 0;
  _readIDtoPartitionIdx   = NULL;
  _readIDtoPartitionID    = NULL;
  _readsPerPartition      = NULL;
  //_readsInThisPartition   = NULL;

  //
  //  READ ONLY
  //

  if ((mode == gkStore_readOnly) &&
      (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-only access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, false) == false) {
      fprintf(stderr, "gkStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", _storePath);
      exit(1);
    }

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, NO APPEND (also for building a partitioned store)
  //

  else if ((mode == gkStore_modify) &&
           (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-write access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, false) == false) {
      fprintf(stderr, "gkStore()--  failed to open '%s' for read-write access: store doesn't exist.\n", _storePath);
      exit(1);
    }

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readWrite);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, APPEND, open mmap'd files, but copy them entirely to local memory
  //

  else if ((mode == gkStore_extend) &&
           (partID == UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' for read-write and append access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, true) == false)
      AS_UTL_mkdir(_storePath);

    _librariesAlloc = MAX(64, 2 * _info.numLibraries);
    _libraries      = new gkLibrary [_librariesAlloc];

    sprintf(name, "%s/libraries", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _librariesMMap  = new memoryMappedFile (name, memoryMappedFile_readOnly);

      memcpy(_libraries, _librariesMMap->get(0), sizeof(gkLibrary) * (_info.numLibraries + 1));

      delete _librariesMMap;
      _librariesMMap = NULL;;
    }

    _readsAlloc     = MAX(128, 2 * _info.numReads);
    _reads          = new gkRead [_readsAlloc];

    sprintf(name, "%s/reads", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _readsMMap      = new memoryMappedFile (name, memoryMappedFile_readOnly);

      memcpy(_reads, _readsMMap->get(0), sizeof(gkRead) * (_info.numReads + 1));

      delete _readsMMap;
      _readsMMap = NULL;
    }

    sprintf(name, "%s/blobs", _storePath);

    _blobsMMap     = NULL;
    _blobs         = NULL;

    errno = 0;
    _blobsFile     = fopen(name, "a+");
    if (errno)
      fprintf(stderr, "gkStore()--  Failed to open blobs file '%s' for appending: %s\n",
              name, strerror(errno)), exit(1);
  }

  //
  //  PARTITIONED, no modifications, no appends
  //
  //  BIG QUESTION: do we want to partition the read metadata too, or is it small enough
  //  to load in every job?  For now, we load all the metadata.

  else if ((mode == gkStore_readOnly) &&
           (partID != UINT32_MAX)) {
    //fprintf(stderr, "gkStore()--  opening '%s' partition '%u' for read-only access.\n", _storePath, partID);

    //  For partitioned reads, we need to have a uint32 map of readID to partitionReadID so we can
    //  lookup the metadata in the partitoned _reads data.  This is 4 bytes per read, compared to 24
    //  bytes for the full meta data.  Assuming 100x of 3kb read coverage on human, that's 100
    //  million reads, so 0.400 GB vs 2.4 GB.

    sprintf(name, "%s/partitions/map", _storePath);

    errno = 0;
    FILE *F = fopen(name, "r");
    if (errno)
      fprintf(stderr, "gkStore::gkStore()-- failed to open '%s' for reading: %s\n",
              name, strerror(errno)), exit(1);

    AS_UTL_safeRead(F, &_numberOfPartitions, "gkStore::_numberOfPartitions", sizeof(uint32), 1);

    _partitionID            = partID;
    _readsPerPartition      = new uint32 [_numberOfPartitions   + 1];  //  No zeroth element in any of these
    _readIDtoPartitionID    = new uint32 [gkStore_getNumReads() + 1];
    _readIDtoPartitionIdx   = new uint32 [gkStore_getNumReads() + 1];

    AS_UTL_safeRead(F, _readsPerPartition,    "gkStore::_readsPerPartition",    sizeof(uint32), _numberOfPartitions   + 1);
    AS_UTL_safeRead(F, _readIDtoPartitionID,  "gkStore::_readIDtoPartitionID",  sizeof(uint32), gkStore_getNumReads() + 1);
    AS_UTL_safeRead(F, _readIDtoPartitionIdx, "gkStore::_readIDtoPartitionIdx", sizeof(uint32), gkStore_getNumReads() + 1);

    fclose(F);

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _libraries);

    sprintf(name, "%s/partitions/reads.%04"F_U32P"", _storePath, partID);
    _readsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _reads         = (gkRead *)_readsMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _reads);

    sprintf(name, "%s/partitions/blobs.%04"F_U32P"", _storePath, partID);
    _blobsMMap     = new memoryMappedFile (name, memoryMappedFile_readOnly);
    _blobs         = (void *)_blobsMMap->get(0);
    //fprintf(stderr, " -- openend '%s' at "F_X64"\n", name, _blobs);
  }

  //  Info only, no access to reads or libraries.

  else if (mode == gkStore_infoOnly) {
    //fprintf(stderr, "gkStore()--  opening '%s' for info-only access.\n", _storePath);
  }

  else {
    fprintf(stderr, "gkStore::gkStore()-- invalid mode '%s' with partition ID %u.\n",
            toString(mode), partID);
    assert(0);
  }
}






gkStore::~gkStore() {
  char   N[FILENAME_MAX];
  FILE  *F;

  //  Should check that inf on disk is the same as inf in memory, and update if needed.

  bool   needsInfoUpdate = false;

  //  Write N+1 because we write, but don't count, the [0] element.
  
  if        (_librariesMMap) {
    delete _librariesMMap;

  } else if (_libraries) {
    sprintf(N, "%s/libraries", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, _libraries, "libraries", sizeof(gkLibrary), gkStore_getNumLibraries() + 1);
    fclose(F);

    delete [] _libraries;

    needsInfoUpdate = true;
  }


  if        (_readsMMap) {
    delete _readsMMap;

  } else if (_reads) {
    sprintf(N, "%s/reads", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, _reads, "reads", sizeof(gkRead), gkStore_getNumReads() + 1);
    fclose(F);

    delete [] _reads;

    needsInfoUpdate = true;
  }


  if (needsInfoUpdate) {
    sprintf(N, "%s/info", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    AS_UTL_safeWrite(F, &_info, "info", sizeof(gkStoreInfo), 1);
    fclose(F);


    sprintf(N, "%s/info.txt", gkStore_path());
    errno = 0;
    F = fopen(N, "w");
    if (errno)
      fprintf(stderr, "gkStore::~gkStore()-- failed to open '%s' for writing: %s\n",
              N, strerror(errno)), exit(1);

    _info.writeInfoAsText(F);

    fclose(F);
  }


  if (_blobsMMap)
    delete _blobsMMap;

  if (_blobsFile)
    fclose(_blobsFile);

  delete [] _readIDtoPartitionIdx;
  delete [] _readIDtoPartitionID;
  delete [] _readsPerPartition;
};


gkLibrary *
gkStore::gkStore_addEmptyLibrary(char const *name) {

  assert(_librariesMMap == NULL);
  assert(_info.numLibraries <= _librariesAlloc);

  //  Library[0] doesn't exist, see comments in addEmptyRead below.

  _info.numLibraries++;

  if (_info.numLibraries == _librariesAlloc)
    increaseArray(_libraries, _info.numLibraries, _librariesAlloc, 128);

  //  Bullet proof the library name - so we can make files with this prefix.

  char    libname[LIBRARY_NAME_SIZE];
  uint32  libnamepos = 0;
  bool    modified   = false;
  bool    truncated  = false;

  for (char const *orig=name; *orig; orig++) {
    if        (*orig == '/') {
      libname[libnamepos++] = '_';
      libname[libnamepos++] = '-';
      libname[libnamepos++] = '_';
      modified = true;

    } else if (isspace(*orig) == 0) {
      libname[libnamepos++] = *orig;

    } else {
      libname[libnamepos++] = '_';
      modified = true;
    }

    if (libnamepos >= LIBRARY_NAME_SIZE) {
      libname[LIBRARY_NAME_SIZE-1] = 0;
      truncated = true;
      break;
    }
  }

  libname[libnamepos] = 0;

#if 0
  if (modified || truncated)
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s' (original name '%s')\n",
            libname, name);
  else
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s'\n",
            libname);
#endif

  _libraries[_info.numLibraries] = gkLibrary();

  strcpy(_libraries[_info.numLibraries]._libraryName, libname);

  _libraries[_info.numLibraries]._libraryID = _info.numLibraries;

  return(_libraries + _info.numLibraries);
}


gkRead *
gkStore::gkStore_addEmptyRead(gkLibrary *lib) {

  assert(_readsMMap == NULL);
  assert(_info.numReads <= _readsAlloc);
  assert(_mode != gkStore_readOnly);

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.numReads.

  _info.numReads++;

  if (_info.numReads == _readsAlloc)
    increaseArray(_reads, _info.numReads, _readsAlloc, _info.numReads/2);

  _reads[_info.numReads]            = gkRead();
  _reads[_info.numReads]._readID    = _info.numReads;
  _reads[_info.numReads]._libraryID = lib->gkLibrary_libraryID();

  //fprintf(stderr, "ADD READ %u = %u alloc = %u\n", _info.numReads, _reads[_info.numReads]._readID, _readsAlloc);

  return(_reads + _info.numReads);
}





void
gkRead::gkRead_copyDataToPartition(void *blobs, FILE **partfiles, uint64 *partfileslen, uint32 partID) {

  //  Stash away the location of the partitioned data

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

  //  Figure out where the blob actually is, and make sure that it really is a blob

  uint8  *blob    = (uint8 *)blobs + _mPtr * BLOB_BLOCK_SIZE;
  uint32  blobLen = 8 + *((uint32 *)blob + 1);

#if 0
  fprintf(stderr, "blobs "F_X64"\n", blobs);
  fprintf(stderr, "blobs "F_X64"\n", blobs + 1);
  fprintf(stderr, "blobs "F_X64"\n", (uint32 *)blobs);
  fprintf(stderr, "blobs "F_X64"\n", (uint32 *)blobs + 1);

  fprintf(stderr, "blob "F_X64"\n", blob);
  fprintf(stderr, "blob "F_X64"\n", (uint32 *)blob + 1);

  fprintf(stderr, "copy to partition blob="F_X64" len="F_U32" mPtr="F_U64"\n", blob, blobLen, _mPtr);
#endif

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the blob to the partition, update the length of the partition

  AS_UTL_safeWrite(partfiles[partID], blob, "gkRead::gkRead_copyDataToPartition::blob", sizeof(char), blobLen);

  //  Update the read to the new location of the blob in the partitioned data.

  _mPtr = partfileslen[partID];
  _pID  = partID;

  //  And finalize by remembering the length.

  partfileslen[partID] += blobLen;

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

}


void
gkStore::gkStore_buildPartitions(uint32 *partitionMap) {
  char              name[FILENAME_MAX];

  //  Store cannot be partitioned already, and it must be readOnly (for safety) as we don't need to
  //  be changing any of the normal store data.

  assert(_numberOfPartitions == 0);
  assert(_mode               == gkStore_readOnly);

  //  Figure out what the last partition is

  uint32  maxPartition = 0;
  uint32  unPartitioned = 0;

  assert(partitionMap[0] == UINT32_MAX);

  for (uint32 fi=1; fi<=gkStore_getNumReads(); fi++) {
    if (partitionMap[fi] == UINT32_MAX)
      unPartitioned++;

    else if (maxPartition < partitionMap[fi])
      maxPartition = partitionMap[fi];
  }

  fprintf(stderr, "Found "F_U32" unpartitioned reads and maximum partition of "F_U32"\n",
          unPartitioned, maxPartition);

  //  Create the partitions by opening N copies of the data stores,
  //  and writing data to each.

  FILE         **blobfiles    = new FILE * [maxPartition + 1];
  uint64        *blobfileslen = new uint64 [maxPartition + 1];            //  Offset, in bytes, into the blobs file
  FILE         **readfiles    = new FILE * [maxPartition + 1];
  uint32        *readfileslen = new uint32 [maxPartition + 1];            //  aka _readsPerPartition
  uint32        *readIDmap    = new uint32 [gkStore_getNumReads() + 1];   //  aka _readIDtoPartitionIdx

  //  Be nice and put all the partitions in a subdirectory.

  sprintf(name,"%s/partitions", _storePath);

  if (AS_UTL_fileExists(name, true, true) == false)
    AS_UTL_mkdir(name);

  //  Open all the output files -- fail early if we can't open that many files.

  blobfiles[0]    = NULL;
  blobfileslen[0] = UINT64_MAX;
  readfiles[0]    = NULL;
  readfileslen[0] = UINT32_MAX;

  for (uint32 i=1; i<=maxPartition; i++) {
    sprintf(name,"%s/partitions/blobs.%04d", _storePath, i);

    errno = 0;
    blobfiles[i]    = fopen(name, "w");
    blobfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);

    sprintf(name,"%s/partitions/reads.%04d", _storePath, i);

    errno = 0;
    readfiles[i]    = fopen(name, "w");
    readfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);
  }

  //  Open the output partition map file -- we might as well fail early if we can't make it also.

  sprintf(name,"%s/partitions/map", _storePath);

  errno = 0;
  FILE *rIDmF = fopen(name, "w");
  if (errno)
    fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition map file '%s': %s\n",
            name, strerror(errno)), exit(1);

  //  Copy the blob from the master file to the partitioned file, update pointers.

  readIDmap[0] = UINT32_MAX;    //  There isn't a zeroth read, make it bogus.

  for (uint32 fi=1; fi<=gkStore_getNumReads(); fi++) {
    uint32  pi = partitionMap[fi];

    assert(pi != 0);  //  No zeroth partition, right?

    if (pi == UINT32_MAX)
      //  Deleted reads are not assigned a partition; skip them
      continue;

    //  Make a copy of the read, then modify it for the partition, then write it to the partition.
    //  Without the copy, we'd need to update the master record too.

    gkRead  partRead = _reads[fi];  //*gkStore_getRead(fi);

    partRead.gkRead_copyDataToPartition(_blobs, blobfiles, blobfileslen, pi);

#if 1
    fprintf(stderr, "read "F_U32"="F_U32" len "F_U32" -- blob master "F_U64" -- to part "F_U32" new read id "F_U32" blob "F_U64"/"F_U64" -- at readIdx "F_U32"\n",
            fi, _reads[fi].gkRead_readID(), _reads[fi].gkRead_sequenceLength(),
            _reads[fi]._mPtr,
            pi,
            partRead.gkRead_readID(), partRead._pID, partRead._mPtr,
            readfileslen[pi]);
#endif

    AS_UTL_safeWrite(readfiles[pi], &partRead, "gkStore::gkStore_buildPartitions::read", sizeof(gkRead), 1);

    readIDmap[fi] = readfileslen[pi]++;
  }

  //  There isn't a zeroth read.

  AS_UTL_safeWrite(rIDmF, &maxPartition,  "gkStore::gkStore_buildPartitions::maxPartition", sizeof(uint32), 1);
  AS_UTL_safeWrite(rIDmF,  readfileslen,  "gkStore::gkStore_buildPartitions::readfileslen", sizeof(uint32), maxPartition + 1);
  AS_UTL_safeWrite(rIDmF,  partitionMap,  "gkStore::gkStore_buildPartitions::partitionMap", sizeof(uint32), gkStore_getNumReads() + 1);
  AS_UTL_safeWrite(rIDmF,  readIDmap,     "gkStore::gkStore_buildPartitions::readIDmap",    sizeof(uint32), gkStore_getNumReads() + 1);

  //  cleanup -- close all the files, delete storage

  fclose(rIDmF);

  for (uint32 i=1; i<=maxPartition; i++) {
    fprintf(stderr, "partition "F_U32" has "F_U32" reads\n", i, readfileslen[i]);

    errno = 0;

    fclose(blobfiles[i]);
    fclose(readfiles[i]);

    if (errno)
      fprintf(stderr, "  warning: %s\n", strerror(errno));
  }

  delete [] readIDmap;
  delete [] readfileslen;
  delete [] readfiles;
  delete [] blobfileslen;
  delete [] blobfiles;
}



void
gkStore::gkStore_delete(void) {
  char path[FILENAME_MAX];

  delete [] _libraries;
  delete [] _reads;

  _libraries = NULL;
  _reads     = NULL;

  gkStore_deletePartitions();

  sprintf(path, "%s/info",      gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/libraries", gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/reads",     gkStore_path());  AS_UTL_unlink(path);
  sprintf(path, "%s/blobs",     gkStore_path());  AS_UTL_unlink(path);

  AS_UTL_unlink(path);
}



void
gkStore::gkStore_deletePartitions(void) {
  char path[FILENAME_MAX];

  sprintf(path, "%s/partitions/map", gkStore_path());

  if (AS_UTL_fileExists(path, false, false) == false)
    return;

  //  How many partitions?

  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "ERROR: failed to open partition meta data '%s': %s\n", path, strerror(errno)), exit(1);

  fread(&_numberOfPartitions, sizeof(uint32), 1, F);

  fclose(F);

  //  Yay!  Delete!

  AS_UTL_unlink(path);

  for (uint32 ii=0; ii<_numberOfPartitions; ii++) {
    sprintf(path, "%s/partitions/reads.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
    sprintf(path, "%s/partitions/blobs.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
  }
}






void
gkStoreStats::init(gkStore *gkp) {

#if 0
  gkFragment    fr;
  gkStream     *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  numActiveFrag     = 0;
  numMatedFrag      = 0;
  readLength        = 0;
  clearLength       = 0;

  lowestID          = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  highestID         = new uint32 [gkp->gkStore_getNumLibraries() + 1];

  numActivePerLib   = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numMatedPerLib    = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  readLengthPerLib  = new uint64 [gkp->gkStore_getNumLibraries() + 1];
  clearLengthPerLib = new uint64 [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkp->gkStore_getNumLibraries() + 1; i++) {
    lowestID[i]          = 0;
    highestID[i]         = 0;

    numActivePerLib[i]   = 0;
    numMatedPerLib[i]    = 0;
    readLengthPerLib[i]  = 0;
    clearLengthPerLib[i] = 0;
  }

  while (fs->next(&fr)) {
    uint32     lid = fr.gkFragment_getLibraryID();
    uint32     rid = fr.gkFragment_getReadID();

    if (lowestID[lid] == 0) {
      lowestID[lid]  = rid;
      highestID[lid] = rid;
    }
    if (highestID[lid] < rid) {
      highestID[lid] = rid;
    }

    numActiveFrag++;
    numActivePerLib[lid]++;

    readLength             += fr.gkFragment_getSequenceLength();
    readLengthPerLib[lid]  += fr.gkFragment_getSequenceLength();

    clearLength            += fr.gkFragment_getClearRegionLength();
    clearLengthPerLib[lid] += fr.gkFragment_getClearRegionLength();
  }

  delete fs;
#endif

}

