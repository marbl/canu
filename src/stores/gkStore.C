
#include "gkStore.H"

#include "AS_UTL_fileIO.H"


bool
gkStore::gkStore_loadReadData(uint32  readID, gkReadData *readData) {

  return(_reads[readID].gkRead_loadData(readData,
                                        _blobs,
                                        (_numberOfPartitions == 0) ? false : true));
};


bool
gkRead::gkRead_loadData(gkReadData *readData, void *blobs, bool partitioned) {

  readData->_read = this;

  if (readData->_seq == NULL)
    readData->_seq  = new char [_seqLen + 1];

  if (readData->_qlt == NULL)
    readData->_qlt  = new char [_seqLen + 1];

  //  Where, or where!, is the data?

  uint64  offset = (partitioned == false) ? (_mPtr * BLOB_BLOCK_SIZE) : (_pPtr * BLOB_BLOCK_SIZE);

  //  One might be tempted to set the readData blob to point to the blob data in the mmap,
  //  but doing so will cause it to be written out again.

  readData->_blobLen = 0;
  readData->_blobMax = 0;
  readData->_blob    = NULL;

  //  Instead, we'll use someting horribly similar.

  uint8  *blob    = (uint8 *)blobs + offset;
  char    chunk[5];

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  uint32  blobLen = *((uint32 *)blob + 1);

  blob += 8;

  while ((blob[0] != 'S') &&
         (blob[1] != 'T') &&
         (blob[2] != 'O') &&
         (blob[3] != 'P')) {
    chunk[0] = blob[0];
    chunk[1] = blob[1];
    chunk[2] = blob[2];
    chunk[3] = blob[3];
    chunk[4] = 0;

    uint32   chunkLen = *((uint32 *)blob + 1);

    if      (strncmp(chunk, "TYPE", 4) == 0) {
      readData->_type = *((uint32 *)blob + 2);
    }

    else if (strncmp(chunk, "VERS", 4) == 0) {
    }

    else if (strncmp(chunk, "QSEQ", 4) == 0) {
      fprintf(stderr, "QSEQ not supported.\n");
    }

    else if (strncmp(chunk, "USEQ", 4) == 0) {
      assert(_seqLen <= chunkLen);
      memcpy(readData->_seq, blob + 8, _seqLen);
      readData->_seq[_seqLen] = 0;
      //fprintf(stderr, "READ USEQ chunklen "F_U32" len "F_U64" %c%c%c\n", chunkLen, _seqLen, chunk[8], chunk[9], chunk[10]);
    }

    else if (strncmp(chunk, "UQLT", 4) == 0) {
      assert(_seqLen <= chunkLen);
      memcpy(readData->_qlt, blob + 8, _seqLen);
      readData->_qlt[_seqLen] = 0;
    }

    else if (strncmp(chunk, "2SEQ", 4) == 0) {
      fprintf(stderr, "2SEQ not supported.\n");
    }

    else if (strncmp(chunk, "3SEQ", 4) == 0) {
      fprintf(stderr, "3SEQ not supported.\n");
    }

    else if (strncmp(chunk, "QVAL", 4) == 0) {
      uint32  qval = *((uint32 *)blob + 2);

      for (uint32 ii=0; ii<_seqLen; ii++)
        readData->_qlt[ii] = qval;
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

  fprintf(stderr, "STASH read %u at position "F_SIZE_T"\n", read->gkRead_readID(), AS_UTL_ftell(_blobsFile));

  if (_numberOfPartitions == 0) {
    read->_mPtr = AS_UTL_ftell(_blobsFile);

  } else {

    read->_pID  = _partitionID;
    read->_pPtr = AS_UTL_ftell(_blobsFile);
  }

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

  //  Finally, update the blob length.

  _blobLen -= 8;

  memcpy(_blob + 4, &_blobLen, sizeof(uint32));

  _blobLen += 8;
}


gkReadData *
gkRead::gkRead_encodeSeqQlt(char *H, char *S, char *Q) {
  gkReadData *rd = new gkReadData;

  uint32  Slen = strlen(S);
  uint32  Qlen = strlen(Q);

  assert(Slen == Qlen);

  uint32  blobType = GKREAD_TYPE_SEQ_QLT;
  uint32  blobVers = 0x00000001;

  uint8  *qseq = NULL;

  _seqLen    = Slen;

  _clrBgn    = 0;
  _clrEnd    = Slen;

  rd->gkReadData_encodeBlobChunk("BLOB",    0,  NULL);
  rd->gkReadData_encodeBlobChunk("TYPE",    4, &blobType);
  rd->gkReadData_encodeBlobChunk("VERS",    4, &blobVers);
  //rd->gkReadData_encodeBlobChunk("QSEQ",    0,  qseq);      //  Encoded sequence and quality
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

  _clrBgn    = 0;
  _clrEnd    = Slen;

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
gkLibrary::gkLibrary_setInitialTrim(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _initialTrim = INITIALTRIM_NONE;

  else if (strcasecmp(t, "mer") == 0)
    _initialTrim = INITIALTRIM_MER_BASED;

  else if (strcasecmp(t, "flow") == 0)
    _initialTrim = INITIALTRIM_FLOW_BASED;

  else if (strcasecmp(t, "quality") == 0)
    _initialTrim = INITIALTRIM_QUALITY_BASED;

  else
    fprintf(stderr, "gkLibrary::gkLibrary_setInitialTrim()--  ERROR: unknown initial trim '%s'\n",
            t), exit(1);
}


void
gkLibrary::gkLibrary_setFinalTrim(char *t) {

  if      (strcasecmp(t, "none") == 0)
    _finalTrim = FINALTRIM_NONE;

  else if (strcasecmp(t, "largest") == 0)
    _finalTrim = FINALTRIM_LARGEST_COVERED;

  else if (strcasecmp(t, "evidence") == 0)
    _finalTrim = FINALTRIM_EVIDENCE_BASED;

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

  //  If the info file exists, load it, otherwise, initialize a new empty store.

  if (AS_UTL_fileExists(name, false, false) == true) {
    errno = 0;
    FILE *I = fopen(name, "r");
    AS_UTL_safeRead(I, &_info, "gkStore::_info", sizeof(gkStoreInfo), 1);
    fclose(I);
  } else {
    //  already initialized.
    assert(_info.gkLibrarySize == sizeof(gkLibrary));
    assert(_info.gkReadSize    == sizeof(gkRead));
  }

  //  Clear ourself, to make valgrind happier.

  _librariesMMap   = NULL;
  _librariesAlloc  = 0;
  _libraries       = NULL;

  _readsMMap       = NULL;
  _readsAlloc      = 0;
  _reads           = NULL;

  _blobsMMap       = NULL;
  _blobs           = NULL;
  _blobsFile       = NULL;

  //
  //  READ ONLY
  //

  if (mode == gkStore_readOnly) {
    fprintf(stderr, "gkStore()--  opening '%s' for read-only access.\n", _storePath);

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, false);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, false);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, false);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, NO APPEND
  //

  else if (mode == gkStore_modify) {
    fprintf(stderr, "gkStore()--  opening '%s' for read-write access.\n", _storePath);

    sprintf(name, "%s/libraries", _storePath);
    _librariesMMap = new memoryMappedFile (name, true);
    _libraries     = (gkLibrary *)_librariesMMap->get(0);

    sprintf(name, "%s/reads", _storePath);
    _readsMMap     = new memoryMappedFile (name, true);
    _reads         = (gkRead *)_readsMMap->get(0);

    sprintf(name, "%s/blobs", _storePath);
    _blobsMMap     = new memoryMappedFile (name, true);
    _blobs         = (void *)_blobsMMap->get(0);
  }

  //
  //  MODIFY, APPEND
  //

  else if (mode == gkStore_extend) {
    fprintf(stderr, "gkStore()--  opening '%s' for read-write and append access.\n", _storePath);

    if (AS_UTL_fileExists(_storePath, true, true) == false)
      AS_UTL_mkdir(_storePath);

    _librariesAlloc = MAX(64, 2 * _info.numLibraries);
    _libraries      = new gkLibrary [_librariesAlloc];

    sprintf(name, "%s/libraries", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _librariesMMap  = new memoryMappedFile (name, true);

      memcpy(_libraries, _librariesMMap->get(0), sizeof(gkLibrary) * (_info.numLibraries + 1));

      delete _librariesMMap;
      _librariesMMap = NULL;;
    }


    _readsAlloc     = MAX(1048576, 2 * _info.numReads);
    _reads          = new gkRead [_readsAlloc];

    sprintf(name, "%s/reads", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _readsMMap      = new memoryMappedFile (name, true);

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

  else if (mode == gkStore_partitioned) {
  }

  else {
    fprintf(stderr, "gkStore::gkStore()-- unknown mode 0x%02x.\n", mode);
    exit(1);
  }

  //  Continue clearing ourself.

  _numberOfPartitions     = 0;
  _partitionID            = 0;
  _partitionIDmap         = NULL;
  _readsPerPartition      = NULL;
  _readsInThisPartition   = NULL;
}






gkStore::~gkStore() {
  char   N[FILENAME_MAX];
  FILE  *F;

  //  Should check that inf on disk is the same as inf in memory, and update if needed.

  bool   needsInfoUpdate = false;

  //  Write N+1 because we write, but don't count, the [0] element.
  
  if (_librariesMMap) {
    delete _librariesMMap;

  } else {
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


  if (_readsMMap) {
    delete _readsMMap;

  } else {
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
  }


  if (_blobsMMap)
    delete _blobsMMap;

  if (_blobsFile)
    fclose(_blobsFile);

  delete [] _partitionIDmap;
  delete [] _readsPerPartition;
  delete [] _readsInThisPartition;
};


gkLibrary *
gkStore::gkStore_addEmptyLibrary(char const *name) {

  assert(_librariesMMap == NULL);
  assert(_info.numLibraries <= _librariesAlloc);

  if (_librariesAlloc == _info.numLibraries) {
    fprintf(stderr, "NEED TO REALLOC LIBRARIES.\n");
    assert(0);
  }

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

  if (modified || truncated)
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s' (original name '%s')\n",
            libname, name);
  else
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s'\n",
            libname);

  //  Library[0] doesn't exist, see comments in addEmptyRead below.
  _info.numLibraries++;

  _libraries[_info.numLibraries] = gkLibrary();

  strcpy(_libraries[_info.numLibraries]._libraryName, libname);

  _libraries[_info.numLibraries]._libraryID = _info.numLibraries;

  return(_libraries + _info.numLibraries);
}


gkRead *
gkStore::gkStore_addEmptyRead(gkLibrary *lib) {

  assert(_readsMMap == NULL);
  assert(_info.numReads <= _readsAlloc);

  if (_readsAlloc == _info.numReads) {
    fprintf(stderr, "NEED TO REALLOC READS.\n");
    assert(0);
  }

  //  We reserve the zeroth read for "null".  This is easy to accomplish
  //  here, just pre-increment the number of reads.  However, we need to be sure
  //  to iterate up to and including _info.numReads.

  _info.numReads++;

  _reads[_info.numReads] = gkRead();
  _reads[_info.numReads]._readID    = _info.numReads;
  _reads[_info.numReads]._libraryID = lib->gkLibrary_libraryID();

  return(_reads + _info.numReads);
}



bool
gkStore::gkStore_loadPartition(uint32 partID) {
  char     path[FILENAME_MAX];
  uint32   nReads = 0;

  sprintf(path, "%s/part", gkStore_path());

  if (AS_UTL_fileExists(path, false, false) == false)
    //  Not partitioned.
    return(false);

  //  Count the number of reads in this partition, and generate a list of them (later).

  for (uint32 ii=1; ii<gkStore_getNumReads() + 1; ii++)
    if (_reads[ii]._pID == partID)
      nReads++;

  //  Open the partition data

  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "ERROR: failed to open partition meta data '%s': %s\n", path, strerror(errno)), exit(1);

  //  Load it.

  fread(&_numberOfPartitions, sizeof(uint32), 1, F);

  _partitionID          = partID;
  _partitionIDmap       = new uint32 [gkStore_getNumReads() + 1];
  _readsPerPartition    = new uint32 [_numberOfPartitions];
  _readsInThisPartition = new uint32 [nReads];

  fread(_partitionIDmap, sizeof(uint32), gkStore_getNumReads() + 1, F);

  fclose(F);

  //  Now, generate a list of the reads in this partition (not dependent on the loaded
  //  data, just on the allocation done with the rest of them).

  nReads = 0;

  for (uint32 ii=1; ii<gkStore_getNumReads() + 1; ii++)
    if (_reads[ii]._pID == _partitionID)
      _readsInThisPartition[nReads++] = ii;

  //  We probably want to load the sequence data too, or at least mmap it.

  sprintf(path, "%s/%s.%04u", gkStore_path(), "seqs", _partitionID);

  _blobsMMap = new memoryMappedFile(path);
  _blobs     = (void *)_blobsMMap->get(0);

  return(true);
}


void
gkStore::gkStore_buildPartitions(uint32 *partitionMap, uint32 maxPartition) {
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

  //  Probably need to load the 'part' info.

  sprintf(path, "%s/partitions", gkStore_path());

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

  for (uint32 ii=0; ii<_numberOfPartitions; ii++)
    sprintf(path, "%s/blobs.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
}






void
gkStoreStats::init(gkStore *gkp) {

#if 0
  gkFragment    fr;
  gkStream     *fs = new gkStream(gkp, 0, 0, GKFRAGMENT_INF);

  numActiveFrag     = 0;
  numDeletedFrag    = 0;
  numMatedFrag      = 0;
  readLength        = 0;
  clearLength       = 0;

  lowestID          = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  highestID         = new uint32 [gkp->gkStore_getNumLibraries() + 1];

  numActivePerLib   = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numDeletedPerLib  = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  numMatedPerLib    = new uint32 [gkp->gkStore_getNumLibraries() + 1];
  readLengthPerLib  = new uint64 [gkp->gkStore_getNumLibraries() + 1];
  clearLengthPerLib = new uint64 [gkp->gkStore_getNumLibraries() + 1];

  for (uint32 i=0; i<gkp->gkStore_getNumLibraries() + 1; i++) {
    lowestID[i]          = 0;
    highestID[i]         = 0;

    numActivePerLib[i]   = 0;
    numDeletedPerLib[i]  = 0;
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

    if (fr.gkFragment_getIsDeleted()) {
      numDeletedFrag++;
      numDeletedPerLib[lid]++;
    } else {
      numActiveFrag++;
      numActivePerLib[lid]++;

      //if (fr.gkFragment_getMateID() > 0) {
      //  numMatedFrag++;
      //  numMatedPerLib[lid]++;
      //}

      readLength             += fr.gkFragment_getSequenceLength();
      readLengthPerLib[lid]  += fr.gkFragment_getSequenceLength();

      clearLength            += fr.gkFragment_getClearRegionLength();
      clearLengthPerLib[lid] += fr.gkFragment_getClearRegionLength();
    }
  }

  delete fs;
#endif

}

