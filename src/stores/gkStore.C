
#include "gkStore.H"

#include "AS_UTL_fileIO.H"


bool
gkStore::gkStore_loadData(uint32  readID, gkReadData *readData) {

  return(_reads[readID].gkRead_loadData(readData,
                                        _blobs,
                                        (_numberOfPartitions == 0) ? false : true));
};


bool
gkRead::gkRead_loadData(gkReadData *readData, void *blobs, bool partitioned) {

  readData->_read = this;

  readData->_seq  = new char [_seqLen];
  readData->_qlt  = new char [_seqLen];

  uint64  offset = (partitioned == false) ? (_mPtr * BLOB_BLOCK_SIZE) : (_pPtr * BLOB_BLOCK_SIZE);

  switch (_technology) {
  case GKREAD_TECH_FASTA:
    return (gkRead_loadFASTAData(readData, (char *)blobs + offset));
    break;
  case GKREAD_TECH_FASTQ:
    return (gkRead_loadFASTQData(readData, (char *)blobs + offset));
    break;
  case GKREAD_TECH_PACBIO:
    return (gkRead_loadPacBioData(readData, (char *)blobs + offset));
    break;
  case GKREAD_TECH_NANOPORE:
    return (gkRead_loadNanoporeData(readData, (char *)blobs + offset));
    break;
  default:
    fprintf(stderr, "gkRead()-- ERROR: no known technology for type 0x%02lx.\n", _technology);
    return(false);
    exit(1);
  }

  return(false);
};


bool
gkRead::gkRead_loadFASTAData(gkReadData *readData, void *blob) {
  assert(_technology == GKREAD_TECH_FASTA);

  return(true);
};

bool
gkRead::gkRead_loadFASTQData(gkReadData *readData, void *blob) {
  assert(_technology == GKREAD_TECH_FASTQ);

  return(true);
};

bool
gkRead::gkRead_loadPacBioData(gkReadData *readData, void *blob) {
  assert(_technology == GKREAD_TECH_PACBIO);

  return(true);
};

bool
gkRead::gkRead_loadNanoporeData(gkReadData *readData, void *blob) {
  assert(_technology == GKREAD_TECH_NANOPORE);

  return(true);
};




////////////////////////////////////////

//  The N valid modes for a 'new gkpStore' call:
//
//  1)  Add new reads/libraries, modify old ones.  gkStore(path, true, true)
//  2)  No addition, but can modify old ones.      gkStore(path, true)
//  3)  No addition, no modification.              gkStore(path);
//
gkStore::gkStore(char const *path, gkStore_mode mode, uint32 partID) {
  char    name[FILENAME_MAX];

  strcpy(_storePath, path);
  strcpy(_storeName, path);  //  Broken.

  sprintf(name, "%s/info", _storePath);

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

      memcpy(_libraries, _librariesMMap->get(0), sizeof(gkLibrary) * _info.numLibraries);

      delete _librariesMMap;
      _librariesMMap = NULL;;
    }


    _readsAlloc     = MAX(1048576, 2 * _info.numLibraries);
    _reads          = new gkRead [_readsAlloc];

    sprintf(name, "%s/reads", _storePath);
    if (AS_UTL_fileExists(name, false, false) == true) {
      _readsMMap      = new memoryMappedFile (name, true);

      memcpy(_reads, _readsMMap->get(0), sizeof(gkRead) * _info.numReads);

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
}






gkStore::~gkStore() {

  //  Should check that inf on disk is the same as inf in memory, and update if needed.

  bool   needsInfoUpdate = false;

  if (_librariesMMap) {
    delete _librariesMMap;

  } else {
    //  dump the new data
    needsInfoUpdate = true;
  }

  if (_blobsFile)
    fclose(_blobsFile);
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

  for (char const *orig=name; *orig; orig++) {
    if        (*orig == '/') {
      libname[libnamepos++] = '_';
      libname[libnamepos++] = '-';
      libname[libnamepos++] = '_';
    } else if (isspace(*orig) == 0) {
      libname[libnamepos++] = *orig;
    } else {
      libname[libnamepos++] = '_';
    }

    if (libnamepos >= LIBRARY_NAME_SIZE) {
      libname[LIBRARY_NAME_SIZE-1] = 0;
      fprintf(stderr, "gkStore_addEmptyLibrary()--  WARNING: library name '%s' truncated to '%s'\n",
              name, libname);
      break;
    }
  }

  libname[libnamepos] = 0;

  if (strcmp(libname, name) != 0)
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s' (original name '%s')\n",
            libname, name);
  else
    fprintf(stderr, "gkStore_addEmptyLibrary()--  added library '%s'\n",
            libname);

  _libraries[_info.numLibraries] = gkLibrary();
  _libraries[_info.numLibraries].gkLibrary_setLibraryName(libname);

  _info.numLibraries++;

  return(_libraries + _info.numLibraries);
}


gkRead *
gkStore::gkStore_addEmptyRead(void) {
  gkRead     empty;

  assert(_readsMMap == NULL);
  assert(_info.numReads <= _readsAlloc);

  if (_readsAlloc == _info.numReads) {
    fprintf(stderr, "NEED TO REALLOC READS.\n");
    assert(0);
  }

  _reads[_info.numReads] = gkRead();

  _info.numReads++;

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

  for (uint32 ii=0; ii<gkStore_getNumReads(); ii++)
    if (_reads[ii]._pID == partID)
      nReads++;

  //  Open the partition data

  FILE *F = fopen(path, "r");
  if (errno)
    fprintf(stderr, "ERROR: failed to open partition meta data '%s': %s\n", path, strerror(errno)), exit(1);

  //  Load it.

  fread(&_numberOfPartitions, sizeof(uint32), 1, F);

  _partitionID          = partID;
  _partitionIDmap       = new uint32 [gkStore_getNumReads()];
  _readsPerPartition    = new uint32 [_numberOfPartitions];
  _readsInThisPartition = new uint32 [nReads];

  fread(_partitionIDmap, sizeof(uint32), gkStore_getNumReads(), F);

  fclose(F);

  //  Now, generate a list of the reads in this partition (not dependent on the loaded
  //  data, just on the allocation done with the rest of them).

  nReads = 0;

  for (uint32 ii=0; ii<gkStore_getNumReads(); ii++)
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

