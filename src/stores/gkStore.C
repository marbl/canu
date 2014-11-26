
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

  uint64  offset = (partitioned == false) ? (_mPtr * DATA_BLOCK_SIZE) : (_pPtr * DATA_BLOCK_SIZE);

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



gkStore::gkStore(char const *path, bool readOnly, bool createOnly) {
  _librariesFile = NULL;
  _libraries     = (gkLibrary *)_librariesFile->get(0);

  _readsFile     = NULL;
  _reads         = (gkRead *)_readsFile->get(0);

  _blobsFile     = NULL;
  _blobs         = (void *)_blobsFile->get(0);
}



//  Open a partitioned store, only reading is supported.
gkStore::gkStore(char const *path, uint32 partID) {
  _librariesFile = NULL;
  _libraries     = (gkLibrary *)_librariesFile->get(0);

  _readsFile     = NULL;
  _reads         = (gkRead *)_readsFile->get(0);

  _blobsFile     = NULL;
  _blobs         = (void *)_blobsFile->get(0);
}



gkStore::~gkStore() {

  //  Should check that inf on disk is the same as inf in memory, and update if needed.

  delete _librariesFile;   _libraries = NULL;
  delete _readsFile;       _reads     = NULL;
  delete _blobsFile;       _blobs     = NULL;
};




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

  _blobsFile = new memoryMappedFile(path);
  _blobs     = (void *)_blobsFile->get(0);

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

