
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
 *    src/stores/gkStoreConstructor.C
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
#include "ovOverlap.H"

#include "files.H"



void
sqStore::sqStore_loadMetadata(void) {
  char    name[FILENAME_MAX+1];

  _librariesAlloc = _info.sqInfo_numLibraries() + 1;
  _readsAlloc     = _info.sqInfo_numReads()     + 1;

  _libraries      = new sqLibrary [_librariesAlloc];
  _reads          = new sqRead    [_readsAlloc];

  AS_UTL_loadFile(_storePath, '/', "libraries", _libraries, _librariesAlloc);
  AS_UTL_loadFile(_storePath, '/', "reads",     _reads,     _readsAlloc);
}










sqStore::sqStore(char const    *storePath_,
                 char const    *clonePath_,
                 sqStore_mode   mode_,
                 uint32         partID_) {
  char    nameI[FILENAME_MAX+1];
  char    nameL[FILENAME_MAX+1];
  char    nameR[FILENAME_MAX+1];
  char    nameB[FILENAME_MAX+1];

  //  Set the sqStore pointer in an overlap.

  ovOverlap::sqStoreAttach(this);

  //  Clear ourself, to make valgrind happier.

  memset(_storePath, 0, sizeof(char) * (FILENAME_MAX + 1));
  memset(_clonePath, 0, sizeof(char) * (FILENAME_MAX + 1));

  _mode                   = mode_;

  _librariesAlloc         = 0;
  _libraries              = NULL;

  _readsAlloc             = 0;
  _reads                  = NULL;

  _blobsData              = NULL;

  _blobsFilesMax          = 0;
  _blobsFiles             = NULL;

  _blobsWriter            = NULL;

  _numberOfPartitions     = 0;
  _partitionID            = 0;
  _readIDtoPartitionIdx   = NULL;
  _readIDtoPartitionID    = NULL;
  _readsPerPartition      = NULL;

  //  Save the path and name.

  if (storePath_)   strncpy(_storePath, storePath_, FILENAME_MAX);   //  storePath must always exist though.
  if (clonePath_)   strncpy(_clonePath, clonePath_, FILENAME_MAX);   //  clonePath is definitely optional.

  //  If the info file exists, load it.

  snprintf(nameI, FILENAME_MAX, "%s/info", _storePath);

  if (fileExists(nameI) == true)
    AS_UTL_loadFile(nameI, &_info, 1);

  //  Check sizes are correct.

  if (_info.checkInfo() == false) {
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR:  Can't open store '%s': parameters in sqStore.H and sqRead.H are incompatible with the store.\n", _storePath);
    exit(1);
  }

  //
  //  CREATE - allocate some memory for saving libraries and reads, and create a file to dump the data into.
  //

  if (_mode == sqStore_create) {
    if (partID_ != UINT32_MAX)
      fprintf(stderr, "sqStore()-- Illegal combination of sqStore_create with defined partID.\n"), exit(1);

    if (directoryExists(_storePath) == true)
      fprintf(stderr, "ERROR:  Can't create store '%s': store already exists.\n", _storePath), exit(1);

    AS_UTL_mkdir(_storePath);

    _librariesAlloc = 32;           //  _libraries and
    _readsAlloc     = 32768;        //  _reads MUST be preallocated.

    _libraries      = new sqLibrary [_librariesAlloc];
    _reads          = new sqRead    [_readsAlloc];

    _blobsWriter    = new sqStoreBlobWriter(_storePath, 0);

    return;
  }

  //
  //  Not creating, so the store MUST exist.  Check some other conditions too.
  //

  if (directoryExists(_storePath) == false)
    fprintf(stderr, "sqStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", _storePath), exit(1);

  if ((_mode == sqStore_extend) &&
      (partID_ != UINT32_MAX))
    fprintf(stderr, "sqStore()-- Illegal combination of sqStore_extend with defined partID.\n"), exit(1);

  //
  //  EXTEND - just load the metadata, allocate some stuff, and return.
  //

  if (_mode == sqStore_extend) {
    sqStore_loadMetadata();

    _blobsFilesMax = omp_get_max_threads();
    _blobsFiles    = new sqStoreBlobReader [_blobsFilesMax];

    _blobsWriter   = new sqStoreBlobWriter(_storePath, _info.sqInfo_numBlobs());

    return;
  }

  //
  //  BUILDING PARTITIONS - load metadata and return.
  //

  if (_mode == sqStore_buildPart) {
    sqStore_loadMetadata();

    _blobsFilesMax = omp_get_max_threads();
    _blobsFiles    = new sqStoreBlobReader [_blobsFilesMax];

    return;
  }

  //
  //  READ ONLY non-partitioned - just load the metadata and return.
  //

  if (partID_ == UINT32_MAX) {       //  READ ONLY, non-partitioned (also for creating partitions)
    sqStore_loadMetadata();

    _blobsFilesMax = omp_get_max_threads();
    _blobsFiles    = new sqStoreBlobReader [_blobsFilesMax];

    return;
  }

  //
  //  READ ONLY partitioned.  A whole lotta work to do.
  //

  snprintf(nameI, FILENAME_MAX, "%s/partitions/map", _storePath);

  FILE *F = AS_UTL_openInputFile(nameI);

  loadFromFile(_numberOfPartitions, "sqStore::_numberOfPartitions", F);

  _partitionID            = partID_;
  _readsPerPartition      = new uint32 [_numberOfPartitions   + 1];  //  No zeroth element in any of these
  _readIDtoPartitionID    = new uint32 [sqStore_getNumReads() + 1];
  _readIDtoPartitionIdx   = new uint32 [sqStore_getNumReads() + 1];

  loadFromFile(_readsPerPartition,    "sqStore::_readsPerPartition",    _numberOfPartitions   + 1, F);
  loadFromFile(_readIDtoPartitionID,  "sqStore::_readIDtoPartitionID",  sqStore_getNumReads() + 1, F);
  loadFromFile(_readIDtoPartitionIdx, "sqStore::_readIDtoPartitionIdx", sqStore_getNumReads() + 1, F);

  AS_UTL_closeFile(F, nameI);

  //  Load the rest of the data, just suck in entire files.

  snprintf(nameL, FILENAME_MAX, "%s/libraries", _storePath);
  snprintf(nameR, FILENAME_MAX, "%s/partitions/reads.%04" F_U32P, _storePath, _partitionID);
  snprintf(nameB, FILENAME_MAX, "%s/partitions/blobs.%04" F_U32P, _storePath, _partitionID);

  _librariesAlloc = _info.sqInfo_numLibraries() + 1;
  _readsAlloc     = _readsPerPartition[_partitionID];

  uint64 bs       = AS_UTL_sizeOfFile(nameB);

  _libraries = new sqLibrary [_librariesAlloc];
  _reads     = new sqRead    [_readsAlloc];
  _blobsData = new uint8     [bs];

  AS_UTL_loadFile(nameL, _libraries, _librariesAlloc);
  AS_UTL_loadFile(nameR, _reads,     _readsAlloc);
  AS_UTL_loadFile(nameB, _blobsData,  bs);
}






sqStore::~sqStore() {
  char    No[FILENAME_MAX+1];
  char    Nn[FILENAME_MAX+1];
  uint32  V = 1;

  //  Save original metadata.

  if (_mode == sqStore_extend) {
    snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, _storePath, V);
    while (directoryExists(No) == true) {
      V++;
      snprintf(No, FILENAME_MAX, "%s/version.%03" F_U32P, _storePath, V);
    }

    AS_UTL_mkdir(No);

    snprintf(No, FILENAME_MAX, "%s/libraries", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/libraries", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/reads", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/reads", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info", _storePath, V);
    AS_UTL_rename(No, Nn);

    snprintf(No, FILENAME_MAX, "%s/info.txt", _storePath);
    snprintf(Nn, FILENAME_MAX, "%s/version.%03" F_U32P "/info.txt", _storePath, V);
    AS_UTL_rename(No, Nn);
  }

  //  Recount.

  if ((_mode == sqStore_create) ||
      (_mode == sqStore_extend)) {
    _info.recountReads(_reads);
    _info.setLastBlob(_blobsWriter);
  }

  //  Write updated metadata.

  if ((_mode == sqStore_create) ||
      (_mode == sqStore_extend)) {
    AS_UTL_saveFile(_storePath, '/', "libraries", _libraries, sqStore_getNumLibraries() + 1);
    AS_UTL_saveFile(_storePath, '/', "reads",     _reads,     sqStore_getNumReads()     + 1);
    AS_UTL_saveFile(_storePath, '/', "info",     &_info,                                  1);

    FILE *F = AS_UTL_openOutputFile(_storePath, '/', "info.txt");   //  Used by Canu/Gatekeeper.pm
    _info.writeInfoAsText(F);                                       //  Do not remove!
    AS_UTL_closeFile(F, _storePath, '/', "info.txt");
  }

  //  Write original metadata to the clone.

  if (_mode == sqStore_buildPart) {
    AS_UTL_saveFile(_clonePath, '/', "libraries", _libraries, sqStore_getNumLibraries() + 1);
    AS_UTL_saveFile(_clonePath, '/', "info",     &_info,                                  1);

    FILE *F = AS_UTL_openOutputFile(_clonePath, '/', "info.txt");   //  Used by Canu/Gatekeeper.pm
    _info.writeInfoAsText(F);                                       //  Do not remove!
    AS_UTL_closeFile(F, _clonePath, '/', "info.txt");
  }

  //  Clean up.

  delete [] _libraries;
  delete [] _reads;
  delete [] _blobsData;
  delete [] _blobsFiles;

  delete    _blobsWriter;

  delete [] _readIDtoPartitionIdx;
  delete [] _readIDtoPartitionID;
  delete [] _readsPerPartition;
};



sqStore *
sqStore::sqStore_open(char const *path, sqStore_mode mode, uint32 partID) {

  //  If an instance exists, return it, otherwise, make a new one.

#pragma omp critical
  {
    if (_instance != NULL) {
      _instanceCount++;
    } else {
      _instance      = new sqStore(path, NULL, mode, partID);
      _instanceCount = 1;
    }
  }

  return(_instance);
}



sqStore *
sqStore::sqStore_open(char const *storePath, char const *clonePath) {

  //  Only one instance can be opened at a time.

#pragma omp critical
  {
    assert(_instance == NULL);

    _instance      = new sqStore(storePath, clonePath, sqStore_buildPart, UINT32_MAX);
    _instanceCount = 1;
  }

  return(_instance);
}



void
sqStore::sqStore_close(void) {

#pragma omp critical
  {
    _instanceCount--;

    if (_instanceCount == 0) {
      delete _instance;
      _instance = NULL;
    }
  }
}



void
sqStore::sqStore_delete(void) {
  char path[FILENAME_MAX+1];

  sqStore_deletePartitions();

  snprintf(path, FILENAME_MAX, "%s/info",      _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/libraries", _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads",     _storePath);  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/blobs",     _storePath);  AS_UTL_unlink(path);

  AS_UTL_rmdir(_storePath);
}



void
sqStore::sqStore_deletePartitions(void) {
  char path[FILENAME_MAX+1];

  snprintf(path, FILENAME_MAX, "%s/partitions/map", _storePath);

  if (fileExists(path) == false)
    return;

  //  How many partitions?

  FILE *F = AS_UTL_openInputFile(path);

  loadFromFile(_numberOfPartitions, "sqStore_deletePartitions::numberOfPartitions", F);

  AS_UTL_closeFile(F, path);

  //  Yay!  Delete!

  AS_UTL_unlink(path);

  for (uint32 ii=0; ii<_numberOfPartitions; ii++) {
    snprintf(path, FILENAME_MAX, "%s/partitions/reads.%04u", _storePath, ii+1);  AS_UTL_unlink(path);
    snprintf(path, FILENAME_MAX, "%s/partitions/blobs.%04u", _storePath, ii+1);  AS_UTL_unlink(path);
  }

  //  And the directory.

  snprintf(path, FILENAME_MAX, "%s/partitions", _storePath);

  AS_UTL_rmdir(path);
}


