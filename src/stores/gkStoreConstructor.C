
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



//  The valid modes for a 'new gkpStore' call:
//
//  1)  Add new reads/libraries, modify old ones.  gkStore(path, true, true)
//  2)  No addition, but can modify old ones.      gkStore(path, true)
//  3)  No addition, no modification.              gkStore(path);




void
gkStore::gkStore_loadMetadata(void) {
  char    name[FILENAME_MAX];


  _librariesAlloc = _info.numLibraries + 1;
  _libraries      = new gkLibrary [_librariesAlloc];

  snprintf(name, FILENAME_MAX, "%s/libraries", _storePath);
  AS_UTL_loadFile(name, _libraries, _librariesAlloc);


  _readsAlloc = _info.numReads + 1;
  _reads      = new gkRead [_readsAlloc];

  snprintf(name, FILENAME_MAX, "%s/reads", _storePath);
  AS_UTL_loadFile(name, _reads, _readsAlloc);
}



void
gkStore::gkStore_openBlobs(void) {
  char    name[FILENAME_MAX];

  snprintf(name, FILENAME_MAX, "%s/blobs", _storePath);

  _blobsFiles    = new FILE * [omp_get_max_threads()];
  for (uint32 ii=0; ii<omp_get_max_threads(); ii++)
    _blobsFiles[ii] = AS_UTL_openInputFile(_storePath, '/', "blobs");
}



gkStore::gkStore(char const *path, gkStore_mode mode, uint32 partID) {
  char    name[FILENAME_MAX];

  memset(_storePath, 0, sizeof(char) * FILENAME_MAX);
  memset(_storeName, 0, sizeof(char) * FILENAME_MAX);

  strncpy(_storePath, path, FILENAME_MAX-1);
  strncpy(_storeName, path, FILENAME_MAX-1);  //  Broken.

  //  If the info file exists, load it.

  snprintf(name, FILENAME_MAX, "%s/info", _storePath);

  if (AS_UTL_fileExists(name, false, false) == true) {
    FILE *I = AS_UTL_openInputFile(name);
    AS_UTL_safeRead(I, &_info, "gkStore::_info", sizeof(gkStoreInfo), 1);
    AS_UTL_closeFile(I, name);
  }

  //  Check sizes are correct.

  uint32  failed = 0;

  if (_info.gkLibrarySize      != sizeof(gkLibrary))
    failed += fprintf(stderr, "ERROR:  gkLibrary size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _info.gkLibrarySize, sizeof(gkLibrary));

  if (_info.gkReadSize         != sizeof(gkRead))
    failed += fprintf(stderr, "ERROR:  gkRead size in store = " F_U32 ", differs from executable = " F_SIZE_T "\n",
                      _info.gkReadSize, sizeof(gkRead));

  if (_info.gkMaxLibrariesBits != AS_MAX_LIBRARIES_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_LIBRARIES_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _info.gkMaxLibrariesBits, AS_MAX_LIBRARIES_BITS);

  if (_info.gkLibraryNameSize  != LIBRARY_NAME_SIZE)
    failed += fprintf(stderr, "ERROR:  LIBRARY_NAME_SIZE in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _info.gkLibraryNameSize, LIBRARY_NAME_SIZE);

  if (_info.gkMaxReadBits      != AS_MAX_READS_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READS_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _info.gkMaxReadBits, AS_MAX_READS_BITS);

  if (_info.gkMaxReadLenBits   != AS_MAX_READLEN_BITS)
    failed += fprintf(stderr, "ERROR:  AS_MAX_READLEN_BITS in store = " F_U32 ", differs from executable = " F_U32 "\n",
                      _info.gkMaxReadLenBits, AS_MAX_READLEN_BITS);

  if (failed)
    fprintf(stderr, "ERROR:\nERROR:  Can't open store '%s': parameters in src/AS_global.H are incompatible with the store.\n", _storePath), exit(1);

  assert(_info.gkLibrarySize      == sizeof(gkLibrary));
  assert(_info.gkReadSize         == sizeof(gkRead));

  assert(_info.gkMaxLibrariesBits == AS_MAX_LIBRARIES_BITS);
  assert(_info.gkLibraryNameSize  == LIBRARY_NAME_SIZE);
  assert(_info.gkMaxReadBits      == AS_MAX_READS_BITS);
  assert(_info.gkMaxReadLenBits   == AS_MAX_READLEN_BITS);

  //  Clear ourself, to make valgrind happier.

  _librariesAlloc         = 0;
  _libraries              = NULL;

  _readsAlloc             = 0;
  _reads                  = NULL;

  _blobs                  = NULL;
  _blobsWriter            = NULL;
  _blobsFiles             = NULL;

  _mode                   = mode;

  _numberOfPartitions     = 0;
  _partitionID            = 0;
  _readIDtoPartitionIdx   = NULL;
  _readIDtoPartitionID    = NULL;
  _readsPerPartition      = NULL;
  //_readsInThisPartition   = NULL;


  //
  //  CREATE - allocate some memory for saving libraries and reads, and create a file to dump the data into.
  //

  if (mode == gkStore_create) {
    if (partID != UINT32_MAX)
      fprintf(stderr, "gkStore()-- Illegal combination of gkStore_create with defined partID.\n"), exit(1);

    if (AS_UTL_fileExists(_storePath, true, true) == true)
      fprintf(stderr, "ERROR:  Can't create store '%s': store already exists.\n", _storePath), exit(1);

    AS_UTL_mkdir(_storePath);

    _librariesAlloc = 32;
    _libraries      = new gkLibrary [_librariesAlloc];

    _readsAlloc     = 32768;
    _reads          = new gkRead [_readsAlloc];

    snprintf(name, FILENAME_MAX, "%s/blobs", _storePath);
    _blobsWriter   = new writeBuffer(name, "a");

    return;
  }


  //
  //  EXTEND - load libraries and reads into core, open blobs for reading and appending.
  //

  if ((mode == gkStore_extend) &&
      (partID != UINT32_MAX))
    fprintf(stderr, "gkStore()-- Illegal combination of gkStore_extend with defined partID.\n"), exit(1);

  if (mode == gkStore_extend) {
    gkStore_loadMetadata();
    gkStore_openBlobs();

    snprintf(name, FILENAME_MAX, "%s/blobs", _storePath);
    _blobsWriter   = new writeBuffer(name, "a");

    return;
  }




  //
  //  READ ONLY - Two options: normal and partitioned.
  //
  //  If normal, 


  assert(mode == gkStore_readOnly);

  if (AS_UTL_fileExists(_storePath, true, false) == false)
    fprintf(stderr, "gkStore()--  failed to open '%s' for read-only access: store doesn't exist.\n", _storePath), exit(1);


  //  If normal, nothing special; load the metadata and open the blob files, one file per thread.

  if (partID == UINT32_MAX) {
    gkStore_loadMetadata();
    gkStore_openBlobs();
  }


  //  If partitioned, we need to have a uint32 map of readID to partitionReadID so we can
  //  lookup the metadata in the partitoned _reads data.  This is 4 bytes per read, compared to 24
  //  bytes for the full meta data.  Assuming 100x of 3kb read coverage on human, that's 100
  //  million reads, so 0.400 GB vs 2.4 GB.
  

  if (partID != UINT32_MAX) {
    snprintf(name, FILENAME_MAX, "%s/partitions/map", _storePath);

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

    AS_UTL_closeFile(F, name);

    //  Libraries are easy, just the normal file.

    snprintf(name, FILENAME_MAX, "%s/libraries", _storePath);

    _librariesAlloc = _info.numLibraries + 1;
    _libraries      = new gkLibrary [_librariesAlloc];

    AS_UTL_loadFile(name, _libraries, _librariesAlloc);

    //  Reads are partitioned.

    snprintf(name, FILENAME_MAX, "%s/partitions/reads.%04" F_U32P, _storePath, partID);

    _readsAlloc = _readsPerPartition[partID];
    _reads      = new gkRead [_readsAlloc];

    AS_UTL_loadFile(name, _reads, _readsAlloc);

    //  Blobs are partitioned too.  Sadly, we don't know the size.

    snprintf(name, FILENAME_MAX, "%s/partitions/blobs.%04" F_U32P, _storePath, partID);

    uint64 bs      = AS_UTL_sizeOfFile(name);
    _blobs         = new uint8 [bs];

    AS_UTL_loadFile(name, _blobs, bs);
  }
}






gkStore::~gkStore() {

  //  If we've potentially changed the store, recount the number of raw, corrected and trimmed
  //  reads.  Then dump library and read meta data, then the info files.

  //  Canu uses info.txt (and libraries.txt, but it makes that one) to determine the number of
  //  reads/bases in the store.

  if ((_mode == gkStore_create) ||
      (_mode == gkStore_extend)) {
    FILE  *F;

    _info.numRawReads = _info.numCorrectedReads = _info.numTrimmedReads = 0;
    _info.numRawBases = _info.numCorrectedBases = _info.numTrimmedBases = 0;

    for (uint32 ii=0; ii<_info.numReads + 1; ii++) {
      if (_reads[ii]._rseqLen > 0) {
        _info.numRawReads++;
        _info.numRawBases += _reads[ii]._rseqLen;
      }

      if (_reads[ii]._cseqLen > 0) {
        _info.numCorrectedReads++;
        _info.numCorrectedBases += _reads[ii]._cseqLen;
      }

      if (_reads[ii]._clearBgn < _reads[ii]._clearEnd) {
        _info.numTrimmedReads++;
        _info.numTrimmedBases += _reads[ii]._clearEnd - _reads[ii]._clearBgn;
      }
    }

    F = AS_UTL_openOutputFile(gkStore_path(), '/', "libraries");
    AS_UTL_safeWrite(F, _libraries, "libraries", sizeof(gkLibrary), gkStore_getNumLibraries() + 1);
    AS_UTL_closeFile(F);

    F = AS_UTL_openOutputFile(gkStore_path(), '/', "reads");
    AS_UTL_safeWrite(F, _reads, "reads", sizeof(gkRead), gkStore_getNumReads() + 1);
    AS_UTL_closeFile(F);

    F = AS_UTL_openOutputFile(gkStore_path(), '/', "info");
    AS_UTL_safeWrite(F, &_info, "info", sizeof(gkStoreInfo), 1);
    AS_UTL_closeFile(F);

    F = AS_UTL_openOutputFile(gkStore_path(), '/', "info.txt");
    _info.writeInfoAsText(F);
    AS_UTL_closeFile(F);
  }

  //  Clean up.

  delete [] _libraries;
  delete [] _reads;
  delete [] _blobs;
  delete    _blobsWriter;

  for (uint32 ii=0; ii<omp_get_max_threads(); ii++)
    if ((_blobsFiles) && (_blobsFiles[ii]))
      AS_UTL_closeFile(_blobsFiles[ii]);

  delete [] _blobsFiles;

  delete [] _readIDtoPartitionIdx;
  delete [] _readIDtoPartitionID;
  delete [] _readsPerPartition;
};



gkStore *
gkStore::gkStore_open(char const *path, gkStore_mode mode, uint32 partID) {

  //  If an instance exists, return it, otherwise, make a new one.

#pragma omp critical
  {
    if (_instance != NULL) {
      _instanceCount++;
      //fprintf(stderr, "gkStore_open(%s) from thread %d, %u instances now\n", path, omp_get_thread_num(), _instanceCount);
    } else {
      _instance      = new gkStore(path, mode, partID);
      _instanceCount = 1;
      //fprintf(stderr, "gkStore_open(%s) form thread %d, first instance, create store\n", path, omp_get_thread_num());
    }
  }

  return(_instance);
}



void
gkStore::gkStore_close(void) {

#pragma omp critical
  {
    _instanceCount--;

    if (_instanceCount == 0) {
      delete _instance;
      _instance = NULL;
      //fprintf(stderr, "gkStore_close(%s) from thread %d, no instances remain, delete store\n",
      //        _storeName, omp_get_thread_num());
    }

    else {
      //fprintf(stderr, "gkStore_close(%s) from thread %d, %u instances remain\n",
      //        _storeName, omp_get_thread_num(), _instanceCount);
    }
  }
}



void
gkStore::gkStore_clone(char *originalPath, char *clonePath) {
  char cPath[FILENAME_MAX];
  char sPath[FILENAME_MAX];

  getcwd(cPath, FILENAME_MAX);

  AS_UTL_mkdir(clonePath);

  chdir(clonePath);

  snprintf(sPath, FILENAME_MAX, "%s/info",      originalPath);
  AS_UTL_symlink(sPath, "info");

  snprintf(sPath, FILENAME_MAX, "%s/libraries", originalPath);
  AS_UTL_symlink(sPath, "libraries");

  snprintf(sPath, FILENAME_MAX, "%s/reads",     originalPath);
  AS_UTL_symlink(sPath, "reads");

  snprintf(sPath, FILENAME_MAX, "%s/blobs",     originalPath);
  AS_UTL_symlink(sPath, "blobs");

  chdir(cPath);
}



void
gkStore::gkStore_delete(void) {
  char path[FILENAME_MAX];

  gkStore_deletePartitions();

  snprintf(path, FILENAME_MAX, "%s/info",      gkStore_path());  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/libraries", gkStore_path());  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/reads",     gkStore_path());  AS_UTL_unlink(path);
  snprintf(path, FILENAME_MAX, "%s/blobs",     gkStore_path());  AS_UTL_unlink(path);

  AS_UTL_rmdir(gkStore_path());
}



void
gkStore::gkStore_deletePartitions(void) {
  char path[FILENAME_MAX];

  snprintf(path, FILENAME_MAX, "%s/partitions/map", gkStore_path());

  if (AS_UTL_fileExists(path, false, false) == false)
    return;

  //  How many partitions?

  FILE *F = AS_UTL_openInputFile(path);

  AS_UTL_safeRead(F, &_numberOfPartitions, "gkStore_deletePartitions::numberOfPartitions", sizeof(uint32), 1);

  AS_UTL_closeFile(F, path);

  //  Yay!  Delete!

  AS_UTL_unlink(path);

  for (uint32 ii=0; ii<_numberOfPartitions; ii++) {
    snprintf(path, FILENAME_MAX, "%s/partitions/reads.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
    snprintf(path, FILENAME_MAX, "%s/partitions/blobs.%04u", gkStore_path(), ii+1);  AS_UTL_unlink(path);
  }

  //  And the directory.

  snprintf(path, FILENAME_MAX, "%s/partitions", gkStore_path());

  AS_UTL_rmdir(path);
}


