
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


void
gkRead::gkRead_copyDataToPartition(void     *blobs,
                                   FILE    **partfiles,
                                   uint64   *partfileslen,
                                   uint32    partID) {

  if (partID == UINT32_MAX)  //  If an invalid partition, don't do anything.
    return;

  //  Figure out where the blob actually is, and make sure that it really is a blob

  uint8  *blob    = (uint8 *)blobs + _mPtr;
  uint32  blobLen = 8 + *((uint32 *)blob + 1);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  The partfile should be at what we think is the end.

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

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
gkRead::gkRead_copyDataToPartition(FILE    **blobsFiles,
                                   FILE    **partfiles,
                                   uint64   *partfileslen,
                                   uint32    partID) {

  //  Load the blob from disk.  We must always read the data, even if we don't want
  //  to write it.  Or, I suppose, we could skip and seek.

  uint8   tag[4]  = {0};
  uint32  blobLen = 0;
  FILE   *file    = blobsFiles[omp_get_thread_num()];

  //  Ideally, we'd do one read to get the whole blob.  Without knowing
  //  the length, we're forced to do two.  Or maybe three.

  AS_UTL_fseek(file, _mPtr, SEEK_SET);

  AS_UTL_safeRead(file,  tag,     "gkStore::gkStore_loadDataFromFile::tag",     sizeof(int8),   4);
  AS_UTL_safeRead(file, &blobLen, "gkStore::gkStore_loadDataFromFile::blobLen", sizeof(uint32), 1);

  uint8 *blob     = new uint8 [8 + blobLen];

  memcpy(blob,    tag,     sizeof(uint8)  * 4);
  memcpy(blob+4, &blobLen, sizeof(uint32) * 1);

  AS_UTL_safeRead(file, blob+8, "gkStore::gkStore_loadDataFromFile::blob", sizeof(char), blobLen);

  assert(blob[0] == 'B');
  assert(blob[1] == 'L');
  assert(blob[2] == 'O');
  assert(blob[3] == 'B');

  //  Write the data and update pointers and lengths.

  //fprintf(stderr, "COPY blob of length %8u from position _mPtr %8lu to position %8lu in partition %4u\n",
  //        blobLen, _mPtr, partfileslen[partID], partID);

  AS_UTL_safeWrite(partfiles[partID], blob, "gkRead::gkRead_copyDataToPartition::blob", sizeof(char), blobLen + 8);

  _mPtr = partfileslen[partID];   //  Update the read to point to this data
  _pID  = partID;                 //  in the new blob and partition.

  partfileslen[partID] += blobLen + 8;

  assert(partfileslen[partID] == AS_UTL_ftell(partfiles[partID]));

  delete [] blob;
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

  fprintf(stderr, "Found " F_U32 " unpartitioned reads and maximum partition of " F_U32 "\n",
          unPartitioned, maxPartition);

  //  Create the partitions by opening N copies of the data stores,
  //  and writing data to each.

  FILE         **blobfiles    = new FILE * [maxPartition + 1];
  uint64        *blobfileslen = new uint64 [maxPartition + 1];            //  Offset, in bytes, into the blobs file
  FILE         **readfiles    = new FILE * [maxPartition + 1];
  uint32        *readfileslen = new uint32 [maxPartition + 1];            //  aka _readsPerPartition
  uint32        *readIDmap    = new uint32 [gkStore_getNumReads() + 1];   //  aka _readIDtoPartitionIdx

  //  Be nice and put all the partitions in a subdirectory.

  snprintf(name, FILENAME_MAX, "%s/partitions", _storePath);

  if (AS_UTL_fileExists(name, true, true) == false)
    AS_UTL_mkdir(name);

  //  Open all the output files -- fail early if we can't open that many files.

  blobfiles[0]    = NULL;
  blobfileslen[0] = UINT64_MAX;
  readfiles[0]    = NULL;
  readfileslen[0] = UINT32_MAX;

  for (uint32 i=1; i<=maxPartition; i++) {
    snprintf(name, FILENAME_MAX, "%s/partitions/blobs.%04d", _storePath, i);

    errno = 0;
    blobfiles[i]    = fopen(name, "w");
    blobfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);

    snprintf(name, FILENAME_MAX, "%s/partitions/reads.%04d", _storePath, i);

    errno = 0;
    readfiles[i]    = fopen(name, "w");
    readfileslen[i] = 0;

    if (errno)
      fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition %u file '%s' for write: %s\n",
              i, name, strerror(errno)), exit(1);
  }

  //  Open the output partition map file -- we might as well fail early if we can't make it also.

  snprintf(name, FILENAME_MAX, "%s/partitions/map", _storePath);

  errno = 0;
  FILE *rIDmF = fopen(name, "w");
  if (errno)
    fprintf(stderr, "gkStore::gkStore_buildPartitions()-- ERROR: failed to open partition map file '%s': %s\n",
            name, strerror(errno)), exit(1);

  //  Copy the blob from the master file to the partitioned file, update pointers.

  readIDmap[0] = UINT32_MAX;    //  There isn't a zeroth read, make it bogus.

  for (uint32 fi=1; fi<=gkStore_getNumReads(); fi++) {
    uint32  pi = partitionMap[fi];

    if (pi == UINT32_MAX)
      continue;

    assert(pi != 0);  //  No zeroth partition, right?

    //  Make a copy of the read, then modify it for the partition, then write it to the partition.
    //  Without the copy, we'd need to update the master record too.

    gkRead  partRead = _reads[fi];

    if (_blobs)
      partRead.gkRead_copyDataToPartition(_blobs, blobfiles, blobfileslen, pi);
    if (_blobsFiles)
      partRead.gkRead_copyDataToPartition(_blobsFiles, blobfiles, blobfileslen, pi);

    //  Because the blobsFiles copyDataToPartition variant is streaming through the file,
    //  we need to let it load (and ignore) deleted reads.  After they're loaded (and ignored)
    //  we can then skip it.

    if (pi < UINT32_MAX) {
      fprintf(stderr, "read " F_U32 "=" F_U32 " len " F_U32 " -- blob master " F_U64 " -- to part " F_U32 " new read id " F_U32 " part " F_U64 " blob " F_U64 " -- at readIdx " F_U32 "\n",
              fi, _reads[fi].gkRead_readID(), _reads[fi].gkRead_sequenceLength(),
              _reads[fi]._mPtr,
              pi,
              partRead.gkRead_readID(), partRead._pID, partRead._mPtr,
              readfileslen[pi]);

      AS_UTL_safeWrite(readfiles[pi], &partRead, "gkStore::gkStore_buildPartitions::read", sizeof(gkRead), 1);

      readIDmap[fi] = readfileslen[pi]++;
    }

    else {
      fprintf(stderr, "read " F_U32 "=" F_U32 " len " F_U32 " -- blob master " F_U64 " -- DELETED\n",
              fi, _reads[fi].gkRead_readID(), _reads[fi].gkRead_sequenceLength(),
              _reads[fi]._mPtr);
    }
  }

  //  There isn't a zeroth read.

  AS_UTL_safeWrite(rIDmF, &maxPartition,  "gkStore::gkStore_buildPartitions::maxPartition", sizeof(uint32), 1);
  AS_UTL_safeWrite(rIDmF,  readfileslen,  "gkStore::gkStore_buildPartitions::readfileslen", sizeof(uint32), maxPartition + 1);
  AS_UTL_safeWrite(rIDmF,  partitionMap,  "gkStore::gkStore_buildPartitions::partitionMap", sizeof(uint32), gkStore_getNumReads() + 1);
  AS_UTL_safeWrite(rIDmF,  readIDmap,     "gkStore::gkStore_buildPartitions::readIDmap",    sizeof(uint32), gkStore_getNumReads() + 1);

  //  cleanup -- close all the files, delete storage

  AS_UTL_closeFile(rIDmF, name);

  for (uint32 i=1; i<=maxPartition; i++) {
    fprintf(stderr, "partition " F_U32 " has " F_U32 " reads\n", i, readfileslen[i]);

    AS_UTL_closeFile(blobfiles[i]);
    AS_UTL_closeFile(readfiles[i]);
  }

  delete [] readIDmap;
  delete [] readfileslen;
  delete [] readfiles;
  delete [] blobfileslen;
  delete [] blobfiles;
}


