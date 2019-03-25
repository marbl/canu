
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
 *    src/stores/gkStorePartition.C
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



void
sqStore::sqStore_buildPartitions(uint32 *partitionMap) {
  char              name[FILENAME_MAX];

  //  Store cannot be partitioned already, and it must be readOnly (for safety) as we don't need to
  //  be changing any of the normal store data.

  assert(_clonePath[0]       != 0);     //  Clone path must be set.

  assert(_blobsData          == NULL);  //  Can't already be partitioned!
  assert(_numberOfPartitions == 0);

  assert(_mode               == sqStore_buildPart);

  //  Figure out what the last partition is

  uint32  maxPartition       = 0;
  uint32  readsPartitioned   = 0;
  uint32  readsUnPartitioned = 0;

  assert(partitionMap[0] == UINT32_MAX);

  for (uint32 fi=1; fi<=sqStore_getNumReads(); fi++) {
    if (partitionMap[fi] == UINT32_MAX) {
      readsUnPartitioned++;
      continue;
    }

    readsPartitioned++;

    if (maxPartition < partitionMap[fi])
      maxPartition = partitionMap[fi];
  }

  fprintf(stderr, "Creating " F_U32 " partitions with " F_U32 " reads.  Ignoring " F_U32 " reads.\n",
          maxPartition, readsPartitioned, readsUnPartitioned);

  //  Create the partitions by opening N copies of the data stores,
  //  and writing data to each.

  FILE         **partfiles    = new FILE * [maxPartition + 1];
  uint64        *partfileslen = new uint64 [maxPartition + 1];            //  Offset, in bytes, into the blobs file
  FILE         **readfiles    = new FILE * [maxPartition + 1];
  uint32        *readfileslen = new uint32 [maxPartition + 1];            //  aka _readsPerPartition
  uint32        *readIDmap    = new uint32 [sqStore_getNumReads() + 1];   //  aka _readIDtoPartitionIdx

  //  Be nice and put all the partitions in a subdirectory.

  if (directoryExists(_clonePath) == false)
    AS_UTL_mkdir(_clonePath);

  snprintf(name, FILENAME_MAX, "%s/partitions", _clonePath);

  if (directoryExists(name) == false)
    AS_UTL_mkdir(name);

  //  Open all the output files -- fail early if we can't open that many files.

  partfiles[0]    = NULL;
  partfileslen[0] = UINT64_MAX;
  readfiles[0]    = NULL;
  readfileslen[0] = UINT32_MAX;

  for (uint32 i=1; i<=maxPartition; i++) {
    snprintf(name, FILENAME_MAX, "%s/partitions/blobs.%04d", _clonePath, i);
    partfiles[i]    = AS_UTL_openOutputFile(name);
    partfileslen[i] = 0;

    snprintf(name, FILENAME_MAX, "%s/partitions/reads.%04d", _clonePath, i);
    readfiles[i]    = AS_UTL_openOutputFile(name);
    readfileslen[i] = 0;
  }

  FILE *mapFile = AS_UTL_openOutputFile(_clonePath, '/', "partitions/map");

  //  Copy the blob from the master file to the partitioned file, update pointers.

  readIDmap[0] = UINT32_MAX;    //  There isn't a zeroth read, make it bogus.

  for (uint32 fi=1; fi<=sqStore_getNumReads(); fi++) {
    uint32  pi = partitionMap[fi];

    //  Skip reads not in a partition.

    if (pi == UINT32_MAX)
      continue;

    assert(pi != 0);  //  No zeroth partition, right?

    //  Load the blob from disk.  We must always read the data, even if we don't want
    //  to write it.  Or, I suppose, we could skip and seek.
    //
    //  Ideally, we'd do one read to get the whole blob.  Without knowing
    //  the length, we're forced to do two.  Or maybe three.

    uint8   tag[4]  = {0};
    uint32  blobLen = 0;

    FILE *blobFile = _blobsFiles[omp_get_thread_num()].getFile(_storePath, &_reads[fi]);  //  NOTE!  _storePath for original data!

    loadFromFile(tag,     "sqStore::sqStore_buildPartitions::tag",     4, blobFile);
    loadFromFile(blobLen, "sqStore::sqStore_buildPartitions::blobLen",    blobFile);

    uint8 *blob     = new uint8 [8 + blobLen];

    memcpy(blob,    tag,     sizeof(uint8)  * 4);
    memcpy(blob+4, &blobLen, sizeof(uint32) * 1);

    loadFromFile(blob+8, "sqStore::sqStore_buildPartitions::blob", blobLen, blobFile);

    assert(blob[0] == 'B');
    assert(blob[1] == 'L');
    assert(blob[2] == 'O');
    assert(blob[3] == 'B');

    //  Write the data and update pointers and lengths.

    //fprintf(stderr, "READ %7u copy blob of length %8u from position mPart %4lu mByte %8lu ---> partition %4u position %8lu\n",
    //        fi, blobLen, _reads[fi].sqRead_mPart(), _reads[fi].sqRead_mByte(), pi, partfileslen[pi]);

    //  Make a copy of the read, then modify it to point to the new blob, then write it to the partition.

    sqRead  partRead = _reads[fi];

    partRead._mSegm     = 0;
    partRead._mByteHigh = partfileslen[pi] >> 32;
    partRead._mByteLow  = partfileslen[pi]  & 0xffffffffllu;
    partRead._mPart     = pi;

    //  Write the data.

    writeToFile(blob,     "sqRead::sqRead_buildPartitions::blob",   blobLen + 8, partfiles[pi]);
    writeToFile(partRead, "sqStore::sqStore_buildPartitions::read",              readfiles[pi]);

    delete [] blob;

    //  Update position pointers.

    readIDmap[fi]     = readfileslen[pi];

    partfileslen[pi] += blobLen + 8;               assert(partfileslen[pi] == AS_UTL_ftell(partfiles[pi]));
    readfileslen[pi] += 1;
  }

  //  There isn't a zeroth read.

  writeToFile(maxPartition,  "sqStore::sqStore_buildPartitions::maxPartition",                            mapFile);
  writeToFile(readfileslen,  "sqStore::sqStore_buildPartitions::readfileslen", maxPartition + 1,          mapFile);
  writeToFile(partitionMap,  "sqStore::sqStore_buildPartitions::partitionMap", sqStore_getNumReads() + 1, mapFile);
  writeToFile(readIDmap,     "sqStore::sqStore_buildPartitions::readIDmap",    sqStore_getNumReads() + 1, mapFile);

  //  cleanup -- close all the files, delete storage

  AS_UTL_closeFile(mapFile, _clonePath, '/', "partitions/map");

  for (uint32 i=1; i<=maxPartition; i++) {
    snprintf(name, FILENAME_MAX, "%s/partitions/blobs.%04d", _clonePath, i);
    AS_UTL_closeFile(partfiles[i], name);

    snprintf(name, FILENAME_MAX, "%s/partitions/reads.%04d", _clonePath, i);
    AS_UTL_closeFile(readfiles[i], name);
  }

  delete [] readIDmap;
  delete [] readfileslen;
  delete [] readfiles;
  delete [] partfileslen;
  delete [] partfiles;

  fprintf(stderr, "Partitions created.  Bye.\n");
}


