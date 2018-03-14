
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
 *    Brian P. Walenz from 2014-NOV-26 to 2015-AUG-14
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-09
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef GKSTORE_H
#define GKSTORE_H

#include "AS_global.H"
#include "writeBuffer.H"

#include <vector>

using namespace std;


#define GK_MAGIC   0x504b473a756e6162lu      //  canu:GKP
#define GK_VERSION 0x0000000000000004lu


//  The number of library IIDs we can handle.
//
#define AS_MAX_LIBRARIES_BITS      6
#define AS_MAX_LIBRARIES           (((uint32)1 << AS_MAX_LIBRARIES_BITS) - 1)

#define LIBRARY_NAME_SIZE          128

//  Maximum length of reads.
//
//  If 16, an overlap is only 20 bytes.  (5x 32 bit words)
//  If 17-21, an overlap is 24 bytes.    (3x 64 bit words)
//  If 22-32, an overlap is 32 bytes.    (4x 64 bit words)
//
//  if 26, bogart has issues with storing the error rate
//  If 28, alignment/alignment-drivers.C won't link
//  If 29, alignment/alignment-drivers.C won't link
//  If 30, alignment/alignment-drivers.C won't link
//  If 31, alignment/alignment-drivers.C won't compile, len+len+2 == 0
//  If 32, it won't compile because of shifting (uint32)1 << 32 == 0.
//
#define AS_MAX_READLEN_BITS        21
#define AS_MAX_READLEN             (((uint32)1 << AS_MAX_READLEN_BITS) - 1)

//  The number of read IDs we can handle.  Longer reads implies fewer reads.
//    readLen 32 + numLibs 6 -> numReads 26 ( 64  million)
//    readLen 30 + numLibs 6 -> numReads 28 (256  million)
//    readLen 28 + numLibs 6 -> numReads 30 (1024 million)
//    readLen 26 + numLibs 6 -> numReads 32 (4096 million)  //  limited elsewhere!
//    readLen 24 + numLibs 6 -> numReads 34 (4096 million)  //  limited elsewhere!
//    readLen 22 + numLibs 6 -> numReads 36 (4096 million)  //  limited elsewhere!
//    readLen 21 + numLibs 6 -> numReads 37 (4096 million)  //  limited elsewhere!
//    readLen 20 + numLibs 6 -> numReads 38 (4096 million)  //  limited elsewhere!
//
#define AS_MAX_READS_BITS          64 - AS_MAX_READLEN_BITS - AS_MAX_LIBRARIES_BITS
#define AS_MAX_READS               (((uint64)1 << AS_MAX_READS_BITS) - 1)

#define AS_BLOBFILE_MAX_SIZE       1024 * 1024 * 1024


#include "gkLibrary.H"
#include "gkRead.H"
#include "gkStoreBlobReader.H"
#include "gkStoreBlobWriter.H"


//  The default behavior is to open the store for read only, and to load
//  all the metadata into memory.

typedef enum {
  gkStore_create      = 0x00,  //  Open for creating, will fail if files exist already
  gkStore_extend      = 0x01,  //  Open for modification and appending new reads/libraries
  gkStore_readOnly    = 0x02,  //  Open read only
  gkStore_buildPart   = 0x03   //  For building the partitions
} gkStore_mode;


static
const
char *
toString(gkStore_mode m) {
  switch (m) {
    case gkStore_create:       return("gkStore_create");       break;
    case gkStore_extend:       return("gkStore_extend");       break;
    case gkStore_readOnly:     return("gkStore_readOnly");     break;
    case gkStore_buildPart:    return("gkStore_buildPart");    break;
  }

  return("undefined-mode");
}



class gkStoreInfo {
public:
  gkStoreInfo();
  ~gkStoreInfo();

  bool      checkInfo(void);
  void      recountReads(gkRead *reads);
  void      writeInfoAsText(FILE *F);

  uint64    gkInfo_numLibraries(void)        { return(_numLibraries);      };
  uint64    gkInfo_numReads(void)            { return(_numReads);          };

  uint64    gkInfo_numRawReads(void)         { return(_numRawReads);       };
  uint64    gkInfo_numCorrectedReads(void)   { return(_numCorrectedReads); };
  uint64    gkInfo_numTrimmedReads(void)     { return(_numTrimmedReads);   };

  void      gkInfo_addLibrary(void)          { _numLibraries++; };
  void      gkInfo_addRead(void)             { _numReads++;     };

private:
  uint64    _gkMagic;
  uint64    _gkVersion;

  uint32    _gkLibrarySize;      //  Sanity checks that this code can load the data properly.
  uint32    _gkReadSize;
  uint32    _gkMaxLibrariesBits;
  uint32    _gkLibraryNameSize;
  uint32    _gkMaxReadBits;
  uint32    _gkMaxReadLenBits;

  uint32    _numLibraries;       //  Counts of types of things we have loaded (next
  uint32    _numReads;           //  available index into _libraries and _reads in gkStore)

  uint32    _numRawReads;        //  If any corrected reads exist, the store will return
  uint32    _numCorrectedReads;  //  only corrected reads by default.
  uint32    _numTrimmedReads;

  uint64    _numRawBases;
  uint64    _numCorrectedBases;
  uint64    _numTrimmedBases;
};






class gkStore {

private:
  gkStore(char const *storePath, char const *clonePath, gkStore_mode mode, uint32 partID);
  ~gkStore();

  void         gkStore_loadMetadata(void);
  void         gkStore_checkInfo(void);

public:
  static
  gkStore     *gkStore_open(char const *path, gkStore_mode mode=gkStore_readOnly, uint32 partID=UINT32_MAX);

  static
  gkStore     *gkStore_open(char const *storePath, char const *clonePath);

  void         gkStore_close(void);

public:
  const char  *gkStore_path(void) { return(_storePath); };  //  Returns the path to the store

  void         gkStore_buildPartitions(uint32 *partitionMap);

  void         gkStore_delete(void);             //  Deletes the files in the store.
  void         gkStore_deletePartitions(void);   //  Deletes the files for a partition.

  uint32       gkStore_getNumLibraries(void)       { return(_info.gkInfo_numLibraries()); };

  uint32       gkStore_getNumReads(void)           { return(_info.gkInfo_numReads()); };
  uint32       gkStore_getNumRawReads(void)        { return(_info.gkInfo_numRawReads()); };
  uint32       gkStore_getNumCorrectedReads(void)  { return(_info.gkInfo_numCorrectedReads()); };
  uint32       gkStore_getNumTrimmedReads(void)    { return(_info.gkInfo_numTrimmedReads()); };

  gkLibrary   *gkStore_getLibrary(uint32 id)       { return(&_libraries[id]); };

  //  The ONLY two approved methods for getting a read and data are:
  //    gkStore_getRead(uint32 id)
  //    gkStore_loadReadData(gkRead *read)  -- implies gkStore_getRead() was called already.
  //    gkStore_loadReadData(uint32  id)    -- calls gkStore_getRead(), then loadReadData(gkRead).

  gkRead      *gkStore_getRead(uint32 id);
  void         gkStore_loadReadData(gkRead *read,   gkReadData *readData);
  void         gkStore_loadReadData(uint32  readID, gkReadData *readData);

  void         gkStore_stashReadData(gkReadData *data);

  bool         gkStore_readInPartition(uint32 id) {        //  True if read is in this partition.
    return((_readIDtoPartitionID     == NULL) ||           //    Not partitioned, read in partition!
           (_readIDtoPartitionID[id] == _partitionID));    //    Partitioned, and in this one!
  };

  gkLibrary   *gkStore_addEmptyLibrary(char const *name);
  gkReadData  *gkStore_addEmptyRead(gkLibrary *lib);

  void         gkStore_setClearRange(uint32 id, uint32 bgn, uint32 end);

  //  Used in utgcns, for the package format.
  static
  void         gkStore_loadReadFromStream(FILE *S, gkRead *read, gkReadData *readData);
  void         gkStore_saveReadToStream(FILE *S, uint32 id);

private:
  static gkStore      *_instance;
  static uint32        _instanceCount;

  gkStoreInfo          _info;  //  All the stuff stored on disk.

  char                 _storePath[FILENAME_MAX+1];    //  Needed to create files
  char                 _clonePath[FILENAME_MAX+1];    //  Path to the partitioned store data, for when building the partitions

  gkStore_mode         _mode;  //  What mode this store is opened as, sanity checking


  uint32               _librariesAlloc;  //  Size of allocation
  gkLibrary           *_libraries;       //  In core data

  uint32               _readsAlloc;      //  Size of allocation
  gkRead              *_reads;           //  In core data

  uint8               *_blobsData;       //  For partitioned data, in-core data.

  uint32               _blobsFilesMax;   //  For normal store, loading reads
  gkStoreBlobReader   *_blobsFiles;      //  directly, one per thread.

  gkStoreBlobWriter   *_blobsWriter;

  //  If the store is openend partitioned, this data is loaded from disk

  uint32               _numberOfPartitions;     //  Total number of partitions that exist
  uint32               _partitionID;            //  Which partition this is
  uint32              *_readsPerPartition;      //  Number of reads in each partition, mostly sanity checking
  uint32              *_readIDtoPartitionIdx;   //  Map from global ID to local partition index
  uint32              *_readIDtoPartitionID;    //  Map from global ID to partition ID
};



#endif
