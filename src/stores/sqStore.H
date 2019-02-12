
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
 *    src/stores/gkStore.H
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

#ifndef SQSTORE_H
#define SQSTORE_H

#include "AS_global.H"
#include "files.H"

#include <vector>

using namespace std;


#define SQ_MAGIC   0x504b473a756e6162lu      //  canu:GKP
#define SQ_VERSION 0x0000000000000006lu


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


#include "sqLibrary.H"
#include "sqRead.H"
#include "sqStoreBlobReader.H"
#include "sqStoreBlobWriter.H"


//  The default behavior is to open the store for read only, and to load
//  all the metadata into memory.

typedef enum {
  sqStore_create      = 0x00,  //  Open for creating, will fail if files exist already
  sqStore_extend      = 0x01,  //  Open for modification and appending new reads/libraries
  sqStore_readOnly    = 0x02,  //  Open read only
  sqStore_buildPart   = 0x03   //  For building the partitions
} sqStore_mode;


static
const
char *
toString(sqStore_mode m) {
  switch (m) {
    case sqStore_create:       return("sqStore_create");       break;
    case sqStore_extend:       return("sqStore_extend");       break;
    case sqStore_readOnly:     return("sqStore_readOnly");     break;
    case sqStore_buildPart:    return("sqStore_buildPart");    break;
  }

  return("undefined-mode");
}



class sqStoreInfo {
public:
  sqStoreInfo();
  ~sqStoreInfo();

  bool      checkInfo(void);

  void      recountReads(sqRead *reads);
  void      setLastBlob(sqStoreBlobWriter *writer);

  void      writeInfoAsText(FILE *F);

  uint64    sqInfo_numLibraries(void)        { return(_numLibraries);      };
  uint64    sqInfo_numReads(void)            { return(_numReads);          };
  uint32    sqInfo_numBlobs(void)            { return(_numBlobs);          };

  uint64    sqInfo_numRawReads(void)         { return(_numRawReads);       };
  uint64    sqInfo_numCorrectedReads(void)   { return(_numCorrectedReads); };
  uint64    sqInfo_numTrimmedReads(void)     { return(_numTrimmedReads);   };

  void      sqInfo_addLibrary(void)          { _numLibraries++; };
  void      sqInfo_addRead(void)             { _numReads++;     };
  void      sqInfo_addBlob(void)             { _numBlobs++;     };

private:
  uint64    _sqMagic;
  uint64    _sqVersion;

  uint32    _sqLibrarySize;      //  Sanity checks that this code can load the data properly.
  uint32    _sqReadSize;
  uint32    _sqMaxLibrariesBits;
  uint32    _sqLibraryNameSize;
  uint32    _sqMaxReadBits;
  uint32    _sqMaxReadLenBits;

  uint32    _numLibraries;       //  Counts of types of things we have loaded (next
  uint32    _numReads;           //  available index into _libraries and _reads in sqStore)
  uint32    _numBlobs;

  uint32    _numRawReads;        //  If any corrected reads exist, the store will return
  uint32    _numCorrectedReads;  //  only corrected reads by default.
  uint32    _numTrimmedReads;

  uint64    _numRawBases;
  uint64    _numCorrectedBases;
  uint64    _numTrimmedBases;
};






class sqStore {

private:
  sqStore(char const *storePath, char const *clonePath, sqStore_mode mode, uint32 partID);
  ~sqStore();

  void         sqStore_loadMetadata(void);
  void         sqStore_checkInfo(void);

public:
  static
  sqStore     *sqStore_open(char const *path, sqStore_mode mode=sqStore_readOnly, uint32 partID=UINT32_MAX);

  static
  sqStore     *sqStore_open(char const *storePath, char const *clonePath);

  void         sqStore_close(void);

public:
  const char  *sqStore_path(void) { return(_storePath); };  //  Returns the path to the store

  void         sqStore_buildPartitions(uint32 *partitionMap);

  void         sqStore_delete(void);             //  Deletes the files in the store.
  void         sqStore_deletePartitions(void);   //  Deletes the files for a partition.

  uint32       sqStore_getNumLibraries(void)       { return(_info.sqInfo_numLibraries()); };

  uint32       sqStore_getNumReads(void)           { return(_info.sqInfo_numReads()); };
  uint32       sqStore_getNumRawReads(void)        { return(_info.sqInfo_numRawReads()); };
  uint32       sqStore_getNumCorrectedReads(void)  { return(_info.sqInfo_numCorrectedReads()); };
  uint32       sqStore_getNumTrimmedReads(void)    { return(_info.sqInfo_numTrimmedReads()); };

  sqLibrary   *sqStore_getLibrary(uint32 id)       { return(&_libraries[id]); };

  //  The ONLY two approved methods for getting a read and data are:
  //    sqStore_getRead(uint32 id)
  //    sqStore_loadReadData(sqRead *read)  -- implies sqStore_getRead() was called already.
  //    sqStore_loadReadData(uint32  id)    -- calls sqStore_getRead(), then loadReadData(sqRead).

  uint8       *sqStore_loadReadBlob(uint32 readID);  //  Returns encoded data.

  sqRead      *sqStore_getRead(uint32 id);
  void         sqStore_loadReadData(sqRead *read,   sqReadData *readData);
  void         sqStore_loadReadData(uint32  readID, sqReadData *readData);

  void         sqStore_stashReadData(sqReadData *data);

  bool         sqStore_readInPartition(uint32 id) {        //  True if read is in this partition.
    return((_readIDtoPartitionID     == NULL) ||           //    Not partitioned, read in partition!
           (_readIDtoPartitionID[id] == _partitionID));    //    Partitioned, and in this one!
  };

  sqLibrary   *sqStore_addEmptyLibrary(char const *name);
  sqReadData  *sqStore_addEmptyRead(sqLibrary *lib);

  void         sqStore_setClearRange(uint32 id, uint32 bgn, uint32 end);
  void         sqStore_setIgnore(uint32 id);

  //  Used in utgcns, for the package format.  Needs to be static for use in tgTig::importData().
  static
  void         sqStore_loadReadFromStream(FILE *S, sqRead *read, sqReadData *readData);
  void         sqStore_saveReadToStream(FILE *S, uint32 id);

private:
  static sqStore      *_instance;
  static uint32        _instanceCount;

  sqStoreInfo          _info;  //  All the stuff stored on disk.

  char                 _storePath[FILENAME_MAX+1];    //  Needed to create files
  char                 _clonePath[FILENAME_MAX+1];    //  Path to the partitioned store data, for when building the partitions

  sqStore_mode         _mode;  //  What mode this store is opened as, sanity checking


  uint32               _librariesAlloc;  //  Size of allocation
  sqLibrary           *_libraries;       //  In core data

  uint32               _readsAlloc;      //  Size of allocation
  sqRead              *_reads;           //  In core data

  uint8               *_blobsData;       //  For partitioned data, in-core data.

  uint32               _blobsFilesMax;   //  For normal store, loading reads
  sqStoreBlobReader   *_blobsFiles;      //  directly, one per thread.

  sqStoreBlobWriter   *_blobsWriter;

  //  If the store is openend partitioned, this data is loaded from disk

  uint32               _numberOfPartitions;     //  Total number of partitions that exist
  uint32               _partitionID;            //  Which partition this is
  uint32              *_readsPerPartition;      //  Number of reads in each partition, mostly sanity checking
  uint32              *_readIDtoPartitionIdx;   //  Map from global ID to local partition index
  uint32              *_readIDtoPartitionID;    //  Map from global ID to partition ID
};

#endif  //  SQSTORE_H
