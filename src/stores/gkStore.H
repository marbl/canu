
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
#include "memoryMappedFile.H"

#include <vector>

using namespace std;

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


//  Per-library options.

//  Read type

#define GK_READTYPE_GENERIC            0x0000
#define GK_READTYPE_CONTIG             0x0001
#define GK_READTYPE_PACBIO_RAW         0x0002
#define GK_READTYPE_PACBIO_CORRECTED   0x0003
#define GK_READTYPE_NANOPORE_RAW       0x0004
#define GK_READTYPE_NANOPORE_CORRECTED 0x0005

//  Correction algorithm

#define GK_CORRECTION_NONE             0x0000
#define GK_CORRECTION_CONSENSUS        0x0001
#define GK_CORRECTION_MER              0x0002

//  Trimming algorithm

#define GK_FINALTRIM_NONE              0x0000
#define GK_FINALTRIM_LARGEST_COVERED   0x0001  //  largest region covered by good overlaps
#define GK_FINALTRIM_BEST_EDGE         0x0002  //  largest region covered by best overlaps (broken)



class gkLibrary {
public:
  gkLibrary() {
    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);
    strncpy(_libraryName, "UNDEFINED", LIBRARY_NAME_SIZE-1);

    gkLibrary_clearFeatures();
  };

  gkLibrary(char const *name) {
    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);
    strncpy(_libraryName, name, LIBRARY_NAME_SIZE-1);

    gkLibrary_clearFeatures();
  };

  ~gkLibrary() {
  };

  void gkLibrary_clearFeatures(void) {
    //  DO NOT change defaults without updating gkLibrary_parsePreset().
    _isNonRandom                = 0;

    _readCorrection             = GK_CORRECTION_NONE;
    _readType                   = GK_READTYPE_GENERIC;

    _finalTrim                  = GK_FINALTRIM_LARGEST_COVERED;
    _removeDuplicateReads       = 1;
    _removeSpurReads            = 1;
    _removeChimericReads        = 1;
    _checkForSubReads           = 1;

    _defaultQV                  = 20;
  };

public:

  char const    *gkLibrary_libraryName(void)           { return(_libraryName);          };
  uint32         gkLibrary_libraryID(void)             { return(_libraryID);            };

  uint32         gkLibrary_isNonRandom(void)           { return(_isNonRandom);          };
  uint32         gkLibrary_readType(void)              { return(_readType);             };
  char const    *gkLibrary_readTypeString(void);
  uint32         gkLibrary_readCorrection(void)        { return(_readCorrection);       };
  char const    *gkLibrary_readCorrectionString(void);
  uint32         gkLibrary_finalTrim(void)             { return(_finalTrim);            };
  char const    *gkLibrary_finalTrimString(void);
  uint32         gkLibrary_removeDuplicateReads(void)  { return(_removeDuplicateReads); };
  uint32         gkLibrary_removeSpurReads(void)       { return(_removeSpurReads);      };
  uint32         gkLibrary_removeChimericReads(void)   { return(_removeChimericReads);  };
  uint32         gkLibrary_checkForSubReads(void)      { return(_checkForSubReads);     };

  uint32         gkLibrary_defaultQV(void)             { return(_defaultQV);            };


  void           gkLibrary_setIsNonRandom(bool f)           { _isNonRandom = f;          };
  void           gkLibrary_setReadType(char *f);
  void           gkLibrary_setReadCorrection(char *t);
  void           gkLibrary_setFinalTrim(char *t);
  void           gkLibrary_setRemoveDuplicateReads(bool f)  { _removeDuplicateReads = f; };
  void           gkLibrary_setRemoveSpurReads(bool f)       { _removeSpurReads = f;      };
  void           gkLibrary_setRemoveChimericReads(bool f)   { _removeChimericReads = f;  };
  void           gkLibrary_setCheckForSubReads(bool f)      { _checkForSubReads = f;     };

  void           gkLibrary_setDefaultQV(double qv)          { _defaultQV = qv;           };

  void           gkLibrary_parsePreset(char *t);


private:
  char           _libraryName[LIBRARY_NAME_SIZE];
  uint32         _libraryID;

  //  If set, reads are from a non-random library and shouldn't count toward coverage stats.
  uint32         _isNonRandom;

  //  What generated these reads?
  uint32         _readType;

  //  Should reads be corrected?  How?
  uint32         _readCorrection;

  //  Should reads be trimmed based on overlaps?  How?
  uint32         _finalTrim;

  //  Should duplicate reads (based on overlaps) be removed?
  uint32         _removeDuplicateReads;

  //  Should spur reads be cleaned up?  How?
  uint32         _removeSpurReads;

  //  Should chimeric reads be cleaned up?  How?
  uint32         _removeChimericReads;

  //  Should PacBio circular sub-reads be cleaned up?  How?
  uint32         _checkForSubReads;

  //  For reads with no QVs, use this.
  uint32         _defaultQV;

  friend class gkStore;
};



class gkRead;

class gkReadData {
public:
  gkReadData() {
    _read = NULL;

    _seq       = NULL;
    _qlt       = NULL;
    _seqAlloc  = 0;

    _blobLen   = 0;
    _blobMax   = 0;
    _blob      = NULL;

    _pacbio._deletionQV     = NULL;
    _pacbio._deletionTag    = NULL;
    _pacbio._insertionQV    = NULL;
    _pacbio._substitutionQV = NULL;
    _pacbio._mergeQV        = NULL;

    _nanopore._someQV       = NULL;
  };

  ~gkReadData() {
    delete [] _seq;
    delete [] _qlt;

    delete [] _blob;

    delete [] _pacbio._deletionQV;
    delete [] _pacbio._deletionTag;
    delete [] _pacbio._insertionQV;
    delete [] _pacbio._substitutionQV;
    delete [] _pacbio._mergeQV;

    delete [] _nanopore._someQV;
  };

  gkRead  *gkReadData_getRead(void)         { return(_read); };

  bool     gkReadData_hasSequence(void)     { return(_seq != NULL); };
  bool     gkReadData_hasQualities(void)    { return(_qlt != NULL); };

  char    *gkReadData_getSequence(void)     { return(_seq);  };
  char    *gkReadData_getQualities(void)    { return(_qlt);  };

private:
  gkRead            *_read;     //  Pointer to the mmap'd read

  char              *_seq;      //  Everyone has sequence
  char              *_qlt;      //    and quality
  uint32             _seqAlloc;

  uint32             _blobLen;
  uint32             _blobMax;
  uint8             *_blob;     //  And maybe even an encoded blob of data from the store.

  //  Used by the store for adding a read.

  void     gkReadData_encodeBlobChunk(char const *tag, uint32 len, void *dat);

  //  Some reads have more data than others.  This used to be a union, but that
  //  only complicated things with no clear beneift.  Yay, we saved a few bytes on
  //  a structure that we have only a few copies in core at a time.

  struct pacbio {
    char          *_deletionQV;
    char          *_deletionTag;
    char          *_insertionQV;
    char          *_substitutionQV;
    char          *_mergeQV;
  } _pacbio;

  struct nanopore {
    char          *_someQV;
  } _nanopore;


  friend class gkRead;
  friend class gkStore;
};




class gkRead {
public:
  gkRead() {
    _readID     = 0;
    _libraryID  = 0;
    _seqLen     = 0;

    _mPtr       = 0;
    _pID        = 0;
  };

  ~gkRead() {
  };

  uint32      gkRead_readID(void)             { return((uint32)_readID);           };

  uint32      gkRead_libraryID(void)          { return((uint32)_libraryID);        };

  uint32      gkRead_sequenceLength(void)     { return((uint32)_seqLen);           };

  uint64      gkRead_mPtr(void) { return(_mPtr); };  //  For debugging, in gatekeeperDumpMetatData
  uint64      gkRead_pID(void)  { return(_pID);  };

public:
  bool        gkRead_loadData(gkReadData *readData, void *blob);

private:
  uint32      gkRead_encode2bit(uint8  *&chunk, char *seq, uint32 seqLen);
  uint32      gkRead_encode3bit(uint8  *&chunk, char *seq, uint32 seqLen);
  uint32      gkRead_encode4bit(uint8  *&chunk, char *qlt, uint32 seqLen);
  uint32      gkRead_encode5bit(uint8  *&chunk, char *qlt, uint32 seqLen);

  bool        gkRead_decode2bit(uint8  *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
  bool        gkRead_decode3bit(uint8  *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
  bool        gkRead_decode4bit(uint8  *chunk, uint32 chunkLen, char *qlt, uint32 seqLen);
  bool        gkRead_decode5bit(uint8  *chunk, uint32 chunkLen, char *qlt, uint32 seqLen);

public:
  //  These are called by gatekeeperCreate.
  gkReadData *gkRead_encodeSeqQlt(char *H, char *S, char *Q, uint32 qv);
  //gkReadData *gkRead_encodePacBio(char *H, char *S, char *Q);  //  Not implemented.
  //gkReadData *gkRead_encodeMinION(char *H, char *S, char *Q);

private:
  char       *gkRead_encodeSequence(char *sequence, char *encoded);
  char       *gkRead_decodeSequence(char *encoded,  char *sequence);

  char       *gkRead_encodeQuality(char *sequence, char *encoded);
  char       *gkRead_decodeQuality(char *encoded,  char *sequence);

private:
  //  Used by the store to copy data to a partition
  void     gkRead_copyDataToPartition(void *blobs, FILE **partfiles, uint64 *partfileslen, uint32 partID);

private:

  uint64   _readID       : AS_MAX_READS_BITS;
  uint64   _libraryID    : AS_MAX_LIBRARIES_BITS;
  uint64   _seqLen       : AS_MAX_READLEN_BITS;

  uint64   _mPtr         : 48;      //  Pointer to blob of data in blob file, 0..256 TB
  uint64   _pID          : 16;      //  Partition file id, 0...65536

  friend class gkStore;
};





//  gkStoreInfo is saved on disk.
//  gkStore is the in memory structure used to access the data.
//
class gkStoreInfo {
public:
  gkStoreInfo() {
    gkMagic            = 0x504b473a756e6163llu;  //  canu:GKP
    gkVersion          = 0x0000000000000001llu;

    gkLibrarySize      = sizeof(gkLibrary);
    gkReadSize         = sizeof(gkRead);
    gkMaxLibrariesBits = AS_MAX_LIBRARIES_BITS;
    gkLibraryNameSize  = LIBRARY_NAME_SIZE;
    gkMaxReadBits      = AS_MAX_READS_BITS;
    gkMaxReadLenBits   = AS_MAX_READLEN_BITS;
    gkUNUSED           = 0;

    numLibraries       = 0;
    numReads           = 0;
  };
  ~gkStoreInfo() {
  };

  void      writeInfoAsText(FILE *F) {
    fprintf(F, "gkMagic            = 0x"F_X64"\n", gkMagic);
    fprintf(F, "gkVersion          = 0x"F_X64"\n", gkVersion);
    fprintf(F, "\n");
    fprintf(F, "gkLibrarySize      = "F_U32"\n", gkLibrarySize);
    fprintf(F, "gkReadSize         = "F_U32"\n", gkReadSize);
    fprintf(F, "gkMaxLibrariesBits = "F_U32"\n", gkMaxLibrariesBits);
    fprintf(F, "gkLibraryNameSize  = "F_U32"\n", gkLibraryNameSize);
    fprintf(F, "gkMaxReadBits      = "F_U32"\n", gkMaxReadBits);
    fprintf(F, "gkMaxReadLenBits   = "F_U32"\n", gkMaxReadLenBits);
    fprintf(F, "\n");
    fprintf(F, "numLibraries       = "F_U32"\n", numLibraries);
    fprintf(F, "numReads           = "F_U32"\n", numReads);
  };

private:
  uint64    gkMagic;
  uint64    gkVersion;

  uint32    gkLibrarySize;      //  Sanity checks that this code can load the data properly.
  uint32    gkReadSize;
  uint32    gkMaxLibrariesBits;
  uint32    gkLibraryNameSize;
  uint32    gkMaxReadBits;
  uint32    gkMaxReadLenBits;
  uint32    gkUNUSED;           //  Used to hold a blob block size that was never implemented

  uint32    numLibraries;       //  Counts of types of things we have loaded
  uint32    numReads;

  friend class gkStore;
};



//  The default behavior is to open the store for read only, and to load
//  all the metadata into memory.

typedef enum {
  gkStore_readOnly    = 0x00,  //  Open read only
  gkStore_modify      = 0x01,  //  Open for modification
  gkStore_extend      = 0x02,  //  Open for modification and appending new reads/libraries
  gkStore_infoOnly    = 0x03   //  Open read only, but only load the info on the store; no access to reads or libraries
} gkStore_mode;


static
const
char *
toString(gkStore_mode m) {
  switch (m) {
    case gkStore_readOnly:     return("gkStore_readOnly");     break;
    case gkStore_modify:       return("gkStore_modify");       break;
    case gkStore_extend:       return("gkStore_extend");       break;
    case gkStore_infoOnly:     return("gkStore_infoOnly");     break;
  }

  return("undefined-mode");
}



class gkStore {

private:
  gkStore(char const *path, gkStore_mode mode, uint32 partID);
  ~gkStore();

public:
  static
  gkStore     *gkStore_open(char const *path, gkStore_mode mode=gkStore_readOnly, uint32 partID=UINT32_MAX) {

#pragma omp critical
    {
      //  If an instance exists, return it.

      if (_instance != NULL) {
        _instanceCount++;
        //fprintf(stderr, "gkStore_open(%s), %u instances now\n", path, _instanceCount);
      }

      //  Otherwise, make a new one.

      else {
        _instance      = new gkStore(path, mode, partID);
        _instanceCount = 1;
        //fprintf(stderr, "gkStore_open(%s), first instance, create store\n", path);
      }
    }

    return(_instance);
  };


  void         gkStore_close(void) {

#pragma omp critical
    {
      _instanceCount--;

      if (_instanceCount == 0) {
        delete _instance;
        _instance = NULL;
        //fprintf(stderr, "gkStore_close(%s), no instances remain, delete store\n", _storeName, _instanceCount);
      }

      else {
        //fprintf(stderr, "gkStore_close(%s), %u instances remain\n", _storeName, _instanceCount);
      }
    }

  };


public:
  const char  *gkStore_path(void) { return(_storePath); };  //  Returns the path to the store
  const char  *gkStore_name(void) { return(_storeName); };  //  Returns the name, e.g., name.gkpStore

  void         gkStore_buildPartitions(uint32 *partitionMap);

  void         gkStore_delete(void);             //  Deletes the files in the store.
  void         gkStore_deletePartitions(void);   //  Deletes the files for a partition.

  uint32       gkStore_getNumLibraries(void)       { return(_info.numLibraries); };
  uint32       gkStore_getNumReads(void)           { return(_info.numReads); };

  gkLibrary   *gkStore_getLibrary(uint32 id)       { return(&_libraries[id]); };

  //  Returns a read, using the copy in the partition if the partition exists.
  gkRead      *gkStore_getRead(uint32 id)          {
    if ((_readIDtoPartitionID)   &&
        (_readIDtoPartitionID[id] != _partitionID)) {
      fprintf(stderr, "getRead()--  WARNING: access to read %u in partition %u is slow when partition %u is loaded.\n",
              id, _readIDtoPartitionID[id], _partitionID);
      assert(0);
    }

    if ((_readIDtoPartitionID)   &&
        (_readIDtoPartitionID[id] == _partitionID)) {
      //fprintf(stderr, "getRead()-- id=%u mapped=%u mappedid=%u\n",
      //        id, _readIDtoPartitionIdx[id], _reads[_readIDtoPartitionIdx[id]].gkRead_readID());
      return(_reads + _readIDtoPartitionIdx[id]);
    }

    return(_reads + id);
  };

  //  Returns a read, but only if it is in the currently loaded partition.
  gkRead      *gkStore_getReadInPartition(uint32 id) {

    if (_readIDtoPartitionID == NULL)              //  Not partitioned, return regular read.
      return(gkStore_getRead(id));

    if (_readIDtoPartitionID[id] != _partitionID)  //  Patitioned, and not in this partition.
      return(NULL);

    //fprintf(stderr, "getRead()--  SUCCESS: access to read %u in partition %u is FAST when partition %u is loaded.\n",
    //        id, _readIDtoPartitionID[id], _partitionID);

    return(_reads + _readIDtoPartitionIdx[id]);
  }

  gkLibrary   *gkStore_addEmptyLibrary(char const *name);
  gkRead      *gkStore_addEmptyRead(gkLibrary *lib);

  bool         gkStore_loadReadData(gkRead *read,   gkReadData *readData) {
    return(read->gkRead_loadData(readData, _blobs));
  };
  bool         gkStore_loadReadData(uint32  readID, gkReadData *readData) {
    return(gkStore_getRead(readID)->gkRead_loadData(readData, _blobs));
  };

  void         gkStore_stashReadData(gkRead *read, gkReadData *data);

  static
  void         gkStore_loadReadFromStream(FILE *S, gkRead *read, gkReadData *readData);
  void         gkStore_saveReadToStream(FILE *S, uint32 id);

private:
  static gkStore      *_instance;
  static uint32        _instanceCount;

  gkStoreInfo          _info;  //  All the stuff stored on disk.

  char                 _storePath[FILENAME_MAX];  //  Needed to create files
  char                 _storeName[FILENAME_MAX];  //  Useful for log files in other programs

  gkStore_mode         _mode;  //  What mode this store is opened as, sanity checking

  //  If these are memory mapped, then multiple processes on the same host can share the
  //  (read-only) data.

  memoryMappedFile    *_librariesMMap;
  uint32               _librariesAlloc;  //  If zero, the mmap is used.
  gkLibrary           *_libraries;

  memoryMappedFile    *_readsMMap;
  uint32               _readsAlloc;      //  If zero, the mmap is used.
  gkRead              *_reads;

  memoryMappedFile    *_blobsMMap;       //  Either the full blobs, or the partitioned blobs.
  void                *_blobs;
  FILE                *_blobsFile;

  //  If the store is openend partitioned, this data is loaded from disk

  uint32               _numberOfPartitions;     //  Total number of partitions that exist
  uint32               _partitionID;            //  Which partition this is
  uint32              *_readsPerPartition;      //  Number of reads in each partition, mostly sanity checking
  uint32              *_readIDtoPartitionIdx;   //  Map from global ID to local partition index
  uint32              *_readIDtoPartitionID;    //  Map from global ID to partition ID
};




class gkStoreStats {
public:
  gkStoreStats(char const *gkStoreName) {
    gkStore *gkp = gkStore::gkStore_open(gkStoreName);
    init(gkp);
    gkp->gkStore_close();
  };

  gkStoreStats(gkStore *gkp) {
    init(gkp);
  };

  ~gkStoreStats() {
    delete [] lowestID;
    delete [] highestID;

    delete [] numActivePerLib;

    delete [] readLengthPerLib;
    delete [] clearLengthPerLib;
  };

  void    init(gkStore *gkp);

  //  Global stats over the whole store

  uint32  numActiveFrag;
  uint64  readLength;
  uint64  clearLength;

  //  Per library stats

  uint32  *lowestID;
  uint32  *highestID;

  uint32  *numActivePerLib;
  uint64  *readLengthPerLib;
  uint64  *clearLengthPerLib;
};


#endif
