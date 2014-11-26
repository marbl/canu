
#ifndef GKSTORE_H
#define GKSTORE_H

#include "AS_global.H"
#include "memoryMappedFile.H"

//  The number of library IIDs we can handle.
//
#define AS_MAX_LIBRARIES_BITS      8
#define AS_MAX_LIBRARIES           ((1 << AS_READ_MAX_LIBRARIES_BITS) - 1)

#define LIBRARY_NAME_SIZE          128
#define LIBRARY_NAME_DEFAULT      "UNDEFINED"

//  The number of read IDs we can handle.  This has more implications for
//  ovlStore than gkpStore.  The default of 27 = 134.2 million reads.
//
#define AS_MAX_READS_BITS          27
#define AS_MAX_READS               ((1 << AS_READ_MAX_READS_BITS) - 1)

//  Be aware that overlapper does not work well at large sizes.  Consensus allocates excessive
//  amounts of memory too.
//
#define AS_MAX_READLEN_BITS        16
#define AS_MAX_READLEN             ((1 << AS_READ_MAX_LEN_BITS) - 1)

//  This indirectly controls how much sequence can be loaded, in that the store
//  can address 2^32 blocks of this many bytes.  QVs are stored in 5/8 of a byte,
//  so with the default of 64, the store can hold at most
//    2^32 * 64 * 5 / 8 = 171.8 Gbp.
//  With the default 134.2 million reads, that's only 1280 bp/read.
//
#define DATA_BLOCK_SIZE            1024

//  Per-library algorithm options.
//
#define CORRECTION_NONE            0x00
#define CORRECTION_CONSENSUS       0x01
#define CORRECTION_MER             0x02

#define INITIALTRIM_NONE           0x00
#define INITIALTRIM_MER_BASED      0x01
#define INITIALTRIM_FLOW_BASED     0x02
#define INITIALTRIM_QUALITY_BASED  0x03

#define FINALTRIM_NONE             0x00
#define FINALTRIM_LARGEST_COVERED  0x01
#define FINALTRIM_EVIDENCE_BASED   0x02
#define FINALTRIM_BEST_EDGE        0x03

//  Sequence data storage methods, one per technology supported.
//
#define GKREAD_TECH_UNKNOWN   0x00
#define GKREAD_TECH_FASTA     0x01
#define GKREAD_TECH_FASTQ     0x02
#define GKREAD_TECH_PACBIO    0x03
#define GKREAD_TECH_NANOPORE  0x04
#define GKREAD_TECH_UNUSED1   0x05
#define GKREAD_TECH_UNUSED2   0x06
#define GKREAD_TECH_UNUSED3   0x07
#define GKREAD_TECH_NUM       8


class gkLibrary {
public:
  gkLibrary() {
    memset(_libraryName, 0, sizeof(char) * LIBRARY_NAME_SIZE);
    strcpy(_libraryName, LIBRARY_NAME_DEFAULT);

    gkLibrary_clearFeatures();
  };
  ~gkLibrary() {
  };

  void gkLibrary_clearFeatures(void) {
    _isNotRandom                = 0;

    _trustHomopolymerRuns       = 1;
    _doCorrection               = CORRECTION_NONE;

    _initialTrim                = INITIALTRIM_NONE;
    _finalTrim                  = FINALTRIM_NONE;
    _removeDuplicateReads       = 0;
    _removeSpurReads            = 0;
    _removeChimericReads        = 0;
    _checkForSubReads           = 0;
  };

public:
  char const    *gkLibrary_libraryName(void) {
    return(_libraryName);
  };

  uint32         gkLibrary_isNotRandom(void) {
    return(_isNotRandom);
  };

  uint32         gkLibrary_trustHomopolymerRuns(void) {
    return(_trustHomopolymerRuns);
  };

  uint32         gkLibrary_correctBases(void) {
    return(_doCorrection);
  };

  uint32         gkLibrary_initialTrim(void) {
    return(_initialTrim);
  };

  uint32         gkLibrary_finalTrim(void) {
    return(_finalTrim);
  };

  uint32         gkLibrary_removeDuplicateReads(void) {
    return(_removeDuplicateReads);
  };

  uint32         gkLibrary_removeSpurReads(void) {
    return(_removeSpurReads);
  };

  uint32         gkLibrary_removeChimericReads(void) {
    return(_removeChimericReads);
  };

  uint32         gkLibrary_checkForSubReads(void) {
    return(_checkForSubReads);
  };

private:
  char           _libraryName[LIBRARY_NAME_SIZE];

  //  If set, reads are from a non-random library and shouldn't count toward coverage stats.
  uint32         _isNotRandom;

  //  If set, reads have issues with homopolymer runs.
  uint32         _trustHomopolymerRuns;

  //  Should reads be corrected?  How?
  uint32         _doCorrection;

  //  Should reads be trimmed based on read data alone?  How?
  uint32         _initialTrim;

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
};



class gkRead;

class gkReadData {
public:
  gkReadData() {
    _read = NULL;

    _technology = GKREAD_TECH_UNKNOWN;

    _seq = NULL;
    _qlt = NULL;
  };
  ~gkReadData() {
    delete [] _seq;
    delete [] _qlt;

    if ((_technology == GKREAD_TECH_PACBIO)) {
      delete [] _data._pacbio._deletionQV;
      delete [] _data._pacbio._deletionTag;
      delete [] _data._pacbio._insertionQV;
      delete [] _data._pacbio._substitutionQV;
      delete [] _data._pacbio._mergeQV;
    }
    if ((_technology == GKREAD_TECH_NANOPORE)) {
      delete [] _data._nanopore._someQV;
    }
  };

  uint32   technology(void)   { return(_technology); };

  bool     hasSequence(void)  { return(_seq != NULL); };
  bool     hasQualities(void) { return(_qlt != NULL); };

  char    *sequence(void)     { return(_seq);  };
  char    *qualities(void)    { return(_qlt);  };

  //void    initialize(uint32 seqLen) {
  //  _seq = new char [seqLen + 1];
  //  _qlt = new char [seqLen + 1];
  //};

private:
  gkRead            *_read;     //  Pointer to the mmap'd read

  uint32             _technology;

  char              *_seq;      //  Everyone has sequence
  char              *_qlt;      //    and quality

  union {
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
  } _data;

  friend class gkRead;
};




class gkRead {
public:
  gkRead() {
    _readID     = 0;
    _libraryID  = 0;
    _seqLen     = 0;
    _deleted    = true;
    _technology = GKREAD_TECH_UNKNOWN;
    _flagpad    = 0;

    _clrBgn     = 0;
    _clrEnd     = 0;
    _sizepad    = 0;

    _mPtr       = 0;
    _pID        = 0;
    _pPtr       = 0;
  };

  ~gkRead() {
  };

  uint32    gkRead_readID(void) {
    return(_readID);
  };

  uint32    gkRead_libraryID(void) {
    return(_libraryID);
  };

  uint32    gkRead_sequenceLength(void) {
    return(_seqLen);
  };

  bool      gkRead_isDeleted(void) {
    return(_deleted);
  };
  void      gkRead_isDeleted(bool deleted) {
    _deleted = deleted;
  };

  uint32    gkRead_clearRegionBegin(void) {
    return(_clrBgn);
  };
  uint32    gkRead_clearRegionEnd(void) {
    return(_clrEnd);
  };
  uint32    gkRead_clearRegionLength(void) {
    return(_clrEnd - _clrBgn);
  };

  void      gkRead_setClearRegion(uint32 bgn, uint32 end) {
    _clrBgn = bgn;
    _clrEnd = end;
  };

public:
  bool      gkRead_loadData(gkReadData *readData, void *blob, bool partitioned);

private:
  bool      gkRead_loadFASTAData   (gkReadData *readData, void *blob);
  bool      gkRead_loadFASTQData   (gkReadData *readData, void *blob);
  bool      gkRead_loadPacBioData  (gkReadData *readData, void *blob);
  bool      gkRead_loadNanoporeData(gkReadData *readData, void *blob);


private:
  char     *gkRead_encodeSequence(char *sequence, char *encoded);
  char     *gkRead_decodeSequence(char *encoded,  char *sequence);
  
  char     *gkRead_encodeQuality(char *sequence, char *encoded);
  char     *gkRead_decodeQuality(char *encoded,  char *sequence);

private:
  uint64   _readID       : AS_MAX_READS_BITS;       //  Up to 128 million reads
  uint64   _libraryID    : AS_MAX_LIBRARIES_BITS;   //  Up to 256 libraries
  uint64   _seqLen       : AS_MAX_READLEN_BITS;     //  Up to 128 kilobases
  uint64   _deleted      :  1;                      //  Deleted?
  uint64   _technology   :  3;                      //  Type of sequence data this represents
  uint64   _flagpad      : (64 -
                            AS_MAX_READS_BITS -
                            AS_MAX_LIBRARIES_BITS -
                            AS_MAX_READLEN_BITS - 1 - 3);

  uint64   _clrBgn       : AS_MAX_READLEN_BITS;
  uint64   _clrEnd       : AS_MAX_READLEN_BITS;
  uint64   _sizepad      : 64 - 2 * AS_MAX_READLEN_BITS;

  uint64   _mPtr         : 32;      //  Pointer to blob of data in master blob file (4e9)
  uint64   _pID          : 10;      //  Partition file id, 0...1024
  uint64   _pPtr         : 22;      //  Pointer to blob of data in partitioned blob file (4e6)

  friend class gkStore;
};





//  gkStoreInfo is saved on disk.
//  gkStore is the in memory structure used to access the data.
//
class gkStoreInfo {
public:
  gkStoreInfo() {
  };
  ~gkStoreInfo() {
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
  uint32    gkOffsetBlockSize;

  uint32    numLibraries;       //  Counts of types of things we have loaded
  uint32    numReads;

  friend class gkStore;
};



//  The default behavior is to open the store for read only, and to load
//  all the metadata into memory.

const bool gkStore_readOnly = true;
const bool gkStore_writable = false;

const bool gkStore_creatable = false;

class gkStore {
public:

  //  Open a full gkpStore.
  gkStore(char const *path, bool readOnly=true, bool createOnly=true);

  //  Open a partitioned store, only reading is supported.
  gkStore(char const *path, uint32 partID);

  ~gkStore();

  const char  *gkStore_path(void) { return(_storePath); };  //  Returns the path to the store
  const char  *gkStore_name(void) { return(_storeName); };  //  Returns the name, e.g., name.gkpStore

  bool         gkStore_loadPartition(uint32 partID);

  void         gkStore_buildPartitions(uint32 *partitionMap, uint32 maxPartition);

  void         gkStore_delete(void);             //  Deletes the files in the store.
  void         gkStore_deletePartitions(void);   //  Deletes the files for a partition.

  uint32       gkStore_getNumLibraries(void)       { return(_inf.numLibraries); }
  uint32       gkStore_getNumReads(void)           { return(_inf.numReads); }

  gkLibrary   *gkStore_getLibrary(uint32 id)       { return(NULL); }
  gkRead      *gkStore_getRead(uint32 id)          { return(NULL); }

  //void         gkStore_addLibrary(gkLibrary *lb) {
  //  inf.libLoaded++;
  //};
  //void         gkStore_addFragment(gkReadData *fr) {
  //  inf.frgLoaded++;
  //};

  bool         gkStore_loadData(gkRead *read,   gkReadData *readData) {
    return(gkStore_loadData(read->gkRead_readID(), readData));
  };
  bool         gkStore_loadData(uint32  readID, gkReadData *readData);

private:
  gkStoreInfo          _inf;  //  All the stuff stored on disk.

  char                 _storePath[FILENAME_MAX];  //  Needed?
  char                 _storeName[FILENAME_MAX];  //  Needed?

  //  If these are memory mapped, then multiple processes on the same host can share the
  //  (read-only) data.

  memoryMappedFile    *_librariesFile;
  gkLibrary           *_libraries;

  memoryMappedFile    *_readsFile;
  gkRead              *_reads;

  memoryMappedFile    *_blobsFile;  //  Either the full blobs, or the partitioned blobs.
  void                *_blobs;

  //  If the store is openend partitioned, this data is loaded from disk

  uint32               _numberOfPartitions;   //  Total number of partitions that exist
  uint32               _partitionID;          //  Which partition this is
  uint32              *_partitionIDmap;       //  Map from global ID to partition ID
  uint32              *_readsPerPartition;    //  Number of reads per partition, probably not useful
  uint32              *_readsInThisPartition; //  A list of the (global ID) reads in this partition
};




class gkStoreStats {
public:
  gkStoreStats(char const *gkStoreName) {
    gkStore *gkp = new gkStore(gkStoreName);
    init(gkp);
    delete gkp;
  };

  gkStoreStats(gkStore *gkp) {
    init(gkp);
  };

  ~gkStoreStats() {
    delete [] lowestID;
    delete [] highestID;

    delete [] numActivePerLib;
    delete [] numDeletedPerLib;

    for (uint32 i=0; i<GKREAD_TECH_NUM; i++)
      delete [] numByTypePerLib[i];

    delete [] readLengthPerLib;
    delete [] clearLengthPerLib;
  };

  void    init(gkStore *gkp);

  //  Global stats over the whole store

  uint32  numActiveFrag;
  uint32  numDeletedFrag;
  uint32  numByType[GKREAD_TECH_NUM];
  uint64  readLength;
  uint64  clearLength;

  //  Per library stats

  uint32  *lowestID;
  uint32  *highestID;

  uint32  *numActivePerLib;
  uint32  *numDeletedPerLib;
  uint32 **numByTypePerLib;
  uint64  *readLengthPerLib;
  uint64  *clearLengthPerLib;
};







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



#endif
