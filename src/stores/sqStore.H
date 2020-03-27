
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef SQSTORE_H
#define SQSTORE_H

#include "runtime.H"
#include "files.H"

#include <vector>

using namespace std;

//  Versions:
//
//   4 - first to be called seqStore.               18 APR 2018.  01d182304a951846677a641197bb761674a5e662
//   5 - tracking of the last blob used.            01 MAY 2018.  1e1ca668fbe16763ef22d7d69a04413cc2b95fb1
//   6 - unlimited blobs, segs and parts.           04 MAY 2018.  5654d275dfe096465eedd36d0dc316a6c79c59ec
//   7 - skipped, used during development.
//   8 - allow clear/compr for any sequence type.   05 JUL 2019.  40763d8b21bfc7f65530c8d768abffbc64e777c3
//   9 - add, again, numBlobs.                      18 OCT 2019.  (unknown)

#define SQ_MAGIC8  0x504b473a756e6162lu      //  banu:GKP (for versions before 9)
#define SQ_MAGIC   0x5145533a756e6163lu      //  canu:SEQ
#define SQ_VERSION 0x0000000000000009lu


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


//  The default behavior is to open the store for read only, and to load
//  all the metadata into memory.

typedef enum {
  sqStore_create      = 0x00,  //  Open for creating, will fail if files exist already
  sqStore_extend      = 0x01,  //  Open for modification and appending new reads/libraries
  sqStore_readOnly    = 0x02,  //  Open read only
} sqStore_mode;


static
const
char *
toString(sqStore_mode m) {
  switch (m) {
    case sqStore_create:       return("sqStore_create");       break;
    case sqStore_extend:       return("sqStore_extend");       break;
    case sqStore_readOnly:     return("sqStore_readOnly");     break;
  }

  return("undefined-mode");
}


static
void
makeBlobName(char *storePath, uint32 blobNumber, char *blobName) {
  snprintf(blobName, FILENAME_MAX, "%s/blobs.%04" F_U32P , storePath, blobNumber);
}


class sqStoreInfo {
public:
  sqStoreInfo();
  ~sqStoreInfo();

private:
  bool      readInfo8(char *storePath);
public:
  void      readInfo(char *storePath);
  void      writeInfo(char *storePath);

  void      writeInfoAsText(FILE *F);

  uint32    sqInfo_lastLibraryID(void)          { return(_numLibraries); };
  uint32    sqInfo_lastReadID(void)             { return(_numReads);     };

  uint32    sqInfo_numReads(sqRead_which w)     { return(_reads[w]);     };

public:
  bool      examineRead(uint32 ii, sqReadSeq *seq, sqRead_which w);

private:
  bool      checkInfo(void);
  void      update(sqReadSeq  *rawU,
                   sqReadSeq  *rawC,
                   sqReadSeq  *corU,
                   sqReadSeq  *corC);

public:
  //  For info in the store, these control the number of objects
  //  that are allocated in the metadata arrays.
  //
  //  For info in dumpMetaData, these are used to count the number
  //  of objects seen in the scan.
  //
  void      sqInfo_addLibrary(void)             { _numLibraries++;       };
  void      sqInfo_addRead(void)                { _numReads++;           };

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
  uint32    _numReads;           //    available index into _libraries and _reads in sqStore)
  uint32    _numBlobs;

  uint64    _reads[sqRead_largest];
  uint64    _bases[sqRead_largest];

  friend class sqStore;
  friend class sqStoreBlobWriter;
};





class sqStoreBlobWriter {
public:
  sqStoreBlobWriter(const char *storePath, sqStoreInfo *info);
  ~sqStoreBlobWriter();

  void           writeData(sqReadDataWriter *readData);

private:
  char          _storePath[FILENAME_MAX+1];        //  Path to the seqStore.
  char          _blobName[FILENAME_MAX+1];         //  A temporary to make life easier.

  sqStoreInfo  *_info;
  writeBuffer  *_buffer;
};



//  Manages access to blob data.  You need one of these per thread.
//
class sqStoreBlobReader {
public:
  sqStoreBlobReader(const char *storePath);
  ~sqStoreBlobReader();

  readBuffer    *getBuffer(sqReadMeta *meta);
  readBuffer    *getBuffer(sqReadMeta &meta)   { return(getBuffer(&meta)); };

private:
  char          _storePath[FILENAME_MAX+1];        //  Path to the seqStore.
  char          _blobName[FILENAME_MAX+1];         //  A temporary to make life easier.

  uint32        _buffersMax;
  readBuffer  **_buffers;   //  One per blob file.
};





class sqStore {
public:
  sqStore(char const *storePath, sqStore_mode mode=sqStore_readOnly);
  ~sqStore();
private:
  void         sqStore_loadMetadata(void);

public:
  const char  *sqStore_path(void) { return(_storePath); };  //  Returns the path to the store

  void         sqStore_delete(void);             //  Deletes the files in the store.

  uint32       sqStore_lastLibraryID(void)            { return(_info.sqInfo_lastLibraryID()); };
  uint32       sqStore_lastReadID(void)               { return(_info.sqInfo_lastReadID());    };

  //  If iterating over reads, use sqStore_lastReadID(), not sqStore_getNumReads().
  uint32       sqStore_getNumReads(sqRead_which w)    { return(_info.sqInfo_numReads(w));     };

  sqLibrary   *sqStore_getLibrary(uint32 id)          { return(&_libraries[id]); };

  uint32       sqStore_getLibraryIDForRead(uint32 id) { return(_meta[id].sqRead_libraryID()); };
  sqLibrary   *sqStore_getLibraryForRead(uint32 id)   { return(&_libraries[_meta[id].sqRead_libraryID()]); };

public:
  readBuffer  *sqStore_getReadBuffer(uint32 readID);
  sqRead      *sqStore_getRead(uint32 readID, sqRead *read);

public:
  static
  bool         sqStore_loadReadFromBuffer(readBuffer *B, sqRead *read);
  void         sqStore_saveReadToBuffer(writeBuffer *B, uint32 id, sqRead *rd, sqReadDataWriter *wr);

private:
  sqReadSeq   *sqStore_getSeq(sqRead_which w) {
    bool  isRaw = ((w & sqRead_raw)        == sqRead_raw);
    bool  isCor = ((w & sqRead_corrected)  == sqRead_corrected);
    bool  isTrm = ((w & sqRead_trimmed)    == sqRead_trimmed);
    bool  isCmp = ((w & sqRead_compressed) == sqRead_compressed);

    if ((isRaw == true) && (isCmp == false))          { return(_rawU); };
    if ((isRaw == true) && (isCmp == true))           { return(_rawC); };

    if ((isCor == true) && (isCmp == false))          { return(_corU); };
    if ((isCor == true) && (isCmp == true))           { return(_corC); };

    fprintf(stderr, "sqStore_getSeq()-- Unknown which '%s'\n", toString(w));

    assert(0);
    return(NULL);
  };

  //  Accessors to read data.
public:
  uint64             sqStore_getReadSegm(uint32 id);
  uint64             sqStore_getReadByte(uint32 id);

  sqReadSeq         *sqStore_getReadSeq(uint32 id, sqRead_which w=sqRead_defaultVersion);

  //  Read length is zero if the read isn't present, isn't valid, is ignored,
  //  or is trimmed out.  Ideally, those cases should generate an error -- we
  //  should be testing for those cases before deciding to use a read -- but
  //  too many places are expecting length==0 in thse cases.
  //
  uint32             sqStore_getReadLength(uint32 id, sqRead_which w=sqRead_defaultVersion);

  //  Unlike read length, we'll insist that accessing clear ranges only occur
  //  for trimmed reads (enforced by sqReadSeq).
  //
  uint32             sqStore_getClearBgn  (uint32 id, sqRead_which w=sqRead_defaultVersion);
  uint32             sqStore_getClearEnd  (uint32 id, sqRead_which w=sqRead_defaultVersion);

  bool               sqStore_isValidRead(uint32 id, sqRead_which w=sqRead_defaultVersion);
  bool               sqStore_isIgnoredRead(uint32 id, sqRead_which w=sqRead_defaultVersion);
  bool               sqStore_isTrimmedRead(uint32 id, sqRead_which w=sqRead_defaultVersion);

  //  For use ONLY by sqStoreCreate, to add new libraries and reads to a
  //  store.  The first two allocate a new metadata object in the store,
  //  while the last loads read sequence data.
  //
public:
  sqLibrary         *sqStore_addEmptyLibrary(char const *name, sqLibrary_tech techType);
  sqReadDataWriter  *sqStore_addEmptyRead(sqLibrary *lib, const char *name);

  void               sqStore_addRead(sqReadDataWriter *rdw) {
    _blobWriter->writeData(rdw);
  };

  //  Used when initially loading reads into seqStore, and when loading
  //  trimmed reads.  It sets the ignore flag in both the normal and
  //  compressed metadata for a given read.  Select 'raw' or 'corrected' with
  //  aqRead_which, then select untrimmed or trimmed with the flags.
public:
  void               sqStore_setIgnored(uint32       id,
                                        bool         untrimmed,
                                        bool         trimmed,
                                        sqRead_which w=sqRead_defaultVersion);

  //  This loads a read, computes the compressed sequence, and sets both the
  //  normal and compressed clear ranges.
public:
  void               sqStore_setClearRange(uint32 id,
                                           uint32 bgn, uint32 end, bool bogus,
                                           sqRead_which w=sqRead_defaultVersion);

private:
  sqStoreInfo          _info;  //  All the stuff stored on disk.

  char                 _storePath[FILENAME_MAX+1];    //  Needed to create files

  sqStore_mode         _mode;  //  What mode this store is opened as, sanity checking

  uint32               _librariesAlloc;  //  Size of allocation
  sqLibrary           *_libraries;       //  In core data

  uint32               _readsAlloc;      //  Size of allocation
  sqReadMeta          *_meta;            //  In core data
  sqReadSeq           *_rawU;            //  Metadata for raw sequence
  sqReadSeq           *_rawC;            //  Metadata for raw compressed sequence
  sqReadSeq           *_corU;
  sqReadSeq           *_corC;

  sqStoreBlobReader   *_blobReader;
  sqStoreBlobWriter   *_blobWriter;
};




//  The rest are just accessors to reads.


inline
uint64
sqStore::sqStore_getReadSegm(uint32 id) {
  return(_meta[id].sqRead_mSegm());
}


inline
uint64
sqStore::sqStore_getReadByte(uint32 id) {
  return(_meta[id].sqRead_mByte());
}


inline
sqReadSeq *
sqStore::sqStore_getReadSeq(uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  return((seq == NULL) ? NULL : &seq[id]);
}


//  Read length is zero if the read isn't present, isn't valid, is ignored,
//  or is trimmed out.  Ideally, those cases should generate an error -- we
//  should be testing for those cases before deciding to use a read -- but
//  too many places are expecting length==0 in thse cases.
//
inline
uint32
sqStore::sqStore_getReadLength(uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  if (w & sqRead_trimmed)
    return(((seq == NULL) ||
            (seq[id].sqReadSeq_trimmed() == false) ||
            (seq[id].sqReadSeq_valid()   == false) ||
            (seq[id].sqReadSeq_ignoreT() == true)) ? 0 : (seq[id].sqReadSeq_clearEnd() - seq[id].sqReadSeq_clearBgn()));

  else
    return(((seq == NULL) ||
            (seq[id].sqReadSeq_valid()   == false) ||
            (seq[id].sqReadSeq_ignoreU() == true)) ? 0 : (seq[id].sqReadSeq_length()));
}


//  Unlike read length, we'll insist that accessing clear ranges only occur
//  for trimmed reads (enforced by sqReadSeq).
//
inline
uint32
sqStore::sqStore_getClearBgn  (uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  if (seq == NULL)
    return(0);

  if (w & sqRead_trimmed)
    return(seq[id].sqReadSeq_clearBgn());
  else
    return(0);
}


inline
uint32
sqStore::sqStore_getClearEnd  (uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  if (seq == NULL)
    return(0);

  if (w & sqRead_trimmed)
    return(seq[id].sqReadSeq_clearEnd());
  else
    return(seq[id].sqReadSeq_length());
}


inline
bool
sqStore::sqStore_isValidRead(uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  if (seq == NULL)
    return(false);

  if (w & sqRead_trimmed)
    return(seq[id].sqReadSeq_valid() && seq[id].sqReadSeq_trimmed());
  else
    return(seq[id].sqReadSeq_valid());
}


inline
bool
sqStore::sqStore_isIgnoredRead(uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  if (seq == NULL)
    return(false);

  if (sqStore_isValidRead(id, w) == false)
    return(true);

  if (w & sqRead_trimmed)
    return(seq[id].sqReadSeq_ignoreT());
  else
    return(seq[id].sqReadSeq_ignoreU());
}


inline
bool
sqStore::sqStore_isTrimmedRead(uint32 id, sqRead_which w) {
  sqReadSeq  *seq = sqStore_getSeq(w);

  assert(id > 0);

  return((seq == NULL) ? false : seq[id].sqReadSeq_trimmed());
}







#endif  //  SQSTORE_H
