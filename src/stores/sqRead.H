
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
 *    src/stores/gkRead.H
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

#ifndef SQREAD_H
#define SQREAD_H


//  DO NOT INCLUDE THIS FILE DIRECTLY, include sqStore.H.


typedef enum {
  sqRead_latest    = 0x00,  //  Return the latest version, default mode for most functions
  sqRead_raw       = 0x01,  //  Return the raw read
  sqRead_corrected = 0x02,  //  Return the corrected read
  sqRead_trimmed   = 0x03,  //  Return the trimmed read
} sqRead_version;

extern
sqRead_version  sqRead_defaultVersion;

inline
void
sqRead_setDefaultVersion(sqRead_version v) {
  sqRead_defaultVersion = v;
};

inline
const char *
toString(sqRead_version v) {
  switch (v) {
    case sqRead_latest:     return("latest");     break;
    case sqRead_raw:        return("raw");        break;
    case sqRead_corrected:  return("corrected");  break;
    case sqRead_trimmed:    return("trimmed");    break;
  }

  return("unknown");
}



class sqRead;
class sqLibrary;
class sqCache;

class sqReadData {
public:
  sqReadData() {
    //  Set in sqStore_loadReadData()         (when reading from a store)
    //         sqStore_loadReadFromStream()   (when reading from a file, utgcns import)
    //         sqStore_addEmptyRead()         (when making a new read)
    //  The first two then call sqRead_loadDataFromStream() -> sqReadData_loadFromBlob()
    //
    //  Possibly set in stashReadData() and loadFromBlob().
    _read      = NULL;
    _library   = NULL;

    _name      = NULL;
    _nameAlloc = 0;

    _aseq      = NULL;
    _aqlt      = NULL;

    _rseq      = NULL;   //  Raw read sequence and quality.
    _rqlt      = NULL;
    _rseqAlloc = 0;
    _rqltAlloc = 0;

    _cseq      = NULL;   //  Corrected read sequence and quality.
    _cqlt      = NULL;
    _cseqAlloc = 0;
    _cqltAlloc = 0;

    _tseq      = NULL;   //  Trimmed read sequence and quality.
    _tqlt      = NULL;

    _blobLen   = 0;
    _blobMax   = 0;
    _blob      = NULL;
  };

  ~sqReadData() {
    delete [] _name;

    delete [] _rseq;
    delete [] _rqlt;

    delete [] _cseq;
    delete [] _cqlt;

    //delete [] _tseq;  //  The trimmed read is just a
    //delete [] _tqlt;  //  pointer into the corrected read.

    delete [] _blob;
  };

  sqRead     *sqReadData_getRead(void)                { return(_read); };
  sqLibrary  *sqReadData_getLibrary(void)             { return(_library); };

  char       *sqReadData_getName(void)                { return(_name); };

  char       *sqReadData_getRawSequence(void)         { return(_rseq);  };
  char       *sqReadData_getCorrectedSequence(void)   { return(_cseq);  };
  char       *sqReadData_getTrimmedSequence(void)     { return(_tseq);  };

  char       *sqReadData_getSequence(sqRead_version vers = sqRead_defaultVersion) {
    if      (vers == sqRead_raw)          return(_rseq);
    else if (vers == sqRead_corrected)    return(_cseq);
    else if (vers == sqRead_trimmed)      return(_tseq);
    else                                  return(_aseq);
  };

  uint8      *sqReadData_getRawQualities(void)        { return(_rqlt);  };
  uint8      *sqReadData_getCorrectedQualities(void)  { return(_cqlt);  };
  uint8      *sqReadData_getTrimmedQualities(void)    { return(_tqlt);  };

  uint8      *sqReadData_getQualities(sqRead_version vers = sqRead_defaultVersion) {
    if      (vers == sqRead_raw)          return(_rqlt);
    else if (vers == sqRead_corrected)    return(_cqlt);
    else if (vers == sqRead_trimmed)      return(_tqlt);
    else                                  return(_aqlt);
  };

public:
  //  Set the name of the read.
  //
  //  Set the sequence and quality values of the read.  S and Q must be
  //  the same length.  S is a character string, NUL terminated.  Q is
  //  8-byte integers, also NUL terminated.
  //
  //  If Q[0] == 255, then no QVs exists, and the library
  //  default is used instead.

  void        sqReadData_setName(char *H);
  void        sqReadData_setBasesQuals(char *S, uint8 *Q);

private:
  uint32      sqReadData_encode2bit(uint8  *&chunk, char  *seq, uint32 seqLen);
  uint32      sqReadData_encode3bit(uint8  *&chunk, char  *seq, uint32 seqLen);
  uint32      sqReadData_encode4bit(uint8  *&chunk, uint8 *qlt, uint32 qltLen);
  uint32      sqReadData_encode5bit(uint8  *&chunk, uint8 *qlt, uint32 qltLen);

  uint32      sqReadData_encodeConstantQV(uint8 *qlt, uint32 qltLen, uint32 defaultQV);

  void        sqReadData_encodeBlobChunk(char const *tag, uint32 len, void *dat);
  void        sqReadData_encodeBlob(void);


  bool        sqReadData_decode2bit(uint8  *chunk, uint32 chunkLen, char  *seq, uint32 seqLen);
  bool        sqReadData_decode3bit(uint8  *chunk, uint32 chunkLen, char  *seq, uint32 seqLen);
  bool        sqReadData_decode4bit(uint8  *chunk, uint32 chunkLen, uint8 *qlt, uint32 qltLen);
  bool        sqReadData_decode5bit(uint8  *chunk, uint32 chunkLen, uint8 *qlt, uint32 qltLen);

  void        sqReadData_loadFromBlob(uint8 *blob);

private:
  sqRead            *_read;     //  Pointer to the read         set in sqStore_addEmptyRead() and
  sqLibrary         *_library;  //  Pointer to the library             sqStore_loadReadData().

  char              *_name;
  uint32             _nameAlloc;

  char              *_aseq;       //  Pointer to 'active' data.  If trimmed, this will point
  uint8             *_aqlt;       //  to the middle of _cseq (which will be terminated early).

  char              *_rseq;
  uint8             *_rqlt;
  uint32             _rseqAlloc;
  uint32             _rqltAlloc;

  char              *_cseq;
  uint8             *_cqlt;
  uint32             _cseqAlloc;
  uint32             _cqltAlloc;

  char              *_tseq;
  uint8             *_tqlt;

  uint32             _blobLen;
  uint32             _blobMax;
  uint8             *_blob;     //  And maybe even an encoded blob of data from the store.

  friend class sqRead;
  friend class sqStore;
  friend class sqCache;
};




class sqRead {
public:
  sqRead() {
    _readID       = 0;            //  Must be zero so first read (id 0) is
    _libraryID    = 0;            //  initialized correctly.

    _rseqLen      = 0;
    _cseqLen      = 0;

    _clearBgn     = 0;
    _clearEnd     = 0;

    _blobLen      = 0;
    _mByteHigh    = 0;
    _unusedB      = 0;
    _unusedC      = 0;
    _ignore       = false;
    _unusedF      = 0;
    _cExists      = false;
    _tExists      = false;

    _mSegm        = 0;
    _mByteLow     = 0;
    _mPart        = 0;
  };

  ~sqRead() {
  };

  uint32      sqRead_readID(void)     { return(_readID);    };
  uint32      sqRead_libraryID(void)  { return(_libraryID); };

  uint32      sqRead_sequenceLength(sqRead_version vers = sqRead_defaultVersion) {
    if      (_ignore)
      return(0);

    if      (vers == sqRead_raw)          return(_rseqLen);
    else if (vers == sqRead_corrected)    return(_cseqLen);
    else if (vers == sqRead_trimmed)      return(_clearEnd - _clearBgn);
    else {
      if      (_tExists)   return(_clearEnd - _clearBgn);
      else if (_cExists)   return(_cseqLen);
      else                 return(_rseqLen);
    }
  };

  uint32      sqRead_clearBgn(void)   { return(_clearBgn); };
  uint32      sqRead_clearEnd(void)   { return(_clearEnd); };

  bool        sqRead_ignore(void)     { return(_ignore);   };
  bool        sqRead_cExists(void)    { return(_cExists);  };
  bool        sqRead_tExists(void)    { return(_tExists);  };

  uint64      sqRead_mSegm(void)      { return(_mSegm);    };
  uint64      sqRead_mByte(void)      { return((((uint64)_mByteHigh) << 32) | ((uint64)_mByteLow));    };
  uint64      sqRead_mPart(void)      { return(_mPart);    };

private:
  void        sqRead_loadDataFromStream(sqReadData *readData, FILE *file);  //  'file' MUST be at correct position

private:
  //  Description of the read.  10 32-bit words.

  //  Each sqRead needs to know which read (raw, corrected or trimmed) is to be returned.  In
  //  particular, if corrections are done, but there is no corrected read for this raw read,
  //  _cseqLen will be zero, and we need to return zero.  _cExists and _tExists are used for this.

  uint32   _readID;
  uint32   _libraryID;

  uint32   _rseqLen;        //  Length of raw read stored.
  uint32   _cseqLen;        //  Length of corrected read stored.

  uint32   _clearBgn;       //  Trim points in corrected read.
  uint32   _clearEnd;

  uint32   _blobLen;        //  For easier loading of reads.
  uint32   _mByteHigh : 8;  //  The 8 high bits of a 40-bit mByte field.
  uint32   _unusedB   : 8;  //  Unused, to keep this struct 64-bit aligned.
  uint32   _unusedC   : 8;  //  (It's 5 64-bit words)
  uint32   _ignore    : 1;  //  Ignore this read; it's too short for our needs.
  uint32   _unusedF   : 5;
  uint32   _cExists   : 1;  //  Set by sqStore when the sqRead is returned to the user.
  uint32   _tExists   : 1;

  //  Last but not least, we need to remember where the actual data is stored in the blob.
  //
  //  While blobs in an unpartitioned store are limited to 1 GB (AS_BLOBFILE_MAX_SIZE), a blob in a
  //  partitioned store needs to hold all data for one partition.  Occasionally these have exceeded
  //  4 GB, and then fail to find the blob data.  Thus, the byte address into a blob has been
  //  extended to 48 bits, as [_mByteHigh | _mByteLow].  The keep the struct layout the same as
  //  before, the two values are in different locations.

  uint32   _mByteLow;       //  (very big)      - Pointer to the blob data in the blobs file.
  uint16   _mSegm;          //  65k files       - Pointer to the blobs file we are in.
  uint16   _mPart;          //  65k partitions  - If partitioned, which partition is it in.

  friend class sqReadData;
  friend class sqStore;
};

//  Even though we can store up to 4GB blob files, we artificially limit it to 1 GB
//  for (presumed) better caching on object storage systems.  Though there are
//  65k files allowed, pieces that stream through the store (correction, RED, OEA)
//  run out of file handles well before that.

const uint64 AS_BLOBFILE_MAX_SIZE  = 1024 * 1024 * 1024;

#endif  //  SQREAD_H
