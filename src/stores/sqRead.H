
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

#ifndef SQREAD_H
#define SQREAD_H

//  DO NOT INCLUDE THIS FILE DIRECTLY, include sqStore.H.

#include "arrays.H"
#include "sequence.H"

using namespace merylutil;  //  HUGE cop out.

class sqRead;
class sqLibrary;

class sqStore;
class sqStoreBlobWriter;

class sqCache;

//  Even though we can store up to 4GB blob files, we artificially limit it to 1 GB
//  for (presumed) better caching on object storage systems.  Though there are
//  65k files allowed, pieces that stream through the store (correction, RED, OEA)
//  run out of file handles well before that.

const uint64 AS_BLOBFILE_MAX_SIZE  = 1024 * 1024 * 1024;



//  The default version is set either by the user explicitly, or by the store
//  when it is opened.  It should never be unset.
//
//  If the store is created (new sqStore(...)) with no default set, it will
//  examine what reads are present and return the 'latest' available (raw,
//  corrected, trimmed, in that order).  If the homopolymer compression flag
//  exists, compressed reads will be returned.
//
//  If the default version is set before the store is opened, exactly those
//  reads will be returned.
//
//  The need for sqRead_normal arises when the store is set to return
//  homopoly compressed reads by default (if file 'homopolymerCompression'
//  exists).  sqRead_normal prevents sqStore_loadMetadata from enabling
//  homopoly compression.  It is otherwise not used (and so gets a quite
//  bogus value).
//
//  The 'exclude' function works ONLY AFTER the store is opened - since there
//  is no default set we can't remove a flag yet.
//

typedef uint32  sqRead_which;

const sqRead_which sqRead_unset      = 0x0000;   //  Nothing specified, query the store to decide what to return.
const sqRead_which sqRead_raw        = 0x0001;   //  Use noisy sequence.
const sqRead_which sqRead_corrected  = 0x0002;   //  Use corrected sequence.
const sqRead_which sqRead_normal     = 0x0080;   //  Return normal uncompressed bases.
const sqRead_which sqRead_compressed = 0x0004;   //  Return compressed bases.
const sqRead_which sqRead_trimmed    = 0x0008;   //  Return trimmed bases.
const sqRead_which sqRead_largest    = 0x0010;   //  Used to size an array.

extern
sqRead_which  sqRead_defaultVersion;

static
const
char *
toString(sqRead_which w) {
  const sqRead_which  c  = sqRead_compressed;
  const sqRead_which  t  = sqRead_trimmed;
  const sqRead_which  ct = sqRead_compressed | sqRead_trimmed;

  w &= ~sqRead_normal;  //  Disable flag so we report the normal uncompressed string.

  switch (w) {
    case sqRead_unset:            return("unset");                           break;

    case sqRead_raw:              return("raw");                             break;
    case sqRead_raw | c:          return("raw-compressed");                  break;
    case sqRead_raw | t:          return("raw-trimmed");                     break;
    case sqRead_raw | ct:         return("raw-compressed-trimmed");          break;

    case sqRead_corrected:        return("corrected");                       break;
    case sqRead_corrected | c:    return("corrected-compressed");            break;
    case sqRead_corrected | t:    return("corrected-trimmed");               break;
    case sqRead_corrected | ct:   return("corrected-compressed-trimmed");    break;

    case sqRead_compressed:       return("compressed");                      break;
    case sqRead_trimmed:          return("trimmed");                         break;

    default:                      return("undefined-mode");                  break;
  }

  return("undefined-mode");
}

static
void
sqRead_setDefaultVersion(sqRead_which v) {
  sqRead_defaultVersion = v;
}

static
void
sqRead_setDefaultVersionInclude(sqRead_which v) {
  sqRead_defaultVersion |= v;
}

static
void
sqRead_setDefaultVersionExclude(sqRead_which v) {
  sqRead_defaultVersion &= ~v;
}

static    //  Return true of the default version has mode 'v' set.
bool
sqRead_defaultIs(sqRead_which v) {
  return((sqRead_defaultVersion & v) == v);
}

static    //  Return true of the default version has mode 'v' not set.
bool
sqRead_defaultIsNot(sqRead_which v) {
  return((sqRead_defaultVersion & v) == sqRead_unset);
}

static
const
char *
sqRead_getDefaultVersion(void) {
  return(toString(sqRead_defaultVersion));
}








//  On disk sequence metadata.  Sequence data itself is in the blobs.
//
//  In general, you should not be directly using this class:
//    sqReadSeq_length() ALWAYS returns the untrimmed length of the read.
//
//    sqReadSeq_clearBgn() and sqReadSeq_clearEnd() will both return zero for
//    a read with no clear range set.
//
//  There are too many places where encapsulation needs to be broken to
//  make it truly private.  For example:
//    Marking reads as ignore in sqStoreCreate (to get rid of extra coverage).
//    Loading clear ranges.
//
class sqReadSeq {
public:
  void        sqReadSeq_initialize(void) {
    _seqValid  = 0;
    _unused1   = 0;
    _seqLength = 0;

    _ignoreU   = 0;
    _ignoreT   = 0;
    _clearBgn  = 0;

    _trimmed   = 0;
    _unused2   = 0;
    _clearEnd  = 0;
  };

private:
  uint32      sqReadSeq_length(void)          { return(_seqLength);  };   //  ALWAYS untrimmed length.

  uint32      sqReadSeq_clearBgn(void)        { assert(_trimmed);  return(_clearBgn);   };   //  NOT valid unless trimmed.
  uint32      sqReadSeq_clearEnd(void)        { assert(_trimmed);  return(_clearEnd);   };

private:
  bool        sqReadSeq_valid(void)           { return(_seqValid);   };   //  True if there is data.
  bool        sqReadSeq_trimmed(void)         { return(_trimmed);    };   //  True if the clear range is set.

  bool        sqReadSeq_ignoreU(void)         { return(_ignoreU);               };   //  True if this read should be ignored.
  bool        sqReadSeq_ignoreT(void)         { return(_ignoreU | _ignoreT);    };   //  True if the trimmed version should be ignored.


  //  Call ONLY for initializing with a newly added sequence.
  //  Do NOT call for trimmed reads.
  //  The only caller should be sqReadDataWriter::sqReadDataWriter_writeBlob().
  //  If you're thinking you want to call this, think again.
private:
  void        sqReadSeq_setLength(char *bases, uint32 basesLen, bool doCompress) {
    uint32  sl = 0;

    assert(bases[basesLen] == 0);

    if (doCompress == false)
      sl = basesLen;
    else
      sl = merylutil::homopolyCompress(bases, basesLen);

    assert(_seqValid == 0);

    _seqValid  = 1;
    _unused1   = 0;
    _seqLength = sl;

    _ignoreU   = 0;
    _ignoreT   = 0;
    _clearBgn  = 0;

    _trimmed   = 0;
    _unused2   = 0;
    _clearEnd  = sl;
  };


  //  These are public, but unless the store is opened for appending, changes
  //  will be lost.
public:
  void        sqReadSeq_setAllClear(void) {
    _clearBgn = 0;
    _clearEnd = _seqLength;
    _trimmed  = true;

    _ignoreT |= _ignoreU;   //  If the untrimmed is ignored, ignore the trimmed too.
  };

  void        sqReadSeq_setClearRange(uint32 bgn, uint32 end, bool set=true) {
    _clearBgn = bgn;
    _clearEnd = end;
    _trimmed  = set;

    _ignoreT |= _ignoreU;   //  If the untrimmed is ignored, ignore the trimmed too.
  };

private:
  //  Only access from sqStore_setIgnored().
  void        sqReadSeq_setIgnoreU(void) {
    _ignoreU  = true;
    _ignoreT |= true;   //  Also set ignoreT if untrimmed is ignored.
  };

  void        sqReadSeq_setIgnoreT(void) {
    _ignoreT = true;
  };


private:
  //  The data are logically out of order, so we can fit things into three 32-bit words.

  uint32  _seqValid  : 1;      //  The sequence for this record is present in the blob.
  uint32  _unused1   : 1;      //  (unused)
  uint32  _seqLength : 30;     //  The length of the sequence stored in the blob.

  uint32  _ignoreU   : 1;      //  Ignore both the untrimmed and trimmed versions.
  uint32  _ignoreT   : 1;      //  Ignore only the trimmed version.
  uint32  _clearBgn  : 30;     //

  uint32  _trimmed   : 1;      //  The trim points are valid.
  uint32  _unused2   : 1;      //  (unused)
  uint32  _clearEnd  : 30;     //

  //  Friends are only allowed to access the methods, not directly the data!

  friend class sqRead;
  friend class sqStore;
  friend class sqCache;
  friend class sqStoreInfo;
  friend class sqReadDataWriter;
};



//  On disk read metadata, in particular, the pointer to the blob data.
class sqReadMeta {
public:
  void        sqReadMeta_initialize(uint32 readID    = 0,
                                    uint32 libraryID = 0) {
    _readID            = readID;      assert(_readID    == readID);
    _libraryID         = libraryID;   assert(_libraryID == libraryID);
    _assignment        = 0;
    _assignmentScore   = 0;

    _unused            = 0;
    _mSegm             = 0;
    _mByte             = 0;
  };

  uint32      sqRead_readID(void)           { return(_readID);          };
  uint32      sqRead_libraryID(void)        { return(_libraryID);       };

  uint32      sqRead_assignment(void)       { return(_assignment);      };
  uint32      sqRead_assignmentScore(void)  { return(_assignmentScore); };

  uint64      sqRead_mSegm(void)            { return(_mSegm);           };
  uint64      sqRead_mByte(void)            { return(_mByte);           };

  void        sqRead_setPosition(uint64  mSegm,
                                 uint64  mByte) {
    _mSegm = mSegm;   assert(_mSegm == mSegm);   //  Check that mSegm and mByte are in
    _mByte = mByte;   assert(_mByte == mByte);   //  the range that can be stored.
  };

private:
  uint64           _readID          : 30;   //    1 billion
  uint64           _libraryID       : 12;   //    4 thousand
  uint64           _assignment      : 15;   //   32 thousand
  uint64           _assignmentScore : 7;    //  128
#if 30 + 12 + 15 + 7 != 64
#error sqReadMeta bits #1 are the wrong size.
#endif

  uint64           _unused          :  8;   //   -- FOR FUTURE USE
  uint64           _mSegm           : 16;   //  64 thousand files       - Pointer to the blobs file we are in.
  uint64           _mByte           : 40;   //    1 TB                  - Pointer to the blob data in the blobs file.
#if 8 + 16 + 40 != 64
#error sqReadMeta bits #2 are the wrong size.
#endif
};



//  In core read representation.  Only instantiated as needed, and sequence data is only loaded
//  as requested.
//
//  This is really the old sqReadData structure.  It just got promoted.
//
class sqRead {
public:
  ~sqRead() {
    delete [] _metaA;
    delete [] _rseqA;

    delete [] _blob;

    delete [] _name;
    delete [] _rawBases;
    delete [] _corBases;
    delete [] _retBases;
  };

  uint32      sqRead_readID(void)       { return(_meta->sqRead_readID());    };
  uint32      sqRead_libraryID(void)    { return(_meta->sqRead_libraryID()); };
  sqLibrary  *sqRead_library(void)      { return(_library);                  };
  char       *sqRead_name(void)         { return(_name);                     };

  //  Like sqStore, return 0 for reads we shouldn't use.
  uint32      sqRead_length  (sqRead_which w=sqRead_defaultVersion) {
    sqReadSeq  *seq = sqRead_getSeq(w);

    if (w & sqRead_trimmed)
      return(((seq == NULL) ||
              (seq->sqReadSeq_trimmed() == false) ||
              (seq->sqReadSeq_valid()   == false) ||
              (seq->sqReadSeq_ignoreT() == true)) ? 0 : (seq->sqReadSeq_clearEnd() - seq->sqReadSeq_clearBgn()));

    else
      return(((seq == NULL) ||
              (seq->sqReadSeq_valid()   == false) ||
              (seq->sqReadSeq_ignoreU() == true)) ? 0 : (seq->sqReadSeq_length()));
  };

  uint32      sqRead_clearBgn(sqRead_which w=sqRead_defaultVersion) {
    sqReadSeq  *seq = sqRead_getSeq(w);

    if (seq == NULL)
      return(0);

    if (w & sqRead_trimmed)
      return(seq->sqReadSeq_clearBgn());
    else
      return(0);
  };

  uint32      sqRead_clearEnd(sqRead_which w=sqRead_defaultVersion) {
    sqReadSeq  *seq = sqRead_getSeq(w);

    if (seq == NULL)
      return(0);

    if (w & sqRead_trimmed)
      return(seq->sqReadSeq_clearBgn());
    else
      return(seq->sqReadSeq_length());
  };

  char       *sqRead_sequence(sqRead_which w=sqRead_defaultVersion) {
    bool  comp = ((w & sqRead_compressed) == sqRead_unset) ? false : true;
    bool  trim = ((w & sqRead_trimmed)    == sqRead_unset) ? false : true;

    //  Figure out which bases to return, either raw or corrected, and the
    //  length of that sequence.  'bases' must be a valid pointer (we're
    //  either loading raw or corrected bases, but the string can be length
    //  zero.  In all cases, though, it must be NUL terminated.

    sqReadSeq   *seq      = sqRead_getSeq(w);
    char        *bases    = NULL;
    uint32       basesLen = 0;

    if (w & sqRead_raw) {
      bases    = _rawBases;
      basesLen = _rawU->sqReadSeq_length();
    }

    if (w & sqRead_corrected) {
      bases    = _corBases;
      basesLen = _corU->sqReadSeq_length();
    }

    assert(bases           != NULL);
    assert(bases[basesLen] == 0);

    //  If neither compressed or trimmed, just return the sequence we already have.
    //
    if ((comp == false) &&
        (trim == false)) {
      return(bases);
    }

    //  If not compressed but trimmed, copy the trimmed bases into
    //  temporary space and return that.
    //
    if ((comp == false) &&
        (trim == true)) {
      uint32  bgn = seq->sqReadSeq_clearBgn();   //  Only valid if trimmed, do not make global!
      uint32  end = seq->sqReadSeq_clearEnd();

      merylutil::resizeArray(_retBases, 0, _retBasesAlloc, end - bgn + 1, _raAct::doNothing);

      memmove(_retBases, bases + bgn, end - bgn);
      _retBases[end-bgn] = 0;

      return(_retBases);
    }

    //  Otherwise, we need to homopolymer compress the sequence.  It'll be no
    //  longer than the uncompressed sequence, so we can allocate that much,
    //  instead of tracking down the actual length.

    merylutil::resizeArray(_retBases, 0, _retBasesAlloc, basesLen + 1, _raAct::doNothing);

    homopolyCompress(bases, basesLen, _retBases);

    //  Trim the read if needed.
    //
    if (trim == true) {
      uint32  bgn = seq->sqReadSeq_clearBgn();   //  Only valid if trimmed, do not make global!
      uint32  end = seq->sqReadSeq_clearEnd();

      memmove(_retBases, _retBases + bgn, end - bgn);
      _retBases[end-bgn] = 0;
    }

    //  But either way, we return the same thing.
    return(_retBases);
  };

private:
  void        sqRead_fetchBlob(merylutil::readBuffer *B);
  void        sqRead_decodeBlob(void);

private:
  sqReadSeq  *sqRead_getSeq(sqRead_which w) {

    if (w == sqRead_unset)         //  If w is somehow set to the unset state,
      w = sqRead_defaultVersion;   //  go back to using the global default.

    bool  isRaw = ((w & sqRead_raw)        == sqRead_raw);         //  Return the sqReadSeq object
    bool  isCor = ((w & sqRead_corrected)  == sqRead_corrected);   //  appropriate for the supplied
    bool  isTrm = ((w & sqRead_trimmed)    == sqRead_trimmed);     //  sqRead_which.
    bool  isCmp = ((w & sqRead_compressed) == sqRead_compressed);

    if ((isRaw == true) && (isCmp == false))   { return(_rawU); };
    if ((isRaw == true) && (isCmp == true))    { return(_rawC); };

    if ((isCor == true) && (isCmp == false))   { return(_corU); };
    if ((isCor == true) && (isCmp == true))    { return(_corC); };

    fprintf(stderr, "sqRead_getSeq()-- Unknown which '%s'\n", toString(w));

    assert(0);
    return(NULL);
  };

private:
  sqReadMeta   *_meta = nullptr;                 //  Pointers to arrays owned by sqStore.
  sqReadSeq    *_rawU = nullptr;
  sqReadSeq    *_rawC = nullptr;
  sqReadSeq    *_corU = nullptr;
  sqReadSeq    *_corC = nullptr;

  sqReadMeta   *_metaA = nullptr;                //  When loading without a store, these
  sqReadSeq    *_rseqA = nullptr;                //  hold the data for the above pointers.

  sqLibrary    *_library = nullptr;

  bool          _blobLoaded = false;

  char          _blobName[4] = {0};
  uint32        _blobLen     =  0;
  uint32        _blobMax     =  0;
  uint8        *_blob        =  nullptr;

  uint32        _nameAlloc = 0;
  char         *_name      = nullptr;

  uint32        _rawBasesAlloc = 0;              //  The raw sequence, as loaded from disk.
  char         *_rawBases      = nullptr;

  uint32        _corBasesAlloc = 0;              //  The corrected sequence, as loaded from disk.
  char         *_corBases      = nullptr;

  sqRead_which  _retFlags      = sqRead_unset;   //  Remember what the returned sequence is.
  uint32        _retBasesAlloc = 0;              //  Scratch space for computing trimmed and
  char         *_retBases      = nullptr;        //  compressed sequences to return to the user.

  friend class sqStore;
  friend class sqCache;
  friend class sqReadDataWriter;
};




class sqReadDataWriter {
public:
  sqReadDataWriter(sqReadMeta *meta=NULL,
                   sqReadSeq  *rawu=NULL,
                   sqReadSeq  *rawc=NULL,
                   sqReadSeq  *coru=NULL,
                   sqReadSeq  *corc=NULL) {
    _meta          = meta;
    _rawU          = rawu;
    _rawC          = rawc;
    _corU          = coru;
    _corC          = corc;

    _charMap['A'] = _charMap['a'] = 'A';
    _charMap['C'] = _charMap['c'] = 'C';
    _charMap['G'] = _charMap['g'] = 'G';
    _charMap['T'] = _charMap['t'] = 'T';
    _charMap['U'] = _charMap['u'] = 'T';   //  Map U to T.
    _charMap['N'] = _charMap['n'] = 'N';
  };

  ~sqReadDataWriter() {
    delete [] _name;
    delete [] _rawBases;
    delete [] _corBases;
  };

public:
  void        sqReadDataWriter_importData(sqRead *read);
  void        sqReadDataWriter_importData(char const *readName,
                                          char const *readBases,
                                          uint32      readLength, uint32 bgn, uint32 end,
                                          sqReadSeq  *rawU,
                                          sqReadSeq  *rawC,
                                          sqReadSeq  *corU,
                                          sqReadSeq  *corC);

  void        sqReadDataWriter_setName(const char *N);
  void        sqReadDataWriter_setRawBases(const char *S, uint32 Slen);
  void        sqReadDataWriter_setCorrectedBases(const char *S, uint32 Slen);

  void        sqReadDataWriter_writeBlob(writeBuffer *buffer);

  uint32      sqReadDataWriter_getRawLength(bool compressed) {
    return((compressed == false) ? _rawU->sqReadSeq_length() : _rawC->sqReadSeq_length());
  };

  uint32      sqReadDataWriter_getCorrectedLength(bool compressed) {
    return((compressed == false) ? _corU->sqReadSeq_length() : _corC->sqReadSeq_length());
  };

private:
  sqReadMeta  *_meta = nullptr;           //  Pointer to the read metadata.
  sqReadSeq   *_rawU = nullptr;           //    In sqStoreCreate, these pointers are
  sqReadSeq   *_rawC = nullptr;           //    back to the data in sqStore itself.
  sqReadSeq   *_corU = nullptr;
  sqReadSeq   *_corC = nullptr;

  uint32       _nameAlloc     = 0;
  uint32       _nameLen       = 0;
  char        *_name          = nullptr;

  uint32       _rawBasesAlloc = 0;        //  The raw sequence.  The length is the
  uint32       _rawBasesLen   = 0;        //  length of the array, not of the string.
  char        *_rawBases      = nullptr;

  uint32       _corBasesAlloc = 0;        //  The corrected sequence.  The length is the
  uint32       _corBasesLen   = 0;        //  length of the array, not of the string.
  char        *_corBases      = nullptr;

  char         _charMap[256]  = { 0 };

  friend class sqStore;
  friend class sqStoreBlobWriter;
};



#endif  //  SQREAD_H
