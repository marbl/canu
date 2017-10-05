
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

#ifndef GKREAD_H
#define GKREAD_H


//  DO NOT INCLUDE THIS FILE DIRECTLY, include gkStore.H.


class gkRead;
class gkLibrary;

class gkReadData {
public:
  gkReadData() {
    _read      = NULL;   //  Set in stashReadData() and loadFromBlob()
    _library   = NULL;

    _name      = NULL;
    _nameAlloc = 0;

    _rseq      = NULL;   //  Raw read sequence and quality.
    _rqlt      = NULL;
    _rseqAlloc = 0;
    _rqltAlloc = 0;

    _cseq      = NULL;   //  Corrected read sequence and quality.
    _cqlt      = NULL;
    _cseqAlloc = 0;
    _cqltAlloc = 0;

    _blobLen   = 0;
    _blobMax   = 0;
    _blob      = NULL;
  };

  ~gkReadData() {
    delete [] _name;

    delete [] _rseq;
    delete [] _rqlt;

    delete [] _cseq;
    delete [] _cqlt;

    delete [] _blob;
  };

  gkRead     *gkReadData_getRead(void)                { return(_read); };
  gkLibrary  *gkReadData_getLibrary(void)             { return(_library); };

  char       *gkReadData_getName(void)                { return(_name); };

  char       *gkReadData_getRawSequence(void)         { return(_rseq);  };
  char       *gkReadData_getCorrectedSequence(void)   { return(_cseq);  };
  char       *gkReadData_getSequence(void)            { return(_aseq);  };

  uint8      *gkReadData_getRawQualities(void)        { return(_rqlt);  };
  uint8      *gkReadData_getCorrectedQualities(void)  { return(_cqlt);  };
  uint8      *gkReadData_getQualities(void)           { return(_aqlt);  };

public:
  //  Set the name of the read.
  //
  //  Set the sequence and quality values of the read.  S and Q must be
  //  the same length.  S is a character string, NUL terminated.  Q is
  //  8-byte integers, also NUL terminated.
  //
  //  If Q[0] == 255, then no QVs exists, and the library
  //  default is used instead.

  void        gkReadData_setName(char *H);
  void        gkReadData_setBasesQuals(char *S, uint8 *Q);

private:
  uint32      gkReadData_encode2bit(uint8  *&chunk, char  *seq, uint32 seqLen);
  uint32      gkReadData_encode3bit(uint8  *&chunk, char  *seq, uint32 seqLen);
  uint32      gkReadData_encode4bit(uint8  *&chunk, uint8 *qlt, uint32 qltLen);
  uint32      gkReadData_encode5bit(uint8  *&chunk, uint8 *qlt, uint32 qltLen);

  void        gkReadData_encodeBlobChunk(char const *tag, uint32 len, void *dat);
  void        gkReadData_encodeBlob(void);


  bool        gkReadData_decode2bit(uint8  *chunk, uint32 chunkLen, char  *seq, uint32 seqLen);
  bool        gkReadData_decode3bit(uint8  *chunk, uint32 chunkLen, char  *seq, uint32 seqLen);
  bool        gkReadData_decode4bit(uint8  *chunk, uint32 chunkLen, uint8 *qlt, uint32 qltLen);
  bool        gkReadData_decode5bit(uint8  *chunk, uint32 chunkLen, uint8 *qlt, uint32 qltLen);

  void        gkReadData_loadFromBlob(uint8 *blob);

private:
  gkRead            *_read;     //  Pointer to the read         set in gkStore_addEmptyRead() and
  gkLibrary         *_library;  //  Pointer to the library             gkStore_loadReadData().

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

  uint32             _blobLen;
  uint32             _blobMax;
  uint8             *_blob;     //  And maybe even an encoded blob of data from the store.

  friend class gkRead;
  friend class gkStore;
};




class gkRead {
public:
  gkRead() {
    _readID       = 0;            //  Must be zero so first read (id 0) is
    _libraryID    = 0;            //  initialized correctly.

    _rseqLen      = 0;
    _cseqLen      = 0;

    _clearBgn     = 0;
    _clearEnd     = 0;

    _cExists      = false;
    _tExists      = false;

    _mPtr         = 0;
    _pID          = 0;
  };

  ~gkRead() {
  };

  uint32      gkRead_readID(void)             { return(_readID);               };

  uint32      gkRead_libraryID(void)          { return(_libraryID);            };

  //  Move these to read data?  It's nice not having to load anything to get lengths.
  //  But then we need someway of telling if corrected/trimmed reads exist.

  uint32      gkRead_rawLength(void)          { return(_rseqLen);              };
  uint32      gkRead_correctedLength(void)    { return(_cseqLen);              };
  uint32      gkRead_trimmedLength(void)      { return(_clearEnd - _clearBgn); };

  uint32      gkRead_sequenceLength(void)     {
    if      (_tExists)   return(gkRead_trimmedLength());
    else if (_cExists)   return(gkRead_correctedLength());
    else                 return(gkRead_rawLength());
  };

  uint32      gkRead_clearBgn(void)           { return(_clearBgn); };
  uint32      gkRead_clearEnd(void)           { return(_clearEnd); };

  bool        gkRead_cExists(void)            { return(_cExists); };
  bool        gkRead_tExists(void)            { return(_tExists); };

  uint64      gkRead_mPtr(void) { return(_mPtr); };  //  For debugging, in gatekeeperDumpMetatData
  uint64      gkRead_pID(void)  { return(_pID);  };

private:
  void        gkRead_loadDataFromStream(gkReadData *readData, FILE *file);  //  file at position
  void        gkRead_loadDataFromFile  (gkReadData *readData, FILE *file);  //  move file pointer before reading
  void        gkRead_loadDataFromCore  (gkReadData *readData, void *blob);  //  read from blob in core

private:
  void        gkRead_copyDataToPartition(void  *blobs,      FILE **partfiles, uint64 *partfileslen, uint32 partID);
  void        gkRead_copyDataToPartition(FILE **blobsFiles, FILE **partfiles, uint64 *partfileslen, uint32 partID);

private:
  //  Description of the read.

  uint32   _readID;
  uint32   _libraryID;

  uint32   _rseqLen;        //  Length of raw read stored.
  uint32   _cseqLen;        //  Length of corrected read stored.

  uint32   _clearBgn;       //  Trim points in corrected read.
  uint32   _clearEnd;

  //  Each gkRead needs to know which read (raw, corrected or trimmed) is to be returned.
  //  In particular, if corrections are done, but there is no corrected read for this raw
  //  read, _cseqLen will be zero, and we need to return zero.

  uint64   _cExists : 1;    //  Set by gkStore when the gkRead is returned to the user.
  uint64   _tExists : 1;

  //  Last but not least, we need to remember where the actual data is stored in the blob.

  uint64   _mPtr    : 46;
  uint64   _pID     : 18;


  friend class gkReadData;
  friend class gkStore;
};


#endif  //  GKREAD_H
