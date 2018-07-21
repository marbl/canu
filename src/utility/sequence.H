
/******************************************************************************
 *
 *  This file is part of 'sequence' and/or 'meryl', software programs for
 *  working with DNA sequence files and k-mers contained in them.
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2018-FEB-26
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.license' in the root directory of this distribution contains
 *  full conditions and disclaimers.
 */

#ifndef LIBSEQUENCE_H
#define LIBSEQUENCE_H

#include "AS_global.H"
#include "types.H"
#include "files.H"


class dnaSeqIndexEntry;   //  Internal use only, sorry.



class dnaSeq {
public:
  dnaSeq() {
    _nameMax = 0;
    _name    = NULL;
    _seqMax  = 0;
    _seq     = NULL;
    _qlt     = NULL;
    _seqLen  = 0;
  };

  ~dnaSeq() {
    delete [] _name;
    delete [] _seq;
    delete [] _qlt;
  };


  char             *name(void)         { return(_name);   };
  char             *bases(void)        { return(_seq);    };
  uint8            *quals(void)        { return(_qlt);    };

  uint64            length(void)       { return(_seqLen); };

private:
  uint32            _nameMax = 0;
  char             *_name    = NULL;
  uint64            _seqMax  = 0;
  char             *_seq     = NULL;
  uint8            *_qlt     = NULL;
  uint64            _seqLen  = 0;

  friend class dnaSeqFile;
};




class dnaSeqFile {
public:
  dnaSeqFile(const char *filename, bool indexed=false);
  ~dnaSeqFile();

  compressedFileReader  *_file;
  readBuffer            *_buffer;

  dnaSeqIndexEntry      *_index;
  uint64                 _indexLen;
  uint64                 _indexMax;

private:
  bool     loadIndex(void);
  void     saveIndex(void);

public:
  void     generateIndex(void);

  //  If indexed, searches the index for the proper sequence.
  //
  //  If not indexed, searches forward in the file for the sequence.  If not found,
  //  the file will be at the end.
  //
  //  In both cases, the file is left positioned at the start of the sequence header.
  //
  //  Returns true if found, false if not.
  //
  bool     findSequence(uint64 i);
  bool     findSequence(const char *name);

  //  Returns the number of sequences in the file.
  uint64   numberOfSequences(void) {
    return(_indexLen);
  };

  //  Returns the length of sequence i.  If no such sequence, returns UINT64_MAX.
  uint64   sequenceLength(uint64 i);

  char    *filename(void) {
    return(_file->filename());
  }

private:
  uint64
  loadFASTA(char   *&name,     uint32   nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64   seqMax);

  uint64
  loadFASTQ(char   *&name,     uint32   nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64   seqMax);


public:
  //  Return the next sequence in the file.
  //  Returns the number of bases loaded; zero if EOF.
  //
  bool   loadSequence(char   *&name,     uint32   nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64   seqMax,
                      uint64  &seqLen);

  bool   loadSequence(dnaSeq &seq) {
    return(loadSequence(seq._name,    seq._nameMax,
                        seq._seq,
                        seq._qlt,     seq._seqMax,
                        seq._seqLen));
  };

  //  Returns a chunk of sequence from the file, up to 'maxLength' bases or
  //  the end of the current sequence.
  //
  //  Returns false if EOF is hit and no bases were loaded.
  //
  bool   loadBases(char    *seq,
                   uint64   maxLength,
                   uint64  &seqLength,
                   bool    &endOfSequence);
};


#endif  //  LIBSEQUENCE_H
