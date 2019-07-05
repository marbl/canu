
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef LIBSEQUENCE_H
#define LIBSEQUENCE_H

#include "AS_global.H"
#include "types.H"
#include "files.H"


void   reverseComplementSequence(char *seq, int len);
char  *reverseComplementCopy(char *seq, int len);

template<typename qvType>
void   reverseComplement(char *seq, qvType *qlt, int len);

uint32 homopolyCompress(char *bases, char *compr=NULL);


//  Encode a sequence into chunk.  The length of the chunk in bytes is returned.
//  If the chunk is NULL, it is allocated.  Otherwise, it must be
//  at least seqLen bytes in length.
uint32 encode2bitSequence(uint8 *&chunk, char *seq, uint32 seqLen);
uint32 encode3bitSequence(uint8 *&chunk, char *seq, uint32 seqLen);
uint32 encode8bitSequence(uint8 *&chunk, char *seq, uint32 seqLen);


//  Decode an encoded sequence (in chunk) of length chunkLen.
//  If seq is NULL, space will be allocated, otherwise, seqLen+1 bytes
//  must be allocated.
//  seqLen must be the length of the sequence to decode.
void   decode2bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
void   decode3bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen);
void   decode8bitSequence(uint8 *chunk, uint32 chunkLen, char *seq, uint32 seqLen);




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
  uint32            _nameMax;
  char             *_name;
  uint64            _seqMax;
  char             *_seq;
  uint8            *_qlt;
  uint64            _seqLen;

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
  loadFASTA(char   *&name,     uint32  &nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64  &seqMax);

  uint64
  loadFASTQ(char   *&name,     uint32  &nameMax,
            char   *&seq,
            uint8  *&qlt,      uint64  &seqMax);


public:
  //  Return the next sequence in the file.
  //  Returns false if EOF, true otherwise, even if the sequence is length zero.
  //
  bool   loadSequence(char   *&name,     uint32  &nameMax,
                      char   *&seq,
                      uint8  *&qlt,      uint64  &seqMax,
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
