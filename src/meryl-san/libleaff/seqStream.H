
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
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef SEQSTREAM_H
#define SEQSTREAM_H

#include "seqFile.H"

struct seqStreamIndex {
  uint32   _iid;  //  seqFile IID
  uint32   _len;  //  length of the sequence
  uint64   _bgn;  //  begin position in the stream
};


class seqStream {
public:
  seqStream(const char *filename);
  seqStream(const char *sequence, uint32 length);
  ~seqStream();

  //  Separate sequences with this letter.  Non-ACGT is always
  //  returned as 'N'.  Changing the length of the separator AFTER
  //  setting the range will result in the wrong range being used.
  //
  void              setSeparator(char sep, uint32 len);

  //  get() returns one letter per input letter -- a gap of size n
  //  will return n gap symbols.
  //
  unsigned char     get(void);
  bool              eof(void)        { return(_eof); };

  //  Returns to the start of the range.
  //
  void              rewind(void);

  //  Set the range of ACGT sequence we will return.  Coordinates are
  //  space-based.  Example:
  //
  //  >one
  //  AAA
  //  >two
  //  C
  //  >three
  //  GGG
  //
  //  We separate these sequences with three '-' letters.
  //
  //  strPos    012...3...456
  //            AAA---C---GGG
  //
  //  range(0,0) -> nothing
  //  range(0,1) -> A
  //  range(0,3) -> AAA
  //  range(0,4) -> AAAnnnC
  //  range(0,5) -> AAAnnnCnnnG
  //
  void              setRange(uint64 bgn, uint64 end);
  void              setPosition(uint64 pos);

  //  seqPos() is the position we are at in the current sequence;
  //  seqIID() is the iid of that sequence;
  //  strPos() is the position we are at in the chained sequence
  //
  //  Values are not defined if the letter is a separator.
  //
  uint32            seqPos(void)     { return(_currentPos);            };
  uint32            seqIID(void)     { return(_idx[_currentIdx]._iid); };
  uint64            strPos(void)     { return(_streamPos);             };

  uint32            numberOfSequences(void) { return(_idxLen); };

  //  Return the length of, position of (in the chain) and IID of the
  //  (s)th sequence in the chain.
  //
  uint32            lengthOf(uint32 s) { return((s >= _idxLen) ? ~uint32ZERO : _idx[s]._len); };
  uint32            IIDOf(uint32 s)    { return((s >= _idxLen) ? ~uint32ZERO : _idx[s]._iid); };
  uint64            startOf(uint32 s)  { return((s >= _idxLen) ? ~uint64ZERO : _idx[s]._bgn); };

  //  For a chain position p, returns the s (above) for that position.
  //
  uint32            sequenceNumberOfPosition(uint64 p);
  void              tradeSpaceForTime(void);

private:
  void              fillBuffer(void);

  seqFile         *_file;        //  Backed by a seqFile.
  char            *_string;      //  Backed by a character string.

  uint64           _bgn;         //  Begin/End position in chained sequence
  uint64           _end;

  uint32           _currentIdx;  //  index into _idx of the current sequence
  uint32           _currentPos;  //  position in the current sequence
  uint64           _streamPos;   //  position in the chained sequence

  // Buffer for holding sequence from the seqFile.

  uint32           _bufferMax;   //  Max size of the buffer
  uint32           _bufferLen;   //  Actual size of the buffer
  uint32           _bufferPos;   //  Position we are at in the buffer
  uint32           _bufferSep;   //  How much of the buffer is separator
  char            *_buffer;

  //  Info about the raw sequences

  uint32           _idxLen;
  seqStreamIndex  *_idx;

  uint32          *_seqNumOfPos;

  uint64           _lengthOfSequences;

  bool             _eof;

  char             _separator;
  uint32           _separatorLength;
};



#endif  //  SEQSTREAM_H
