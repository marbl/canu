
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

#ifndef LIBKMER_ITERATOR_H
#define LIBKMER_ITERATOR_H

#include "kmers.H"


//  Converts a buffer of characters (or a file of characters) into kmers, one
//  kmer at a time.


class kmerIterator {
public:
  kmerIterator(void) {
    reset();
    addSequence(NULL, 0);
  };
  kmerIterator(FILE *input);
  kmerIterator(char *buffer, uint64 bufferLen) {
    reset();
    addSequence(buffer, bufferLen);
  };

  void       reset(void) {
    _kmerLoad   = 0;
    _kmerValid  = _fmer.merSize() - 1;
  };

  void       addSequence(char *buffer, uint64 bufferLen) {
    _buffer    = buffer;
    _bufferLen = bufferLen;
    _bufferPos = 0;
  };

  bool       nextMer(void) {
  nextMer_anotherBase:
    if (_bufferPos >= _bufferLen)      //  No more sequence, and not a valid kmer.
      return(false);

    if ((_buffer[_bufferPos] != 'A') && (_buffer[_bufferPos] != 'a') &&
        (_buffer[_bufferPos] != 'C') && (_buffer[_bufferPos] != 'c') &&
        (_buffer[_bufferPos] != 'G') && (_buffer[_bufferPos] != 'g') &&
        (_buffer[_bufferPos] != 'T') && (_buffer[_bufferPos] != 't')) {
      _kmerLoad = 0;                   //  Not a valis base.  Clear the current
      _bufferPos++;                    //  kmer and move to the next base.
      goto nextMer_anotherBase;
    }

    _fmer.addR(_buffer[_bufferPos]);   //  A valid base, so push it onto
    _rmer.addL(_buffer[_bufferPos]);   //  the kmer.

    _bufferPos++;

    if (_kmerLoad < _kmerValid) {      //  Not a full kmer, keep loading
      _kmerLoad++;                     //  bases from the buffer.
      goto nextMer_anotherBase;
    }

    return(true);                      //  Valid kmer!
  };

  kmerTiny   fmer(void)      { return(_fmer);                        };
  kmerTiny   rmer(void)      { return(_rmer);                        };
  uint64     position(void)  { return(_bufferPos - _fmer.merSize()); };

private:
  uint32    _kmerLoad;
  uint32    _kmerValid;

  char     *_buffer;
  uint64    _bufferLen;
  uint64    _bufferPos;

  kmerTiny  _fmer;
  kmerTiny  _rmer;
};


#endif  //  LIBKMER_ITERATOR_H
