
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
 *    kmer/libutil/readBuffer.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2003-APR-15 to 2004-JAN-06
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAR-17 to 2004-MAY-24
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2006-JUN-24 to 2014-APR-11
 *      are Copyright 2006-2008,2010,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-08
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef READ_BUFFER_H
#define READ_BUFFER_H

#include "AS_global.H"
#include "memoryMappedFile.H"

class readBuffer {
public:
  readBuffer(const char *filename, uint64 bufferMax = 32 * 1024);
  readBuffer(FILE *F, uint64 bufferMax = 32 * 1024);
  ~readBuffer();

  bool                 eof(void) { return(_eof); };

  char                 peek(void);

  char                 read(void);
  uint64               read(void *buf, uint64 len);
  uint64               read(void *buf, uint64 maxlen, char stop);

  void                 seek(uint64 pos);
  uint64               tell(void) { return(_filePos); };

  const char          *filename(void) { return(_filename); };

private:
  void                 fillBuffer(void);
  void                 init(int fileptr, const char *filename, uint64 bufferMax);

  char               *_filename;

  int                 _file;
  uint64              _filePos;

  memoryMappedFile   *_mmap;
  bool                _stdin;

  bool                _eof;

  //  If bufferMax is zero, then we are using the mmapped interface, otherwise,
  //  we are using a open()/read() and a small buffer.

  uint64              _bufferPos;
  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;
};


//  Returns the next letter in the buffer, but DOES NOT advance past
//  it.  Might have some wierd interaction with EOF -- if you peek()
//  and the next thing is eof , the _eof flag might get set.
//
inline
char
readBuffer::peek(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  return(_buffer[_bufferPos]);
}


//  Returns the next letter in the buffer.  Returns EOF (0) if there
//  is no next letter.
//
inline
char
readBuffer::read(void) {

  if ((_eof == false) && (_bufferPos >= _bufferLen))
    fillBuffer();

  if (_eof)
    return(0);

  _bufferPos++;
  _filePos++;

  return(_buffer[_bufferPos-1]);
}


#endif  //  READ_BUFFER_H
