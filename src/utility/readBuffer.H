
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
  char                 readuntil(char stop);
  uint64               read(void *buf, uint64 len);
  uint64               read(void *buf, uint64 maxlen, char stop);

  void                 skipLine(void);
  void                 skipAhead(char stop);
  uint64               copyUntil(char stop, char *dest, uint64 destLen);

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


#include "readBuffer-implementation.H"


#endif  //  READ_BUFFER_H
