
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

#ifndef WRITE_BUFFER_H
#define WRITE_BUFFER_H

#include "AS_global.H"
#include "files.H"

class writeBuffer {
public:
  writeBuffer(const char *filename, const char *filemode, uint64 bufferMax = 1024 * 1024);
  ~writeBuffer();

  const char          *filename(void) { return(_filename); };
  uint64               tell(void)     { return(_filePos);  };

  void                 write(void *data, uint64 length);

private:
  void                 open(void);
  void                 writeToDisk(void *data, uint64 length);
  void                 flush(void);

  char                _filename[FILENAME_MAX+1];
  char                _filemode[17];

  FILE               *_file;
  uint64              _filePos;

  uint64              _bufferLen;
  uint64              _bufferMax;
  char               *_buffer;
};

#endif  //  WRITE_BUFFER_H
