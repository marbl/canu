
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
 *    src/utility/memoryMappedFile.C
 *
 *  Modifications by:
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "files.H"

#include <fcntl.h>
#include <sys/mman.h>



memoryMappedFile::memoryMappedFile(const char           *name,
                                   memoryMappedFileType  type) {

  strncpy(_name, name, FILENAME_MAX-1);

  _type = type;

  errno = 0;
  _fd = (_type == memoryMappedFile_readOnly) ? open(_name, O_RDONLY | O_LARGEFILE)
                                             : open(_name, O_RDWR   | O_LARGEFILE);
  if (errno)
    fprintf(stderr, "memoryMappedFile()-- Couldn't open '%s' for mmap: %s\n", _name, strerror(errno)), exit(1);

  struct stat  sb;

  fstat(_fd, &sb);
  if (errno)
    fprintf(stderr, "memoryMappedFile()-- Couldn't stat '%s' for mmap: %s\n", _name, strerror(errno)), exit(1);

  _length = sb.st_size;
  _offset = 0;

  if (_length == 0)
    fprintf(stderr, "memoryMappedFile()-- File '%s' is empty, can't mmap.\n", _name), exit(1);

  //  Map the file to memory, or grab some anonymous space for the file to be copied to.

  if (_type == memoryMappedFile_readOnly)
    _data = mmap(0L, _length, PROT_READ,              MAP_FILE | MAP_PRIVATE, _fd, 0);

  if (_type == memoryMappedFile_readOnlyInCore)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_ANON | MAP_PRIVATE, -1, 0);

  if (_type == memoryMappedFile_readWrite)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, _fd, 0);

  if (_type == memoryMappedFile_readWriteInCore)
    _data = mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);

  //  If loading into core, read the file into core.

  if ((_type == memoryMappedFile_readOnlyInCore) ||
      (_type == memoryMappedFile_readWriteInCore))
    read(_fd, _data, _length);

  //  Close the file if we're done with it.

  if (_type != memoryMappedFile_readWriteInCore)
    close(_fd), _fd = -1;

  //  Catch any and all errors.

  if (errno)
    fprintf(stderr, "memoryMappedFile()-- Couldn't mmap '%s' of length " F_SIZE_T ": %s\n", _name, _length, strerror(errno)), exit(1);


  //fprintf(stderr, "memoryMappedFile()-- File '%s' of length %lu is mapped.\n", _name, _length);
};



memoryMappedFile::~memoryMappedFile() {

  errno = 0;

  if (_type == memoryMappedFile_readWrite)
    msync(_data, _length, MS_SYNC);

  if (_type == memoryMappedFile_readWriteInCore)
    write(_fd, _data, _length), close(_fd);

  if (errno)
    fprintf(stderr, "memoryMappedFile()-- Failed to close mmap '%s' of length " F_SIZE_T ": %s\n", _name, _length, strerror(errno)), exit(1);

  //  Destroy the mapping.

  munmap(_data, _length);
};


