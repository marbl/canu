
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
 *    src/AS_BAT/memoryMappedFile.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2012-FEB-16 to 2013-AUG-01
 *      are Copyright 2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-NOV-26 to 2015-APR-21
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Sergey Koren beginning on 2015-DEC-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef AS_UTL_MEMORYMAPPEDFILE_H
#define AS_UTL_MEMORYMAPPEDFILE_H

#include "AS_global.H"
#include "AS_UTL_fileIO.H"

#include <fcntl.h>
#include <sys/mman.h>

#include <set>
#include <map>
using namespace std;

//  The BSD's are able to map to an arbitrary position in the file, but the Linux's can only map to
//  multiples of pagesize.  Thus, this class maps the whole file into address space, then passes out
//  pointers to pieces in it.  This is slightly unfortunate, because array out-of-bounds will not be
//  caught.  To be fair, on the BSD's the file is mapped to a length that is a multiple of pagesize,
//  so it would take a big out-of-bounds to fail.

enum memoryMappedFileType {
  memoryMappedFile_readOnly   = 0x00,
  memoryMappedFile_readWrite  = 0x01
};


#ifndef MAP_POPULATE
#define MAP_POPULATE 0
#endif

class memoryMappedFile {
public:
  memoryMappedFile(const char           *name,
                   memoryMappedFileType  type = memoryMappedFile_readOnly) {

    strcpy(_name, name);

    _type = type;

    errno = 0;
    int fd = (_type == memoryMappedFile_readOnly) ? open(_name, O_RDONLY | O_LARGEFILE)
                                                  : open(_name, O_RDWR   | O_LARGEFILE);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't open '%s' for mmap: %s\n", _name, strerror(errno)), exit(1);

    struct stat  sb;

    fstat(fd, &sb);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't stat '%s' for mmap: %s\n", _name, strerror(errno)), exit(1);

    _length = sb.st_size;
    _offset = 0;

    if (_length == 0)
      fprintf(stderr, "memoryMappedFile()-- File '%s' is empty, can't mmap.\n", _name), exit(1);

    //  Map a region that allows reading, or reading and shared writing.  One could add PROT_WRITE
    //  to the readOnly, and modifications will be kept private (if MAP_SHARED is switched to
    //  MAP_PRIVATE) to the process (and discarded at the end).
    //
    //  FreeBSD supports MAP_NOCORE which will exclude the region from any core files generated.  Linux does not support it.
    //
    //  Linux supports MAP_NORESERVE which will not reserve swap space for the file.  When reserved, a write is guaranteed to succeed.
    //
    //  NOTA BENE!!  Even though it is writable, it CANNOT be extended.

    _data = (_type == memoryMappedFile_readOnly) ? mmap(0L, _length, PROT_READ,              MAP_FILE | MAP_PRIVATE | MAP_POPULATE, fd, 0)
                                                 : mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, fd, 0);

    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't mmap '%s' of length "F_SIZE_T": %s\n", _name, _length, strerror(errno)), exit(1);

    close(fd);

    //fprintf(stderr, "memoryMappedFile()-- File '%s' of length %lu is mapped.\n", _name, _length);
  };

  ~memoryMappedFile() {
    if (_type == memoryMappedFile_readWrite)
      msync(_data, _length, MS_SYNC);

    munmap(_data, _length);
  };

  //  Return a pointer to position 'offset' in the file, and set the current position to 'offset +
  //  length'.
  //
  //  The length parameter is checked against the length of the file, and a fatal error occurs if
  //  'offset + length' exceeds the bounds of the file.
  //
  void  *get(size_t offset, size_t length) {

    if (length == 0)
      length = _length - offset;

    if (offset + length > _length)
      fprintf(stderr, "memoryMappedFile()-- Requested "F_SIZE_T" bytes at position "F_SIZE_T" in file '%s', but only "F_SIZE_T" bytes in file.\n",
              length, offset, _name, _length), exit(1);

    _offset = offset + length;

    return((uint8 *)_data + offset);
  };

  //  Return a pointer to the current position in the file, and move the current position ahead by
  //  'length' bytes.
  //
  //  get() (or get(0)) returns the current position without moving it.
  //
  void *get(size_t length=0) {
    return(get(_offset, length));
  };


  size_t  length(void) {
    return(_length);
  };

  memoryMappedFileType   type(void) {
    return(_type);
  };


private:
  char                    _name[FILENAME_MAX];

  memoryMappedFileType    _type;

  size_t                  _length;  //  Length of the mapped file
  size_t                  _offset;  //  File pointer for reading

  void                   *_data;
};

#endif  //  AS_UTL_MEMORYMAPPEDFILE_H
