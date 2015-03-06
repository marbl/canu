
/**************************************************************************
 * Copyright (C) 2011, J Craig Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

#ifndef AS_UTL_MEMORYMAPPEDFILE_H
#define AS_UTL_MEMORYMAPPEDFILE_H

//static const char *rcsid = "$Id$";

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
    //  NOTA BENE!!  Even though it is writable, it CANNOT be extended.

    _data = (_type == memoryMappedFile_readOnly) ? mmap(0L, _length, PROT_READ,              MAP_FILE | MAP_SHARED, fd, 0)
                                                 : mmap(0L, _length, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, fd, 0);

    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't mmap '%s': %s\n", _name, strerror(errno)), exit(1);

    close(fd);
  };

  ~memoryMappedFile() {
    if (_type == memoryMappedFile_readWrite)
      msync(_data, _length, MS_ASYNC);

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
