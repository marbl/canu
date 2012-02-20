
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

//static const char *rcsid = "$Id: memoryMappedFile.H,v 1.3 2012-02-20 01:53:38 brianwalenz Exp $";

#include <sys/types.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/sysctl.h>

#include <set>
#include <map>
using namespace std;

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

//  Note that the BSD's are able to map to an arbitrary position in the file,
//  but the Linux's can only map to multiples of pagesize.  Thus, this class
//  maps the whole file into address space, then passes out pointers to
//  pieces in it.  This is slightly unfortunate, because array out-of-bounds
//  will not be caught.  To be fair, on the BSD's the file is mapped to a length
//  that is a multiple of pagesize, so it would take a big out-of-bounds to fail.

class memoryMappedFile {
public:
  memoryMappedFile(const char *name_, size_t offset_=0) {

    strcpy(nm, name_);

    errno = 0;
    fd = open(nm, O_RDONLY | O_LARGEFILE);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't open '%s' for mapping: %s", nm, strerror(errno)), exit(1);

    struct stat  sb;

    fstat(fd, &sb);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't stat '%s': %s", nm, strerror(errno)), exit(1);

    ln = sb.st_size;
    co = offset_;

    pt = mmap(0L, ln, PROT_READ | PROT_WRITE, MAP_FILE | MAP_PRIVATE, fd, 0);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't mmap '%s': %s", nm, strerror(errno)), exit(1);

    close(fd);
  };


  ~memoryMappedFile() {
    munmap(pt, ln);
  };


  void  *get(size_t offset_, size_t length_) {

    if (length_ == 0)
      length_ = ln - offset_;

    if (offset_ + length_ > ln)
      fprintf(stderr, "memoryMappedFile()-- Requested "F_SIZE_T" bytes at position "F_SIZE_T" in file '%s', but only "F_SIZE_T" bytes in file.\n",
              length_, offset_, nm, ln), exit(1);

    co = offset_ + length_;

    return((char *)pt + offset_);
  };


  void *get(size_t length_) {
    return(get(co, length_));
  };


private:
  char                    nm[FILENAME_MAX];

  int                     fd;

  size_t                  ln;
  size_t                  co;

  void                   *pt;
};

#endif  //  AS_UTL_MEMORYMAPPEDFILE_H
