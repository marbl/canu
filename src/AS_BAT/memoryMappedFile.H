
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

//static const char *rcsid = "$Id: memoryMappedFile.H,v 1.1 2012-02-16 20:14:26 brianwalenz Exp $";

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

class memoryMappedFile {
public:
  memoryMappedFile(const char *name_, size_t offset_=0) {

    strcpy(name, name_);

    errno = 0;
    fd = open(name, O_RDONLY | O_LARGEFILE);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't open '%s' for mapping: %s", name, strerror(errno)), exit(1);

    struct stat  sb;

    fstat(fd, &sb);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Couldn't stat '%s': %s", name, strerror(errno)), exit(1);

    totalLength = sb.st_size;

    currentOffset = offset_;
  };


  ~memoryMappedFile() {
    map<void *, size_t>::iterator  it;

    for (it=addressToLength.begin(); it != addressToLength.end(); it++)
      if (it->second > 0)
        munmap(it->first, it->second);

    close(fd);
  };


  void  *get(size_t offset_, size_t length_) {

    if (length_ == 0)
      length_ = totalLength - offset_;

    if (offset_ + length_ > totalLength)
      fprintf(stderr, "memoryMappedFile()-- Requested "F_SIZE_T" bytes at position "F_SIZE_T" in file '%s', but only "F_SIZE_T" bytes in file.\n",
              length_, offset_, name, totalLength), exit(1);

    void *ptr = mmap(0L, length_, PROT_READ, MAP_FILE | MAP_SHARED, fd, offset_);
    if (errno)
      fprintf(stderr, "memoryMappedFile()-- Requested "F_SIZE_T" bytes at position "F_SIZE_T" in file '%s' (length "F_SIZE_T" bytes) failed: %s\n",
              length_, offset_, name, totalLength, strerror(errno)), exit(1);

    addressToLength[ptr] = length_;

    currentOffset = offset_ + length_;

    return(ptr);
  };


  void *get(size_t length_) {
    return(get(currentOffset, length_));
  };


  bool   release(void *ptr) {
    size_t   length = addressToLength[ptr];

    if (length == 0)
      return(false);

    munmap(ptr, length);

    addressToLength[ptr] = 0;
    return(true);
  };

private:
  char                    name[FILENAME_MAX];

  int                     fd;
  size_t                  totalLength;

  size_t                  currentOffset;

  map<void *, size_t>     addressToLength;
};

#endif  //  AS_UTL_MEMORYMAPPEDFILE_H
