
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
 *    Brian P. Walenz from 2003-MAY-06 to 2004-AUG-02
 *      are Copyright 2003-2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-JAN-12 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Clark Mobarry on 2004-FEB-18
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-MAR-09 to 2014-APR-11
 *      are Copyright 2005-2008,2011-2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <sys/statvfs.h>

#include "util.h"


int
isHuman(FILE *F) {
  return(isatty(fileno(F)));
}


#ifdef __alpha
unsigned long __sbrk_override = 1;  //  See malloc(3) for details.

#define MMAPFLAGS    (MAP_FILE | MAP_VARIABLE | MAP_SHARED)
#endif

#ifdef _AIX
#define MMAPFLAGS    (MAP_FILE | MAP_VARIABLE | MAP_SHARED)
#endif

#ifdef __CYGWIN__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __linux
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __FreeBSD__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif

#ifdef __sun
#define MMAPFLAGS    (MAP_SHARED)
#endif

#ifdef __APPLE__
#define MMAPFLAGS    (MAP_FILE | MAP_SHARED)
#endif




FILE*
makeTempFile(char *path) {
  char   template[PATH_MAX + 1];
  int    fildes;
  FILE  *F;

  if (path) {
    strcpy(template, path);
    strcat(template, "/XXXXXX");
  } else {
    strcpy(template, "XXXXXX");
  }

  errno = 0;
  fildes = mkstemp(template);
  if (errno) {
    fprintf(stderr, "Failed to create temporary file '%s': %s\n", template, strerror(errno));
    exit(1);
  }

  errno = 0;
  F = fdopen(fildes, "w+");
  if (errno) {
    fprintf(stderr, "Failed to open temporary file '%s': %s\n", template, strerror(errno));
    exit(1);
  }

  errno = 0;
  unlink(template);
  if (errno) {
    fprintf(stderr, "Failed to hide temporary file '%s': %s\n", template, strerror(errno));
    exit(1);
  }

  return(F);
}










void*
mapFile(const char *filename, uint64 *length, char mode) {
  void        *ptr = 0L;
  struct stat  sb;
  int          f;
  int          openMode = O_RDONLY | O_LARGEFILE;
  int          mapMode  = O_RDWR   | O_LARGEFILE;

  switch (mode) {
    case 'r':
      openMode = O_RDONLY | O_LARGEFILE;
      mapMode  = PROT_READ;
      break;
    case 'w':
      openMode = O_RDWR   | O_LARGEFILE;
      mapMode  = PROT_READ | PROT_WRITE;
      break;
    default:
      fprintf(stderr, "Invalid mode to mapFile; must be 'r' or 'w'\n");
      exit(1);
      break;
  }

  errno = 0;
  f = open(filename, openMode);
  if (errno) {
    fprintf(stderr, "Couldn't open() '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  fstat(f, &sb);
  if (errno) {
    fprintf(stderr, "Couldn't fstat() '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  *length = sb.st_size;

  ptr = mmap(0L, *length, mapMode, MMAPFLAGS, f, (off_t)0);
  if (errno) {
    fprintf(stderr, "Couldn't mmap() '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  close(f);

  return(ptr);
}



void
unmapFile(void *addr, uint64 length) {
#ifdef __sun
  //  This might work in general, but sun definitely needs the cast.
  //
  (void)munmap((caddr_t)addr, length);
#else
  (void)munmap(addr, length);
#endif
}





//  Copies all of srcFile to dstFile, returns the number of bytes written
//
off_t
copyFile(char *srcName, FILE *dstFile) {
  off_t  srcSize     = 0;
  off_t  bytesRemain = 0;
  off_t  bytesRead   = 0;
  int    bufferSize  = 1024 * 1024;
  char  *buffer      = 0L;
  FILE  *srcFile     = 0L;

  buffer     = (char *)malloc(sizeof(char) * bufferSize);
  if (buffer == 0L) {
    fprintf(stderr, "copyFile()-- Can't allocate buffer.\n");
    exit(1);
  }

  srcSize     = sizeOfFile(srcName);
  bytesRemain = srcSize;

  errno = 0;
  srcFile = fopen(srcName, "r");
  if (errno) {
    fprintf(stderr, "copyFile()-- failed to open the '%s' during merge: %s\n", srcName, strerror(errno));
    exit(1);
  }

  while (bytesRemain > 0) {

    errno = 0;

    if (bytesRemain > bufferSize)
      bytesRead = fread(buffer, sizeof(char), (size_t)bufferSize, srcFile);
    else
      bytesRead = fread(buffer, sizeof(char), (size_t)bytesRemain, srcFile);

    if (errno) {
      fprintf(stderr, "copyFile()-- Error reading source: %s\n", strerror(errno));
      exit(1);
    }

    if (bytesRead == 0) {
      fprintf(stderr, "copyFile()-- Short read (%d bytes) on source: %s\n", (int)bytesRead, strerror(errno));
      exit(1);
    }

    if (bytesRead > 0) {
      fwrite(buffer, sizeof(char), (size_t)bytesRead, dstFile);

      if (errno) {
        fprintf(stderr, "copyFile()-- Error writing %d bytes to destination: %s\n", (int)bytesRead, strerror(errno));
        exit(1);
      }
    }

    bytesRemain -= bytesRead;
  }

  fclose(srcFile);
  free(buffer);

  return(srcSize);
}





//  Takes a path to a file (that possibly doesn't exist) and returns
//  the number of MB (1048576 bytes) free in the directory of that
//  file.
//
uint32
freeDiskSpace(char *path) {
  char          *p, *t;
  struct statvfs dst;
  struct stat    fst;
  uint64         ret = 0;

  //  Stat the path; if it exists, we're golden.
  //
  if (stat(path, &fst) == 0) {
    if (statvfs(path, &dst) == -1) {
      perror("statvfs");
      exit(1);
    }
  } else {
    //  Doesn't exist.  Try to find the directory that the file goes into.
    //
    //  Copy the input path to a temporary string.  Strip off
    //  the last component (probably a file prefix, but it could also
    //  be a directory -- see below) and return the free space on
    //  that device.
    //
    p = (char *)malloc(sizeof(char) * (strlen(path) + 1));
    strcpy(p, path);
    t = strrchr(p, '/');

    if (t) {
      *t = 0;
    } else {
      p[0] = '.';
      p[1] = 0;
    }

    if (statvfs(p, &dst) == -1) {
      perror("statvfs");
      exit(1);
    }

    free(p);
  }

  ret   = dst.f_frsize;
  ret  *= dst.f_bavail;
  ret >>= 20;

  return((uint32)ret);
}






//  Split writes/reads into smaller pieces, check the result of each
//  piece.  Really needed by OSF1 (V5.1).
//
void
safeWrite(int filedes, const void *buffer, const char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  towrite  = 0;
  size_t  written  = 0;

  while (position < nbytes) {
    towrite = length;
    if (position + towrite > nbytes)
      towrite = nbytes - position;

    errno = 0;
    written = write(filedes, ((char *)buffer) + position, towrite);

    if ((errno) || (towrite != written)) {
      fprintf(stderr, "safeWrite()-- Write failure on %s: %s\n", desc, strerror(errno));
      fprintf(stderr, "safeWrite()-- Wanted to write "int64FMT" bytes, wrote "int64FMT".\n", (int64)towrite, (int64)written);
      exit(1);
    }

    position += written;
  }
}

int
safeRead(int filedes, const void *buffer, const char *desc, size_t nbytes) {
  size_t  position = 0;
  size_t  length   = 32 * 1024 * 1024;
  size_t  toread   = 0;
  size_t  written  = 0;  //  readen?
  int     failed   = 0;

  while (position < nbytes) {
    toread = length;
    if (position + toread > nbytes)
      toread = nbytes - position;

    errno = 0;
    written = read(filedes, ((char *)buffer) + position, toread);

    failed = errno;
#ifdef VERY_SAFE
    if (toread != written)
      failed = 1;
#endif

    if ((failed) && (errno != EINTR)) {
      fprintf(stderr, "safeRead()-- Read failure on %s: %s.\n", desc, strerror(errno));
      fprintf(stderr, "safeRead()-- Wanted to read "int64FMT" bytes, read "int64FMT".\n", (int64)toread, (int64)written);
      exit(1);
    }

    if (written == 0)
      break;

    position += written;
  }

  return(position);
}



void
closeFile(FILE *F, const char *path) {

  //  If we're given the path name, see if we need to pclose(),
  //  otherwise just fclose() the file.

  if ((path) &&
      ((strcmp(path + strlen(path) - 4, ".bz2") == 0) ||
       (strcmp(path + strlen(path) - 3, ".gz") == 0))) {
    pclose(F);
  } else {
    fclose(F);
  }
}

FILE*
openFile(const char *path, const char *mode) {
  FILE *F         = 0L;
  int   isBz      = 0;
  int   isGz      = 0;
  int   isRead    = 0;
  int   isWrite   = 0;
  int   isRW      = 1;
  char  cmd[1024] = { 0 };;

  //  Yes, one could make this significantly simpler by saving the
  //  compression command into a variable, instead of the isBz and
  //  isGz flags.  Maybe instead we should find a compression command
  //  that uses different flags.

  if (strcmp(path + strlen(path) - 4, ".bz2") == 0)
    isBz = 1;
  if (strcmp(path + strlen(path) - 3, ".gz") == 0)
    isGz = 1;

  if (strcmp(mode, "w") == 0) {
    isRead   = 0;
    isWrite  = 1;
    isRW     = 0;
  }
  if (strcmp(mode, "r") == 0) {
    isRead   = 1;
    isWrite  = 0;
    isRW     = 0;
  }

  if (isBz) {
    if        (isRead) {
      sprintf(cmd, "bzip2 -dc %s", path);
    } else if (isWrite) {
      sprintf(cmd, "bzip2 -9c > %s", path);
    } else {
      fprintf(stderr, "openFile()-- Error!  Requested mode '%s' unavailable for bzip2 file '%s'\n", mode, path);
      exit(1);
    }
  } else if (isGz) {
    if        (isRead) {
      sprintf(cmd, "gzip -dc %s", path);
    } else if (isWrite) {
      sprintf(cmd, "gzip -9c > %s", path);
    } else {
      fprintf(stderr, "openFile()-- Error!  Requested mode '%s' unavailable for gzip file '%s'\n", mode, path);
      exit(1);
    }
  } else {
    //  Must be a normal file!
  }


  if (cmd[0]) {
    errno = 0;
    F = popen(cmd, mode);
    //  popen doesn't reliably set errnoman
    //if (errno)
    //  fprintf(stderr, "openFile()--  Failed to open pipe '%s': %s\n", cmd, strerror(errno)), exit(1);
    if (F == 0L)
      fprintf(stderr, "openFile()--  Failed to open pipe '%s'\n", cmd), exit(1);
  } else {
    errno = 0;
    F = fopen(path, mode);
    if (errno)
      fprintf(stderr, "openFile()--  Failed to open '%s': %s\n", path, strerror(errno)), exit(1);
  }

  return(F);
}

