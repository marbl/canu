#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#ifndef _AIX
#include <sys/param.h>
#include <sys/mount.h>
#endif

#ifdef _AIX
#include <sys/statfs.h>
#endif

#ifdef linux
#include <sys/vfs.h>
#endif

#include "bri.h"


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
mapFile(const char *filename, size_t *length, char mode) {
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

  ptr = mmap(0L, *length, mapMode, MMAPFLAGS, f, 0);
  if (errno) {
    fprintf(stderr, "Couldn't mmap() '%s'\n%s\n", filename, strerror(errno));
    exit(1);
  }

  close(f);

  return(ptr);
}



void
unmapFile(void *addr, size_t length) {
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
    fprintf(stderr, "merStreamFileBuilder::build()-- failed to open the '%s' during merge: %s\n", srcName, strerror(errno));
    exit(1);
  }

  while (bytesRemain > 0) {

    errno = 0;

    if (bytesRemain > bufferSize)
      bytesRead = fread(buffer, sizeof(char), bufferSize, srcFile);
    else
      bytesRead = fread(buffer, sizeof(char), bytesRemain, srcFile);

    if (errno) {
      fprintf(stderr, "copyFile()-- Error reading source: %s\n", strerror(errno));
      exit(1);
    }

    if (bytesRead == 0) {
      fprintf(stderr, "copyFile()-- Short read (%d bytes) on source: %s\n", (int)bytesRead, strerror(errno));
      exit(1);
    }

    if (bytesRead > 0) {
      fwrite(buffer, sizeof(char), bytesRead, dstFile);

      if (errno) {
        fprintf(stderr, "copyFile()-- Error writing %d bytes to destination: %s\n", (int)bytesRead, strerror(errno));
        exit(1);
      }
    }

    bytesRemain -= bytesRead;
  }

  free(buffer);

  return(srcSize);
}





//  Takes a path to a file (that possibly doesn't exist) and returns
//  the number of MB (1048576 bytes) free in the directory of that
//  file.
//
u32bit
freeDiskSpace(char *path) {
  char          *p, *t;
  struct statfs  dst;
  struct stat    fst;
  u64bit         ret = 0;

  //  Stat the path; if it exists, we're golden.
  //
  if (stat(path, &fst) == 0) {
    if (statfs(path, &dst) == -1) {
      perror("statfs");
      exit(-1);
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

    if (statfs(p, &dst) == -1) {
      perror("statfs");
      exit(-1);
    }

    free(p);
  }

  ret   = dst.f_bsize;
  ret  *= dst.f_bavail;
  ret >>= 20;

  return((u32bit)ret);
}
