#include "stat.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

stat_s*
stat_onPath(const char *path, stat_s *sb) {
  errno = 0;

  if (sb == NULL) {
    sb = (stat_s *)malloc(sizeof(stat_s));
    if (errno) {
      fprintf(stderr, "stat_onPath()-- Can't allocate stat_s\n%s\n", strerror(errno));
      exit(1);
    }
  }

  if (stat(path, sb) != 0) {
    fprintf(stderr, "Couldn't stat() '%s'\n%s\n", path, strerror(errno));
    exit(1);
  }

  return(sb);
}


stat_s*
stat_onLink(const char *path, stat_s *sb) {
  errno = 0;

  if (sb == NULL) {
    sb = (stat_s *)malloc(sizeof(stat_s));
    if (errno) {
      fprintf(stderr, "stat_onLink()-- Can't allocate stat_s\n%s\n", strerror(errno));
      exit(1);
    }
  }

  if (lstat(path, sb) != 0) {
    fprintf(stderr, "Couldn't lstat() '%s'\n%s\n", path, strerror(errno));
    exit(1);
  }

  return(sb);
}


stat_s*
stat_onDescriptor(int file, stat_s *sb) {
  errno = 0;

  if (sb == NULL) {
    sb = (stat_s *)malloc(sizeof(stat_s));
    if (errno) {
      fprintf(stderr, "stat_onDescriptor()-- Can't allocate stat_s\n%s\n", strerror(errno));
      exit(1);
    }
  }

  if (fstat(file, sb) != 0) {
    fprintf(stderr, "Couldn't fstat() descriptor %d\n%s\n", file, strerror(errno));
    exit(1);
  }

  return(sb);
}


stat_s*
stat_onFile(FILE *F, stat_s *sb) {
  return(stat_onDescriptor(fileno(F), sb));
}


void      stat_free(stat_s *sb) {
  free(sb);
}




int       stat_fileIsPipe(stat_s *sb) {
  return(sb->st_mode & S_IFIFO);
}
int       stat_fileIsCharacterSpecial(stat_s *sb) {
  return(sb->st_mode & S_IFCHR);
}
int       stat_fileIsDirectory(stat_s *sb) {
  return(sb->st_mode & S_IFDIR);
}
int       stat_fileIsBlockSpecial(stat_s *sb) {
  return(sb->st_mode & S_IFBLK);
}
int       stat_fileIsRegular(stat_s *sb) {
  return(sb->st_mode & S_IFREG);
}
int       stat_fileIsSymbolicLink(stat_s *sb) {
  return(sb->st_mode & S_IFLNK);
}
int       stat_fileIsSocket(stat_s *sb) {
  return(sb->st_mode & S_IFSOCK);
}

#ifdef __FreeBSD__
int       stat_fileIsWhiteout(stat_s *sb) {
  return(sb->st_mode & S_IFWHT);
}
#endif

uid_t     stat_getUID(stat_s *sb) {
  return(sb->st_uid);
}
gid_t     stat_getGID(stat_s *sb) {
  return(sb->st_gid);
}


double    stat_getAccessTime(stat_s *sb) {
#ifdef __FreeBSD__
#ifndef _POSIX_SOURCE
  return(sb->st_atimespec.tv_sec + (double)sb->st_atimespec.tv_nsec * 1.0e-9);
#else
  return(sb->st_atime + (double)sb->st_atimensec * 1.0e-9);
#endif
#else
  return(sb->st_atime);
#endif
}
double    stat_getModificationTime(stat_s *sb) {
#ifdef __FreeBSD__
#ifndef _POSIX_SOURCE
  return(sb->st_mtimespec.tv_sec + (double)sb->st_mtimespec.tv_nsec * 1.0e-9);
#else
  return(sb->st_mtime + (double)sb->st_mtimensec * 1.0e-9);
#endif
#else
  return(sb->st_mtime);
#endif
}
double    stat_getStatusTime(stat_s *sb) {
#ifdef __FreeBSD__
#ifndef _POSIX_SOURCE
  return(sb->st_ctimespec.tv_sec + (double)sb->st_ctimespec.tv_nsec * 1.0e-9);
#else
  return(sb->st_ctime + (double)sb->st_ctimensec * 1.0e-9);
#endif
#else
  return(sb->st_ctime);
#endif
}


off_t     stat_getSize(stat_s *sb) {
  return(sb->st_size);
}





int
fileExists(char *path) {
  stat_s  s;
  int     r;

  r = stat(path, &s);

  //  XXX:  We should be more careful on handling errors here.

  return(r == 0);
}


off_t
sizeOfFile(char *path) {
  struct stat s;

  stat_onPath(path, &s);
  return(s.st_size);
}
