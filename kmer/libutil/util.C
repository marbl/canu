#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>


double
getTime(void) {
  struct timeval  tp;
  gettimeofday(&tp, NULL);
  return(tp.tv_sec + (double)tp.tv_usec / 1000000.0);
}


void
write_rusage(FILE *F) {
    struct utsname un;
    struct rusage  ru;

    fprintf(F, "uname()---------------------------------\n");
    errno = 0;
    if (uname(&un) == -1) {
      fprintf(F, "uname() call failed: %s\n", strerror(errno));
    } else {
      fprintf(F, "sysname:     %s\n", un.sysname);
      fprintf(F, "nodename:    %s\n", un.nodename);
      fprintf(F, "release:     %s\n", un.release);
      fprintf(F, "version:     %s\n", un.version);
      fprintf(F, "machine:     %s\n", un.machine);
    }

    fprintf(F, "\n");

    fprintf(F, "rusage()--------------------------------\n");
    errno = 0;
    if (getrusage(RUSAGE_SELF, &ru) == -1) {
      fprintf(F, "getrusage() call failed: %s\n", strerror(errno));
    } else {
      fprintf(F, "userTime:    %f\n", ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000.0);
      fprintf(F, "systemTime:  %f\n", ru.ru_stime.tv_sec + (double)ru.ru_stime.tv_usec / 1000000.0);
      fprintf(F, "maxrss:      %ld\n", ru.ru_maxrss);
      fprintf(F, "ixrss:       %ld\n", ru.ru_ixrss);
      fprintf(F, "idrss:       %ld\n", ru.ru_idrss);
      fprintf(F, "isrss:       %ld\n", ru.ru_isrss);
      fprintf(F, "minflt:      %ld\n", ru.ru_minflt);
      fprintf(F, "majflt:      %ld\n", ru.ru_majflt);
      fprintf(F, "nswap:       %ld\n", ru.ru_nswap);
      fprintf(F, "inblock:     %ld\n", ru.ru_inblock);
      fprintf(F, "oublock:     %ld\n", ru.ru_oublock);
      fprintf(F, "msgsnd:      %ld\n", ru.ru_msgsnd);
      fprintf(F, "msgrcv:      %ld\n", ru.ru_msgrcv);
      fprintf(F, "nsignals:    %ld\n", ru.ru_nsignals);
      fprintf(F, "nvcsw:       %ld\n", ru.ru_nvcsw);
      fprintf(F, "nivcsw:      %ld\n", ru.ru_nivcsw);
    }

    fprintf(F, "\n");
}

#if 0

#ifdef __alpha
#define MMAPFLAGS    (MAP_FILE | MAP_VARIABLE | MAP_SHARED)
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

//  See malloc(3) for details.
//
unsigned long __sbrk_override = 1;


void*
mapFile(char *filename, size_t &length, bool verbose) {
  void   *ptr;

  errno = 0;
  int f = open(filename, O_RDONLY);
  if ((f < 0) || (errno)) {
    if (verbose) {
      fprintf(stderr, "Couldn't open '%s'\n", filename);
      perror("open");
    }
    return(0L);
  }

  errno = 0;
  struct stat  sb;
  fstat(f, &sb);
  if (errno) {
    close(f);
    if (verbose) {
      fprintf(stderr, "Couldn't stat '%s'\n", filename);
      perror("fstat\n");
    }
    return(0L);
  }

  length = sb.st_size;

  errno = 0;
  ptr = mmap(0L, length, PROT_READ, MMAPFLAGS, f, 0);
  if (errno) {
    close(f);
    if (verbose) {
      fprintf(stderr, "Couldn't map '%s'\n", filename);
      perror("mmap");
    }
    return(0L);
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


void*
mapFileForWrite(char *filename, size_t &length, bool verbose) {
  void   *ptr;

  errno = 0;
  int f = open(filename, O_RDWR);
  if ((f < 0) || (errno)) {
    if (verbose) {
      fprintf(stderr, "Couldn't open '%s'\n", filename);
      perror("open");
    }
    return(0L);
  }

  struct stat  sb;
  errno = 0;
  fstat(f, &sb);
  if (errno) {
    close(f);
    if (verbose) {
      fprintf(stderr, "Couldn't stat '%s'\n", filename);
      perror("fstat\n");
    }
    return(0L);
  }

  length = sb.st_size;

  errno = 0;
  ptr = mmap(0L, length, PROT_READ | PROT_WRITE, MMAPFLAGS, f, 0);
  if (errno) {
    close(f);
    if (verbose) {
      fprintf(stderr, "Couldn't map '%s'\n", filename);
      perror("mmap");
    }
    return(0L);
  }

  close(f);

  return(ptr);
}

#endif
