#ifndef FILE_H
#define FILE_H

#include <sys/types.h>

#ifdef __cplusplus
extern "C" {
#endif

//  Useful routines for dealing with the existence of files

int   fileExists(char *path);
off_t sizeOfFile(char *path);
int   isHuman(FILE *F);

//  Handles mmap() of files.  Write is not tested -- in particluar,
//  the test main() in mmap.c fails.
//
void*
mapFile(char     *filename,
        size_t   *length,
        char      mode);

void
unmapFile(void   *addr,
          size_t  length);

#ifdef __cplusplus
}
#endif

#endif  //  FILE_H

