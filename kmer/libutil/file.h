#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <sys/types.h>

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__CYGWIN__)
#include <sys/syslimits.h>
#endif

//  Create the O_LARGEFILE type for open(), if it doesn't already
//  exist (FreeBSD, Tru64).  We assume that by including the stuff
//  needed for open(2) we'll get any definition of O_LARGEFILE.
//
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifndef O_LARGEFILE
#define O_LARGEFILE    0
#endif


#ifdef __cplusplus
extern "C" {
#endif

//  Useful routines for dealing with the existence of files

int   isHuman(FILE *F);

//  Handles mmap() of files.  Write is not tested -- in particluar,
//  the test main() in mmap.c fails.
//
void*
mapFile(const char *filename,
        size_t     *length,
        char        mode);

void
unmapFile(void     *addr,
          size_t    length);



//  Creates a hidden temporary file.  If path is given, the temporary
//  file is created in that directory.  The temoprary file is unlinked
//  after it is created, so once you close the file, it's gone.
//
FILE *makeTempFile(char *path);


//  Copies all of srcFile to dstFile, returns the number of bytes written
//
off_t copyFile(char *srcName, FILE *dstFile);




#ifdef __cplusplus
}
#endif

#endif  //  FILE_H

