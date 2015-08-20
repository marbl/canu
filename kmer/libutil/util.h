
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
 *    Brian P. Walenz from 2004-APR-21 to 2004-OCT-10
 *      are Copyright 2004 Brian P. Walenz, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2004-MAY-06 to 2004-AUG-02
 *      are Copyright 2004 Applera Corporation, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2005-FEB-07 to 2014-APR-11
 *      are Copyright 2005-2008,2010-2012,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef UTIL_H
#define UTIL_H

//  ISO C99 says that to get INT32_MAX et al, these must be defined. (7.18.2, 7.18.4, 7.8.1)
#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>

#include <inttypes.h>

//  Useful types.
//
//  *MASK(x) is only defined for unsigned types, with x != 0 and less
//  than the datawidth.

typedef uint64_t       uint64;
typedef uint32_t       uint32;
typedef uint16_t       uint16;
typedef uint8_t        uint8;

typedef int64_t         int64;
typedef int32_t         int32;
typedef int16_t         int16;
typedef int8_t          int8;


#if defined(__alpha) || defined(_AIX) || defined(__LP64__) || defined(_LP64)
#define TRUE64BIT
#define  uint64NUMBER(X) X ## LU
#define  uint32NUMBER(X) X ## U
#else
#define  uint64NUMBER(X) X ## LLU
#define  uint32NUMBER(X) X ## LU
#endif


#define  sizetFMT        "%zd"

#define  uint64ZERO      uint64NUMBER(0x0000000000000000)
#define  uint64ONE       uint64NUMBER(0x0000000000000001)
#define  uint64MAX       uint64NUMBER(0xffffffffffffffff)
#define  uint64MASK(X)   ((~uint64ZERO) >> (64 - (X)))
#define  uint64FMTW(X)   "%" #X PRIu64
#define  uint64FMT       "%"PRIu64
#define  uint64HEX       "0x%016"PRIx64
#define  int64FMTW(X)    "%" #X PRId64
#define  int64FMT        "%"PRId64

#define  uint32ZERO      uint32NUMBER(0x00000000)
#define  uint32ONE       uint32NUMBER(0x00000001)
#define  uint32MAX       uint32NUMBER(0xffffffff)
#define  uint32MASK(X)   ((~uint32ZERO) >> (32 - (X)))
#define  uint32FMTW(X)   "%" #X PRIu32
#define  uint32FMT       "%"PRIu32
#define  uint32HEX       "0x%08"PRIx32
#define  int32FMTW(X)    "%" #X PRId32
#define  int32FMT        "%"PRId32

#define  uint16ZERO      (0x0000)
#define  uint16ONE       (0x0001)
#define  uint16MAX       (0xffff)
#define  uint16MASK(X)   ((~uint16ZERO) >> (16 - (X)))
#define  uint16FMTW(X)   "%" #X PRIu16
#define  uint16FMT       "%"PRIu16

#define  uint8ZERO       (0x00)
#define  uint8ONE        (0x01)
#define  uint8MAX        (0xff)
#define  uint8MASK(X)    ((~uint8ZERO) >> (8 - (X)))

#define  strtouint32(N,O) (uint32)strtoul(N, O, 10)
#define  strtouint64(N,O) (uint64)strtoul(N, O, 10)




#ifdef __cplusplus
extern "C" {
#endif




////////////////////////////////////////
//
//  time
//
double  getTime(void);



////////////////////////////////////////
//
//  file
//

//  Create the O_LARGEFILE type for open(), if it doesn't already
//  exist (FreeBSD, Tru64).  We assume that by including the stuff
//  needed for open(2) we'll get any definition of O_LARGEFILE.
//
#ifndef O_LARGEFILE
#define O_LARGEFILE    0
#endif


uint64   getProcessSizeCurrent(void);
uint64   getProcessSizeLimit(void);


//  Useful routines for dealing with the existence of files

int   isHuman(FILE *F);

//  Handles mmap() of files.  Write is not tested -- in particluar,
//  the test main() in mmap.c fails.
//
void*
mapFile(const char *filename,
        uint64     *length,
        char        mode);

void
unmapFile(void     *addr,
          uint64    length);



//  Creates a hidden temporary file.  If path is given, the temporary
//  file is created in that directory.  The temoprary file is unlinked
//  after it is created, so once you close the file, it's gone.
//
FILE *makeTempFile(char *path);


//  Copies all of srcFile to dstFile, returns the number of bytes written
//
off_t copyFile(char *srcName, FILE *dstFile);


//  Takes a path to a file (that possibly doesn't exist) and returns
//  the number of MB (1048576 bytes) free in the directory of that
//  file.
//
uint32 freeDiskSpace(char *path);

//  Safer read(2) and write(2).
//
void   safeWrite(int filedes, const void *buffer, const char *desc, size_t nbytes);
int    safeRead(int filedes, const void *buffer, const char *desc, size_t nbytes);



////////////////////////////////////////
//
int       fileExists(const char *path);
off_t     sizeOfFile(const char *path);
uint64    timeOfFile(const char *path);

//  Open a file, read/write, using compression based on the file name
//
FILE *openFile(const char *path, const char *mode);
void  closeFile(FILE *F, const char *path);

////////////////////////////////////////
//
void    *memdup(const void *orig, size_t size);


////////////////////////////////////////
//
//  Pac-Man's memory allocator.
//
//  Grabs big chunks of memory, then gives out little pieces.  You can
//  only free ALL memory, not single blocks.
//
//  This is useful when one needs to malloc() tens of millions of
//  things, at which point the overhead of finding a free block is
//  large.
//
void   *palloc(size_t size);
void    pfree(void);

//  A thread-safe(r) implementation just forces the user to use a
//  handle.  This also lets us use palloc() for collections of things
//  -- e.g., twice in a program.  If you don't give a handle, the
//  default one is used.
//
void   *palloc2(size_t size, void *handle);
void    pfree2(void *handle);

//  Get a new handle, release a used one.  The size is the same
//  as for psetblocksize().
//
void   *pallochandle(size_t size);
void    pfreehandle(void *handle);

//  The block size can only be changed before the first call to
//  palloc().  Calling psetblocksize() after that has no effect.
//
void    psetblocksize(size_t size);
size_t  pgetblocksize(void);

//  Not generally useful - just dumps the allocated blocks to stdout.
//  Uses internal structures, and used in the test routine.
//
//  psetdebug() enables reporting of allocations.
//
void    pdumppalloc(void *handle);
void    psetdebug(int on);


////////////////////////////////////////
//
//  md5
//


typedef struct {
  uint64  a;
  uint64  b;
  uint32  i;    //  the iid, used in leaff
  uint32  pad;  //  keep us size compatible between 32- and 64-bit machines.
} md5_s;

#define MD5_BUFFER_SIZE   32*1024

typedef struct {
  uint64           a;
  uint64           b;
  void            *context;
  int              bufferPos;
  unsigned char    buffer[MD5_BUFFER_SIZE];
} md5_increment_s;


//  Returns -1, 0, 1 depending on if a <, ==, > b.  Suitable for
//  qsort().
//
int     md5_compare(void const *a, void const *b);


//  Converts an md5_s into a character string.  s must be at least
//  33 bytes long.
//
char   *md5_toascii(md5_s *m, char *s);


//  Computes the md5 checksum on the string s.
//
md5_s  *md5_string(md5_s *m, char *s, uint32 l);


//  Computes an md5 checksum piece by piece.
//
//  If m is NULL, a new md5_increment_s is allocated and returned.
//
md5_increment_s  *md5_increment_char(md5_increment_s *m, char s);
md5_increment_s  *md5_increment_block(md5_increment_s *m, char *s, uint32 l);
void              md5_increment_finalize(md5_increment_s *m);
void              md5_increment_destroy(md5_increment_s *m);


////////////////////////////////////////
//
//  Matsumoto and Nichimura's Mersenne Twister pseudo random number
//  generator.  The struct and functions are defined in external/mt19937ar.[ch]
//
typedef struct mtctx mt_s;

mt_s          *mtInit(uint32 s);
mt_s          *mtInitArray(uint32 *init_key, uint32 key_length);
uint32         mtRandom32(mt_s *mt);

//  A uint64 random number
//
#define        mtRandom64(MT) ( (((uint64)mtRandom32(MT)) << 32) | (uint64)mtRandom32(MT) )

//  Real valued randomness
//    mtRandomRealOpen()    -- on [0,1) real interval
//    mtRandomRealClosed()  -- on [0,1] real interval
//    mrRandomRealOpen53()  -- on [0,1) real interval, using 53 bits
//
//  "These real versions are due to Isaku Wada, 2002/01/09 added" and were taken from
//  the mt19937ar.c distribution (but they had actual functions, not macros)
//
//  They also had
//    random number in (0,1) as (mtRandom32() + 0.5) * (1.0 / 4294967296.0)
//
#define        mtRandomRealOpen(MT)   ( (double)mtRandom32(MT) * (1.0 / 4294967296.0) )
#define        mtRandomRealClosed(MT) ( (double)mtRandom32(MT) * (1.0 / 4294967295.0) )
#define        mtRandomRealOpen53(MT) ( ((mtRandom32(MT) >> 5) * 67108864.0 + (mtRandom32(MT) >> 6)) * (1.0 / 9007199254740992.0) )

//  returns a random number with gaussian distribution, mean of zero and std.dev. of 1
//
double  mtRandomGaussian(mt_s *mt);


////////////////////////////////////////
//
//  FreeBSD's multithreaded qsort.
//
void
qsort_mt(void *a,
         size_t n,
         size_t es,
         int (*cmp)(const void *, const void *),
         int maxthreads,
         int forkelem);

//#define qsort(A, N, ES, CMP)  qsort_mt((A), (N), (ES), (CMP), 4, 64 * 1024)



////////////////////////////////////////
//
//  perl's chomp is pretty nice
//
#ifndef chomp
#define chomp(S)    { char *t=S; while (*t) t++; t--; while (isspace(*t)) { *t--=0; } }
#define chompL(S,L) { char *t=S; while (*t) t++; t--; while (isspace(*t)) { *t--=0; L--; } }
#endif

#ifndef munch
#define munch(S)    { while (*(S) &&  isspace(*(S))) (S)++; }
#endif

#ifndef crunch
#define crunch(S)   { while (*(S) && !isspace(*(S))) (S)++; }
#endif


#ifndef MIN
#define  MIN(x,y)        (((x) > (y)) ? (y) : (x))
#endif

#ifndef MAX
#define  MAX(x,y)        (((x) < (y)) ? (y) : (x))
#endif


#ifdef __cplusplus
}
#endif

#endif  //  UTIL_H
