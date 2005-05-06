#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__CYGWIN__)
#include <sys/syslimits.h>
#endif




//  Useful types.
//
//  *MASK(x) is only defined for unsigned types, with x != 0 and less
//  than the datawidth.

#if defined(__alpha) || defined(_AIX)
#define TRUE64BIT
#endif

#ifdef TRUE64BIT

//  Generic 64-bit platform

typedef unsigned long       u64bit;
typedef unsigned int        u32bit;
typedef signed long         s64bit;
typedef signed int          s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitNUMBER(X) X ## LU
#define  u64bitZERO      (0x0000000000000000LU)
#define  u64bitONE       (0x0000000000000001LU)
#define  u64bitMAX       (0xffffffffffffffffLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMTW(X)   "%" #X "lu"
#define  u64bitFMT       "%lu"
#define  u64bitHEX       "0x%016lx"
#define  s64bitFMTW(X)   "%" #X "ld"
#define  s64bitFMT       "%ld"

#define  u32bitNUMBER(X) X ## U
#define  u32bitZERO      (0x00000000U)
#define  u32bitONE       (0x00000001U)
#define  u32bitMAX       (0xffffffffU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMTW(X)   "%" #X "u"
#define  u32bitFMT       "%u"
#define  u32bitHEX       "0x%08x"
#define  s32bitFMTW(X)   "%" #X "d"
#define  s32bitFMT       "%d"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))
#define  u16bitFMT       "%hd"

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#define  strtou32bit(N,O) (u32bit)strtoul(N, O, 0)
#define  strtou64bit(N,O) (u64bit)strtoul(N, O, 0)

#else

//  Generic 32-bit platform

typedef unsigned long long  u64bit;
typedef unsigned long       u32bit;
typedef signed long long    s64bit;
typedef signed long         s32bit;
typedef unsigned short      u16bit;
typedef signed short        s16bit;
typedef unsigned char       u8bit;
typedef signed char         s8bit;

#define  u64bitNUMBER(X) X ## LLU
#define  u64bitZERO      (0x0000000000000000LLU)
#define  u64bitONE       (0x0000000000000001LLU)
#define  u64bitMAX       (0xffffffffffffffffLLU)
#define  u64bitMASK(X)   ((~u64bitZERO) >> (64 - (X)))
#define  u64bitFMTW(X)   "%" #X "llu"
#define  u64bitFMT       "%llu"
#define  u64bitHEX       "0x%016llx"
#define  s64bitFMTW(X)   "%" #X "lld"
#define  s64bitFMT       "%lld"

#define  u32bitNUMBER(X) X ## LU
#define  u32bitZERO      (0x00000000LU)
#define  u32bitONE       (0x00000001LU)
#define  u32bitMAX       (0xffffffffLU)
#define  u32bitMASK(X)   ((~u32bitZERO) >> (32 - (X)))
#define  u32bitFMTW(X)   "%" #X "lu"
#define  u32bitFMT       "%lu"
#define  u32bitHEX       "0x%08lx"
#define  s32bitFMTW(X)   "%" #X "ld"
#define  s32bitFMT       "%ld"

#define  u16bitZERO      (0x0000)
#define  u16bitONE       (0x0001)
#define  u16bitMAX       (0xffff)
#define  u16bitMASK(X)   ((~u16bitZERO) >> (16 - (X)))
#define  u16bitFMT       "%hd"

#define  u8bitZERO       (0x00)
#define  u8bitONE        (0x01)
#define  u8bitMAX        (0xff)
#define  u8bitMASK(X)    ((~u8bitZERO) >> (8 - (X)))

#define  strtou32bit(N,O) (u32bit)strtoul(N, O, 0)
#define  strtou64bit(N,O) (u64bit)strtoull(N, O, 0)

#endif  //  TRUE64BIT



#ifdef __cplusplus
extern "C" {
#endif




////////////////////////////////////////
//
//  time
//
double  getTime(void);
void    write_rusage(FILE *F);



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


//  Takes a path to a file (that possibly doesn't exist) and returns
//  the number of MB (1048576 bytes) free in the directory of that
//  file.
//
u32bit freeDiskSpace(char *path);




////////////////////////////////////////
//
//  stat
//


//
//  Simple wrapper around stat(), lstat() and fstat().
//

typedef struct stat stat_s;

stat_s   *stat_onPath(const char *path, stat_s *sb);
stat_s   *stat_onLink(const char *path, stat_s *sb);
stat_s   *stat_onDescriptor(int file, stat_s *sb);
stat_s   *stat_onFile(FILE *F, stat_s *sb);

void      stat_free(stat_s *sb);

int       stat_fileIsPipe(stat_s *sb);
int       stat_fileIsCharacterSpecial(stat_s *sb);
int       stat_fileIsDirectory(stat_s *sb);
int       stat_fileIsBlockSpecial(stat_s *sb);
int       stat_fileIsRegular(stat_s *sb);
int       stat_fileIsSymbolicLink(stat_s *sb);
int       stat_fileIsSocket(stat_s *sb);
int       stat_fileIsWhiteout(stat_s *sb);

uid_t     stat_getUID(stat_s *sb);
gid_t     stat_getGID(stat_s *sb);

double    stat_getAccessTime(stat_s *sb);
double    stat_getModificationTime(stat_s *sb);
double    stat_getStatusTime(stat_s *sb);

off_t     stat_getSize(stat_s *sb);


//  Convenience functions
//
int   fileExists(const char *path);
off_t sizeOfFile(const char *path);


////////////////////////////////////////
//
//  memory
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

//  Get a new handle, release a used one.
//
void   *pallochandle(void);
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
void    pdumppalloc(void);
void    psetdebug(int on);


////////////////////////////////////////
//
//  md5
//


typedef struct {
  u64bit  a;
  u64bit  b;
  u32bit  i;  //  the iid, used in leaff
#ifdef __APPLE__
  u32bit  pad;
#endif
} md5_s;

#define MD5_BUFFER_SIZE   32*1024

typedef struct {
  u64bit           a;
  u64bit           b;
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
md5_s  *md5_string(md5_s *m, char *s, u32bit l);


//  Computes an md5 checksum piece by piece.
//
//  If m is NULL, a new md5_increment_s is allocated and returned.
//
md5_increment_s  *md5_increment_char(md5_increment_s *m, char s);
md5_increment_s  *md5_increment_block(md5_increment_s *m, char *s, u32bit l);
void              md5_increment_finalize(md5_increment_s *m);
void              md5_increment_destroy(md5_increment_s *m);


////////////////////////////////////////
//
//  Matsumoto and Nichimura's Mersenne Twister pseudo random number
//  generator.  The struct and functions are defined in external/mt19937ar.[ch]
//
typedef struct mt mt_s;

mt_s          *mtInit(u32bit s);
mt_s          *mtInitArray(u32bit *init_key, u32bit key_length);
u32bit         mtRandom32(mt_s *mt);

//  A u64bit random number
//
#define        mtRandom64(MT) ( (((u64bit)mtRandom32(MT)) << 32) | (u64bit)mtRandom32(MT) )

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
//  Kaz Kylheku <kaz@ashi.footprints.net> library.
//
#include "kazlib/dict.h"
#include "kazlib/except.h"
#include "kazlib/hash.h"
#include "kazlib/list.h"
#include "kazlib/sfx.h"


////////////////////////////////////////
//
//  perl's chomp is pretty nice
//
#define chomp(S) { char *t=S; while (*t) t++; t--; while (isspace(*t)) *t--=0; }

#ifdef __cplusplus
}
#endif

#endif  //  UTIL_H
