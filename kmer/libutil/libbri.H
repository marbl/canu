#ifndef LIBBRI_H
#define LIBBRI_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////
//
//  Global DNA encoding arrays.
//
extern const unsigned char  compressSymbol[256];
extern const unsigned char  validSymbol[256];
extern const unsigned char  decompressSymbol[256];
extern const unsigned char  complementSymbol[256];
extern const unsigned char  validCompressedSymbol[256];


//  Create the O_LARGEFILE type for open(), if it doesn't already
//  exist (FreeBSD)
//
#ifdef __FreeBSD__
#define O_LARGEFILE    0
#endif

#ifdef __alpha
#define O_LARGEFILE    0
#endif

#ifdef __apple
void
nanosleep(struct timespec ts, struct timespec tr) {
  usleep(ts.tv_nsec >> 10);
}
#endif







#include "libbritypes.H"


//  XXX:  Put this someplace better
//
//  Inplace reverse-complement the sequence.  A pointer
//  the the character string is returned.
//
inline
unsigned char *
reverseComplementSequence(unsigned char *seq, u32bit seqlen) {
  unsigned char  *s = seq;
  unsigned char  *e = seq + seqlen - 1;
  unsigned char   t;
  u32bit          c = seqlen / 2;

  while (c--) {
    t = complementSymbol[*s];
    *(s++) = complementSymbol[*e];
    *(e--) = t;
  }

  if (s == e)
    *s = complementSymbol[*s];

  return(seq);
}




////////////////////////////////////////////////////////////////////////////////
//
//  Miscelaneous system utility functions
//
double getTime(void);
void   write_rusage(FILE *F);

#if 0
void  *mapFile(char *filename, size_t &length, bool verbose=false);
void   unmapFile(void *addr, size_t length);
void  *mapFileForWrite(char *filename, size_t &length, bool verbose=false);
#endif

inline
bool
fileExists(char *path) {
  struct stat s;

  return(stat(path, &s) == 0);
}

inline
off_t
sizeOfFile(char *path) {
  struct stat s;

  errno = 0;
  if (stat(path, &s) != 0) {
    fprintf(stderr, "Couldn't stat '%s': %s\n", path, strerror(errno));
    exit(1);
  }

  return(s.st_size);
}

#if 0
//  Linux uses _fileno
inline
bool
isHuman(FILE *F) {
  return(isatty(fileno(F)));
}
#endif


////////////////////////////////////////////////////////////////////////////////
//
//  Bit operations
//
//
//  countNumberOfSetBits:  returns the number of set bits.
//  The super-slow way to do this is:
//    u32bit r = u32bitZERO;
//    while (x) {
//      r++;
//      x >>= 1;
//    }
//    return(r);
//  But we have a super-fast way.
//  Credit:  http://remus.rutgers.edu/~rhoads/Code/bitcount3.c
//
//
//  logBaseTwo:  Returns the 'index' of the highest bit set.
//
//
//  reverseBits:  Reverse the bits in a word.
//  Credit:  The BSD fortunes file.
//
inline
u32bit
countNumberOfSetBits(u32bit x) {
  x -= (0xaaaaaaaa & x) >> 1;
  x  = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x += x >> 4;
  x &= 0x0f0f0f0f;
  x += x >> 8;
  x += x >> 16;
  x &= 0x000000ff;

  return(x);
}

inline
u32bit
logBaseTwo(u32bit x) {

  //  Step one, set all the bits below the highest bit set.
  //  Credit:  Ross Lippert
  //  
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;

  //  Step two, count the number of bits set.
  //
  return(countNumberOfSetBits(x));
}

inline
u32bit
logBaseTwo(u64bit x) {
  u32bit r = 64;
  u64bit m = u64bitONE << 63;

  while ((x & m) == 0) {
    r--;
    m = u64bitONE << (r-1);
  }

  return(r);
}

inline
u32bit
reverseBits(u32bit x) {
  x = ((x >>  1) & 0x55555555) | ((x <<  1) & 0xaaaaaaaa);
  x = ((x >>  2) & 0x33333333) | ((x <<  2) & 0xcccccccc);
  x = ((x >>  4) & 0x0f0f0f0f) | ((x <<  4) & 0xf0f0f0f0);
  x = ((x >>  8) & 0x00ff00ff) | ((x <<  8) & 0xff00ff00);
  x = ((x >> 16) & 0x0000ffff) | ((x << 16) & 0xffff0000);

  return(x);
}


////////////////////////////////////////////////////////////////////////////////
//
//  Bit packing operations
//
#include "bit-packing.H"



////////////////////////////////////////////////////////////////////////////////
//
//  All the rest
//
#include "bitPackedFile.H"
#include "fastastream.H"
#include "merstream.H"
#include "speedCounter.H"
#include "splittowords.H"

#endif  //  LIBBRI_H
