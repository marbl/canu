
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

//  This is the global include file that all C files in the AS
//  subsystem should include.

#ifndef AS_GLOBAL_H
#define AS_GLOBAL_H

static const char *rcsid_AS_GLOBAL_H = "$Id: AS_global.h,v 1.50 2011-09-06 01:11:56 mkotelbajcvi Exp $";

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <exception>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

#include "AS_UTL_alloc.h"

#ifndef TRUE
  #define TRUE true
#endif
#ifndef FALSE
  #define FALSE false
#endif

#ifndef MIN
  #define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
  #define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef _AIX
  typedef int8_t  int8;
  typedef int16_t int16;
  typedef int32_t int32;
  typedef int64_t int64;
#endif

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef void* PtrT;

// Format for wildcard in *scanf functions
#define F_I           "%*"

// Pointers
#define F_P           "%p"
#define F_PP           "p"

// Characters
#define F_C           "%c"
#define F_CP           "c"
#define F_CI     F_I""F_CP

// Strings
#define F_STR         "%s"
#define F_STRP         "s"
#define F_STRI F_I""F_STRP

#if ULONG_MAX == 0xffffffff
  // 32-bit architecture

  #define TRUE32BIT

  typedef uintptr_t INTPTR;

  #ifndef UINT64_MAX
    #define INT64_MAX  LLONG_MAX
    #define UINT64_MAX ULLONG_MAX
  #endif

  // Integers
  #define F_S16         "%d"
  #define F_S16P         "d"
  #define F_S16I F_I""F_S16P
  #define F_U16         "%u"
  #define F_U16P         "u"
  #define F_U16I F_I""F_U16P
  #define F_S32         "%d"
  #define F_S32P         "d"
  #define F_S32I F_I""F_S32P
  #define F_U32         "%u"
  #define F_U32P         "u"
  #define F_U32I F_I""F_U32P
  #define F_S64       "%lld"
  #define F_S64P       "lld"
  #define F_S64I F_I""F_S64P
  #define F_U64       "%llu"
  #define F_U64P       "llu"
  #define F_U64I F_I""F_U64P
  #define F_X64       "%llx"
  #define F_X64P       "llx"
  #define F_X64I F_I""F_X64P

  // Floating points
  #define F_F32        "%f"
  #define F_F32P        "f"
  #define F_F32I F_I""F_32P
  #define F_F64       F_F32
  #define F_F64P     F_F32P
  #define F_F64I     F_F32I
  
  // Standard typedefs
  #define F_SIZE_T     F_U32
  #define F_SIZE_TP   F_U32P
  #define F_SIZE_TI   F_U32I

  #define F_TIME_T           "%ld"
  #define F_TIME_TP           "ld"
  #define F_TIME_TI F_I""F_TIME_TP

  #define F_PID_T      F_S32
  #define F_PID_TP    F_S32P
  #define F_PID_TI    F_S32I

  #if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 32)
    // off_t is 32-bit
    #error I do not support 32-bit off_t.
    #define F_OFF_T          "%ld"
    #define F_OFF_TP          "ld"
    #define F_OFF_TI F_I""F_OFF_TP
  #endif

  // off_t is 64-bit

  #define F_OFF_T   F_S64
  #define F_OFF_TP F_S64P
  #define F_OFF_TI F_S64I

  #define FILEID_MASK       0xffff000000000000ULL
  #define FILEOFFSET_MASK   0x0000ffffffffffffULL
  #define LOCALE_OFFSET            10000000000ULL  // 10^10 > 2^32

  //========== SIMULATOR-specific
  // Fix initial creation time & uid in celsim to support regression testing
  // time stamp of batch message in first dros file...
  #define CREATION_TIME      915170460L
  // UID of first frag in first dros input file...
  #define FIRST_UID     17000001585819ULL

  //========== CGB-specific
  #define CGB_MULTIPLIER               1000

  //========== CGW-specific
  // Used in cgw: 1024^3
  #define MAX_SEQUENCEDB_SIZE        1073741824UL
  // Used in cgw: 256 * 1024^2
  #define MAX_SEQUENCEDB_CACHE_SIZE   268435456UL

#endif  // 32-bit architecture




#if ULONG_MAX == 0xffffffffffffffff
  // 64-bit architecture

  #define TRUE64BIT

  typedef uintptr_t INTPTR;

  #ifndef UINT64_MAX
    #define INT64_MAX  LONG_MAX
    #define UINT64_MAX ULONG_MAX
  #endif

  // Integers
  #define F_S16        "%hd"
  #define F_S16P        "hd"
  #define F_S16I F_I""F_S16P
  #define F_U16        "%hu"
  #define F_U16P        "hu"
  #define F_U16I F_I""F_U16P
  #define F_S32         "%d"
  #define F_S32P         "d"
  #define F_S32I F_I""F_S32P
  #define F_U32         "%u"
  #define F_U32P         "u"
  #define F_U32I F_I""F_U32P
  #define F_S64        "%ld"
  #define F_S64P        "ld"
  #define F_S64I  F_I""F_64P
  #define F_U64        "%lu"
  #define F_U64P        "lu"
  #define F_U64I F_I""F_U64P
  #define F_X64        "%lx"
  #define F_X64P        "lx"
  #define F_X64I F_I""F_X64P
  
  // Floating points
  #define F_F32        "%f"
  #define F_F32P        "f"
  #define F_F32I F_I""F_32P
  #define F_F64       "%lf"
  #define F_F64P       "lf"
  #define F_F64I  F_I""F64P
  
  // Standard typedefs
  #define F_SIZE_T           F_U64
  #define F_SIZE_TP         F_U64P
  #define F_SIZE_TI F_I""F_SIZE_TP

  #ifdef _AIX
    #define F_TIME_T   F_S64
    #define F_TIME_TP F_S64P
    #define F_TIME_TI F_S64I

    #define F_PID_T    F_S64
    #define F_PID_TP  F_S64P
    #define F_PID_TI  F_S64I

    #ifdef _LARGE_FILES
      #define F_OFF_T         "%lld"
      #define F_OFF_TP         "lld"
      #define F_OFF_TI F_I""F_OFF_TP
    #else
      #define F_OFF_T   F_S64
      #define F_OFF_TP F_S64P
      #define F_OFF_TI F_S64I
    #endif
  #else
    // these are valid for __alpha, perhaps not for others...
    #define F_TIME_T   F_S64
    #define F_TIME_TP F_S64P
    #define F_TIME_TI F_S64I

    #define F_PID_T    F_S32
    #define F_PID_TP  F_S32P
    #define F_PID_TI  F_S32I

    #ifdef _KERNEL
      #define F_OFF_T   F_U64
      #define F_OFF_TP F_U64P
      #define F_OFF_TI F_U64I
    #else
      #define F_OFF_T   F_S64
      #define F_OFF_TP F_S64P
      #define F_OFF_TI F_S64I
    #endif
  #endif

  #define FILEID_MASK       0xffff000000000000UL
  #define FILEOFFSET_MASK   0x0000ffffffffffffUL
  #define LOCALE_OFFSET            10000000000UL  // 10^10 > 2^32

  //========== SIMULATOR-specific
  // Fix initial creation time & uid in celsim to support regression testing
  // time stamp of batch message in first dros file...
  #define CREATION_TIME      915170460
  // UID of first frag in first dros input file...
  #define FIRST_UID     17000001585819UL

  //========== CGB-specific
  #define CGB_MULTIPLIER            1000

  //========== CGW-specific
  // Used in cgw: 2 * 1024^3
  #define MAX_SEQUENCEDB_SIZE           2147483648ul

  // Used in cgw: 2 * 256 * 1024^2
  #define MAX_SEQUENCEDB_CACHE_SIZE      536870912ul

#endif  // 64-bit architecture

// Common to 32-bit and 64-bit architechture


#ifdef __alpha
//  BPW hates doing this, but, well, its needed on our OSF1 V5.1 box.
int   fseeko(FILE *stream, off_t offset, int whence );
off_t ftello(FILE *stream );
#endif

//  OSF1 doesn't know UINT32_MAX or INT32_MAX or INT32_MIN.
#ifndef UINT32_MAX
  #define UINT32_MAX UINT_MAX
#endif
#ifndef INT32_MAX
  #define INT32_MAX INT_MAX
#endif
#ifndef INT32_MIN
  #define INT32_MIN INT_MIN
#endif


//  perl's chomp is pretty nice
//  Not a great place to put this, but it's getting used all over.
#ifndef chomp
#define chomp(S)  { char *t=(S); while (*t) t++; t--; while (t >= S && isspace(*t)) *t--=0; }
#endif

#ifndef munch
#define munch(S)  { while (*(S) &&  isspace(*(S))) (S)++; }
#endif

#ifndef crunch
#define crunch(S) { while (*(S) && !isspace(*(S))) (S)++; }
#endif


#define CGB_INVALID_CUTOFF           -12.0f
// A threshold value for Gene^s coverage statistic. Values BELOW this value
// have never been known to associated with unitigs with fragments that ARE
// not contiguous in the genome. They are guaranteed REPEATS.

#define CGB_UNIQUE_CUTOFF            10.0f
//#define CGB_UNIQUE_CUTOFF            12.0f
// A threshold value for Gene^s coverage statistic. Values above this value
// have never been known to associated with unitigs with fragments that are
// not contiguous in the genome.

#define CGB_TANDEM_REPEAT_THRESHOLD  50
// A threshold distance in base pairs
// for the allowed slop between the minimum overlap and maximum overlap before
// calling an overlap definately a tandem repeat.


//  AS_READ_MAX_NORMAL_LEN_BITS must be at least 11.  10 might work, 9 will not compile.
//  AS_READ_MAX_NORMAL_LEN_BITS must be at most 20.  21 will not compile.
//
//  Be aware that overlapper does not work well at large sizes.  Consensus allocates excessive
//  amounts of memory too.
//
#define AS_READ_MAX_PACKED_LEN_BITS       8
#define AS_READ_MAX_NORMAL_LEN_BITS       11

//  AS_READ_MAX_NORMAL_LEN should be a multiple of 8, only to keep things aligned.
//  The total size for storing a single fragment is '24 + n' (as of gkFragment.H r1.8).
//    128 bytes of data allows for length of  104 bases
//    160 bytes of data allows for length of  136 bases
//    192 bytes of data allows for length of  168 bases
//
//  Compare against a gkNormalFragment of fixed size 48 bytes + 10 bits per base.
//    128 bytes of data allows for length of  64 bases
//    160 bytes of data allows for length of  89 bases
//    192 bytes of data allows for length of 115 bases
//  The formula (gkFragment.H r1.8) is 48 + 10/8 * n.
//
//  The catch is that gkPackedFragment always allocates this much space, regardless of the actual
//  length of a fragment, where gkNormalFragment only allocates as much space as used by the
//  sequence.
//
#define AS_READ_MAX_PACKED_LEN            136
#define AS_READ_MAX_NORMAL_LEN            ((1 << AS_READ_MAX_NORMAL_LEN_BITS) - 1)

//  AS_OVL controls both overlapper and bubble popping

extern double AS_OVL_ERROR_RATE;
extern double AS_CGW_ERROR_RATE;  //  former CGW_DP_ERATE
extern double AS_CNS_ERROR_RATE;  //  former CNS_DP_ERATE
extern double AS_MAX_ERROR_RATE;

#define AS_READ_ERROR_RATE         AS_OVL_ERROR_RATE
#define AS_GUIDE_ERROR_RATE        AS_OVL_ERROR_RATE

extern uint32 AS_READ_MIN_LEN;
extern uint32 AS_OVERLAP_MIN_LEN;


int AS_configure(int argc, char **argv);



// These macros are use to eliminate inter-platform differnces between
// calculated results
#define DBL_TO_INT(X)   ((int)((1.0+16.0*DBL_EPSILON)*(X)))
#define ROUNDPOS(X)     (DBL_TO_INT((X)+0.5) )
#define ROUND(X)        (((X)>0.0) ? ROUNDPOS(X) : -ROUNDPOS(-(X)) )
#define ZERO_PLUS       ( 16.0*DBL_EPSILON)
#define ZERO_MINUS      (-16.0*DBL_EPSILON)
#define ONE_PLUS        (1.0+ZERO_PLUS)
#define ONE_MINUS       (1.0+ZERO_MINUS)
#define INT_EQ_DBL(I,D) (fabs((double)(I)-(D)) < 16.0*DBL_EPSILON  )
#define DBL_EQ_DBL(A,B) (fabs((A)-(B))<16.0*DBL_EPSILON)

// cgw and cns use NULLINDEX for a NULL index value
#define NULLINDEX (-1)

// A convenient assert for testing whether ptrs are null
// without bothering lint
#define AssertPtr(ptr) (assert((ptr) != NULL))

#define exitFailure() \
	exit(EXIT_FAILURE)

#define exitSuccess() \
	exit(EXIT_SUCCESS)

static const char LINE_FEED = '\n';
static const char* LINE_FEED_STR = "\n";

static const char CARRIAGE_RETURN = '\r';
static const char* CARRIAGE_RETURN_STR = "\r";

static const char NEWLINE = LINE_FEED;
static const char* NEWLINE_STR = LINE_FEED_STR;

static const char NULL_TERMINATOR = '\0';
static const char* NULL_TERMINATOR_STR = "\0";

static const char PATH_DELIMITER = '/';
static const char* PATH_DELIMITER_STR = "/";

#endif
