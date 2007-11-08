
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

//  $Id: AS_global.h,v 1.24 2007-11-08 12:38:11 brianwalenz Exp $

//  This is the global include file that all C files in the AS
//  subsystem should include.

#ifndef AS_GLOBAL_H
#define AS_GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include <assert.h>

#include <limits.h>
#include <float.h>
#include <inttypes.h>
#include <time.h>

#include "AS_UTL_alloc.h"


#ifdef __alpha
// used in AS_CNS/Array_CNS.h and AS_CNS/MultiAlignment_CNS.h
#include <random.h>
#endif

#ifndef TRUE
  #define TRUE (1)
#endif
#ifndef FALSE
  #define FALSE (0)
#endif


#ifndef MIN
  #define MIN(a,b)		( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
  #define MAX(a,b)		( ((a) > (b)) ? (a) : (b) )
#endif


#ifdef __cplusplus
  #include <cstdlib>
#else
  #include <stdlib.h>
#endif

#ifndef __cplusplus
  #ifdef bool
    #undef bool
  #endif // bool
  typedef int bool;
#endif // cplusplus


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

typedef void *PtrT;

#define F_L     "%ld"
#define F_UL    "%lu"
#define F_LL   "%lld"
#define F_ULL  "%llu"

#if ULONG_MAX == 0xffffffff
  // 32-bit architecture

  #define TRUE32BIT

  typedef uint32   INTPTR;

  #ifndef UINT64_MAX
    #define INT64_MAX  LLONG_MAX
    #define UINT64_MAX ULLONG_MAX
  #endif

  #define F_S16     "%d"
  #define F_U16     "%u"
  #define F_S32     "%d"
  #define F_S32P     "d"
  #define F_U32     "%u"
  #define F_U32P     "u"
  #define F_S64   "%lld"
  #define F_S64P   "lld"
  #define F_U64   "%llu"
  #define F_U64P   "llu"
  #define F_X64   "%llx"
  #define F_X64P   "llx"

  #define F_SIZE_T   "%u"
  #define F_SIZE_TP   "u"

  #define F_TIME_T  "%ld"
  #define F_TIME_TP  "ld"

  #define F_PID_T    "%d"
  #define F_PID_TP    "d"

  #if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 32)
    // off_t is 32-bit
    #error I do not support 32-bit off_t.
    #define F_OFF_T   "%ld"
    #define F_OFF_TP   "ld"
  #endif

  // off_t is 64-bit

  #define F_OFF_T   "%lld"
  #define F_OFF_TP   "lld"

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

  typedef uint64   INTPTR;

  #ifndef UINT64_MAX
    #define INT64_MAX  LONG_MAX
    #define UINT64_MAX ULONG_MAX
  #endif
 
  // included here for downstream .c files
  // I don't know if this is needed anywhere else but on the alphas (MP)
  #ifdef _OSF_SOURCE
   #include <sys/mode.h>
  #endif
    
  #define F_S16    "%d"
  #define F_U16    "%u"
  #define F_S32    "%d"
  #define F_S32P    "d"
  #define F_U32    "%u"
  #define F_U32P    "u"
  #define F_S64   "%ld"
  #define F_S64P   "ld"
  #define F_U64   "%lu"
  #define F_U64P   "lu"
  #define F_X64   "%lx"
  #define F_X64P   "lx"
    
  #define F_SIZE_T  "%lu"
  #define F_SIZE_TP  "lu"

  #ifdef _AIX
    #define F_TIME_T  "%ld"
    #define F_TIME_TP  "ld"

    #define F_PID_T   "%ld"
    #define F_PID_TP   "ld"

    #ifdef _LARGE_FILES
      #define F_OFF_T   "%lld"
      #define F_OFF_TP   "lld"
    #else
      #define F_OFF_T   "%ld"
      #define F_OFF_TP   "ld"
    #endif
  #else
    // these are valid for __alpha, perhaps not for others...
    #define F_TIME_T  "%d"
    #define F_TIME_TP  "d"

    #define F_PID_T   "%d"
    #define F_PID_TP   "d"

    #ifdef _KERNEL
      #define F_OFF_T   "%lu"
      #define F_OFF_TP   "lu"
    #else
      #define F_OFF_T   "%ld"
      #define F_OFF_TP   "ld"
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



#include "AS_UTL_IID.h"
#include "AS_UTL_UID.h"





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


#define AS_FRAG_MAX_LEN 2048
#define AS_FRAG_MIN_LEN 64

#define AS_READ_MAX_LEN AS_FRAG_MAX_LEN
#define AS_READ_MIN_LEN AS_FRAG_MIN_LEN

//  AS_OVL controls both overlapper and bubble popping

extern double AS_OVL_ERROR_RATE;
extern double AS_CGW_ERROR_RATE;  //  former CGW_DP_ERATE
extern double AS_CNS_ERROR_RATE;  //  former CNS_DP_ERATE
extern double AS_MAX_ERROR_RATE;

#define AS_READ_ERROR_RATE         AS_OVL_ERROR_RATE
#define AS_GUIDE_ERROR_RATE        AS_OVL_ERROR_RATE

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


#endif
