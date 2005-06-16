
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
#ifndef CDSH
#define CDSH

#include <limits.h>
#include <float.h>
#include <inttypes.h>
#include <time.h>

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


#ifndef min
  #define min(a,b)		( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
  #define max(a,b)		( ((a) > (b)) ? (a) : (b) )
#endif

#ifdef __cplusplus
  #include <cstdlib>
#else
  #include <stdlib.h>
#endif

#ifdef __cplusplus
#else
  #ifdef bool
    #undef bool
  #endif // bool
  typedef int bool;
#endif // cplusplus

#define CDS_INT8_MAX   CHAR_MAX
#define CDS_INT8_MIN   CHAR_MIN
#define CDS_INT16_MAX  SHRT_MAX
#define CDS_INT16_MIN  SHRT_MIN
#define CDS_INT32_MAX  INT_MAX
#define CDS_INT32_MIN  INT_MIN

#define CDS_UINT8_MAX  UCHAR_MAX
#define CDS_UINT16_MAX USHRT_MAX
#define CDS_UINT32_MAX UINT_MAX

typedef int8_t  cds_int8;
typedef int16_t cds_int16;
typedef int32_t cds_int32;
typedef int64_t cds_int64;
typedef uint8_t  cds_uint8;
typedef uint16_t cds_uint16;
typedef uint32_t cds_uint32;
typedef uint64_t cds_uint64;

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

typedef float  cds_float32;
typedef double cds_float64;
typedef float  float32;
typedef double float64;

#define F_L     "%ld"
#define F_UL    "%lu"
#define F_LL   "%lld"
#define F_ULL  "%llu"

#if ULONG_MAX == 0xffffffff
  /***********************************************************************
   32-bit architecture
   ***********************************************************************/

  #define CDS_INT64_MAX INT64_MAX
  #define CDS_UINT64_MAX UINT64_MAX

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

  #define F_SIZE_T   "%d"
  #define F_SIZE_TP   "d"

  #define F_TIME_T  "%ld"
  #define F_TIME_TP  "ld"

  #define F_PID_T    "%d"
  #define F_PID_TP    "d"

  #if _FILE_OFFSET_BITS==64
    // off_t is 64-bit
    #define F_OFF_T   "%lld"
    #define F_OFF_TP   "lld"
  #else
    // off_t is 32-bin
    #define F_OFF_T   "%ld"
    #define F_OFF_TP   "ld"
  #endif

  #define STR_TO_UID     strtoull
  #define STR_TO_UINT64  strtoull
  #define STR_TO_INT64    strtoll

  #define FILEID_MASK       0xffff000000000000ULL
  #define FILEOFFSET_MASK   0x0000ffffffffffffULL
  #define LOCALE_OFFSET            10000000000ULL  // 10^10 > 2^32

  #ifdef _LARGEFILE_SOURCE
    // ftello & fseeko available & preferable
    #define CDS_FTELL ftello
    #define CDS_FSEEK fseeko
  #else
    #ifdef _LARGEFILE64_SOURCE
      // ftello64 & fseeko64 available & needed?
      #define CDS_FTELL ftello64
      #define CDS_FSEEK fseeko64
    #else
      // use standard ftell & fseek
      #define CDS_FTELL ftell
      #define CDS_FSEEK fseek
    #endif
  #endif

  //========== SIMULATOR-specific
  // Fix initial creation time & uid in celsim to support regression testing 
  // time stamp of batch message in first dros file...
  #define CREATION_TIME      915170460L
  // UID of first frag in first dros input file...
  #define FIRST_UID     17000001585819ULL

  //========== MERYL-specific
  #define CDS_UINT64_ZERO   0x0000000000000000ULL
  #define CDS_UINT64_ONE    0x0000000000000001ULL

  //========== CGB-specific
  #define CGB_MULTIPLIER               1000

  //========== CGW-specific
  // Used in cgw: 1024^3
  #define MAX_SEQUENCEDB_SIZE        1073741824UL
  // Used in cgw: 256 * 1024^2
  #define MAX_SEQUENCEDB_CACHE_SIZE   268435456UL

#else
  /***********************************************************************
   Not 32-bit architecture
   ***********************************************************************/


  #if ULONG_MAX == 0xffffffffffffffff
    /***********************************************************************
     64-bit architecture
     ***********************************************************************/

    #define CDS_INT64_MAX LONG_MAX
    #define CDS_UINT64_MAX ULONG_MAX
    
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
    
    #define STR_TO_UID     strtoul
    #define STR_TO_UINT64  strtoul
    #define STR_TO_INT64    strtol
    
    #define FILEID_MASK       0xffff000000000000UL
    #define FILEOFFSET_MASK   0x0000ffffffffffffUL
    #define LOCALE_OFFSET            10000000000UL  // 10^10 > 2^32
    
    #define CDS_FTELL ftell
    #define CDS_FSEEK fseek
    
    //========== SIMULATOR-specific
    // Fix initial creation time & uid in celsim to support regression testing 
    // time stamp of batch message in first dros file...
    #define CREATION_TIME      915170460
    // UID of first frag in first dros input file...
    #define FIRST_UID     17000001585819UL
    
    //========== MERYL-specific
    #define CDS_UINT64_ZERO   0x0000000000000000UL
    #define CDS_UINT64_ONE    0x0000000000000001UL
    
    //========== CGB-specific
    #ifndef __alpha
      #define CGB_MULTIPLIER            1000
    #else
      #define CGB_MULTIPLIER               1000
    #endif
    
    //========== CGW-specific
    // Used in cgw: 2 * 1024^3
    #define MAX_SEQUENCEDB_SIZE           2147483648ul

    // Used in cgw: 2 * 256 * 1024^2
    #define MAX_SEQUENCEDB_CACHE_SIZE      536870912ul
  #else
    /***********************************************************************
     not 32-bit nor 64-bit architecture
     ***********************************************************************/
    CANNOT SET REQUIRED ASSEMBLER MACROS in cds.h. UNIDENTIFIED ARCHITECTURE.
  #endif
#endif

typedef cds_uint64 CDS_UID_t;
typedef cds_uint32 CDS_IID_t;
typedef cds_int32  CDS_CID_t;
typedef cds_int32  CDS_COORD_t;
#define CDS_IID_MAX     CDS_UINT32_MAX
#define CDS_CID_MAX     CDS_INT32_MAX
#define CDS_COORD_MIN   CDS_INT32_MIN
#define CDS_COORD_MAX   CDS_INT32_MAX

// used in AS_MER
#define CDS_UINT64_MASK(X)   ((~CDS_UINT64_ZERO) >> (64 - (X)))

#define F_UID    F_U64
#define F_UIDP   F_U64P
#define F_IID    F_U32
#define F_IIDP   F_U32P
#define F_CID    F_S32
#define F_CIDP   F_S32P
#define F_COORD  F_S32
#define F_COORDP F_S32P

#endif
