
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

static const char *rcsid_AS_GLOBAL_H = "$Id: AS_global.h,v 1.57 2013-02-22 12:18:52 brianwalenz Exp $";

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

//  Check for GCC parallel STL support.  This appeared in version 4.3, but the lowest we've tested
//  with is 4.6.
//
#if ((__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ > 2)))
//  Not by default!
//#ifndef _GLIBCXX_PARALLEL
//#define _GLIBCXX_PARALLEL
//#endif
#else
#undef  _GLIBCXX_PARALLEL
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>

#include <float.h>
#include <math.h>

#include <assert.h>
#include <errno.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>

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

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;

typedef void* PtrT;

#ifndef UINT32_MAX
#error UINT32_MAX not defined; if new code, AS_global.h must be the first include
#endif

// Pointers
#define F_P           "%p"
#define F_PP           "p"

// Characters
#define F_C           "%c"
#define F_CP           "c"
#define F_CI         "%*c"

// Strings
#define F_STR         "%s"
#define F_STRP         "s"
#define F_STRI       "%*s"

// Integers
#define F_S16    "%"PRId16
#define F_S16P      PRId16
#define F_S16I  "%*"PRId16
#define F_U16    "%"PRIu16
#define F_U16P      PRIu16
#define F_U16I  "%*"PRIu16
#define F_S32    "%"PRId32
#define F_S32P      PRId32
#define F_S32I  "%*"PRId32
#define F_U32    "%"PRIu32
#define F_U32P      PRIu32
#define F_U32I  "%*"PRIu32
#define F_S64    "%"PRId64
#define F_S64P      PRId64
#define F_S64I  "%*"PRId64
#define F_U64    "%"PRIu64
#define F_U64P      PRIu64
#define F_U64I  "%*"PRIu64
#define F_X64    "%"PRIx64
#define F_X64P      PRIx64
#define F_X64I  "%*"PRIx64

// Floating points
#define F_F32         "%f"
#define F_F32P         "f"
#define F_F32I       "%*f"
#define F_F64        "%lf"
#define F_F64P        "lf"
#define F_F64I      "%*lf"

// Standard typedefs
#define F_SIZE_T     "%zd"
#define F_SIZE_TP     "zd"
#define F_SIZE_TI   "%*zd"

#define F_OFF_T     F_S64
#define F_OFF_TP    F_S64P
#define F_OFF_TI    F_S64I

typedef uintptr_t INTPTR;


//  These are used to pad various structs to specific sizes
#if ULONG_MAX == 0xffffffff
#define TRUE32BIT
#else
#define TRUE64BIT
#endif


#if defined(_FILE_OFFSET_BITS) && (_FILE_OFFSET_BITS == 32)
#error I do not support 32-bit off_t.
#endif


//  Enable troublesome asserts.  These typically have work arounds, and occasionally trigger.
//  They're really only useful if the assembly can be debugged.
#undef AGGRESSIVE_ASSERT


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


#define CGB_MULTIPLIER            1000

#define CGB_INVALID_CUTOFF           -12.0f
// A threshold value for Gene's coverage statistic. Values BELOW this value
// have never been known to associated with unitigs with fragments that ARE
// not contiguous in the genome. They are guaranteed REPEATS.

#define CGB_UNIQUE_CUTOFF            10.0f
//#define CGB_UNIQUE_CUTOFF            12.0f
// A threshold value for Gene's coverage statistic. Values above this value
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

//  Do not define this.
#undef OVS_DANGEROUSLY_OVERSIZE

#define AS_READ_MAX_PACKED_LEN            ((1 << AS_READ_MAX_PACKED_LEN_BITS) - 1)
#define AS_READ_MAX_NORMAL_LEN            ((1 << AS_READ_MAX_NORMAL_LEN_BITS) - 1)

//  AS_OVL controls both overlapper and bubble popping

extern double AS_OVL_ERROR_RATE;
extern double AS_CGW_ERROR_RATE;
extern double AS_CNS_ERROR_RATE;
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

#endif
