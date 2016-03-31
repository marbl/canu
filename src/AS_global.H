
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
 *    Brian P. Walenz from 2006-OCT-08 to 2014-MAR-05
 *      are Copyright 2006-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Eli Venter on 2006-NOV-14
 *      are Copyright 2006 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren from 2011-MAY-02 to 2011-OCT-11
 *      are Copyright 2011 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz from 2014-AUG-11 to 2015-MAR-31
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *    Sergey Koren beginning on 2016-FEB-15
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

//  This is the global include file that all C files in the AS
//  subsystem should include.

#ifndef AS_GLOBAL_H
#define AS_GLOBAL_H

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

//typedef void* PtrT;

//  Correct?  Or should it be LLU and LU?
#define  uint64NUMBER(X) X ## LU
#define  uint32NUMBER(X) X ## U

#define  uint64ZERO      uint64NUMBER(0x0000000000000000)
#define  uint64ONE       uint64NUMBER(0x0000000000000001)
#define  uint64MAX       uint64NUMBER(0xffffffffffffffff)
#define  uint64MASK(X)   ((~uint64ZERO) >> (64 - (X)))

#define  uint32ZERO      uint32NUMBER(0x00000000)
#define  uint32ONE       uint32NUMBER(0x00000001)
#define  uint32MAX       uint32NUMBER(0xffffffff)
#define  uint32MASK(X)   ((~uint32ZERO) >> (32 - (X)))

#define  uint16ZERO      (0x0000)
#define  uint16ONE       (0x0001)
#define  uint16MAX       (0xffff)
#define  uint16MASK(X)   ((~uint16ZERO) >> (16 - (X)))

#define  uint8ZERO       (0x00)
#define  uint8ONE        (0x01)
#define  uint8MAX        (0xff)
#define  uint8MASK(X)    ((~uint8ZERO) >> (8 - (X)))


#define  strtouint32(N)   (uint32)strtoul (N, NULL, 10)
#define  strtouint64(N)   (uint64)strtoull(N, NULL, 10)
#define  strtodouble(N)   (double)strtod  (N, NULL)

//inline uint32  strtouint32(char *N, O) (uint32)strtoul (N, O, 10)
//inline uint64  strtouint64(char *N, O) (uint64)strtoull(N, O, 10)
//inline double  strtodouble(char *N, O) (double)strtod  (N, O)


// Pointers
#define F_PTR    "0x%016p"

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
#define F_X64 "%016"PRIx64
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
#define F_SIZE_T     "%zu"
#define F_SIZE_TP     "zu"
#define F_SIZE_TI   "%*zu"

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



//  Tell gcc (and others, maybe) about unused parameters.  This is important for gcc (especially
//  newer ones) that complain about unused parameters.  Thanks to ideasman42 at
//  http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code.
#ifdef __GNUC__
#define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#define UNUSED(x) UNUSED_ ## x
#endif


#include "AS_UTL_alloc.H"



//  Stubs for ignoring OpenMP on clang.
#ifdef BROKEN_CLANG_OpenMP

static
void
omp_set_dynamic(int x) {
    #pragma unused(x)
}

static
int
omp_get_max_threads(void) {
  return(1);
}

static
void
omp_set_num_threads(int x) {
    #pragma unused(x)
}

static
int
omp_get_num_threads(void) {
  return(1);
}

static
int
omp_get_thread_num(void) {
  return(0);
}

typedef int omp_lock_t;

static
void
omp_init_lock(omp_lock_t *a) {
    #pragma unused(a)
}

static
void
omp_set_lock(omp_lock_t *a) {
    #pragma unused(a)
}

static
void
omp_unset_lock(omp_lock_t *a) {
    #pragma unused(a)
}

static
void
omp_destroy_lock(omp_lock_t *a) {
    #pragma unused(a)
}

#endif



#endif
