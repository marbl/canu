/**************************************************************************
 * This file is part of A2Amapper.
 * Copyright (c) 2004 Applera Corporation
 * Author: Clark Mobarry
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
**************************************************************************/


/*************************************************
 Module: cds.h

 Description:
 This is a stub file for the system wide "cds.h" file.
 
 Assumptions:
 
*************************************************/

#ifndef CDSH
#define CDSH
#include <limits.h>

/*
# include <macros.h>
*/
#ifndef min
#	define min(a,b)		( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
#	define max(a,b)		( ((a) > (b)) ? (a) : (b) )
#endif

#ifdef __cplusplus
	/* we collide with LEDA 4.3, which uses abs as a
	 * function parameter somewhere, so we can't use 
	 * the macro definition. Macros are evil anyway!
	 * We can't use a template definition either, sine
	 * this header file is included somewhere in an 
	 * extern "C" block...
	 * OK / 10-2001
	 */

	/* inline int abs(int x) { return ((x >= 0) ? x : -x); } */

	/* even worse, this explicit C++ definition conflicts with the 
	 * extern "C" declaration in stdlib.h and definition in libc.a. 
	 * sigh. get abs(int x) from stdlib.h instead.
	 * Nathan Edwards 11-2001 
	 */

#       include <stdlib.h>

#else

#	ifndef i386
#		ifndef abs
#			define abs(x)			(x>= 0 ? (x) : -(x))
#		endif
#	endif

#endif

/* 
   The following works for:
   (1) Digital cc the alphas,
   (2) GNUC the alphas, and
   (3) GNUC on a Pentium.
   Note that the alphas and Pentiums differ on the sizeof(long int).
   _LONGLONG
*/

/* NEW STYLE - TO BE PHASED IN */
typedef              char      cds_int8;
typedef              short int cds_int16;
typedef                    int cds_int32;
typedef unsigned     char      cds_uint8;
typedef unsigned     short int cds_uint16;
typedef unsigned           int cds_uint32;
typedef                  float cds_float32;
typedef                 double cds_float64;
#define CDS_INT8_MAX   CHAR_MAX
#define CDS_INT8_MIN   CHAR_MIN
#define CDS_INT16_MAX  SHRT_MAX
#define CDS_INT16_MIN  SHRT_MIN
#define CDS_INT32_MAX  INT_MAX
#define CDS_INT32_MIN  INT_MIN

#define CDS_UINT8_MAX  UCHAR_MAX
#define CDS_UINT16_MAX USHRT_MAX
#define CDS_UINT32_MAX UINT_MAX

#include <inttypes.h>

typedef int8_t  int8;
typedef int16_t int16;
typedef int32_t int32;
typedef int64_t int64;
typedef uint8_t  uint8;
typedef uint16_t uint16;
typedef uint32_t uint32;
typedef uint64_t uint64;
typedef                  float float32;
typedef                 double float64;

#if ! (defined( _AIX) || defined(__CYGWIN__))
#define INT8_MAX       CHAR_MAX
#define INT8_MIN       CHAR_MIN
#define INT16_MAX      SHRT_MAX
#define INT16_MIN      SHRT_MIN
#define INT32_MAX      INT_MAX
#define INT32_MIN      INT_MIN

#define UINT8_MAX      UCHAR_MAX
#define UINT16_MAX     USHRT_MAX
#define UINT32_MAX     UINT_MAX
#endif // _AIX

#ifndef  i386
/* OLD STYLE - TO BE PHASED OUT */
#ifndef _AIX
#define UINT64_MAX     ULONG_MAX
#define INT64_MAX      LONG_MAX
#define INT64_MIN      LONG_MIN
typedef          long int int64;
#endif // _AIX
typedef unsigned long int uint64;

/* NEW STYLE - TO BE PHASED IN */
#define CDS_UINT64_MAX ULONG_MAX
#define CDS_INT64_MAX  LONG_MAX
#define CDS_INT64_MIN  LONG_MIN
typedef          long int cds_int64;
typedef unsigned long int cds_uint64;
#else  /*i386*/
/* OLD STYLE - TO BE PHASED OUT */
#if ! (defined(__CYGWIN__))
#define UINT64_MAX  18446744073709551615UL
#define INT64_MAX    9223372036854775807L
#endif
typedef          long long int int64;
typedef unsigned long long int uint64;

/* NEW STYLE - TO BE PHASED IN */
#define CDS_UINT64_MAX  18446744073709551615UL
#define CDS_INT64_MAX    9223372036854775807L
typedef          long long int cds_int64;
typedef unsigned long long int cds_uint64;
#endif /*i386*/

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#endif /*CDSH*/



