
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

#include <omp.h>


#include "types.H"


int AS_configure(int argc, char **argv);


//  Tell gcc (and others, maybe) about unused parameters.  This is important for gcc (especially
//  newer ones) that complain about unused parameters.  Thanks to ideasman42 at
//  http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code.
#ifdef __GNUC__
#define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#define UNUSED(x) UNUSED_ ## x
#endif


#include "arrays.H"


#endif
