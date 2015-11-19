
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
 *      are Copyright 2005-2006,2008,2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz on 2014-DEC-05
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef UTIL_PLUS_PLUS_H
#define UTIL_PLUS_PLUS_H

#include "util.h"

//  These are all inlined, and C doesn't want to listen to that, so
//  they're here.
//
#include "bitOperations.h"
#include "bitPacking.h"
#include "endianess.H"

//  Various methods for encoding numbers into a bitstream.
//
//  Still missing:
//    minimal binary
//    golomb (actually rice, since power of two)
//    teuhola exponential golomb
//
//  And a nice way of getting parameters to those (and generalizedUnary)
//
#include "unaryEncoding.h"
#include "generalizedUnaryEncoding.h"
#include "eliasGammaEncoding.h"
#include "eliasDeltaEncoding.h"
#include "fibonacciEncoding.h"

//  Lists?
#include "uint32List.H"

//  Now the good stuff!
//
#include "speedCounter.H"
//#include "bzipBuffer.H"
#include "readBuffer.H"
#include "splitToWords.H"
#include "bitPackedArray.H"
#include "bitPackedFile.H"
#include "recordFile.H"
#include "intervalList.H"
#include "sweatShop.H"
#include "logMsg.H"

#endif  //  UTIL_PLUS_PLUS_H

