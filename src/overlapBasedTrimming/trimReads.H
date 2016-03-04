
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
 *    Brian P. Walenz on 2015-JUN-16
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef TRIM_READS_H
#define TRIM_READS_H

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "intervalList.H"


#define OBT_MODE_WIGGLE      (5)

#define OBT_CQ_LENGTH        (100)
#define OBT_CQO_LENGTH       (200)
#define OBT_CQO_OVERLAP      (100)
#define OBT_CQ_SHORT         (5)

#define OBT_QLT_CLOSE_5      (10)  // 5,6  use 5'mode, use 5'min>1
#define OBT_QLT_FAR_5        (50)  // 11   use min5'
#define OBT_QLT_MODE3        (150) // 9    use 3'mode
#define OBT_QLT_CLOSE_MAXM3  (30)  // 14   use max>1 close to max
#define OBT_QLT_CLOSE_MAX3   (100) // 12   use max3'


bool
largestCovered(ovOverlap   *ovl,
               uint32       ovlLen,
               gkRead      *read,
               uint32       ibgn,
               uint32       iend,
               uint32      &fbgn,
               uint32      &fend,
               char        *logMsg,
               uint32       errorRate,
               uint32       minOverlap,
               uint32       minCoverage,
               uint32       minReadLength);

bool
bestEdge(ovOverlap   *ovl,
         uint32       ovlLen,
         gkRead      *read,
         uint32       ibgn,
         uint32       iend,
         uint32      &fbgn,
         uint32      &fend,
         char        *logMsg,
         uint32       errorRate,
         uint32       minOverlap,
         uint32       minCoverage,
         uint32       minReadLength);

#endif  //  TRIM_READS_H
