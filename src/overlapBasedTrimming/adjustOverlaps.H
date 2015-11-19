
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef ADJUST_OVERLAPS_H
#define ADJUST_OVERLAPS_H

#include "AS_global.H"

#include "gkStore.H"
#include "ovStore.H"

#include "clearRangeFile.H"


bool
adjustNormal(clearRangeFile  *iniClr,
             gkStore         *gkp,
             ovOverlap       *ovl,
             uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
             uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend);

bool
adjustFlipped(clearRangeFile  *iniClr,
              gkStore         *gkp,
              ovOverlap       *ovl,
              uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
              uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend);


#endif
