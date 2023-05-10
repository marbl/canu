
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#ifndef ADJUST_OVERLAPS_H
#define ADJUST_OVERLAPS_H

#include "runtime.H"

#include "sqStore.H"
#include "ovStore.H"

#include "clearRangeFile.H"


bool
adjustNormal(clearRangeFile  *iniClr,
             ovOverlap       *ovl,
             uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
             uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend);

bool
adjustFlipped(clearRangeFile  *iniClr,
              ovOverlap       *ovl,
              uint32 &aovlbgn,  uint32 &aovlend,  uint32 &bovlbgn,  uint32 &bovlend,
              uint32 &aclrbgn,  uint32 &aclrend,  uint32 &bclrbgn,  uint32 &bclrend,
              sqStore         *seq);


#endif
