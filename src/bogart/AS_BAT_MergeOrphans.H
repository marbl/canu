
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

#ifndef INCLUDE_AS_BAT_MERGEORPHANS
#define INCLUDE_AS_BAT_MERGEORPHANS

#include "runtime.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"

void
mergeOrphans(TigVector    &tigs,
             double        deviation,
             double        similarity,
             bool          isBubble);


#endif  //  INCLUDE_AS_BAT_MERGEORPHANS
