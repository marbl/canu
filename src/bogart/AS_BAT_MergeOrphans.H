
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
 *  This file is derived from:
 *
 *    src/bogart/AS_BAT_PopBubbles.H
 *
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-DEC-07
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_MERGEORPHANS
#define INCLUDE_AS_BAT_MERGEORPHANS

#include "AS_global.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_Unitig.H"

void
mergeOrphans(TigVector &tigs,
             double     deviationBubble);


#endif  //  INCLUDE_AS_BAT_MERGEORPHANS
