
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
 *    Brian P. Walenz beginning on 2016-OCT-03
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_TIGGRAPH
#define INCLUDE_AS_BAT_TIGGRAPH

#include "runtime.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"

#include "AS_BAT_CreateUnitigs.H"

void
reportTigGraph(TigVector      &tigs,
               vector<tigLoc> &tigSource,
               const char     *prefix,
               const char     *label);


#endif  //  INCLUDE_AS_BAT_ASSEMBLYGRAPH
