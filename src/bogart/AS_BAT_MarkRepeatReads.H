
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
 *    Brian P. Walenz beginning on 2016-MAR-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_MARKREPEATREADS
#define INCLUDE_AS_BAT_MARKREPEATREADS

#include "AS_BAT_TigVector.H"



class confusedEdge {
public:
  confusedEdge(uint32 aid_, bool a3p_, uint32 bid_) {
    aid = aid_;
    a3p = a3p_;
    bid = bid_;
  };

  uint32  aid;  //  Read that is confused and needs to be split.
  bool    a3p;  //  End of read that is confused.

  uint32  bid;  //  Read that causes confusion.
};



void
markRepeatReads(AssemblyGraph         *AG,
                TigVector             &tigs,
                double                 deviationRepeat,
                uint32                 confusedAbsolute,
                double                 confusedPercent,
                vector<confusedEdge>  &confusedEdges);



#endif  //  INCLUDE_AS_BAT_MARKREPEATREADS
