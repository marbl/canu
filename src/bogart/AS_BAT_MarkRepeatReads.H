
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

#ifndef INCLUDE_AS_BAT_MARKREPEATREADS
#define INCLUDE_AS_BAT_MARKREPEATREADS

#include "AS_BAT_TigVector.H"
#include "AS_BAT_AssemblyGraph.H"


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
